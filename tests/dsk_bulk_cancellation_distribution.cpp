/*
Whether cancelling machine orders in bulk draws the counts the one-machine-at-a-
time loop in LABOR() draws.

model/dsk_sfc/dsk_sfc_bulk_cancellation.h replaces that loop once the number of
machines to cancel passes a threshold no 600-period run reaches. It takes
different random draws, so it cannot reproduce the loop's run byte for byte;
what it has to reproduce is the distribution of how many cancellations each
customer receives. That distribution can be worked out exactly for small cases,
so this compares both samplers against it rather than against each other.

The exact distribution comes from enumerating the loop's own process: every
cancellation falls on a customer chosen uniformly among those whose order is not
used up. The loop is simulated as written in modules/module_macro_sfc.cpp, with
customers placed among N2 = 200 firm slots and draws that land on a
non-customer or a used-up order thrown away, so the enumeration is also checked
against the code it describes.

Each outcome's observed frequency over R replications is compared with its
exact probability p, and the gap must be below 5 standard deviations of a
binomial proportion, sqrt(p (1 - p) / R). An outcome with p below 1e-4 is left
out of the comparison, since at these R it is seen a handful of times at most.
With about 400 comparisons across the cases the chance of a false alarm is
about 2e-4, and the seeds are fixed, so the result is the same on every run.

A test that only ever passes proves nothing, so a deliberately wrong sampler is
put through the same comparison and must fail it: one that shares the
cancellations out in proportion to each customer's capacity, which gives the same
mean total and the wrong split.

Cases:
  small     capacities 1, 2 and 4, 5 to cancel. bnldev's direct method, and
            capping in most replications. Compared on the joint distribution
            of all three counts.
  capped    capacities 30, 45 and 300, 100 to cancel. bnldev's rejection method,
            and the smallest order is used up in a large share of replications.
            Compared on the joint distribution of the two small customers'
            counts.
  many      40 customers of capacity 1000 each, 30 to cancel. Nothing can be
            capped, so each count is binomial with n = 30 and p = 1/40, which
            is the exact distribution used. The range where bnldev would give a
            Poisson draw instead, and the bulk draw inverts the binomial.
            Compared on customer 0's count.
  huge      capacities 2e9, 2e9 and 3e9, 5e9 to cancel, which takes bnldev past
            its int trial count. Too large to enumerate or to run the loop on.
            Every count must be whole, within its capacity and sum to the
            total, and each customer's mean count must be within 5 standard
            deviations of 5e9 / 3: the capacities are more than 10,000 standard
            deviations above that, so they do not bind and the uncapped
            multinomial is exact.

Writes out/dsk_bulk_cancellation_distribution.txt.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "model/dsk_sfc/dsk_sfc_bulk_cancellation.h"
#include "model/dsk_sfc/auxiliary/ran1.h"

#define REPORT "out/dsk_bulk_cancellation_distribution.txt"
#define N1 20
#define N2 200
#define Z_LIMIT 5.0
#define SMALLEST_COMPARED_PROBABILITY 1e-4

typedef std::vector<long> Counts;
typedef std::map<Counts, double> Distribution;

/* The loop in LABOR(), with the customers' orders as capacities in machines. */
static void loop_sampler(const std::vector<double> &capacity, double to_cancel, double *cancelled, long *idum) {
    const int customers = (int)capacity.size();
    std::vector<int> customer_at_slot(N2 + 1, -1);
    for (int k = 0; k < customers; k++) {
        customer_at_slot[1 + (k * 37) % N2] = k;
        cancelled[k] = 0;
    }
    double reduction = to_cancel;
    while (reduction > 0) {
        int slot = int(ran1(idum) * N1 * N2) % N2 + 1;
        int k = customer_at_slot[slot];
        if (k >= 0 && cancelled[k] < capacity[k]) {
            cancelled[k] += 1;
            reduction -= 1;
        }
    }
}

static void bulk_sampler(const std::vector<double> &capacity, double to_cancel, double *cancelled, long *idum) {
    draw_uniform_cancellations(capacity.data(), (int)capacity.size(), to_cancel, cancelled, idum);
}

/* Shares the cancellations out in proportion to capacity, one at a time. Wrong
   on purpose. */
static void wrong_sampler(const std::vector<double> &capacity, double to_cancel, double *cancelled, long *idum) {
    const int customers = (int)capacity.size();
    for (int k = 0; k < customers; k++) cancelled[k] = 0;
    for (double left = to_cancel; left > 0; left -= 1) {
        double open_capacity = 0;
        for (int k = 0; k < customers; k++) open_capacity += capacity[k] - cancelled[k];
        double u = ran1(idum) * open_capacity;
        int k = 0;
        while (k < customers - 1 && u >= capacity[k] - cancelled[k]) {
            u -= capacity[k] - cancelled[k];
            k++;
        }
        cancelled[k] += 1;
    }
}

/* The exact distribution of the counts, by stepping the loop's process one
   cancellation at a time over every reachable state. Feasible only for a few
   customers. */
static Distribution exact_distribution(const std::vector<double> &capacity, long to_cancel) {
    const int customers = (int)capacity.size();
    std::map<Counts, double> current;
    current[Counts(customers, 0)] = 1.0;
    for (long step = 0; step < to_cancel; step++) {
        std::map<Counts, double> next;
        for (const auto &entry : current) {
            const Counts &state = entry.first;
            int open = 0;
            for (int k = 0; k < customers; k++) open += state[k] < (long)capacity[k];
            for (int k = 0; k < customers; k++) {
                if (state[k] >= (long)capacity[k]) continue;
                Counts moved = state;
                moved[k] += 1;
                next[moved] += entry.second / open;
            }
        }
        current.swap(next);
    }
    return current;
}

/* The count of one customer when no capacity can bind: binomial with
   to_cancel trials and probability one over the number of customers. */
static Distribution uncapped_single_marginal(int customers, long to_cancel) {
    Distribution result;
    const double p = 1.0 / customers;
    for (long k = 0; k <= to_cancel; k++) {
        double log_pmf = std::lgamma((double)to_cancel + 1) - std::lgamma((double)k + 1) -
                         std::lgamma((double)(to_cancel - k) + 1) + k * std::log(p) +
                         (to_cancel - k) * std::log1p(-p);
        result[Counts(1, k)] = std::exp(log_pmf);
    }
    return result;
}

static Distribution marginal(const Distribution &joint, const std::vector<int> &kept) {
    Distribution result;
    for (const auto &entry : joint) {
        Counts key;
        for (int k : kept) key.push_back(entry.first[k]);
        result[key] += entry.second;
    }
    return result;
}

typedef void (*Sampler)(const std::vector<double> &, double, double *, long *);

struct Comparison {
    int compared;
    double worst_z;
    Counts worst_outcome;
};

static Comparison compare(Sampler sampler, const std::vector<double> &capacity, long to_cancel,
                          const std::vector<int> &kept, const Distribution &exact, long replications, long seed) {
    long idum = -seed;
    std::vector<double> cancelled(capacity.size());
    std::map<Counts, long> seen;
    for (long r = 0; r < replications; r++) {
        sampler(capacity, (double)to_cancel, cancelled.data(), &idum);
        Counts key;
        for (int k : kept) key.push_back((long)cancelled[k]);
        seen[key]++;
    }

    Comparison result = {0, 0.0, Counts()};
    for (const auto &entry : exact) {
        double p = entry.second;
        if (p < SMALLEST_COMPARED_PROBABILITY) continue;
        double observed = (double)seen[entry.first] / (double)replications;
        double z = std::fabs(observed - p) / std::sqrt(p * (1 - p) / (double)replications);
        result.compared++;
        if (z > result.worst_z) {
            result.worst_z = z;
            result.worst_outcome = entry.first;
        }
    }
    return result;
}

static std::string outcome_text(const Counts &outcome) {
    std::string text = "(";
    for (size_t k = 0; k < outcome.size(); k++) {
        if (k) text += ", ";
        text += std::to_string(outcome[k]);
    }
    return text + ")";
}

struct Case {
    const char *name;
    std::vector<double> capacity;
    long to_cancel;
    std::vector<int> kept;
    long loop_replications;
    long bulk_replications;
    bool wrong_sampler_differs;
};

int main(void) {
    std::vector<Case> cases = {
        {"small", {1, 2, 4}, 5, {0, 1, 2}, 400000, 400000, true},
        {"capped", {30, 45, 300}, 100, {0, 1}, 40000, 400000, true},
        {"many", std::vector<double>(40, 1000), 30, {0}, 100000, 400000, false},
    };

    FILE *report = std::fopen(REPORT, "w");
    if (!report) {
        std::fprintf(stderr, "dsk_bulk_cancellation_distribution: cannot open %s\n", REPORT);
        return EXIT_FAILURE;
    }
    std::fprintf(report, "Counts of cancellations per customer: the loop in LABOR(), the bulk draw\n"
                         "that replaces it past the threshold, and a deliberately wrong sampler,\n"
                         "each against the exact distribution. Pass: every outcome with exact\n"
                         "probability at least %g within %.0f binomial standard deviations.\n\n",
                 SMALLEST_COMPARED_PROBABILITY, Z_LIMIT);

    int failures = 0;
    long seed = 1;
    for (const Case &c : cases) {
        double smallest_capacity = c.capacity[0];
        for (double capacity : c.capacity) smallest_capacity = std::fmin(smallest_capacity, capacity);
        Distribution exact = (double)c.to_cancel <= smallest_capacity && c.kept.size() == 1
                                 ? uncapped_single_marginal((int)c.capacity.size(), c.to_cancel)
                                 : marginal(exact_distribution(c.capacity, c.to_cancel), c.kept);

        Comparison loop = compare(loop_sampler, c.capacity, c.to_cancel, c.kept, exact, c.loop_replications, seed++);
        Comparison bulk = compare(bulk_sampler, c.capacity, c.to_cancel, c.kept, exact, c.bulk_replications, seed++);
        Comparison wrong = compare(wrong_sampler, c.capacity, c.to_cancel, c.kept, exact, c.loop_replications, seed++);

        bool loop_passes = loop.worst_z < Z_LIMIT;
        bool bulk_passes = bulk.worst_z < Z_LIMIT;
        bool wrong_caught = !c.wrong_sampler_differs || wrong.worst_z >= Z_LIMIT;
        failures += !loop_passes + !bulk_passes + !wrong_caught;

        std::fprintf(report, "%s: capacities", c.name);
        for (size_t k = 0; k < c.capacity.size() && k < 6; k++) std::fprintf(report, " %.0f", c.capacity[k]);
        if (c.capacity.size() > 6) std::fprintf(report, " ... (%zu customers)", c.capacity.size());
        std::fprintf(report, ", %ld to cancel, %d outcomes compared\n", c.to_cancel, loop.compared);
        std::fprintf(report, "  loop   %8ld replications  worst z %6.2f at %s  %s\n", c.loop_replications,
                     loop.worst_z, outcome_text(loop.worst_outcome).c_str(), loop_passes ? "pass" : "FAIL");
        std::fprintf(report, "  bulk   %8ld replications  worst z %6.2f at %s  %s\n", c.bulk_replications,
                     bulk.worst_z, outcome_text(bulk.worst_outcome).c_str(), bulk_passes ? "pass" : "FAIL");
        if (c.wrong_sampler_differs)
            std::fprintf(report, "  wrong  %8ld replications  worst z %6.2f at %s  %s\n", c.loop_replications,
                         wrong.worst_z, outcome_text(wrong.worst_outcome).c_str(),
                         wrong_caught ? "rejected, as it must be" : "NOT REJECTED");
        else
            std::fprintf(report, "  wrong  not compared: with nothing capped it has the right distribution\n");
        std::fprintf(report, "\n");
    }

    /* Beyond what can be enumerated, and beyond bnldev's int. */
    {
        std::vector<double> capacity = {2e9, 2e9, 3e9};
        const double to_cancel = 5e9;
        const long replications = 2000;
        long idum = -(seed++);
        std::vector<double> cancelled(capacity.size());
        std::vector<double> sum(capacity.size(), 0.0);
        int inconsistent = 0;
        for (long r = 0; r < replications; r++) {
            draw_uniform_cancellations(capacity.data(), 3, to_cancel, cancelled.data(), &idum);
            double total = 0;
            for (size_t k = 0; k < capacity.size(); k++) {
                inconsistent += cancelled[k] != std::floor(cancelled[k]) || cancelled[k] < 0 ||
                                cancelled[k] > capacity[k];
                total += cancelled[k];
                sum[k] += cancelled[k];
            }
            inconsistent += total != to_cancel;
        }
        /* Uncapped, each count is binomial with n = 5e9 and p = 1/3, and
           neither mean nor standard deviation reaches a capacity. */
        const double p = 1.0 / 3.0;
        const double mean = to_cancel * p;
        const double sd_of_average = std::sqrt(to_cancel * p * (1 - p) / (double)replications);
        double worst_z = 0;
        for (size_t k = 0; k < capacity.size(); k++)
            worst_z = std::fmax(worst_z, std::fabs(sum[k] / (double)replications - mean) / sd_of_average);
        bool passes = inconsistent == 0 && worst_z < Z_LIMIT;
        failures += !passes;
        std::fprintf(report, "huge: capacities 2e9 2e9 3e9, 5e9 to cancel, %ld replications of the bulk draw\n",
                     replications);
        std::fprintf(report, "  inconsistent counts %d, worst z of a mean count %.2f  %s\n\n", inconsistent,
                     worst_z, passes ? "pass" : "FAIL");
    }

    std::fprintf(report, "%s\n", failures == 0 ? "PASSED" : "FAILED");
    std::fclose(report);

    std::printf("bulk cancellation distribution: %d failures\n", failures);
    std::printf("%s\n", failures == 0 ? "PASSED, 0 failures" : "FAILED");
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
