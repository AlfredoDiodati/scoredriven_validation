// Cancelling machine orders by the million.
//
// When labour is rationed, LABOR() cancels a capital-good firm's lost output
// from its customers' orders one machine at a time: draw one of the N2 firms at
// random, and cancel one machine if it is a customer with an order still
// outstanding. A draw that lands elsewhere is thrown away, so each cancellation
// falls on a customer chosen uniformly among those whose order is not yet used
// up. That costs about N2 over the number of customers draws per machine, and
// the number of machines grows with real output, which the model lets grow
// without limit: hundreds of millions of machines for one firm by period 8,889
// of a long run, which docs/DSK_LONG_HORIZON.md records.
//
// This draws the same outcome as counts. Only the number of cancellations each
// customer receives matters to the model, and those counts have the following
// description. Give every cancellation a label drawn uniformly over the
// customers; a label naming a customer whose order is used up is thrown away,
// exactly as the loop throws away a draw. Then the first `remaining` labels are
// all real cancellations except the ones past a customer's capacity, so drawing
// them as a multinomial count per open customer, capping, and drawing again for
// the ones that were capped away gives the counts the loop would have given, in
// distribution. Each round closes at least one customer or finishes, so there
// are at most as many rounds as customers.
//
// The draws come from the model's own stream through bnldev, so a run stays a
// function of its seed. They are not the draws the loop would have taken, so a
// run that goes through here is not the run the loop would have produced, only
// one with the same distribution. tests/dsk_bulk_cancellation_distribution.cpp
// checks that distribution against exact probabilities worked out by
// enumeration.

#ifndef DSK_SFC_BULK_CANCELLATION_H
#define DSK_SFC_BULK_CANCELLATION_H

#include <vector>
#include <cmath>

#include "auxiliary/bnldev.h"
#include "auxiliary/ran1.h"

// With at least 25 trials and a mean below one, bnldev returns a Poisson draw
// rather than a binomial one, which at 30 trials and probability 1/40 puts 12%
// too much weight on a count of four. Inverting the binomial distribution
// directly is exact, and with a mean below one the search stops after a step or
// two.
inline double binomial_count_by_inversion(double trials, double probability, long *idum)
{
   const double u=ran1(idum);
   const double odds=probability/(1.0-probability);
   double k=0;
   double pmf=std::exp(trials*std::log1p(-probability));
   double cdf=pmf;
   while (u>cdf && k<trials)
   {
      pmf*=(trials-k)/(k+1)*odds;
      k+=1;
      cdf+=pmf;
   }
   return k;
}

// One binomial count over at most an int's worth of trials.
inline double binomial_count_piece(double trials, double probability, long *idum)
{
   if (trials>=25 && trials*probability<1.0)
   {
      return binomial_count_by_inversion(trials,probability,idum);
   }
   return bnldev(probability,(int)trials,idum);
}

// bnldev takes an int trial count. A binomial count over more trials is the sum
// of binomial counts over pieces of them.
inline double binomial_count(double trials, double probability, long *idum)
{
   const double largest_piece=1e9;
   double count=0;
   while (trials>largest_piece)
   {
      count+=binomial_count_piece(largest_piece,probability,idum);
      trials-=largest_piece;
   }
   if (trials>0)
   {
      count+=binomial_count_piece(trials,probability,idum);
   }
   return count;
}

// capacity[k] is how many cancellations customer k can take before its order is
// used up. Writes the number each one takes into cancelled[k]. to_cancel must be
// a whole number no larger than the sum of the capacities.
inline void draw_uniform_cancellations(const double *capacity, int customers, double to_cancel,
                                       double *cancelled, long *idum)
{
   std::vector<int> open;
   for (int k=0; k<customers; k++)
   {
      cancelled[k]=0;
      if (capacity[k]>0)
      {
         open.push_back(k);
      }
   }

   double remaining=to_cancel;
   std::vector<int> still_open;
   while (remaining>0 && !open.empty())
   {
      // One multinomial count of `remaining` over the open customers with equal
      // probabilities, taken as successive binomials.
      double unassigned=remaining;
      double capped_away=0;
      still_open.clear();
      const int n_open=(int)open.size();
      for (int position=0; position<n_open; position++)
      {
         const int k=open[position];
         double share;
         if (position==n_open-1)
         {
            share=unassigned;
         }
         else
         {
            share=binomial_count(unassigned,1.0/(double)(n_open-position),idum);
         }
         unassigned-=share;

         cancelled[k]+=share;
         if (cancelled[k]>=capacity[k])
         {
            capped_away+=cancelled[k]-capacity[k];
            cancelled[k]=capacity[k];
         }
         else
         {
            still_open.push_back(k);
         }
      }
      remaining=capped_away;
      open.swap(still_open);
   }
}

#endif
