// Counting machines in lots.
//
// Machines are whole numbers held in doubles, and a double holds every whole
// number exactly only up to 2^53. The number of machines grows with real
// output, about 0.24% a period, so it passes that near period 11,500 and by
// period 12,908 the totals no longer agree: the second-hand pool comes up one
// machine short of what the entering firms are owed and the loop that hands
// them out never ends. docs/DSK_LONG_HORIZON.md has the measurements.
//
// Nothing in the model depends on how big one machine is, so the answer is to
// make machines bigger: double the output a machine produces, halve every count,
// and the economy is the one it was with half as many machines to count. That
// buys back a bit of headroom each time, and it can be done as often as needed.
//
// What moves with the lot, and which way:
//
//   counts of machines                       halved
//   output a machine makes, and the payback   doubled
//     threshold, which is in units of output
//   money per machine: prices of machines,    doubled
//     their production cost, the capital
//     goods price index, the price each
//     machine was bought at
//   machines per worker and per unit of       halved
//     energy, in use and in the innovation
//     and imitation draws
//
// Everything else - output, capital, employment, productivity per worker, every
// money stock - is unchanged, because capacity is a count times the output a
// machine makes and the two moves cancel.
//
// Unlike redenominating money, this is not exact. Halving a count of an odd
// number of machines leaves half a machine, which is rounded, so a lot change
// perturbs each firm-vintage holding by at most half a lot. At the default
// ceiling that is about one part in 1e12 of a holding. The model's stocks are
// then re-derived from the rounded counts rather than left to disagree with
// them: at a period boundary the model keeps a firm's capacity equal to its
// machines times the output a machine makes, its machine count equal to the sum
// of its holdings, and its capital value equal to the holdings priced at what
// each was bought for, all three exactly, which tests/dsk_machine_lot_rebase.c
// checks before and after.
//
// So a run that changes lots is a different path of the same process, not the
// same path in different units. A run that never reaches the ceiling never calls
// this, and is bit for bit the run upstream produces.

#ifndef DSK_SFC_MACHINE_LOTS_H
#define DSK_SFC_MACHINE_LOTS_H

#include <cmath>
#include <cstdio>
#include <cstring>

// How many times the lot has been doubled so far.
int machine_lot_exponent=0;

// Rounds every holding to a whole number of the new lot and puts the stocks the
// model derives from them back in step.
void _restate_machine_stocks(void)
{
   for (int firm=1; firm<=N2; firm++)
   {
      double machines=0, value=0;
      for (int vintage=t0; vintage<=t; vintage++)
      {
         for (int supplier=1; supplier<=N1; supplier++)
         {
            machines+=g[vintage-1][supplier-1][firm-1];
            value+=g[vintage-1][supplier-1][firm-1]*g_price[vintage-1][supplier-1][firm-1];
         }
      }
      n_mach(firm)=machines;
      K(firm)=machines*dim_mach;
      CapitalStock(1,firm)=value;
   }
}

void REBASE_MACHINE_LOT(int step_exponent)
{
   const double lot=ldexp(1.0,step_exponent);
   const double per_lot=1.0/lot;

   // The holdings themselves, rounded to whole machines in the new lot.
   for (std::size_t k=0; k<g.data.size(); k++) g.data[k]=floor(g.data[k]*per_lot+0.5);
   for (std::size_t k=0; k<gtemp.data.size(); k++) gtemp.data[k]=floor(gtemp.data[k]*per_lot+0.5);
   for (std::size_t k=0; k<g_c2.data.size(); k++) g_c2.data[k]=floor(g_c2.data[k]*per_lot+0.5);
   for (std::size_t k=0; k<g_c3.data.size(); k++) g_c3.data[k]=floor(g_c3.data[k]*per_lot+0.5);
   for (std::size_t r=0; r<g_secondhand.size(); r++)
   {
      for (std::size_t c=0; c<g_secondhand[r].size(); c++)
      {
         g_secondhand[r][c]=floor(g_secondhand[r][c]*per_lot+0.5);
      }
   }

   // Counts the model carries alongside the holdings.
   n_mach=n_mach*per_lot;
   n_mach_entry=n_mach_entry*per_lot;
   D1=D1*per_lot;
   Q1=Q1*per_lot;
   Q1tot*=per_lot;
   Qpast*=per_lot;
   reduction*=per_lot;
   nmachprod*=per_lot;
   nmp_temp*=per_lot;
   scrapmax*=per_lot;
   n_mach_exit*=per_lot;
   n_mach_exit2*=per_lot;
   n_mach_needed*=per_lot;
   n_mach_resid*=per_lot;
   n_mach_resid2*=per_lot;

   // Money for one machine.
   p1=p1*lot;
   p1test*=lot;
   c1=c1*lot;
   kpi*=lot;
   for (std::size_t k=0; k<g_price.data.size(); k++) g_price.data[k]*=lot;
   for (std::size_t r=0; r<g_secondhand_p.size(); r++)
   {
      for (std::size_t c=0; c<g_secondhand_p[r].size(); c++)
      {
         g_secondhand_p[r][c]*=lot;
      }
   }

   // Output for one machine, and the payback threshold, which is the output a
   // machine has to produce to earn back what it cost.
   dim_mach*=lot;
   b*=lot;

   // Machines per worker and per unit of energy.
   A1p=A1p*per_lot;
   A1p_en=A1p_en*per_lot;
   A1pinn=A1pinn*per_lot;
   A1pimm=A1pimm*per_lot;
   EEp_inn=EEp_inn*per_lot;
   EEp_imm=EEp_imm*per_lot;
   pm*=per_lot;

   _restate_machine_stocks();
   machine_lot_exponent+=step_exponent;
}

// The largest holding of any one firm is what decides when to rebase: it is the
// count that reaches the limit first, and the sums that broke ran over the firms.
void REBASE_MACHINE_LOT_IF_NEEDED(void)
{
   if (machine_lot_ceiling_exponent<=0 || machine_lot_step_exponent<=0)
   {
      return;
   }

   double largest=0;
   for (int firm=1; firm<=N2; firm++)
   {
      largest=fmax(largest,fabs(n_mach(firm)));
   }
   if (largest<ldexp(1.0,machine_lot_ceiling_exponent))
   {
      return;
   }

   REBASE_MACHINE_LOT(machine_lot_step_exponent);

   char path[PATH_MAX];
   strcpy(path,filename1);
   char *suffix=strstr(path,".txt");
   if (suffix)
   {
      strcpy(suffix,"_machine_lots.txt");
   }
   ofstream log(path,ios::app);
   log << t << " " << machine_lot_step_exponent << " " << machine_lot_exponent << " " << dim_mach << endl;
}

#endif
