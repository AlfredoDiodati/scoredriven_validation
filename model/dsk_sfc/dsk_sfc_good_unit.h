// Changing the unit the consumption good is measured in.
//
// Real output grows about 0.24% a period and nothing rescales it, so real GDP
// and the capital stock pass the largest double near period 291,000 - a wall of
// their own, after the money side is redenominated and the machines are counted
// in lots. docs/DSK_LONG_HORIZON.md has the measurements.
//
// Nothing in the model depends on whether the good is counted in tonnes or in
// hundreds of tonnes, so the answer is the same as for money: once output grows
// large, count it in a bigger unit, a power of two so that only exponents move.
//
// What goes with the unit, and which way:
//
//   quantities of the good and the capacity to make it     divided
//   how much of it one machine makes, and the payback      divided
//     threshold, which is an amount of the good
//   output per worker and per unit of energy               divided
//   money for one unit of it: prices, unit costs, the      multiplied
//     price index and the floor under it
//
// Machines, money stocks, labour, energy and emissions are counted in their own
// units and do not move. Capacity is machines times what a machine makes, so it
// divides once, not twice.
//
// Unlike the machine lots this is exact: nothing here has to stay a whole
// number, because what the model keeps whole is machines, and a machine count
// is a ratio of two quantities that both divide. tests/dsk_good_unit_invariance.c
// requires every column of the results file to come back unchanged, divided or
// multiplied by exactly the factor applied.

#ifndef DSK_SFC_GOOD_UNIT_H
#define DSK_SFC_GOOD_UNIT_H

#include <cmath>
#include <cstdio>
#include <cstring>

// Powers of two the good's unit has grown by so far.
int good_unit_exponent=0;

void RESCALE_GOOD(double unit)
{
   const double per_unit=1.0/unit;

   // Quantities of the consumption good, and the capacity that makes them.
   Q2=Q2*per_unit;
   Q2temp=Q2temp*per_unit;
   Q2tot*=per_unit;
   Qd=Qd*per_unit;
   D2=D2*per_unit;
   De=De*per_unit;
   D_temp2=D_temp2*per_unit;
   D20*=per_unit;
   CurrentDemand*=per_unit;
   l2=l2*per_unit;
   l2m*=per_unit;
   N=N*per_unit;
   Ne=Ne*per_unit;
   dN=dN*per_unit;
   dNtot*=per_unit;
   I=I*per_unit;
   EI=EI*per_unit;
   SI=SI*per_unit;
   Ip=Ip*per_unit;
   EIp=EIp*per_unit;
   SIp=SIp*per_unit;
   EId=EId*per_unit;
   SId=SId*per_unit;
   Kd=Kd*per_unit;
   Ktrig=Ktrig*per_unit;
   K=K*per_unit;
   K_top*=per_unit;
   K_cur=K_cur*per_unit;
   K_temp=K_temp*per_unit;
   K_temp_sum*=per_unit;
   K_loss=K_loss*per_unit;
   scrap_age=scrap_age*per_unit;
   C_loss=C_loss*per_unit;
   I_loss=I_loss*per_unit;
   lossdouble*=per_unit;
   GDP_r=GDP_r*per_unit;
   Consumption_r*=per_unit;
   Investment_r*=per_unit;
   ExpansionInvestment_r*=per_unit;
   ReplacementInvestment_r*=per_unit;

   // How much of the good one machine makes, and the payback threshold, which is an amount of it.
   good_unit_floor*=per_unit;
   dim_mach*=per_unit;
   b*=per_unit;

   // Output per worker and per unit of energy, in use and in the innovation and imitation draws.
   A2=A2*per_unit;
   A2e=A2e*per_unit;
   A2_mprod=A2_mprod*per_unit;
   A2e2=A2e2*per_unit;
   A1=A1*per_unit;
   A1inn=A1inn*per_unit;
   A1imm=A1imm*per_unit;
   A=A*per_unit;
   Am2*=per_unit;
   A2_en=A2_en*per_unit;
   A2e_en=A2e_en*per_unit;
   A2e_en2=A2e_en2*per_unit;
   A1_en=A1_en*per_unit;
   EE_inn=EE_inn*per_unit;
   EE_imm=EE_imm*per_unit;
   A_en=A_en*per_unit;

   // Upstream's mean productivities add the consumption firms' output per worker
   // to the capital firms' machines per worker, and its energy counterpart does
   // the same with output and machines per unit of energy. Those are amounts of
   // different things, so they cannot be divided by the unit: they are set to
   // what the model itself computes in the new unit, which is what leaves the
   // growth rate the wage rule reads where it was. Am1, the capital firms' own
   // mean, is machines per worker throughout and does not move.
   {
      const double machines_per_worker=A1p.Sum();
      const double mixed=Am(1)*(N1r+N2r);
      if (mixed!=0)
      {
         Am=Am*(((mixed-machines_per_worker)*per_unit+machines_per_worker)/mixed);
      }
   }
   {
      double machines_per_worker=0;
      for (int firm=1; firm<=N1; firm++)
      {
         machines_per_worker+=Ld1(firm)*A1p(firm);
      }
      if (LD2!=0 && Am_a!=0)
      {
         machines_per_worker/=LD2;
         Am_a=(Am_a-machines_per_worker)*per_unit+machines_per_worker;
      }
   }
   {
      const double energy=D2_en_TOT+D1_en_TOT;
      double machines_per_energy=0;
      for (int firm=1; firm<=N1; firm++)
      {
         machines_per_energy+=D1_en(firm)*A1p_en(firm);
      }
      if (energy!=0 && Am_en(1)!=0)
      {
         machines_per_energy/=energy;
         Am_en=Am_en*(((Am_en(1)-machines_per_energy)*per_unit+machines_per_energy)/Am_en(1));
      }
   }

   // Money for one unit of the good: prices, unit costs and the price index.
   p2=p2*unit;
   p2m*=unit;
   p2_entry*=unit;
   ptemp*=unit;
   c2=c2*unit;
   c2e=c2e*unit;
   C=C*unit;
   C_secondhand=C_secondhand*unit;
   cmax*=unit;
   cpi=cpi*unit;
   cpi_temp*=unit;
   cpi_init*=unit;
   pmin*=unit;
   cpi_floor*=unit;
}

// Real GDP is the largest quantity of the good in the model, so it is what
// decides when to change the unit.
void RESCALE_GOOD_IF_NEEDED(void)
{
   if (good_unit_ceiling_exponent<=0 || good_unit_step_exponent<=0)
   {
      return;
   }
   if (fabs(GDP_r(1))<ldexp(1.0,good_unit_ceiling_exponent))
   {
      return;
   }

   RESCALE_GOOD(ldexp(1.0,good_unit_step_exponent));
   good_unit_exponent+=good_unit_step_exponent;

   char path[PATH_MAX];
   strcpy(path,filename1);
   char *suffix=strstr(path,".txt");
   if (suffix)
   {
      strcpy(suffix,"_good_units.txt");
   }
   ofstream log(path,ios::app);
   log << t << " " << good_unit_step_exponent << " " << good_unit_exponent << endl;
}

#endif
