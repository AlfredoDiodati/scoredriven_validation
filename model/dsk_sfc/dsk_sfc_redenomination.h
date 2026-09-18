// Redenominating the money side of the model.
//
// Prices, wages and every money stock in this model grow without limit: over
// periods 250 to 12,750 of one long run the wage grows 0.74% a period, which is
// a bit of precision every 94 periods, and it passes the largest double near
// period 96,000. docs/DSK_LONG_HORIZON.md has the measurements.
//
// Nothing in the model's behaviour depends on the unit money is counted in, so
// the answer is to change the unit: divide every money quantity by a power of
// two once they grow large, the way a currency reform drops zeros. A power of
// two is exact in binary floating point and exactly reversible, so the only
// thing that changes is the exponent of every money number, not its digits.
// A run whose money quantities never reach the ceiling never calls this at all
// and is bit for bit the run upstream produces.
//
// What counts as money is written out below, one list per kind. Machine counts,
// productivities, employment, emissions and every rate or share are not money
// and are left alone. A variable missing from these lists, or wrongly included,
// makes the model inconsistent from the first redenomination onwards, which is
// what tests/dsk_redenomination_invariance.c looks for: it runs the model with
// the ceiling low enough to fire, and requires every column of the results file
// to come back either unchanged or scaled by exactly the factor applied.

#ifndef DSK_SFC_REDENOMINATION_H
#define DSK_SFC_REDENOMINATION_H

#include <cmath>
#include <cstdio>
#include <cstring>

// Powers of two the money side has been divided by so far, as a negative
// exponent. Written to the log below so a reader of the results file can put
// the numbers back in the units the run started in.
int redenomination_exponent=0;

void REDENOMINATE(double factor)
{
   // Prices, wages and unit costs.
   w=w*factor;
   p1=p1*factor;
   p2=p2*factor;
   p1test*=factor;
   p2m*=factor;
   p2_entry*=factor;
   c1=c1*factor;
   c2=c2*factor;
   c2e=c2e*factor;
   c_en=c_en*factor;
   c_de_min*=factor;
   cf_min_ge*=factor;
   cmax*=factor;
   C=C*factor;
   C_secondhand=C_secondhand*factor;
   C_de=C_de*factor;
   c_infra*=factor;
   pf*=factor;
   cpi=cpi*factor;
   kpi*=factor;
   t_CO2*=factor;
   t_CO2_en*=factor;
   initial_price*=factor;

   // Deposits, loans, bonds, reserves and advances.
   Deposits=Deposits*factor;
   Deposits_1=Deposits_1*factor;
   Deposits_2=Deposits_2*factor;
   Deposits_h=Deposits_h*factor;
   Deposits_e=Deposits_e*factor;
   Deposits_hb=Deposits_hb*factor;
   Deposits_eb=Deposits_eb*factor;
   Deposits_fuel=Deposits_fuel*factor;
   Deposits_fuel_cb=Deposits_fuel_cb*factor;
   Loans_2=Loans_2*factor;
   Loans_b=Loans_b*factor;
   Loans_e=Loans_e*factor;
   Loans_preDefault_e*=factor;
   GB=GB*factor;
   GB_b=GB_b*factor;
   GB_cb=GB_cb*factor;
   Advances=Advances*factor;
   Advances_b=Advances_b*factor;
   Reserves=Reserves*factor;
   Reserves_b=Reserves_b*factor;
   bonds_dem=bonds_dem*factor;
   bonds_dem_tot*=factor;
   bonds_purchased=bonds_purchased*factor;
   NewBonds*=factor;
   BankCredit=BankCredit*factor;
   BaselBankCredit=BaselBankCredit*factor;
   CreditDemand=CreditDemand*factor;
   CreditDemand_e=CreditDemand_e*factor;
   riskWeightedAssets=riskWeightedAssets*factor;

   // Net worth, capital and inventories at their money value.
   NW_1=NW_1*factor;
   NW_2=NW_2*factor;
   NW_b=NW_b*factor;
   NW_h=NW_h*factor;
   NW_gov=NW_gov*factor;
   NW_cb=NW_cb*factor;
   NW_e=NW_e*factor;
   NW_f=NW_f*factor;
   NW_1_c=NW_1_c*factor;
   NW_2_c=NW_2_c*factor;
   NW_b_c=NW_b_c*factor;
   NW_h_c*=factor;
   NW_gov_c*=factor;
   NW_cb_c*=factor;
   NW_e_c*=factor;
   NW_f_c*=factor;
   NWSum*=factor;
   RealAssets*=factor;
   CapitalStock=CapitalStock*factor;
   deltaCapitalStock=deltaCapitalStock*factor;
   CapitalStock_e=CapitalStock_e*factor;
   old_capitalStock*=factor;
   Inventories=Inventories*factor;
   dNm=dNm*factor;
   dNmtot*=factor;
   G_ge_n=G_ge_n*factor;
   G_ge_n_0*=factor;
   max_equity*=factor;
   BankEquity_temp=BankEquity_temp*factor;

   // Incomes, payments and transfers.
   Wages*=factor;
   Wages_1=Wages_1*factor;
   Wages_2=Wages_2*factor;
   Wages_en*=factor;
   Dividends=Dividends*factor;
   Dividends_1=Dividends_1*factor;
   Dividends_2=Dividends_2*factor;
   Dividends_b=Dividends_b*factor;
   Dividends_e*=factor;
   Divtot_1*=factor;
   Divtot_2*=factor;
   Benefits*=factor;
   Consumption*=factor;
   Cons*=factor;
   Cres*=factor;
   Cresb*=factor;
   S1=S1*factor;
   S2=S2*factor;
   S1_temp=S1_temp*factor;
   S2_temp=S2_temp*factor;
   Sales1=Sales1*factor;
   Sales2=Sales2*factor;
   mol=mol*factor;
   Pi1=Pi1*factor;
   Pi2=Pi2*factor;
   Pitot1*=factor;
   Pitot2*=factor;
   GDP_n=GDP_n*factor;
   Investment_2=Investment_2*factor;
   Investment_n*=factor;
   ExpansionInvestment_n*=factor;
   ReplacementInvestment_n*=factor;
   EI_n=EI_n*factor;
   SI_n=SI_n*factor;
   Cmach=Cmach*factor;
   CmachEI=CmachEI*factor;
   CmachSI=CmachSI*factor;
   scrap_n*=factor;
   EnergyPayments*=factor;
   EnergyPayments_1=EnergyPayments_1*factor;
   EnergyPayments_2=EnergyPayments_2*factor;
   FuelCost*=factor;
   TransferFuel*=factor;
   TransferCB*=factor;
   govTranfers*=factor;
   FirmTransfers*=factor;
   FirmTransfers_1*=factor;
   FirmTransfers_2*=factor;
   Injection_1=Injection_1*factor;
   Injection_2=Injection_2*factor;
   injection*=factor;
   injection2*=factor;
   EntryCosts*=factor;
   BankTransfer*=factor;
   Bailout*=factor;
   Bailout_b=Bailout_b*factor;
   multip_bailout*=factor;
   capitalRecovered=capitalRecovered*factor;
   capitalRecovered2=capitalRecovered2*factor;
   G*=factor;
   PSBR*=factor;
   Taxes*=factor;
   Taxes_1=Taxes_1*factor;
   Taxes_2=Taxes_2*factor;
   Taxes_b=Taxes_b*factor;
   Taxes_h*=factor;
   Taxes_CO2=Taxes_CO2*factor;
   Taxes_CO2_1=Taxes_CO2_1*factor;
   Taxes_CO2_2=Taxes_CO2_2*factor;
   Taxes_CO2_e*=factor;

   // Interest, debt service and bad debt.
   InterestDeposits=InterestDeposits*factor;
   InterestDeposits_1=InterestDeposits_1*factor;
   InterestDeposits_2=InterestDeposits_2*factor;
   InterestDeposits_h*=factor;
   InterestDeposits_e*=factor;
   InterestBonds*=factor;
   InterestBonds_b=InterestBonds_b*factor;
   InterestBonds_cb*=factor;
   InterestReserves*=factor;
   InterestReserves_b=InterestReserves_b*factor;
   InterestAdvances*=factor;
   InterestAdvances_b=InterestAdvances_b*factor;
   BondRepayments_b=BondRepayments_b*factor;
   BondRepayments_cb*=factor;
   LoanInterest=LoanInterest*factor;
   LoanInterest_2=LoanInterest_2*factor;
   Loan_interest_e*=factor;
   DebtService_2=DebtService_2*factor;
   DebtService_e*=factor;
   DebtRemittances2=DebtRemittances2*factor;
   DebtRemittances_e*=factor;
   DebtWrittenOff_e*=factor;
   DeafaultedDebtRecovered_e*=factor;
   DefaultedDeposits_e*=factor;
   BadDebt_e*=factor;
   baddebt_1=baddebt_1*factor;
   baddebt_2=baddebt_2*factor;
   baddebt_b=baddebt_b*factor;
   baddebt_2_temp*=factor;
   Deposits_recovered_1*=factor;
   Deposits_recovered_2*=factor;
   DepositsCheck_1*=factor;
   DepositsCheck_2*=factor;
   sumDepEn*=factor;
   MaxFunds*=factor;
   mD1*=factor;
   mD2*=factor;

   // The energy sector's own money quantities.
   CF_ge=CF_ge*factor;
   CF_ge_inn*=factor;
   IC_en*=factor;
   IC_en_quota=IC_en_quota*factor;
   PC_en*=factor;
   Rev_en*=factor;
   RD_en_de*=factor;
   RD_en_ge*=factor;
   ProfitEnergy*=factor;
   ProfitCB=ProfitCB*factor;
   BankProfits=BankProfits*factor;
   BankProfits_temp=BankProfits_temp*factor;

   // Research, losses and the reserve accounting residuals.
   RD=RD*factor;
   Loss_Capital=Loss_Capital*factor;
   Loss_Inventories=Loss_Inventories*factor;
   Outflows=Outflows*factor;
   Inflows=Inflows*factor;
   ReserveBalance=ReserveBalance*factor;
   Adjustment=Adjustment*factor;
   Adjustment_cb*=factor;
   cpi_init*=factor;
   GDP_init*=factor;

   // One coefficient is per unit of money rather than money: the chance the
   // energy sector's research succeeds is 1-exp(-o1_en*RD_en) with RD_en a share
   // of revenue, so the product is a pure number only if o1_en moves the other
   // way. The two manufacturing counterparts, o1 and o2, multiply Ld1rd, which
   // is R&D staff rather than money, and are left alone.
   o1_en/=factor;

   // The money constants the model keeps reading as the run goes on. A floor
   // under the wage or a price, and the base of a carbon tax, are amounts in the
   // old unit, and would start binding at the wrong level if left behind.
   w_min*=factor;
   sales_tolerance*=factor;
   pmin*=factor;
   cpi_floor*=factor;
   consumption_residual_floor*=factor;
   t_CO2_0*=factor;
   t_CO2_en_0*=factor;

   // Sectoral balances, the balance-sheet checks and the rest.
   Balance_h*=factor;
   Balance_1*=factor;
   Balance_2*=factor;
   Balance_e*=factor;
   Balance_b*=factor;
   Balance_g*=factor;
   Balance_cb*=factor;
   Balance_f*=factor;
   BalanceSum*=factor;
   Balances_1=Balances_1*factor;
   prior=prior*factor;
   post*=factor;
   prior_cb*=factor;
   post_cb*=factor;
   Deficit*=factor;
   cpi_temp*=factor;
   mi_en*=factor;

   // The per-machine prices carried with the vintages.
   for (std::size_t k=0; k<g_price.data.size(); k++) g_price.data[k]*=factor;
   for (std::size_t r=0; r<g_secondhand_p.size(); r++) for (std::size_t c=0; c<g_secondhand_p[r].size(); c++) g_secondhand_p[r][c]*=factor;
}

// The wage is the fastest growing money quantity in the model, so it is what
// decides when to redenominate.
void REDENOMINATE_IF_NEEDED(void)
{
   if (redenomination_ceiling_exponent<=0 || redenomination_step_exponent<=0)
   {
      return;
   }
   if (fabs(w(1))<ldexp(1.0,redenomination_ceiling_exponent))
   {
      return;
   }

   REDENOMINATE(ldexp(1.0,-redenomination_step_exponent));
   redenomination_exponent-=redenomination_step_exponent;

   // Named after the results file it belongs to, and written only by a run that
   // actually redenominates.
   char path[PATH_MAX];
   strcpy(path,filename1);
   char *suffix=strstr(path,".txt");
   if (suffix)
   {
      strcpy(suffix,"_redenominations.txt");
   }
   ofstream log(path,ios::app);
   log << t << " " << -redenomination_step_exponent << " " << redenomination_exponent << endl;
}

#endif
