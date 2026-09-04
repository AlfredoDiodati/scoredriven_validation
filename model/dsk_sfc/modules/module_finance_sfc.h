#ifndef MODULE_FINANCE_H
#define MODULE_FINANCE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <fenv.h>

//#include <string>
#include <string.h>
#include <sstream>

// Include Newmat Libraries and random number generators
#include "../newmat10/include.h"
#include "../newmat10/newmat.h"
#include "../newmat10/newmatio.h"
#include "../auxiliary/ran1.h"

// Include functions from other modules
#include "../dsk_sfc_functions.h"


// -- Functions -- //   
void TOTCREDIT(void);                                                       // Determines maximum amount of credit banks will extend
void LOANRATES(void);                                                       // Determines loan rates charged to individual borrowers
void BANKING(void);                                                         // Determines bank profits; Banks receive second-hand capital from failing firms
void BAILOUT(void);                                                         // Failing banks are bailed out by government or bought by other banks depending on scenario
void SETTLEMENT(void);                                                      // Settlement of interbank transactions; granting/repayment of CB Advances

//-- Flags --//
extern int              flagbailout;
extern int              flag_rate_setting_loans;

//-- Pars --//
extern double           floor_default_probability;
extern double           upsilon;
extern double           lambdaB1;
extern double           lambdaB2;
extern double           riskWeightLoans;	
extern double           riskWeightGovBonds;	
extern double           capitalAdequacyRatioTarget;
extern int              NB;
extern int              N2;
extern int              N1;
extern double           k_const;
extern double           db;
extern double           aliqb;
extern long             *p_seed;
extern double           b1sup;
extern double           b1inf;
extern double           b2sup;
extern double           b2inf;

// -- Vars -- //
extern int              i;
extern int              j;
extern double           tolerance;
extern RowVector        BankCredit;
extern RowVector        BaselBankCredit;
extern Matrix           NW_b;
extern Matrix           NW_2;
extern RowVector        capitalAdequacyRatio;	
extern RowVector        riskWeightedAssets;
extern Matrix           S2;
extern Matrix           Loans_2;
extern RowVector        CreditDemand;
extern RowVector        DebtServiceToSales2; 
extern Real             DS2_min;   
extern Matrix           DebtServiceToSales2_bank;
extern Matrix           BankMatch_2;
extern Matrix           BankMatch_1;
extern RowVector        DebtServiceToSales2_temp; 
extern std::vector<std::pair<double,int> > DS2_ranked;
extern std::vector<int> DS2_by_rank;
extern int              DS2_min_index; 
extern Matrix           DS2_rating;
extern Matrix           DebtService_2;
extern RowVector        k;
extern RowVector        NL_2;
extern RowVector        NL_1;
extern RowVector        r_deb; 
extern RowVector        r_deb_h;
extern RowVector        FirmDefaultProbability;
extern RowVector        LoanInterest;
extern RowVector        InterestDeposits;
extern RowVector        baddebt_b;  
extern double           r_bonds;
extern Matrix           GB_b;
extern double           r;
extern double           r_cbreserves;
extern Matrix           Reserves_b;
extern Matrix           Advances_b;
extern RowVector        Dividends_b;
extern double           Taxes;
extern RowVector        Taxes_b;
extern Matrix           Deposits;
extern Matrix           Deposits_hb;
extern Matrix           Deposits_eb;
extern RowVector        Deposits_h;
extern RowVector        Outflows;
extern RowVector        Inflows;
extern RowVector        BankProfits;
extern RowVector        Bank_active;
extern double           Bailout;
extern RowVector        BankEquity_temp; 
extern int              maxbank;
extern double           max_equity;
extern double           multip_bailout;
extern RowVector        Bailout_b;
extern RowVector        LossAbsorbed;
extern Matrix           Loans_b;
extern RowVector        BankingSupplier_2;
extern RowVector        BankingSupplier_1;
extern RowVector        capitalRecovered;
extern RowVector        ReserveBalance;
extern RowVector        Advances;
extern RowVector        Reserves;
extern RowVector        Dividends;
extern double           InterestReserves;
extern double           InterestAdvances;
extern RowVector        InterestReserves_b;
extern RowVector        InterestAdvances_b;
extern RowVector        ProfitCB;
extern double           InterestBonds_cb;
extern RowVector        LossEntry_b;
extern RowVector        DepositShare_h;
extern RowVector        DepositShare_e;
extern double           counter_bankfailure;

#endif