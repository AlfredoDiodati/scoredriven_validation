#include "module_finance_sfc.h"
#include "../dsk_sfc_reductions.h"

///FUNCTIONS///

void TOTCREDIT(void)
{
  //Banks determine the maximum amount of credit they are willing to extend based on a regulatory ratio 
  for(j=1; j<=NB; j++)   
  {
    //Consider the case where RWAs are only made up of loans, and government bonds are assumed to have a zero risk weight	
    //In that case, capitalAdequacyRatio(j) = NW_b(1,j)/riskWeightedAssets(j) = NW_b(1,j)/(riskWeightLoans*Loans_b(1,j))	
    //For a given net wealth, banks can reach the regulatory capitalAdequacyRatioTarget by setting their credit supply such that:	
    //BaselBankCredit(j)=NW_b(1,j)/(capitalAdequacyRatioTarget*riskWeightLoans)	
    //As more assets are added to the model, the credit supply equation becomes more complex	
    //In order to avoid constantly modifying the definition of BaselBankCredit(j), one can opt for a more general form using riskWeightedAssets(j) minus loans and their risk weight	
    //As a result, each time an asset is added to the model, you only have to modify riskWeightedAssets(j) while keeping the following part of the code intact	
	  BaselBankCredit(j) = ((NW_b(1,j)/capitalAdequacyRatioTarget) - (riskWeightedAssets(j)-riskWeightLoans*Loans_b(1,j)))/riskWeightLoans;
    BankCredit(j)=BaselBankCredit(j);
  } 
} 

void LOANRATES(void) 
{  
  //Banks rank their C-firm customers according to their debt service to sales ratio
  for (j=1; j<=N2; ++j) 
  {
    //Compute debt service to sales; Add tolerance (small value) to ensure model does not break if
    //S2(2,j)=0
    DebtServiceToSales2(j)=DebtService_2(2,j)/(S2(2,j)+tolerance);
  }

  //Iterate over all banks
  for (i=1; i<=NB; ++i)            
  {
    for (j=1; j<=N2; ++j)    
    {	
      //For every bank i, add the debt service ratios of all C-Firms to column i of the matrix
      //DebtServiceToSales2_bank; even if j is not a customer of i
      DebtServiceToSales2_bank(j,i)=DebtServiceToSales2(j);
    } 
    //Iterate again over all C-Firms; if j is not a customer of bank i, set the corresponding value of
    //DebtServiceToSales2_bank to the maximum across the column plus one
    //Thereby, all firms which are not customers of j will appear at the end of the ranking
    //
    //Each of those assignments makes the value it wrote the new maximum, so
    //the next one is the one before it plus one. Carrying that value forward
    //gives the same sequence of writes as rescanning the column would, and the
    //column is scanned once for the bank rather than once for each of its
    //roughly 180 non-customers.
    Real column_maximum=ColumnMaximum(DebtServiceToSales2_bank,i);
    for (j=1; j<=N2; ++j)    
    {	
      if (BankMatch_2(j,i)==0)
      {
        column_maximum=column_maximum+1;
        DebtServiceToSales2_bank(j,i)=column_maximum;
      }   
    } 
  } 


  for (i=1; i<=NB; ++i)            
  {
    //Store the ranking of bank i in a temporary storage vector
    for (j=1; j<=N2; ++j)   
    {
      DebtServiceToSales2_temp(j) = DebtServiceToSales2_bank(j,i);
    }

    //Rank the bank's borrowers by debt service ratio, lowest first. Written as
    //a selection sort - N2 passes, each scanning all N2 firms for the minimum
    //and all N2 again for the maximum to push the chosen firm out of the way -
    //which is 80,000 comparisons per bank per period. Sorting the firm indices
    //reaches the same ranking in about 1,500.
    //
    //The order is the one the scan produced: by ratio ascending, and among
    //equal ratios the higher firm index first, because Minimum1 keeps the last
    //of equal minima. A ratio that is not a number would make that comparison
    //no longer a strict ordering, which a sort may not be given, so the
    //original passes are kept for that case.
    int sortable=1;
    for (j=1; j<=N2; ++j)
    {
      const double v=DebtServiceToSales2_temp(j);
      if (std::isnan(v)) { sortable=0; break; }
    }

    if (sortable)
    {
      for (j=1; j<=N2; ++j) DS2_order[j-1]=j;
      std::sort(DS2_order.begin(), DS2_order.begin()+N2,
                [](int a, int b)
                {
                  const double va=DebtServiceToSales2_temp(a), vb=DebtServiceToSales2_temp(b);
                  if (va < vb) return true;
                  if (vb < va) return false;
                  return a > b;
                });
      for (j=1; j<=N2; ++j) DS2_rating(DS2_order[j-1],i)=j;
    }
    else
    {
      for (j=1; j<=N2; ++j)
      {
        DS2_min=DebtServiceToSales2_temp.Minimum1(DS2_min_index);
        DS2_rating(DS2_min_index,i)=j;
        DebtServiceToSales2_temp(DS2_min_index)=DebtServiceToSales2_temp.Maximum()+1;
      }
    }
  }
  
  //Set individual loan interest rates to be offered to each borrower
  for(j=1; j<=N2;j++)
  {
    for (i=1; i<=NB; ++i)
    {
      //Only enter here if C-Firm j is a customer of bank i
      if(BankMatch_2(j,i)==1)
      {
        //Calculate the default probability of j as perceived by bank i
        if(NW_2(1,j) == 0)
        {
          FirmDefaultProbability(j)=1;
        }
        else
        {
          //TODO:
          //FirmDefaultProbability(j)=max(floor_default_probability,1-exp(-upsilon*(Loans_2(2,j)+CreditDemand(j))/NW_2(2,j)));
          //This specification should be used when multi-period loans are introduced in the model
          //For now, loans are rolled over every period. So stock of outstanding loans is already included in CreditDemand(j)
          //Tested specifications:
          //FirmDefaultProbability(j)=max(floor_default_probability,1-exp(-upsilon*(CreditDemand(j))/NW_2(1,j)));
          //FirmDefaultProbability(j)=max(floor_default_probability,1-exp(-upsilon*(Pi2(j)/CapitalStock(1,j)));
          FirmDefaultProbability(j)=max(floor_default_probability,1-exp(-upsilon*(DebtServiceToSales2(j))));
        }

        if (flag_rate_setting_loans==0)
        {
          //Original DSK rule whereby mark-up depends on the quartile of the distribution of debt service ratios
          //among customers of i in which j is located
          //Extract the rank of firm j in the ranking of bank i
          k(j)= DS2_rating(j,i);
          //Determine to which quartile it belongs and set the mark-up accordingly
          if(k(j) <= NL_2(i)*0.25)
          {  
            r_deb_h(j)=r_deb(i);
          }
          else if(k(j) > NL_2(i)*0.25 & k(j) <= NL_2(i)*0.5 )
          {              
            r_deb_h(j)=r_deb(i)+k_const;
          }
          else if(k(j) > NL_2(i)*0.5 & k(j) <= NL_2(i)*0.75 )
          {
            r_deb_h(j)=r_deb(i)+2*k_const;
          }
          else
          {
            r_deb_h(j)=r_deb(i)+3*k_const;
          }
        }
        //Alternative new rules whereby the mark-up is additive and depends on perceived default probability
        //or bank's capital adequacy ratio or both
        if (flag_rate_setting_loans==1)
        {
          r_deb_h(j)= max(r_deb(i),r_deb(i) + lambdaB1 * FirmDefaultProbability(j));
        }
        if (flag_rate_setting_loans==2)	
        {	
          r_deb_h(j)= max(r_deb(i),r_deb(i) - lambdaB2*(capitalAdequacyRatio(i)-capitalAdequacyRatioTarget));	
        }	
        if (flag_rate_setting_loans==3)	
        {	
          r_deb_h(j)= max(floor_default_probability,r_deb(i) + lambdaB1 * FirmDefaultProbability(j) - lambdaB2*(capitalAdequacyRatio(i)-capitalAdequacyRatioTarget));	
        }
      }
    }
  }
} 

void BANKING(void)
{
  //Banks determine their profits and pay dividends and taxes
  for(j=1; j<=NB; j++)
  {
    InterestReserves_b(j)=r_cbreserves*Reserves_b(2,j);
    Inflows(j)+=InterestReserves_b(j);
    InterestReserves+=InterestReserves_b(j);

    InterestAdvances_b(j)=r*Advances_b(2,j);
    Outflows(j)+=InterestAdvances_b(j);
    InterestAdvances+=InterestAdvances_b(j);

    BankProfits(j)= LoanInterest(j)+r_bonds*GB_b(2,j)+r_cbreserves*Reserves_b(2,j)-r*Advances_b(2,j)-InterestDeposits(j)-baddebt_b(j)+capitalRecovered(j)-LossEntry_b(j);
        
    if (BankProfits(j)>0) 
    {
      Taxes_b(j)=aliqb*BankProfits(j);
      Taxes+=aliqb*BankProfits(j);
      Dividends_b(j)=db*(BankProfits(j)-Taxes_b(j));
      Deposits(1,j)+=Dividends_b(j);
      Deposits_hb(1,j)+=Dividends_b(j);
      Deposits_h(1)+=Dividends_b(j);
      Dividends(1)+=Dividends_b(j);
      Outflows(j)+=Taxes_b(j);
    }

    //Update net worth and "leverage"	
    NW_b(1,j)=NW_b(2,j)+BankProfits(j)-Dividends_b(j)-Taxes_b(j);	
    //TODO: after introducing more assets held by banks, modify riskWeightedAssets accordingly and create new risk weights	
    //NB: loans' weights vary between 20% and 150% depending on the risk associated, here there are no colaterals, so weight = 100%	
    riskWeightedAssets(j)=riskWeightLoans*Loans_b(1,j) + riskWeightGovBonds*GB_b(1,j);	
    if (NW_b(1,j)>0)	
    {	
      capitalAdequacyRatio(j)=NW_b(1,j)/riskWeightedAssets(j);	
    }	
    else	
    {	
      capitalAdequacyRatio(j)=0;	
    }
  } 

  for(i=1; i<=NB; i++)
  {
    if(Deposits_hb.Row(1).Sum()>0)
    {
      DepositShare_h(i)=Deposits_hb(1,i)/Deposits_hb.Row(1).Sum();
    }
    else
    {
      DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }
  }
}

void BAILOUT(void)
{
  BankEquity_temp=0;

  for (i=1; i<=NB; ++i)
  {
    BankEquity_temp(i)=NW_b(1,i)/(NL_1(i)+NL_2(i));
  }
  
  max_equity=BankEquity_temp.Maximum();
  
  //Failing banks are rescued
  for (j=1; j<=NB; j++)
  {
    //If a bank has previously been bought by another bank, it will be inactive
    if(Bank_active(j)==1)
    {
      //Banks fail if their net worth is negative
      if(NW_b(1,j) < 0)
      { 
        counter_bankfailure+=1;
        //When flagbailout==0, all failing banks are rescued by the government
        if (flagbailout==0) 
        { 
          if(max_equity>0)
          {
            multip_bailout=ran1(p_seed);
            multip_bailout=b1inf+multip_bailout*(b1sup-b1inf);

            if((multip_bailout*max_equity*(NL_1(j)+NL_2(j)))<(capitalAdequacyRatioTarget*Loans_b(1,j)))
            {
              Bailout_b(j)=-NW_b(1,j)+capitalAdequacyRatioTarget*Loans_b(1,j);
            }
            else
            {
              Bailout_b(j)=-NW_b(1,j)+multip_bailout*max_equity*(NL_1(j)+NL_2(j));
            }
            Bailout+=Bailout_b(j);
          }
          else
          {
            multip_bailout=ran1(p_seed);
            multip_bailout=b2inf+multip_bailout*(b2sup-b2inf);
            
            if((multip_bailout*NW_b(2,j))<(capitalAdequacyRatioTarget*Loans_b(1,j)))
            {
              Bailout_b(j)=-NW_b(1,j)+capitalAdequacyRatioTarget*Loans_b(1,j);
            }
            else
            {
              Bailout_b(j)=-NW_b(1,j)+multip_bailout*NW_b(2,j);
            }

            Bailout+=Bailout_b(j);
          }

          Inflows(j)+=Bailout_b(j);
          NW_b(1,j)+=Bailout_b(j);
          capitalAdequacyRatio(j)=NW_b(1,j)/riskWeightedAssets(j);
        }	
        //When flagbailout==1, we check first whether the largest surviving bank can purchase the failing bank
        else
        {
          LossAbsorbed(j)=-NW_b(1,j);

          BankEquity_temp=0;
          for (i=1; i<=NB; ++i)            
          {
            BankEquity_temp(i)=NW_b(1,i); 
          }
          
          maxbank=0;
          max_equity=BankEquity_temp.Maximum1(maxbank);
          // If the failing bank can be purchased by the largest surviving bank, the surviving bank absorbs the negative net worth of the failing bank
		      if (max_equity>0 && (NW_b(1,maxbank) - LossAbsorbed(j)> 0))
          {
			      //All relevant values are added to those of the purchasing bank
            Inflows(maxbank)+=Inflows(j);
            Outflows(maxbank)+=Outflows(j);
            Deposits(1,maxbank)+=Deposits(1,j);
            Deposits(2,maxbank)+=Deposits(2,j);
            Advances_b(1,maxbank)+=Advances_b(1,j);
            Advances_b(2,maxbank)+=Advances_b(2,j);
            Reserves_b(1,maxbank)+=Reserves_b(1,j);
            Reserves_b(2,maxbank)+=Reserves_b(2,j);
            GB_b(1,maxbank)+=GB_b(1,j);
            GB_b(2,maxbank)+=GB_b(2,j);
            Deposits_hb(1,maxbank)+=Deposits_hb(1,j);
            Deposits_hb(2,maxbank)+=Deposits_hb(2,j);
            Deposits_eb(1,maxbank)+=Deposits_eb(1,j);
            Deposits_eb(2,maxbank)+=Deposits_eb(2,j);
            DepositShare_h(maxbank)+=DepositShare_h(j);
            DepositShare_e(maxbank)+=DepositShare_e(j);
            Loans_b(1,maxbank)+=Loans_b(1,j);
            Loans_b(2,maxbank)+=Loans_b(2,j);
            capitalRecovered(maxbank)+=capitalRecovered(j);
            LossEntry_b(maxbank)+=LossEntry_b(j);
            NW_b(1,maxbank)-=LossAbsorbed(j);
            NW_b(2,maxbank)+=NW_b(2,j);
            riskWeightedAssets(maxbank)=riskWeightLoans*Loans_b(1,maxbank) + riskWeightGovBonds*GB_b(1,maxbank);	
            capitalAdequacyRatio(maxbank)=NW_b(1,maxbank)/riskWeightedAssets(maxbank);
            NL_2(maxbank)+=NL_2(j);
            NL_1(maxbank)+=NL_1(j);
            //Failing bank becomes inactive
            Bank_active(j)=0;
            //Purchasing bank receives failing bank's customers
            for (i=1; i<=N2; i++)
            {
              if(BankMatch_2(i,j)==1)
              {
                BankMatch_2(i,maxbank)=1;   
                BankMatch_2(i,j)=0;
                BankingSupplier_2(i)=maxbank;
              }
            }

            for (i=1; i<=N1; i++)
            {
              if(BankMatch_1(i,j)==1)
              {
                BankMatch_1(i,maxbank)=1;   
                BankMatch_1(i,j)=0;
                BankingSupplier_1(i)=maxbank;
              }
            }

            //All relevant variables of failing bank are set to 0
					  Inflows(j)=0;
            Outflows(j)=0;
            Deposits(1,j)=0;
            Advances_b(1,j)=0;
            Reserves_b(1,j)=0;
            GB_b(1,j)=0;
            Deposits_hb(1,j)=0;
            Deposits_eb(1,j)=0;
            DepositShare_h(j)=0;
            DepositShare_e(j)=0;
            Loans_b(1,j)=0;
            capitalRecovered(j)=0;
            LossEntry_b(j)=0;
            NW_b(1,j)=0;
            NL_2(j)=0;
            NL_1(j)=0;
            Deposits(2,j)=0;
            Advances_b(2,j)=0;
            Reserves_b(2,j)=0;
            GB_b(2,j)=0;
            Deposits_hb(2,j)=0;
            Deposits_eb(2,j)=0;
            Loans_b(2,j)=0;
            NW_b(2,j)=0;
            capitalAdequacyRatio(j)=0;
            riskWeightedAssets(j)=0;
          }
          else
          {
            //If the largest bank is unable to save the failing bank, the government steps in
            multip_bailout=0;
            multip_bailout=ran1(p_seed);
            multip_bailout=b2inf+multip_bailout*(b2sup-b2inf);

            BankEquity_temp=0;

            for (i=1; i<=NB; ++i)
            {
              BankEquity_temp(i)=NW_b(1,i)/(NL_1(i)+NL_2(i));
            }
            
            max_equity=BankEquity_temp.Maximum();

            if(max_equity>0)
            {
              multip_bailout=ran1(p_seed);
              multip_bailout=b1inf+multip_bailout*(b1sup-b1inf);

              if((multip_bailout*max_equity*(NL_1(j)+NL_2(j)))<(capitalAdequacyRatioTarget*Loans_b(1,j)))
              {
                Bailout_b(j)=-NW_b(1,j)+capitalAdequacyRatioTarget*Loans_b(1,j);
              }
              else
              {
                Bailout_b(j)=-NW_b(1,j)+multip_bailout*max_equity*(NL_1(j)+NL_2(j));
              }
              Bailout+=Bailout_b(j);
            }
            else
            {
              multip_bailout=ran1(p_seed);
              multip_bailout=b2inf+multip_bailout*(b2sup-b2inf);
              
              if((multip_bailout*NW_b(2,j))<(capitalAdequacyRatioTarget*Loans_b(1,j)))
              {
                Bailout_b(j)=-NW_b(1,j)+capitalAdequacyRatioTarget*Loans_b(1,j);
              }
              else
              {
                Bailout_b(j)=-NW_b(1,j)+multip_bailout*NW_b(2,j);
              }

              Bailout+=Bailout_b(j);
            }

            Inflows(j)+=Bailout_b(j);
            NW_b(1,j)+=Bailout_b(j);
            Bailout+=Bailout_b(j);
            capitalAdequacyRatio(j)=NW_b(1,j)/riskWeightedAssets(j);
          }
        }
      }
    }
  } 

}

void SETTLEMENT(void)
{
  //Calculate central bank profit
  ProfitCB(1)=InterestBonds_cb+InterestAdvances-InterestReserves;
  
  //End of period reserve balances are calculated for each bank by comparing sum of transactions implying outflows of reserves to those implying inflows of reserves
  for (j=1; j<=NB; j++)
  {
    if(Bank_active(j)==1)
    {
      ReserveBalance(j)=Inflows(j)-Outflows(j);
      //If the bank experienced a net inflow of reserves over the period, it repays any CB advances it may have and adds the rest to stock of reserves
      if(ReserveBalance(j)>=0)
      {
        if(Advances_b(1,j)>0)
        {
          if(Advances_b(1,j)>=ReserveBalance(j))
          {
            Advances_b(1,j)-=ReserveBalance(j);
            Advances(1)-=ReserveBalance(j);
            ReserveBalance(j)=0;
          }
          else
          {
            Reserves_b(1,j)+=ReserveBalance(j)-Advances_b(1,j);
            Reserves(1)+=ReserveBalance(j)-Advances_b(1,j);
            Advances(1)-=Advances_b(1,j);
            Advances_b(1,j)=0;
            ReserveBalance(j)=0;
          }
        }
        else
        {
          Reserves_b(1,j)+=ReserveBalance(j);
          Reserves(1)+=ReserveBalance(j);
          ReserveBalance(j)=0;
        }
      }
      //If the bank experienced a net outflow of reserves, it first draws down its stock of reserves and if necessary takes advances from the CB
      else
      {
        if(Reserves_b(1,j)>0)
        {
          if(Reserves_b(1,j)>=(-ReserveBalance(j)))
          {
            Reserves_b(1,j)-=(-ReserveBalance(j));
            Reserves(1)-=(-ReserveBalance(j));
            ReserveBalance(j)=0;
          }
          else
          {
            Advances_b(1,j)+=(-ReserveBalance(j))-Reserves_b(1,j);
            Advances(1)+=(-ReserveBalance(j))-Reserves_b(1,j);
            Reserves(1)-=Reserves_b(1,j);
            Reserves_b(1,j)=0;
            ReserveBalance(j)=0;
          }
        }
        else
        {
          Advances_b(1,j)+=(-ReserveBalance(j));
          Advances(1)+=(-ReserveBalance(j));
          ReserveBalance(j)=0;
        }
      }
    }
  }
}

