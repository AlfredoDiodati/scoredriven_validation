#include "module_energy_sfc.h"
#include "WITCH_input/WITCH_input.h"

void EN_DEM(void)
{
    //Sum up energy demands coming from C-firms and K-firms
    D1_en_TOT=D1_en.Sum();
    D2_en_TOT=D2_en.Sum();
    
    D_en_TOT=ROUND(D1_en_TOT+D2_en_TOT);
} 


void ENERGY(void)
{
    //Determine existing productive capacity for energy
    K_ge=G_ge.Sum(); 
    K_de=G_de.Sum();

    //If existing capacity is insufficient to satisfy demand, expansion investment takes place
    if ((K_ge+K_de)<=D_en_TOT)
    {
        EI_en=D_en_TOT-K_ge-K_de;   
        
        c_de_min=C_de(1)*100000;
        cf_min_ge=CF_ge(1)*100000;

        //Determine best dirty & green technology
        for (tt=1; tt<=t; tt++)
        {
            C_de(tt)=pf/A_de(tt)+t_CO2_en*EM_de(tt);
            if (C_de(tt)<c_de_min)
            {
                c_de_min=C_de(tt);
            }
            
            if (CF_ge(tt)<cf_min_ge) 
            {
                cf_min_ge=CF_ge(tt);
            }
        }

        if(flag_endogenous_exp_quota==1)
        {
            exp_quota=max(0.0, tanh(((c_de_min*payback_en-cf_min_ge)/cf_min_ge)/exp_quota_param));
        }
        
        //If green investment is not constrained, make all investment green if green is superior
        if (flag_energy_exp==0)
        {
            if (c_de_min*payback_en<cf_min_ge)
            {
                EI_en_de=EI_en;
                G_de(t)+=EI_en_de; 
            }
            else
            {
                EI_en_ge=EI_en;
                G_ge(t)+=EI_en_ge;
                IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                G_ge_n(t)+=EI_en_ge*cf_min_ge;
            }
        }
        //Otherwise, may be forced to also invest in dirty even when green is better
        else if(flag_energy_exp==1)
        {
            if (c_de_min*payback_en>cf_min_ge)
            {    
                if (EI_en<(exp_quota*K_gelag+(K_gelag-K_ge)))
                {
                    EI_en_ge=EI_en;
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                }
                else
                {
                    EI_en_ge=(exp_quota*K_gelag+(K_gelag-K_ge));
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                    EI_en_de=EI_en-EI_en_ge;
                    G_de(t)+=EI_en_de; 
                }
            }
            else 
            {
                EI_en_de=EI_en;
                G_de(t)+=EI_en_de; 
            }
        }
        //Case in which share of investment in green is fixed exogenously
        else if (flag_energy_exp==2)
        {
            EI_en_ge=K_ge0_perc*EI_en;
            G_ge(t)+=EI_en_ge;
            IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
            G_ge_n(t)+=EI_en_ge*cf_min_ge;
            EI_en_de=EI_en-EI_en_ge;
            G_de(t)+=EI_en_de; 
        } 
        else
        {
            if (c_de_min*payback_en>cf_min_ge)
            {    
                if (EI_en<(exp_quota*K_gelag+(K_gelag-K_ge)) || (K_delag/(K_gelag+K_delag))<=(exp_quota))
                {
                    EI_en_ge=EI_en;
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                }
                else
                {
                    EI_en_ge=exp_quota*K_gelag+(K_gelag-K_ge);
                    G_ge(t)+=EI_en_ge;
                    IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                    G_ge_n(t)+=EI_en_ge*cf_min_ge;
                    EI_en_de=EI_en-EI_en_ge;
                    G_de(t)+=EI_en_de; 
                }
            }
            else 
            {
                EI_en_ge=K_ge0_perc*EI_en;
                G_ge(t)+=EI_en_ge;
                IC_en_quota(t)=EI_en_ge*cf_min_ge/payback_en;
                G_ge_n(t)+=EI_en_ge*cf_min_ge;
                EI_en_de=EI_en-EI_en_ge;
                G_de(t)+=EI_en_de; 
            }
        }

        // New capacity
        K_ge=G_ge.Sum(); 
        K_de=G_de.Sum();
    }
    
    //Calculate ammortised investment cost from expansion of green energy
    for (tt=1; tt<=t; tt++)
    {
        if (t-tt<=payback_en && G_ge(tt)>0)
        {
            IC_en+=IC_en_quota(tt);
        } 
    }
    
    //Calculate associated labour demand
    LDexp_en=IC_en/w(2);
    
    //Quantity of dirty energy is the residual
    Q_de=D_en_TOT-K_ge; 

    //Only produce dirty energy if green is insufficient
    if (Q_de<=0)
    {
        Q_de=0;
        Q_ge=D_en_TOT;
        //constant marginal cost for green energy
        c_infra=0;
        c_en(1)=mi_en; 
        Emiss_en=0;
    }
    else 
    {
        Q_ge=K_ge;
        for (tt=1; tt<=t; tt++)
        {   
            C_de(tt)=pf/A_de(tt)+t_CO2_en*EM_de(tt);
        }
        G_de_temp=G_de;
        Q_de_temp=Q_de;
        
        //If dirty energy is needed, successively activate dirty plants, starting with the most efficient
        while (Q_de_temp>0)
        {
            c_de_min=C_de(1)*10;
            idmin=1;
            for (tt=1; tt<=t; tt++)
            {
                if (G_de_temp(tt)>0)
                {
                    if (C_de(tt)<=c_de_min)
                    {
                    idmin=tt;
                    c_de_min=C_de(idmin);
                    }
                }
            }
            
            if (Q_de_temp>G_de_temp(idmin))
            {
                PC_en+=G_de(idmin)*C_de(idmin);
                Emiss_en+=G_de(idmin)*EM_de(idmin);
                Q_de_temp-=G_de_temp(idmin);
                FuelCost+=G_de(idmin)*pf/A_de(idmin);
            }
            else
            {
                PC_en+=Q_de_temp*C_de(idmin);
                Emiss_en+=Q_de_temp*EM_de(idmin);
                FuelCost+=Q_de_temp*pf/A_de(idmin);
                Q_de_temp=0;
            }
            G_de_temp(idmin)=0;
            if(Q_de_temp<tolerance)
            {
                Q_de_temp=0;
            }
        }

        //Increasing marginal cost in dirty energy
        c_infra=c_de_min;
        c_en(1)=c_infra+mi_en;    

    }
    
    
    if (flag_share_END==1)
    {
        share_de=K_de/(K_de+K_ge);
    }
    else if (flag_share_END==2)
    {
        share_de=1-Q_ge/D_en_TOT;
    }
    else
    {
        share_de=share_de_0;
    }
    
    //Determine energy sector revenue & R&D expenditures
    Rev_en=c_en(1)*D_en_TOT;

    //Divide between clean and dirty
    if ((Rev_en*share_RD_en)<(Rev_en-IC_en-PC_en))
    {
        RD_en_de=max(0.0,share_RD_en*share_de*Rev_en);
        RD_en_ge=max(0.0,share_RD_en*(1-share_de)*Rev_en);
    }
    else
    {
        RD_en_de=max(0.0,share_de*(Rev_en-PC_en-IC_en));
        RD_en_ge=max(0.0,(1-share_de)*(Rev_en-PC_en-IC_en));
    }
    
    //Compute labour demand for R&D
    LDrd_de=RD_en_de/w(2);     
    LDrd_ge=RD_en_ge/w(2);
    
    parber_en_de=1-exp(-o1_en*RD_en_de);     
    Inn_en_de=bnldev(parber_en_de,1,p_seed);      
    
    parber_en_ge=1-exp(-o1_en*RD_en_ge);    
    Inn_en_ge=bnldev(parber_en_ge,1,p_seed);
   
    //If innovation takes place, update technologies
    if (Inn_en_de == 1)     
    {
        rnd=betadev(b_a11,b_b11,p_seed);          
        rnd=uu1_en+rnd*(uu2_en-uu1_en);                 
        A_de_inn=A_de(t)*(1+rnd);
        
        if (A_de_inn>1)
        {
            A_de_inn=1;
        }
        
        rnd=betadev(b_a11,b_b11,p_seed);          
        rnd=uu1_en+rnd*(uu2_en-uu1_en);                
        EM_de_inn=EM_de(t)*(1-rnd);
        
        if (EM_de_inn<0)
        {
            EM_de_inn=0;
        }

    }
    else
    {
        A_de_inn=A_de(t);
        EM_de_inn=EM_de(t);
    }
    
    if (Inn_en_ge == 1 && CF_ge(t)>0)
    {
        rnd=betadev(b_a11,b_b11,p_seed);         
        rnd=uu1_en+rnd*(uu2_en-uu1_en);
        CF_ge_inn=CF_ge(t)*(1-rnd);
        
        if (CF_ge_inn<0)
        {
            CF_ge_inn=0;
        }
    }
    else
    {
        CF_ge_inn=CF_ge(t);
    }
    
    K_gelag=max(K_ge,(K_ge+K_de)*K_ge0_perc);
    K_delag=K_de;

    if(t<=life_plant)
    {
        G_de(1)-=G_de_0/life_plant;
        G_ge(1)-=G_ge_0/life_plant;
        G_ge_n(1)-=G_ge_n_0/life_plant;
    }

    //Update technology vectors
    if (t < T)
    {
        if (pf/A_de_inn+t_CO2_en*EM_de_inn < pf/A_de(t)+t_CO2_en*EM_de(t))
        {
            A_de(t+1)=A_de_inn;        
            EM_de(t+1)=EM_de_inn;
            C_de(t+1)= pf/A_de_inn+t_CO2_en*EM_de_inn;
        }
        else
        {
            A_de(t+1)=A_de(t);         
            EM_de(t+1)=EM_de(t);
            C_de(t+1)= pf/A_de(t)+t_CO2_en*EM_de(t);
        }

        if (CF_ge_inn < CF_ge(t))
        {
            CF_ge(t+1)=CF_ge_inn;
        }
        else
        {
            CF_ge(t+1)=CF_ge(t);
        }
    
        for (tt=1; tt<=t; tt++)
        {
            if (t-tt>life_plant)
            {
                G_de(tt)=0;
                G_ge(tt)=0;
                G_ge_n(tt)=0;
            }
        }
    }

    //Update nominal value of Energy capital stock; Energy capital valued at production cost --> 0 for dirty
    CapitalStock_e(1)=G_ge_n.Sum();

    //Calculate total emissions
    Emiss_TOT(1)= Emiss_en+Emiss2_TOT+Emiss1_TOT;

    //Calculate total labour demand & wages to be paid by energy sector (paid in next period)
    LDentot=LDexp_en+LDrd_de+LDrd_ge;
    Wages_en=w(2)*LDentot;
}

void ENERGY1()
{
    unit_cost_inv_e_multi_tech=FETCH_ELE_UC_INV(WITCH_scenario);
    uC_prod_mult_tech=FETCH_VAR_COST(WITCH_scenario);
    uC_cap_mult_tech=FETCH_FIX_COST(WITCH_scenario);
    if(t>1)
    {
        for(i=0; i<ene_tecs.size(); i++)
        {
           tech=ene_tecs[i];
           // Capital Depreciates
           capacity_mult_tech[tech]=capacity_mult_tech[tech]-depreciatedCapital_e_multi_tech[ene_tecs[i]](t);
           if(capacity_mult_tech[tech]<-tolerance)
           {
               std::cerr << "\n\n ERROR: Physical capital << tech << energy is negative in "<< t << endl;
                exit(EXIT_FAILURE);
           }

           if((d_en_mult_tech[tech]/en_target_capaciyU)>capacity_mult_tech[tech])
           {
               realInv_e_multi_tech[tech]=(d_en_mult_tech[tech]/en_target_capaciyU)-capacity_mult_tech[tech];
               capacity_mult_tech[tech]+=realInv_e_multi_tech[tech];
               //Set nominal investment across energy technologies
               nom_inv_mult_tech[tech]=unit_cost_inv_e_multi_tech[tech]*realInv_e_multi_tech[tech];
               if(t+e_plant_lifetime[tech]<=T)
               {
                   depreciatedCapital_e_multi_tech[tech](t+e_plant_lifetime[tech])=realInv_e_multi_tech[tech];
               }
           }
           else
           {
               realInv_e_multi_tech[tech]=0;
               nom_inv_mult_tech[tech]=0;
           }
           maintCost_mult_tech[tech]=capacity_mult_tech[tech]*uC_cap_mult_tech[tech];
        }
    }
}

void ENERGY2(void)
{
    LDexp_en=0;
    LDprod_en=0;
    LDmaint_en=0;
    // Update carbon price
    t_CO2_en=FETCH_CARBON_PRICE_WITCH(WITCH_scenario);
    // Ditribute total energy demand across energy technologies
    d_en_mult_tech=DISTRIBUTE_ELE_DEMAND(WITCH_scenario, D_en_TOT);  
    // Compute emissions across energy technologies
    emiss_e_mult_tech=COMPUTE_ENE_EMISSIONS(WITCH_scenario, d_en_mult_tech);
    
    for(i=0; i<ene_tecs.size(); i++)
    {
        tech=ene_tecs[i];
        prodCost_mult_tech[tech]=d_en_mult_tech[tech]*uC_prod_mult_tech[tech];
        if(t==1)
        {
           capacity_mult_tech[tech]=d_en_mult_tech[tech]/en_target_capaciyU;
           for(j=2; j<=21; j++)
           {
            depreciatedCapital_e_multi_tech[tech](j)=capacity_mult_tech[tech]/20;
           }
           maintCost_mult_tech[tech]=capacity_mult_tech[tech]*uC_cap_mult_tech[tech];
        }
        // Compute profits due  to nominal investment 
        inv_redistr_profit_e_multi_tech[tech]=profit_share_energy_inv*nom_inv_mult_tech[tech];
        // Compute labor demand across energy technologies
        LDexp_en_mult_tech[tech]=(1-profit_share_energy_inv)*nom_inv_mult_tech[tech]/w(2);
        LDprod_en_mult_tech[tech]=prodCost_mult_tech[tech]/w(2);
        LDmaint_en_mult_tech[tech]=maintCost_mult_tech[tech]/w(2);
        // Update total labor demand for energy sector
        LDexp_en+=LDexp_en_mult_tech[tech]; 
        LDprod_en+=LDprod_en_mult_tech[tech];
        LDmaint_en+=LDmaint_en_mult_tech[tech];
        // Update wage bill across energy technologies (note: this must equal monetary investment)
        wage_e_mult_tech[tech]=(LDexp_en_mult_tech[tech]+LDprod_en_mult_tech[tech]+LDmaint_en_mult_tech[tech])*w(2);
        
        // Update total emissions in energy sector
        Emiss_en+=emiss_e_mult_tech[tech];
        
        // Update nominal value capital stock across energy technologies. 
        // Unit (in capacity terms) value = production costs of one unit of capacity
        old_capitalStock=capitalStock_e_mult_tech[tech];
        capitalStock_e_mult_tech[tech]=unit_cost_inv_e_multi_tech[tech]*capacity_mult_tech[tech];
        // Compute nominal investment in capacity 
        if(t==1){
            delta_capitalStock_e_mult_tech[tech]=0.0;
        }else{
            delta_capitalStock_e_mult_tech[tech]=capitalStock_e_mult_tech[tech]-old_capitalStock;
        }
    }
    //Calculate total emissions
    Emiss_TOT(1)= Emiss_en+Emiss2_TOT+Emiss1_TOT;
    //Calculate total labour demand in energy sector
    LDentot=LDexp_en+LDprod_en+LDmaint_en;
}

void EMISS_IND(void)
{
    //Calculate total emissions from C-firms and K-firms
    for (i=1; i<=N1; i++)
    {
        Emiss1(i)=A1p_ef(i)/A1p_en(i)*Q1(i);
    }
    
    for (j=1; j<=N2; j++)
    {
        Emiss2(j)= A2e_ef(j)/A2e_en(j)*Q2(j);
    }
    
    Emiss1_TOT=Emiss1.Sum();
    Emiss2_TOT=Emiss2.Sum();
}