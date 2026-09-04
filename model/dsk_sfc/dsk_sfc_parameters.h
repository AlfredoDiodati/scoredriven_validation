#ifndef PARS_H
#define PARS_H

int         N1;			                    // Number of K-Firms
int         N2; 			                // Number of C-firms
int         NB;                             // Numbe of Banks
int         T;			                    // Simulation periods

double      varphi; 		                // Maximum bonds to Loans ratio of banks
double      nu;                             // Fraction of revenue devoted to R&D
double      xi;                             // Divides R&D expenditure between innovation and imitation
double      o1;			                    // Effectiveness of R&D expenditure on innovation
double      o2;			                    // Effectiveness of R&D expenditure on imitation
double      uu11;			                // In DSK Case: Lower bound on machine productivity changes due to R&D; in K+S case: Suppport of uniform distribution for exogenous process innovation
double      uu21;			                // In DSK Case: Upper bound on machine productivity changes due to R&D; in K+S case: Suppport of uniform distribution for exogenous product innovation
double      uu12;			                // In DSK Case: Lower bound on own productivity changes due to R&D
double      uu22;			                // In DSK Case: Upper bound on own productivity changes due to R&D
double      uu31;			                // Lower bound on energy efficiency changes due to R&D; C-firms
double      uu41;			                // Upper bound on energy efficiency changes due to R&D; C-firms
double      uu32;			                // Lower bound on energy efficiency changes due to R&D; K-firms
double      uu42;			                // Upper bound on energy efficiency changes due to R&D; K-firms
double      uu51; 		                    // Lower bound on environmental friendliness changes due to R&D; C-firms
double      uu61; 		                    // Upper bound on environmental friendliness changes due to R&D; C-firms
double      uu52; 		                    // Lower bound on environmental friendliness changes due to R&D; K-firms
double      uu62; 		                    // Upper bound on environmental friendliness changes due to R&D; K-firms
double      b_a11;                          // Parameter alpha for Beta distribution governing machine productivity innovation
double      b_a12;                          // Parameter alpha for Beta distribution governing own productivity innovation
double      b_b11;                          // Parameter beta for Beta distribution governing machine productivity innovation
double      b_b12;                          // Parameter beta for Beta distribution governing own productivity innovation
double      b_a2;		                    // Parameter alpha for Beta distribution governing energy efficiency innovation
double      b_b2;		                    // Parameter beta for Beta distribution governing energy efficiency innovation
double      b_a3;		                    // in DSK: Parameter alpha for Beta distribution governing environmental friendliness innovation; in K+S: Parameter alpha for Beta distribution governing change in technological frontier
double      b_b3;		                    // in DSK: Parameter beta for Beta distribution governing environmental friendliness innovation; in K+S: Parameter beta for Beta distribution governing change in technological frontier
double      mi1;                            // Mark-up K-Firms
double      mi2;                            // Initial Mark-up C-firms
double      Gamma;                          // Determining number of potential new clients contacted by K-firms
double      chi;    	                    // Governing replicator dynamics of C-firm market share
double      omega1;		                    // Weight of relative price in C-firm competitiveness
double      omega2;		                    // Weight of unsatisfied demand in C-firm competitiveness
double      psi1;                           // Wage sensitivity to inflation
double      psi2;                           // Wage sensitivity to productivity
double      psi3;                           // Wage sensitivity to unemployment
double      deltami2;                       // Sensitivity of C-firm mark-up to change in market share
double      w_min;	                        // Subsistence wage
double      pmin;		                    // lower bound on prices
double      u;		                        // C-firms' desired capacity utilisation
double      alfa;                           // Persistence of C-firms' adaptive expectations
double      b;                              // Pay-back period for new technologies
double      dim_mach;                       // Output producible by 1 machine
double      agemax;                         // Maximum lifespan of machine tools
double      floor_default_probability;      // Floor level of the probabiity of default of credit demanders (according to the internal risk-based approach of Basel Accords)
double      upsilon;                        // Sensitivity of default probability to leverage ratio
double      lambdaB1;                       // Weight of the credit risk in the interest rate
double      lambdaB2;                       // Weight of the capital adequacy ratio gap in the interest rate
double      riskWeightLoans;                // Risk weight associated to loans with no colateral	
double      riskWeightGovBonds;             // Risk weight associated to government bonds	
double      capitalAdequacyRatioTarget;     // Regulatory Capital Adequacy Ratio target
double      bankmarkdown;                   // Banks' markdown on deposit interest rate
double      centralbankmarkdown;            // CB's markdown on reserve interest rate
double      d1;                             // K-firm dividend payout rate
double      d2;                             // C-firm dividend payout rate
double      db;                             // Bank dividend payout rate
double      repayment_share;                // Loan repayment share
double      bonds_share;                    // Gov. bond repayment share
double      pareto_a;                       // Alpha parameter for pareto distribution function
double      pareto_k;                       // K parameter for pareto distribution function
double      pareto_p;                       // P parameter for pareto distribution function
double      d_cpi_target;                   // Inflation target
double      ustar;		                    // Target unemployment
double      w1sup;                  		// Upper bound on distribution governing deposits transferred to newly entering K-firms
double      w1inf;                  		// Lower bound on distribution governing deposits transferred to newly entering K-firms
double      w2sup;                          // Upper bound on distribution governing deposits transferred to newly entering C-firms
double      w2inf;                          // Lower bound on distribution governing deposits transferred to newly entering C-firms
double      k_const;                       // Governing interest rates charged to individual bank borrowers (additive markup case)
double      aliqw;                          // Tax rate on wages
double      taylor1;                        // Taylor rule inflation sensitivity
double      taylor2;                        // Taylor rule unemployment sensitivity
double      bondsmarkdown;                  // Markdown on gov. bond interest rate
double      mdw;                            // Maximum variation in factors governing changes in wage
double      phi2;                           // Governing maximum amount of loans demanded by C-firms           
double      b1sup;                          // Upper bound on distribution of net worth of bailed-out banks
double      b1inf;                          // Lower bound on distribution of net worth of bailed-out banks
double      b2sup;                          // Upper bound on distribution of net worth of bailed-out banks (for case when all banks have failed)
double      b2inf;                          // Lower bound on distribution of net worth of bailed-out banks (for case when all banks have failed)
double      aliq;                           // Tax rate on firm profits
double      aliqb;                          // Tax rate on bank profits
double      wu;                             // Unemployment benefit rate
double      r_base;                         // Taylor rule intercept
double      de;                             // Dividend payout rate energy sector
double      a1;                             // Propensity to consume out of wage & benefit income
double      a2;                             // Propensity to consume out of dividend & interest income
double      a3;                             // Propensity to consume out of wealth (deposits)
double      a4;                             // Propensity to consume out gov transfers
double      f2_entry_min;                   // Minimum market share for entering C-firms
double      kappa;                          // Persistence of long-run average mean productivity change
double      taylor;                         // Persistence parameter for CB rate adjustment
double      omicron;                        // Desired inventory to output ratio for C-firms
double      I_max;                          // Maximum desired expansion investment as % of current capital stock
double      omega3;                         // Minimum market share (as % of previous share) remaining post-replicator
double      d_f;                            // Share of fossil fuel cost transfered to households
double      g_ls;                           // Growth rate of labour force
double      redistribute_co2TaxRev;         // Share or carbon tax revenue redistributed to households
double      profit_share_energy_inv;        // Profit share of nominal energy investment 

double      share_RD_en;                    // Share of energy sector revenues devoted to R&D
double      share_de_0;                     // Share of energy R&D devoted to dirty energy
int         payback_en;                     // Payback period for green energy plants
int         life_plant;                     // Maximum lifespan of green energy plants
double      exp_quota;                      // Maximum expansion of green energy capacity
double      o1_en;                          // Effectiveness of R&D expenditure in energy
double      uu1_en;		                    // Lower bound for innovation gain in energy sector
double      uu2_en;		                    // Upper bound for innovation gain in energy sector
double      exp_quota_param;                // Governs the evolution of the constraint on green energy investment when flag_endogenous_exp_quota=1
RowVector   a_0;             		        // Governing influence of temperature anomaly on disaster generating function
RowVector   b_0;                            // Governing influence of temperature variance on disaster generating function
int         nshocks;                        // Number of different climate shock channels
int         t_start_climbox;  	            // Period in which the climate box is activated
double      intercept_temp;		            // Intercept of temperature anomaly function
double      slope_temp;			            // Slope of temperature anomaly function
double      tc1;                            // Scaling parameter for linearly-growing carbon tax
double      tc2;                            // Growth rate of carbon tax for exponential case
int         ndep;                           // Number of ocean layers
RowVector   laydep;                         // Depth of ocean layers (in meters); beginning with top layer
double      fertil;                         // Carbon fertilisation effect of atmospheric carbon on NPP 
double      heatstress;                     // Heat stress effect of warming on NPP
double      humtime;                        // Decaying time of humus (in years)
double      biotime;                        // Decay time of carbon in biosphere (in years)
double      humfrac;                        // Fraction of decaying biomass carbon that ends up in humus
double      eddydif;                        // Eddy diffusion coefficient in m^2/year
double      Conref;                         // Reference amount of carbon in upper ocean layer (= pre-industrial) for atmospheric carbon exchange
double      ConrefT;                        // Influence of temperature on Conref
double      rev0;                           // Standard Revelle factor 
double      revC;                           // Impact of carbon content on revelle factor 
int         niterclim;                      // Number of iterations in the carbon exchange between ocean and atmosphere
double      forCO2;                         // Radiative forcing from e-folding CO2 (W/m^2)
double      otherforcefac;                  // Factor by which CO2-induced forcing is multiplied to account for non-CO2 forcing; if this is changed, initial conditions must also be adjusted
double      outrad;                         // W/m^2/K; outgoing radiation per 1K of surface warming 
double      seasurf;                        // Fraction of planetary surface covered by sea
double      heatcap;                        // J/m^3; heat capacity of water
double      secyr;                          // Seconds per year 
int         freqclim;                       // Climate box is called every time after freqclim time-steps have passed in the economic model (set to 4 if econ model is quarterly)
double      g_emiss_global;                 // Exogenous growth rate of global emissions (exluding endogenous EU emissions)
double      emiss_share;                    // Initial share of EU in global emissions
RowVector   shockexponent1;                 // Exponents used for endogenous change in shape of "disaster generating function"
RowVector   shockexponent2;                 // Exponents used for endogenous change in shape of "disaster generating function"
double      a2_nord;                        // Parameter for DICE-type damage function
double      sd_nord;                        // Governing standard deviation of DICE-based shocks at agent level

int WITCH_scenario;                         // The WITCH scenario we want to run. *counting starts from 0*
                                            // 0=NoPolicy
                                            // 1=Npi2020_1000
                                            // 2=Npi2020_500
                                            // 3=Npi2100
double rescale_ene_price;                   // Rescaling factor used for WITCH energy price 
double rescale_ene_env_filth;               // Rescaling factor used for WITCH energy sectors environmental filthiness
int burnin;                                 // How many periods we use to stabilize the model
int WITCH_t0_year;                          // The first year available in WITCH scenarios
double en_target_capaciyU;                  // Target capcity utilization in energy sector
double iota_e;                              // Allows to set a higher or lower dividend payout rate for some energy technologies
int ntechs;                                 // Number of energy technologies
int lifetime_biomass;                       // Maximum lifetime of biomass plants
int lifetime_coal;                          // Maximum lifetime of coal plants
int lifetime_gas;                           // Maximum lifetime of gas plants
int lifetime_hydro;                         // Maximum lifetime of hydro plants
int lifetime_nuclear;                       // Maximum lifetime of nuclear plants
int lifetime_oil;                           // Maximum lifetime of oil plants
int lifetime_solar;                         // Maximum lifetime of solar plants
int lifetime_wind;                          // Maximum lifetime of wind plants

#endif
