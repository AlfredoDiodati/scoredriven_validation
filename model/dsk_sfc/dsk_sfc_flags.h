#ifndef FLAGS_H
#define FLAGS_H

int flag_cum_emissions;                 // Switches between C-Roads climate box and simple cumulative emissions one
                                        // = 0 [BASELINE] C-ROADS
                                        // = 1  simple cumulative emission linear relation with temp

int flag_tax_CO2;                       // Activates C02 tax 
                                        // = 0 [BASELINE] off
                                        // = 1 on and increasing with inflation
                                        // = 2 on and increasing linearly with time
                                        // = 3 on and increasing exponentially with time + inflation correction
                                        // = 4 on and increasing with nominal GDP

int flag_encapshocks;                   // Shocks to energy sector's productive capacity
                                        // = 0 no shock
                                        // = 1 Energy sector loses some percentage of both brown and green capacity

int flag_popshocks;                     // Shocks to the population/labour force
                                        // = 0 no shock
                                        // = 1 reduce labour force by some percentage

int flag_capshocks;                     // Shocks to C-firms' capital stocks
                                        // = 0 no shock
                                        // = 1 x% shocks to capital stocks of all firms 
                                        // = 2 x% shock to aggregate capital stock, affecting firms with uniform probability

int flag_outputshocks;                  // Shocks to current output of C- and K-firms
                                        // = 0 no shock
                                        // = 1 x% shocks to current output of all C- and K-firms 
                                        // = 2 x% shock to aggregate output of both cons. and cap. goods, affecting firms with uniform probability
                                        // = 3 x% shocks to current output of all C-firms 
                                        // = 4 x% shock to aggregate output of cons. goods, affecting firms with uniform probability

int flag_inventshocks;                  // Shocks to inventories of C-firms
                                        // = 0 no shock
                                        // = 1 x% shocks to inventories of all C-firms
                                        // = 2 x% shock to aggregate inventories of C-firms, affecting firms with uniform probability

int flag_prodshocks1;                   // Shocks to productivity affecting the characteristics of capital vintages (roughly similar to TFP shocks in conventional models)
                                        // = 0 no shock
                                        // = 1 On labour productivity of current vintages
                                        // = 2 On energy efficiency of current vintages
                                        // = 3 Both labour producitivity and energy efficiency of current vintages
                                        // = 4 On labour productivity of all existing vintages
                                        // = 5 On energy efficiency of all existing vintages
                                        // = 6 Both labour producitivity and energy efficiency of all existing vintages

int flag_prodshocks2;                   // Shocks to productivity (not affecting characteristics of capital vintages)
                                        // = 0 no shock
                                        // = 1 On labour productivity of C-firms and K-firms
                                        // = 2 On energy efficiency of C-firms and K-firms
                                        // = 3 Both labour producitivity and energy efficiency

int flag_share_END;                     // Switches on endogenous share of R&D expenditures in dirty vs green energy
                                        // = 0 exogenous, given by share_de_0
                                        // = 1 endogenous, given by share of dirty engergy in productive capacity
                                        // = 2 endogenous, given by share of dirty energy produced in t

int flag_energy_exp;                    //Determines whether maximum expansion of green energy capacity per period is constrained
                                        // = 0 not constrained
                                        // = 1 constrained
                                        // = 2 green energy capacity is expanded in order to keep the green share equal to the initial one
                                        // > 2 expansion is constrained but minimum investment in green to keep share equal to initial one

int flag_endogenous_exp_quota;          //Determines whether the constraint on green investment is exogenous or endogenous
                                        // = 0 exogenous
                                        // = 1 endogenous

int flagbailout;                        // Switches between bailout rules for banks
                                        // = 0 [BASELINE] Banks are always bailed out by government
                                        // = 1 Failing banks are purchased by largest bank; government steps in as last resort
                                
int flag_entry;                         // Determines what happens if households cannot finance K-firm and/or C-firm entry
                                        // = 0 The remaining entry costs are paid by government 
                                        // = 1 The remaining entry costs are booked as a loss for the banks

int flag_nonCO2_force;                  // Determines whether or not to consider non-CO2 radiative forcing in the C-Roads climate box
                                        // = 0 Non-CO2 forcing not included
                                        // = 1 [BASELINE] Non-CO2 forcing included

int flag_validation;                    // Determines whether output files needed to construct validation graphs and statistics should be generated
                                        // = 0 No
                                        // = 1 Yes (need to also have fulloutput=0)

int flag_inventories;                   // Switches C-firm inventories on or off
                                        // = 0 Off
                                        // = 1 On

int flag_rate_setting_loans;            // Determines how the interest rate on loans is set
                                        // = 0 depends on an additive markup over the central bank interest rate, and the ranking of demanders according to their debt service-to-sales ratio
                                        // = 1 depends on an additive markup over the central bank interest rate, and the default probability of the demander
                                        // = 2 depends on an additive markup over the central bank interest rate, and the CAR of the supplier
                                        // = 3 depends on an additive markup over the central bank interest rate, the default probability of the demander, and the CAR of the supplier

int flag_WITCH_on;                      // Determines whether WITCH link is on and replaces much of the endogenous energy sector 
                                        // = 0 WITCH link is off
                                        // = 1 WITCH link is on

int flag_constant_WITCH_input;          // Determines whether WITCH inputs are kept constant at 2005 values
                                        // = 0 changing inputs
                                        // = 1 constant inputs

#endif