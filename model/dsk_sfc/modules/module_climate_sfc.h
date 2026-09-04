#ifndef MODULE_CLIMATE_H
#define MODULE_CLIMATE_H

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <fenv.h>

//#include <string>
#include <string.h>
#include <sstream>

// Include Newmat Libraries and random number generators
#include "../newmat10/include.h"
#include "../newmat10/newmat.h"
#include "../newmat10/newmatio.h"
#include "../auxiliary/betadev.h"
#include "../auxiliary/ran1.h"
#include "../auxiliary/gasdev.h"

// Include functions from other modules
#include "../dsk_sfc_functions.h"


// -- Functions -- //
void CLIMATE_POLICY(void);                          // Sets climate policy (only CO2 tax for now)
void CLIMATEBOX(void);                              // C-Roads climate box
void CLIMATEBOX_CUM_EMISS(void);                    // Simple climate box based on cumulative emissions
void SHOCKS(void);                                  // Imposes repeated endogenous climate shocks
void SHOCKSNORD(void);                              // Imposes repeated endogenous climate shocks following DICE-type damage function
void UPDATECLIMATE(void);                           // Updates climate variables for next period

//-- Flags --//
extern int              flag_tax_CO2;
extern int              flag_nonCO2_force;
extern int              flag_capshocks;
extern int              flag_outputshocks;
extern int              flag_inventshocks;
extern int              flag_encapshocks;
extern int              flag_popshocks;
extern int              flag_prodshocks1;
extern int              flag_prodshocks2;

//-- Inits --//
extern double           A0;
extern double           pm;
extern double           A0_en;

//-- Pars --//
extern int              t_start_climbox;
extern double           intercept_temp;
extern double           slope_temp;
extern int              N1;
extern int              N2;
extern long             *p_seed;
extern RowVector        a_0;
extern RowVector        b_0;
extern double           tc1;
extern double           tc2;
extern int              ndep;
extern RowVector        laydep;
extern double           fertil;
extern double           heatstress;
extern double           humtime;
extern double           biotime;
extern double           humfrac;
extern double           eddydif;
extern double           Conref;
extern double           ConrefT;
extern double           rev0;
extern double           revC;
extern int              niterclim;
extern double           forCO2;
extern double           otherforcefac;
extern double           outrad;
extern double           secyr;
extern double           seasurf;
extern double           heatcap;
extern int              freqclim;
extern double           Tmixedinit1;
extern double           g_emiss_global;
extern double           emiss_share;
extern RowVector        shockexponent1;
extern RowVector        shockexponent2;
extern double           a2_nord;
extern double           sd_nord;

//--Inits--//
extern double           Emiss_yearly_0;
extern double           Cum_emissions_0;
extern double           t_CO2_0;
extern double           t_CO2_en_0;
extern double           NPP0;
extern double           Cat0;

// -- Vars -- //
extern int              t;
extern int              tt;
extern int              i;
extern int              j;
extern int              iterations;
extern int              t0;
extern double           rnd;
extern RowVector        Emiss_yearly;  
extern RowVector        Emiss_TOT;
extern double           g_rate_em_y;
extern RowVector        Emiss_yearly_calib;
extern double           NPP; 
extern RowVector        Tmixed; 
extern double		  	Cum_emissions;
extern RowVector        shocks_machprod;       
extern RowVector        shocks_techprod;       
extern RowVector        shocks_encapstock_de;   
extern RowVector        shocks_encapstock_ge;   
extern RowVector        shocks_capstock;       
extern RowVector        shocks_invent;            
extern RowVector        shocks_labprod1;     
extern RowVector        shocks_labprod2;     
extern RowVector        shocks_eneff1;     
extern RowVector        shocks_eneff2; 
extern RowVector        shocks_output1;     
extern RowVector        shocks_output2; 
extern double           shock_pop;                   
extern RowVector        A1;
extern RowVector        A1p;
extern RowVector        A1_en;   
extern RowVector        A1p_en; 
extern RowVector        G_de;
extern RowVector        G_ge; 
extern RowVector        G_ge_n; 
extern double           LS;
extern RowVector        X_a;                    
extern RowVector        X_b; 
extern RowVector        CapitalStock_e;
extern double           t_CO2;
extern double           t_CO2_en;
extern RowVector        cpi;
extern double           cpi_init;
extern double           GDP_init;
extern RowVector        GDP_n;
extern double           Emiss_gauge;
extern RowVector        Cat;
extern double           humrelease;
extern RowVector        hum;
extern double           biorelease;
extern RowVector        biom;
extern double           Cat1;
extern double           dCat1;
extern RowVector        fluxC;
extern Matrix           Con; 
extern Matrix           Hon;  
extern Matrix           Ton; 
extern double           Con1;
extern double           Ctot1;
extern RowVector        Cax;                                        
extern RowVector        Caxx;                                       
extern RowVector        Cay;                                    
extern RowVector        Cayy;                            
extern RowVector        Caa;           
extern double           FCO2;                                       
extern double           Fin;                                        
extern double           Fout;     
extern RowVector        fluxH;
extern double           Emiss_global;
extern Matrix           A;
extern Matrix           A_en;
extern int              nshocks;
extern double           shock_nord;

#endif