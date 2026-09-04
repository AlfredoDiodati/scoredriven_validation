#ifndef WITCH_INPUT_H
#define WITCH_INPUT_H

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include <map>

using namespace std;

// FUNCTIONS
void FETCH_WITCH_INPUT(void);
void UPDATE_WITCH_YEAR(int time);
void UPDATE_WITCH_PRICE_ADJ(double cpi);
double FETCH_ELEPRICE_WITCH (int scenario);
double FETCH_ENERGY_EFFICIENCY_WITCH(int scenario, double rescaling_factor);
double FETCH_CARBON_PRICE_WITCH (int scenario);
double COMPUTE_EMISSION_PERUNIT_EN(int scenario, double CO2_rescaling, double q_ele_rescaling);
std::map<string, double> DISTRIBUTE_ELE_DEMAND (int scenario, double total_ele_demand);
std::map<string, double> FETCH_INITIAL_CAPACITIES (std::map<string, double> distributed_energy_demand, int scenario);
std::map<string, double> COMPUTE_ENE_EMISSIONS (int scenario, std::map<string, double> ene_production);
std::map<string, double> FETCH_ELE_UC_INV (int scenario);
std::map<string, double> FETCH_FIX_COST (int scenario);
std::map<string, double> FETCH_VAR_COST (int scenario);

// PARAMETERS & FLAGS
extern int burnin; 
extern int WITCH_t0_year; 
extern int WITCH_year;
extern int flag_constant_WITCH_input;
extern double WITCH_price_adj;
extern double rescale_ene_env_filth;
extern double rescale_ene_price; 

// VARIABLES
extern rapidjson::Document  eleprice_WITCH; 
extern rapidjson::Document  carbon_price_WITCH; 
extern rapidjson::Document  share_biomass_WITCH; 
extern rapidjson::Document  share_coal_WITCH; 
extern rapidjson::Document  share_gas_WITCH; 
extern rapidjson::Document  share_hydro_WITCH; 
extern rapidjson::Document  share_nuclear_WITCH; 
extern rapidjson::Document  share_oil_WITCH; 
extern rapidjson::Document  share_solar_WITCH; 
extern rapidjson::Document  share_wind_WITCH;  
extern rapidjson::Document  om_biomass_WITCH; 
extern rapidjson::Document  om_coal_WITCH; 
extern rapidjson::Document  om_gas_WITCH; 
extern rapidjson::Document  om_hydro_WITCH; 
extern rapidjson::Document  om_nuclear_WITCH; 
extern rapidjson::Document  om_oil_WITCH; 
extern rapidjson::Document  om_solar_WITCH; 
extern rapidjson::Document  om_wind_WITCH; 
extern rapidjson::Document  emissEff_biomass_WITCH; 
extern rapidjson::Document  emissEff_coal_WITCH; 
extern rapidjson::Document  emissEff_gas_WITCH; 
extern rapidjson::Document  emissEff_hydro_WITCH; 
extern rapidjson::Document  emissEff_nuclear_WITCH; 
extern rapidjson::Document  emissEff_oil_WITCH; 
extern rapidjson::Document  emissEff_solar_WITCH; 
extern rapidjson::Document  emissEff_wind_WITCH; 
extern rapidjson::Document  varCost_biomass_WITCH; 
extern rapidjson::Document  varCost_coal_WITCH; 
extern rapidjson::Document  varCost_gas_WITCH; 
extern rapidjson::Document  varCost_hydro_WITCH; 
extern rapidjson::Document  varCost_nuclear_WITCH; 
extern rapidjson::Document  varCost_oil_WITCH; 
extern rapidjson::Document  varCost_solar_WITCH; 
extern rapidjson::Document  varCost_wind_WITCH;
extern rapidjson::Document  uc_inv_biomass_WITCH; 
extern rapidjson::Document  uc_inv_coal_WITCH; 
extern rapidjson::Document  uc_inv_gas_WITCH; 
extern rapidjson::Document  uc_inv_hydro_WITCH; 
extern rapidjson::Document  uc_inv_nuclear_WITCH; 
extern rapidjson::Document  uc_inv_oil_WITCH; 
extern rapidjson::Document  uc_inv_solar_WITCH; 
extern rapidjson::Document  uc_inv_wind_WITCH;
extern double initial_price;

#endif
