
#include "WITCH_input.h"

#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;


void FETCH_WITCH_INPUT(void){
    
    char buf_eleprice[0XFFFF];
    char buf_cabon_price[0XFFFF];
    
    char buf_share_biomass[0XFFFF];
    char buf_share_coal[0XFFFF];
    char buf_share_gas[0XFFFF];
    char buf_share_hydro[0XFFFF];
    char buf_share_nuclear[0XFFFF];
    char buf_share_oil[0XFFFF];
    char buf_share_solar[0XFFFF];
    char buf_share_wind[0XFFFF];
    
    char buf_om_biomass[0XFFFF];
    char buf_om_coal[0XFFFF];
    char buf_om_gas[0XFFFF];
    char buf_om_hydro[0XFFFF];
    char buf_om_nuclear[0XFFFF];
    char buf_om_oil[0XFFFF];
    char buf_om_solar[0XFFFF];
    char buf_om_wind[0XFFFF];
    
    char buf_emissEff_biomass[0XFFFF];
    char buf_emissEff_coal[0XFFFF];
    char buf_emissEff_gas[0XFFFF];
    char buf_emissEff_hydro[0XFFFF];
    char buf_emissEff_nuclear[0XFFFF];
    char buf_emissEff_oil[0XFFFF];
    char buf_emissEff_solar[0XFFFF];
    char buf_emissEff_wind[0XFFFF];
    
    char buf_uc_inv_biomass[0XFFFF];
    char buf_uc_inv_coal[0XFFFF];
    char buf_uc_inv_gas[0XFFFF];
    char buf_uc_inv_hydro[0XFFFF];
    char buf_uc_inv_nuclear[0XFFFF];
    char buf_uc_inv_oil[0XFFFF];
    char buf_uc_inv_solar[0XFFFF];
    char buf_uc_inv_wind[0XFFFF];
    
    char buf_varCost_biomass[0XFFFF];
    char buf_varCost_coal[0XFFFF];
    char buf_varCost_gas[0XFFFF];
    char buf_varCost_hydro[0XFFFF];
    char buf_varCost_nuclear[0XFFFF];
    char buf_varCost_oil[0XFFFF];
    char buf_varCost_solar[0XFFFF];
    char buf_varCost_wind[0XFFFF];

    FILE *fp_eleprice = fopen("modules/WITCH_input/eleprice_interpolated.json", "rb");
    FILE *fp_carbon_price = fopen("modules/WITCH_input/cp_interpolated.json", "rb");
    
    FILE *fp_share_biomass = fopen("modules/WITCH_input/share_biomass_interpolated.json", "rb");
    FILE *fp_share_coal = fopen("modules/WITCH_input/share_coal_interpolated.json", "rb");
    FILE *fp_share_gas = fopen("modules/WITCH_input/share_gas_interpolated.json", "rb");
    FILE *fp_share_hydro = fopen("modules/WITCH_input/share_hydro_interpolated.json", "rb");
    FILE *fp_share_nuclear = fopen("modules/WITCH_input/share_nuclear_interpolated.json", "rb");
    FILE *fp_share_oil = fopen("modules/WITCH_input/share_oil_interpolated.json", "rb");
    FILE *fp_share_solar = fopen("modules/WITCH_input/share_solar_interpolated.json", "rb");
    FILE *fp_share_wind = fopen("modules/WITCH_input/share_wind_interpolated.json", "rb");
    
    FILE *fp_om_biomass = fopen("modules/WITCH_input/fixCost_biomass_interpolated.json", "rb");
    FILE *fp_om_coal = fopen("modules/WITCH_input/fixCost_coal_interpolated.json", "rb");
    FILE *fp_om_gas = fopen("modules/WITCH_input/fixCost_gas_interpolated.json", "rb");
    FILE *fp_om_hydro = fopen("modules/WITCH_input/fixCost_hydro_interpolated.json", "rb");
    FILE *fp_om_nuclear = fopen("modules/WITCH_input/fixCost_nuclear_interpolated.json", "rb");
    FILE *fp_om_oil = fopen("modules/WITCH_input/fixCost_oil_interpolated.json", "rb");
    FILE *fp_om_solar = fopen("modules/WITCH_input/fixCost_solar_interpolated.json", "rb");
    FILE *fp_om_wind = fopen("modules/WITCH_input/fixCost_wind_interpolated.json", "rb");
    
    FILE *fp_emissEff_biomass = fopen("modules/WITCH_input/emissEff_biomass_interpolated.json", "rb");
    FILE *fp_emissEff_coal = fopen("modules/WITCH_input/emissEff_coal_interpolated.json", "rb");
    FILE *fp_emissEff_gas = fopen("modules/WITCH_input/emissEff_gas_interpolated.json", "rb");
    FILE *fp_emissEff_hydro = fopen("modules/WITCH_input/emissEff_hydro_interpolated.json", "rb");
    FILE *fp_emissEff_nuclear = fopen("modules/WITCH_input/emissEff_nuclear_interpolated.json", "rb");
    FILE *fp_emissEff_oil = fopen("modules/WITCH_input/emissEff_oil_interpolated.json", "rb");
    FILE *fp_emissEff_solar = fopen("modules/WITCH_input/emissEff_solar_interpolated.json", "rb");
    FILE *fp_emissEff_wind = fopen("modules/WITCH_input/emissEff_wind_interpolated.json", "rb");
    
    FILE *fp_uc_inv_biomass = fopen("modules/WITCH_input/capcost_biomass_interpolated.json", "rb");
    FILE *fp_uc_inv_coal = fopen("modules/WITCH_input/capcost_coal_interpolated.json", "rb");
    FILE *fp_uc_inv_gas = fopen("modules/WITCH_input/capcost_gas_interpolated.json", "rb");
    FILE *fp_uc_inv_hydro = fopen("modules/WITCH_input/capcost_hydro_interpolated.json", "rb");
    FILE *fp_uc_inv_nuclear = fopen("modules/WITCH_input/capcost_nuclear_interpolated.json", "rb");
    FILE *fp_uc_inv_oil = fopen("modules/WITCH_input/capcost_oil_interpolated.json", "rb");
    FILE *fp_uc_inv_solar = fopen("modules/WITCH_input/capcost_solar_interpolated.json", "rb");
    FILE *fp_uc_inv_wind = fopen("modules/WITCH_input/capcost_wind_interpolated.json", "rb");
    
    FILE *fp_varCost_biomass = fopen("modules/WITCH_input/varCost_biomass_interpolated.json", "rb");
    FILE *fp_varCost_coal = fopen("modules/WITCH_input/varCost_coal_interpolated.json", "rb");
    FILE *fp_varCost_gas = fopen("modules/WITCH_input/varCost_gas_interpolated.json", "rb");
    FILE *fp_varCost_hydro = fopen("modules/WITCH_input/varCost_hydro_interpolated.json", "rb");
    FILE *fp_varCost_nuclear = fopen("modules/WITCH_input/varCost_nuclear_interpolated.json", "rb");
    FILE *fp_varCost_oil = fopen("modules/WITCH_input/varCost_oil_interpolated.json", "rb");
    FILE *fp_varCost_solar = fopen("modules/WITCH_input/varCost_solar_interpolated.json", "rb");
    FILE *fp_varCost_wind = fopen("modules/WITCH_input/varCost_wind_interpolated.json", "rb");

    rapidjson::FileReadStream input_eleprice(fp_eleprice, buf_eleprice, sizeof(buf_eleprice));
    rapidjson::FileReadStream input_carbon_price(fp_carbon_price, buf_cabon_price, sizeof(buf_cabon_price));
     
    rapidjson::FileReadStream input_share_biomass(fp_share_biomass, buf_share_biomass, sizeof(buf_share_biomass));
    rapidjson::FileReadStream input_share_coal(fp_share_coal, buf_share_coal, sizeof(buf_share_coal));
    rapidjson::FileReadStream input_share_gas(fp_share_gas, buf_share_gas, sizeof(buf_share_gas));
    rapidjson::FileReadStream input_share_hydro(fp_share_hydro, buf_share_hydro, sizeof(buf_share_hydro));
    rapidjson::FileReadStream input_share_nuclear(fp_share_nuclear, buf_share_nuclear, sizeof(buf_share_nuclear));
    rapidjson::FileReadStream input_share_oil(fp_share_oil, buf_share_oil, sizeof(buf_share_oil));
    rapidjson::FileReadStream input_share_solar(fp_share_solar, buf_share_solar, sizeof(buf_share_solar));
    rapidjson::FileReadStream input_share_wind(fp_share_wind, buf_share_wind, sizeof(buf_share_wind));
    
    rapidjson::FileReadStream input_om_biomass(fp_om_biomass, buf_om_biomass, sizeof(buf_om_biomass));
    rapidjson::FileReadStream input_om_coal(fp_om_coal, buf_om_coal, sizeof(buf_om_coal));
    rapidjson::FileReadStream input_om_gas(fp_om_gas, buf_om_gas, sizeof(buf_om_gas));
    rapidjson::FileReadStream input_om_hydro(fp_om_hydro, buf_om_hydro, sizeof(buf_om_hydro));
    rapidjson::FileReadStream input_om_nuclear(fp_om_nuclear, buf_om_nuclear, sizeof(buf_om_nuclear));
    rapidjson::FileReadStream input_om_oil(fp_om_oil, buf_om_oil, sizeof(buf_om_oil));
    rapidjson::FileReadStream input_om_solar(fp_om_solar, buf_om_solar, sizeof(buf_om_solar));
    rapidjson::FileReadStream input_om_wind(fp_om_wind, buf_om_wind, sizeof(buf_om_wind));
    
    rapidjson::FileReadStream input_emissEff_biomass(fp_emissEff_biomass, buf_emissEff_biomass, sizeof(buf_emissEff_biomass));
    rapidjson::FileReadStream input_emissEff_coal(fp_emissEff_coal, buf_emissEff_coal, sizeof(buf_emissEff_coal));
    rapidjson::FileReadStream input_emissEff_gas(fp_emissEff_gas, buf_emissEff_gas, sizeof(buf_emissEff_gas));
    rapidjson::FileReadStream input_emissEff_hydro(fp_emissEff_hydro, buf_emissEff_hydro, sizeof(buf_emissEff_hydro));
    rapidjson::FileReadStream input_emissEff_nuclear(fp_emissEff_nuclear, buf_emissEff_nuclear, sizeof(buf_emissEff_nuclear));
    rapidjson::FileReadStream input_emissEff_oil(fp_emissEff_oil, buf_emissEff_oil, sizeof(buf_emissEff_oil));
    rapidjson::FileReadStream input_emissEff_solar(fp_emissEff_solar, buf_emissEff_solar, sizeof(buf_emissEff_solar));
    rapidjson::FileReadStream input_emissEff_wind(fp_emissEff_wind, buf_emissEff_wind, sizeof(buf_emissEff_wind));
    
    rapidjson::FileReadStream input_uc_inv_biomass(fp_uc_inv_biomass, buf_uc_inv_biomass, sizeof(buf_uc_inv_biomass));
    rapidjson::FileReadStream input_uc_inv_coal(fp_uc_inv_coal, buf_uc_inv_coal, sizeof(buf_uc_inv_coal));
    rapidjson::FileReadStream input_uc_inv_gas(fp_uc_inv_gas, buf_uc_inv_gas, sizeof(buf_uc_inv_gas));
    rapidjson::FileReadStream input_uc_inv_hydro(fp_uc_inv_hydro, buf_uc_inv_hydro, sizeof(buf_uc_inv_hydro));
    rapidjson::FileReadStream input_uc_inv_nuclear(fp_uc_inv_nuclear, buf_uc_inv_nuclear, sizeof(buf_uc_inv_nuclear));
    rapidjson::FileReadStream input_uc_inv_oil(fp_uc_inv_oil, buf_uc_inv_oil, sizeof(buf_uc_inv_oil));
    rapidjson::FileReadStream input_uc_inv_solar(fp_uc_inv_solar, buf_uc_inv_solar, sizeof(buf_uc_inv_solar));
    rapidjson::FileReadStream input_uc_inv_wind(fp_uc_inv_wind, buf_uc_inv_wind, sizeof(buf_uc_inv_wind));
    
    rapidjson::FileReadStream input_varCost_biomass(fp_varCost_biomass, buf_varCost_biomass, sizeof(buf_varCost_biomass));
    rapidjson::FileReadStream input_varCost_coal(fp_varCost_coal, buf_varCost_coal, sizeof(buf_varCost_coal));
    rapidjson::FileReadStream input_varCost_gas(fp_varCost_gas, buf_varCost_gas, sizeof(buf_varCost_gas));
    rapidjson::FileReadStream input_varCost_hydro(fp_varCost_hydro, buf_varCost_hydro, sizeof(buf_varCost_hydro));
    rapidjson::FileReadStream input_varCost_nuclear(fp_varCost_nuclear, buf_varCost_nuclear, sizeof(buf_varCost_nuclear));
    rapidjson::FileReadStream input_varCost_oil(fp_varCost_oil, buf_varCost_oil, sizeof(buf_varCost_oil));
    rapidjson::FileReadStream input_varCost_solar(fp_varCost_solar, buf_varCost_solar, sizeof(buf_varCost_solar));
    rapidjson::FileReadStream input_varCost_wind(fp_varCost_wind, buf_varCost_wind, sizeof(buf_varCost_wind));
    
    eleprice_WITCH.ParseStream(input_eleprice);
    carbon_price_WITCH.ParseStream(input_carbon_price);
     
    share_biomass_WITCH.ParseStream(input_share_biomass);
    share_coal_WITCH.ParseStream(input_share_coal);
    share_gas_WITCH.ParseStream(input_share_gas);
    share_hydro_WITCH.ParseStream(input_share_hydro);
    share_nuclear_WITCH.ParseStream(input_share_nuclear);
    share_oil_WITCH.ParseStream(input_share_oil);
    share_solar_WITCH.ParseStream(input_share_solar);
    share_wind_WITCH.ParseStream(input_share_wind);
    
    om_biomass_WITCH.ParseStream(input_om_biomass);
    om_coal_WITCH.ParseStream(input_om_coal);
    om_gas_WITCH.ParseStream(input_om_gas);
    om_hydro_WITCH.ParseStream(input_om_hydro);
    om_nuclear_WITCH.ParseStream(input_om_nuclear);
    om_oil_WITCH.ParseStream(input_om_oil);
    om_solar_WITCH.ParseStream(input_om_solar);
    om_wind_WITCH.ParseStream(input_om_wind);
    
    emissEff_biomass_WITCH.ParseStream(input_emissEff_biomass);
    emissEff_coal_WITCH.ParseStream(input_emissEff_coal);
    emissEff_gas_WITCH.ParseStream(input_emissEff_gas);
    emissEff_hydro_WITCH.ParseStream(input_emissEff_hydro);
    emissEff_nuclear_WITCH.ParseStream(input_emissEff_nuclear);
    emissEff_oil_WITCH.ParseStream(input_emissEff_oil);
    emissEff_solar_WITCH.ParseStream(input_emissEff_solar);
    emissEff_wind_WITCH.ParseStream(input_emissEff_wind);
    
    uc_inv_biomass_WITCH.ParseStream(input_uc_inv_biomass);
    uc_inv_coal_WITCH.ParseStream(input_uc_inv_coal);
    uc_inv_gas_WITCH.ParseStream(input_uc_inv_gas);
    uc_inv_hydro_WITCH.ParseStream(input_uc_inv_hydro);
    uc_inv_nuclear_WITCH.ParseStream(input_uc_inv_nuclear);
    uc_inv_oil_WITCH.ParseStream(input_uc_inv_oil);
    uc_inv_solar_WITCH.ParseStream(input_uc_inv_solar);
    uc_inv_wind_WITCH.ParseStream(input_uc_inv_wind);
    
    varCost_biomass_WITCH.ParseStream(input_varCost_biomass);
    varCost_coal_WITCH.ParseStream(input_varCost_coal);
    varCost_gas_WITCH.ParseStream(input_varCost_gas);
    varCost_hydro_WITCH.ParseStream(input_varCost_hydro);
    varCost_nuclear_WITCH.ParseStream(input_varCost_nuclear);
    varCost_oil_WITCH.ParseStream(input_varCost_oil);
    varCost_solar_WITCH.ParseStream(input_varCost_solar);
    varCost_wind_WITCH.ParseStream(input_varCost_wind);
        
    fclose(fp_eleprice);
    fclose(fp_carbon_price);
     
    fclose(fp_share_biomass);
    fclose(fp_share_coal);
    fclose(fp_share_gas);
    fclose(fp_share_hydro);
    fclose(fp_share_nuclear);
    fclose(fp_share_oil);
    fclose(fp_share_solar);
    fclose(fp_share_wind);
    
    fclose(fp_om_biomass);
    fclose(fp_om_coal);
    fclose(fp_om_gas);
    fclose(fp_om_hydro);
    fclose(fp_om_nuclear);
    fclose(fp_om_oil);
    fclose(fp_om_solar);
    fclose(fp_om_wind);
    
    fclose(fp_emissEff_biomass);
    fclose(fp_emissEff_coal);
    fclose(fp_emissEff_gas);
    fclose(fp_emissEff_hydro);
    fclose(fp_emissEff_nuclear);
    fclose(fp_emissEff_oil);
    fclose(fp_emissEff_solar);
    fclose(fp_emissEff_wind);
    
    fclose(fp_uc_inv_biomass);
    fclose(fp_uc_inv_coal);
    fclose(fp_uc_inv_gas);
    fclose(fp_uc_inv_hydro);
    fclose(fp_uc_inv_nuclear);
    fclose(fp_uc_inv_oil);
    fclose(fp_uc_inv_solar);
    fclose(fp_uc_inv_wind);
    
    fclose(fp_varCost_biomass);
    fclose(fp_varCost_coal);
    fclose(fp_varCost_gas);
    fclose(fp_varCost_hydro);
    fclose(fp_varCost_nuclear);
    fclose(fp_varCost_oil);
    fclose(fp_varCost_solar);
    fclose(fp_varCost_wind);
}

void UPDATE_WITCH_YEAR(int time)
{
    if(flag_constant_WITCH_input==1)
    {
        WITCH_year = WITCH_t0_year;
    }
    else
    {
        if(time <= (burnin+1))
        {
            WITCH_year = WITCH_t0_year;
        }
        else if((time-1-burnin)%4==0)
        {
            WITCH_year += 1;
        }
    }
    if(WITCH_year>2100)
    {
        WITCH_year=2100;
    }
}

void UPDATE_WITCH_PRICE_ADJ(double cpi)
{
    WITCH_price_adj=(cpi/initial_price);
}


double FETCH_ELEPRICE_WITCH (int scenario)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }    
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    return max(0.0,eleprice_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_price*WITCH_price_adj);
}

double FETCH_CARBON_PRICE_WITCH (int scenario)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    return carbon_price_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_price*WITCH_price_adj;
}


std::map<string, double> DISTRIBUTE_ELE_DEMAND (int scenario, double total_ele_demand)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();

    double biomass_demand =  share_biomass_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double coal_demand =  share_coal_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double gas_demand =  share_gas_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double hydro_demand =  share_hydro_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double nuclear_demand =  share_nuclear_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double oil_demand =  share_oil_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double solar_demand =  share_solar_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    double wind_demand =  share_wind_WITCH[year_c].GetArray()[scenario].GetDouble()*total_ele_demand;
    
    std::map<string, double> toReturn;
    
    toReturn["biomass"]=biomass_demand;
    toReturn["coal"]=coal_demand;
    toReturn["gas"]=gas_demand;
    toReturn["hydro"]=hydro_demand;
    toReturn["nuclear"]=nuclear_demand;
    toReturn["oil"]=oil_demand;
    toReturn["solar"]=solar_demand;
    toReturn["wind"]=wind_demand;
    
    return toReturn;  
}
  

std::map<string, double> COMPUTE_ENE_EMISSIONS (int scenario, std::map<string, double> ene_production)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    double biomass_env_filth =  emissEff_biomass_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double coal_env_filth =  emissEff_coal_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double gas_env_filth =  emissEff_gas_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double hydro_env_filth =  emissEff_hydro_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double nuclear_env_filth =  emissEff_nuclear_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double oil_env_filth =  emissEff_oil_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double solar_env_filth =  emissEff_solar_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    double wind_env_filth =  emissEff_wind_WITCH[year_c].GetArray()[scenario].GetDouble()*rescale_ene_env_filth;
    
    std::map<string, double> toReturn;
    
    toReturn["biomass"]=(biomass_env_filth*ene_production["biomass"]);
    toReturn["coal"]=(coal_env_filth*ene_production["coal"]);
    toReturn["gas"]=(gas_env_filth*ene_production["gas"]);
    toReturn["hydro"]=(hydro_env_filth*ene_production["hydro"]);
    toReturn["nuclear"]=(nuclear_env_filth*ene_production["nuclear"]);
    toReturn["oil"]=(oil_env_filth*ene_production["oil"]);
    toReturn["solar"]=(solar_env_filth*ene_production["solar"]);
    toReturn["wind"]=(wind_env_filth*ene_production["wind"]);
    
    return toReturn; 
}


std::map<string, double> FETCH_ELE_UC_INV (int scenario)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    double biomass_inv =  uc_inv_biomass_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double coal_inv =  uc_inv_coal_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double gas_inv =  uc_inv_gas_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double hydro_inv =  uc_inv_hydro_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double nuclear_inv =  uc_inv_nuclear_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double oil_inv =  uc_inv_oil_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double solar_inv =  uc_inv_solar_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double wind_inv =  uc_inv_wind_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    
    std::map<string, double> toReturn;
    
    toReturn["biomass"]=biomass_inv;
    toReturn["coal"]=coal_inv;
    toReturn["gas"]=gas_inv;
    toReturn["hydro"]=hydro_inv;
    toReturn["nuclear"]=nuclear_inv;
    toReturn["oil"]=oil_inv;
    toReturn["solar"]=solar_inv;
    toReturn["wind"]=wind_inv;
    
    return toReturn;
}


std::map<string, double> FETCH_FIX_COST (int scenario)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    
    double biomass_fc =  om_biomass_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double coal_fc =  om_coal_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double gas_fc =  om_gas_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double hydro_fc =  om_hydro_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double nuclear_fc =  om_nuclear_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double oil_fc =  om_oil_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double solar_fc =  om_solar_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    double wind_fc =  om_wind_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price/4;
    
    std::map<string, double> toReturn;
    
    toReturn["biomass"]=biomass_fc;
    toReturn["coal"]=coal_fc;
    toReturn["gas"]=gas_fc;
    toReturn["hydro"]=hydro_fc;
    toReturn["nuclear"]=nuclear_fc;
    toReturn["oil"]=oil_fc;
    toReturn["solar"]=solar_fc;
    toReturn["wind"]=wind_fc;
    
    return toReturn;
}


std::map<string, double> FETCH_VAR_COST (int scenario)
{
    int year=WITCH_year;
    if(year==WITCH_t0_year)
    {
      scenario=0;  
    }
    stringstream strs;
    strs << year;
    string temp_str = strs.str();
    char* year_c = (char*) temp_str.c_str();
    
    double biomass_vc =  varCost_biomass_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double coal_vc =  varCost_coal_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double gas_vc =  varCost_gas_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double hydro_vc =  varCost_hydro_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double nuclear_vc =  varCost_nuclear_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double oil_vc =  varCost_oil_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double solar_vc =  varCost_solar_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    double wind_vc =  varCost_wind_WITCH[year_c].GetArray()[scenario].GetDouble()*WITCH_price_adj*rescale_ene_price;
    
    std::map<string, double> toReturn;
    
    toReturn["biomass"]=biomass_vc;
    toReturn["coal"]=coal_vc;
    toReturn["gas"]=gas_vc;
    toReturn["hydro"]=hydro_vc;
    toReturn["nuclear"]=nuclear_vc;
    toReturn["oil"]=oil_vc;
    toReturn["solar"]=solar_vc;
    toReturn["wind"]=wind_vc;
    
    return toReturn;
}