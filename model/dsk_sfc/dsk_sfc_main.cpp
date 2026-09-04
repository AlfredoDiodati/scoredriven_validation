#include "dsk_sfc_include.h"
#include "dsk_sfc_reductions.h"
using namespace std;

int main(int argc, char *argv[])
{
  //Command line parser
  CLI::App app{"DSK_SFC, the Dystopian Schumpeter meeting Keynes Stock Flow Consistent model"};

  //Add command line input: json input file
  string inputstring = "default.json";
  app.add_option("inputfile", inputstring, "A path to the json input file")
    ->required()
    ->check(CLI::ExistingFile);

  //Add command line input: run name
  string str_runname;
  app.add_option("-r,--run", str_runname, "A name for the run")
    ->default_val("test");
  
  //Add command line input: seed
  int exseed{1};
  app.add_option("-s", exseed, "A seed (positive integer)")
    ->default_val(1)
    ->check(CLI::PositiveNumber);

  //Add command line input: full output dummy
  int fullout{0};
  app.add_option("-f,--fulloutput", fullout, "If set to 1, full output will be saved");

  //Add command line input: error printing dummy
  int cerr_enabled{0};
  app.add_option("-c,--cerr", cerr_enabled, "If set to 1, print error messages to the console");
  
  //Add command line input: verbose dummy
  int verbose{0};
  app.add_option("-v,--verbose", verbose, "If set to 1, print simulation progress");

  //Parse command line input
  CLI11_PARSE(app, argc, argv);
  //If run name contains spaces, remove them
  auto noSpace = std::remove(str_runname.begin(), str_runname.end(), ' ');
  str_runname.erase(noSpace, str_runname.end());
  
  //If dummy for error printing to console is set to 0, suppress error messages
  if(cerr_enabled==0){std::cerr.setstate(std::ios_base::failbit);}

  //Initialise variables needed to process console inputs
  std::string seedstring = to_string(exseed);
  char const* seednumber = seedstring.c_str();
  char* exec_dir;
  char const* runname = str_runname.c_str();
  fulloutput=fullout;

  //Path to executable as char; needed to create output folder
  exec_dir = argv[0];
  
  //Seed needs to be converted to negative integer to work with functions generating random draws
  exseed = -exseed;

  //Read JSON file and convert it to a document from which to read parameter, initial and flag values
  using namespace rapidjson;
  std::ifstream ifs { inputstring };
  //Exit if unable to open input file
  if ( !ifs.is_open() )
  {
    cout << "Could not open input file for reading!"<< endl;
    return EXIT_FAILURE;
  }
  //Parse
  IStreamWrapper isw { ifs };
  Document inputs {};
  inputs.ParseStream( isw );
  //Exit if encounter error during parsing
  if (inputs.HasParseError())
  {
    cout << "Parsing of input file failed!"<< endl;
    return EXIT_FAILURE;
  }
  if(verbose){cout << "Finished parsing input file" << endl;}

  //Use the parsed input file to set model parameters, initial values, flags, etc.
  SETPARAMS(inputs);
  if(verbose){cout << "Exiting function SETPARAMS" << endl;}

  //If they do not exist yet, create folders to hold simulation output and error logs
  FOLDERS(exec_dir);
  if(verbose){cout << "Finished creating output folders" << endl;}

  //Generate a string containing the name of the run and the current seed
  char filedesc[32];
  strcpy(filedesc,"_");
  char const* desc=strcat(filedesc,runname);				 
  desc=strcat(filedesc,"_");					
  desc=strcat(filedesc,seednumber);			

  //Path to the output folder
  //char pathname[60];
  std::string outstr(exec_dir);
  for(j=outstr.length(); j>0; j--)
  {
      if(outstr[j-1]=='/')
      {
          break;
      } 
      else
      {
          outstr.pop_back();
      }
  }
  outstr+="output";
  //One past the length, for the terminator strcpy writes. Sized exactly at
  //the length this wrote one byte past the end of the array, which an
  //unoptimised build happened to absorb and an optimised one does not.
  char pathname[outstr.length()+1];
  char* filepath=strcpy(pathname,outstr.c_str());
  //char* filepath=strcpy(pathname, outstr.c_str());

  //Generate the name of the error file
  GENFILEERRORS(filepath,"/errors/Errors",desc);
  //Enter here if dummy is set to save more extensive model output
  if (fulloutput==1)
  {
    //Call functions to generate file names for the various output files
    GENFILEYMC(filepath,"/results", desc);
    GENFILEPRODALL1_en(filepath,"/A1all_en", desc);
    GENFILEPRODALL2_en(filepath,"/A2all_en", desc);
    GENFILEPRODALL1_ef(filepath,"/A1all_ef", desc);
    GENFILEPRODALL2_ef(filepath,"/A2all_ef", desc);
    GENFILEPROD1(filepath,"/A1", desc);
    GENFILEPROD2(filepath,"/A2", desc);
    GENFILEPRODALL1(filepath,"/A1_all", desc);
    GENFILEPRODALL2(filepath,"/A2_all", desc);
    GENFILENWALL1(filepath,"/NW1all",desc);
    GENFILENWALL2(filepath,"/NW2all",desc);
    GENFILENWALL3(filepath,"/NWBall",desc);
    GENFILEDEBALL2(filepath,"/Deb2all",desc);
  }
  else if(flag_validation==1)
  {
    //If performing validation run, generate names for all validation output files
    GENFILEVALIDATION1(filepath,"/validation1", seednumber);
    GENFILEVALIDATION2(filepath,"/validation2", seednumber);
    GENFILEVALIDATION3(filepath,"/validation3", seednumber);
    GENFILEVALIDATION4(filepath,"/validation4", seednumber);
    GENFILEVALIDATION5(filepath,"/validation5", seednumber);
    GENFILEVALIDATION6(filepath,"/validation6", seednumber);
    GENFILEVALIDATION7(filepath,"/validation7", seednumber);
    GENFILEVALIDATION8(filepath,"/validation8", seednumber);
    GENFILEVALIDATION9(filepath,"/validation9", seednumber);
    GENFILEVALIDATION10(filepath,"/validation10", seednumber);
    GENFILEVALIDATION11(filepath,"/validation11", seednumber);
    GENFILEVALIDATION12(filepath,"/validation12", seednumber);
  }
  else
  {
    //Otherwise, create name for standard output file
    GENFILEYMC(filepath,"/results", desc);
  }

  if(verbose){cout << "Finished creating output file names" << endl;}

  //Resize variable arrays to the needed dimensions based on the numbers of agents and periods
  RESIZE();
  if(verbose){cout << "Exiting function RESIZE" << endl;}

  //Only enter here if the soft-link is activated
  if(flag_WITCH_on==1)
  {
    //Load WITCH data
    FETCH_WITCH_INPUT();
  }

  //Initialise endogenous model variables
	INITIALIZE(exseed);
  if(verbose){cout << "Exiting function INITIALIZE; Entering simulation loop" << endl;}

  //This only works on Linux systems (not Windows or Mac)
  //Mostly needed for large batch runs on HPC to avoid getting stuck
  #ifdef __linux__
    //Set the run to time out after T*2 seconds
    signal(SIGALRM, catchAlarm);
    alarm(T*2); 
  #endif

  //Enter loop over simulation periods
  Errors.open(errorfilename,ios::app);
  for (t=1; t<=T; t++)
  {
    
    if(flag_WITCH_on==1)
    {
      UPDATE_WITCH_YEAR(t);
      UPDATE_WITCH_PRICE_ADJ(cpi(2));
      c_en(1)=FETCH_ELEPRICE_WITCH(WITCH_scenario);
    }

    if(verbose){cout << "Entering simulation period " << t << endl;}
    //Resets endogenous variables as needed to begin a new simulation period
    SETVARS();
    if(verbose){cout << "Exiting function SETVARS in period " << t << endl;}
    //Interest on bank deposits is paid
    DEPOSITINTEREST();
    if(verbose){cout << "Exiting function DEPOSITINTEREST in period " << t << endl;}
    //only enter here every freqclim periods
    if(t%freqclim==0 && flag_WITCH_on==0)
    {
      //Update climate policy (only carbon tax for now)
      CLIMATE_POLICY();
      if(verbose){cout << "Exiting function CLIMATE_POLICY in period " << t << endl;}
    }

    //Deliver machines ordered in t-1, calculate unit cost for C and K-Firms, update prices
    MACH();
    if(verbose){cout << "Exiting function MACH in period " << t << endl;}

    //Banks determine maximum amount they are willing to lend
    TOTCREDIT();
    if(verbose){cout << "Exiting function TOTCREDIT in period " << t << endl;}

    //Banks determine the loan interest rate to be charged to each C-Firm customer
    LOANRATES();
    if(verbose){cout << "Exiting function LOANRATES in period " << t << endl;}

    //K-Firms send out brochures to attract customers
    BROCHURE();
    if(verbose){cout << "Exiting function BROCHURE in period " << t << endl;}

    //C-Firms set expected demand & desired production and determine desired investment
    INVEST();
    if(verbose){cout << "Exiting function INVEST in period " << t << endl;}
    
    if(flag_WITCH_on==1)
    {
      //Load WITCH data
      ENERGY1();
      if(verbose){cout << "Exiting function ENERGY1 in period " << t << endl;}
    }

    ALLOCATECREDIT();
    if(verbose){cout << "Exiting function ALLOCATECREDIT in period " << t << endl;}

    PRODMACH();
    if(verbose){cout << "Exiting function PRODMACH in period " << t << endl;}

    EMISS_IND();
    if(verbose){cout << "Exiting function EMISS_IND in period " << t << endl;}

    if(flag_WITCH_on==1)
    {
      //Load WITCH data
      ENERGY2();
      if(verbose){cout << "Exiting function ENERGY2 in period " << t << endl;}
    }
    else
    {
      ENERGY();
      if(verbose){cout << "Exiting function ENERGY in period " << t << endl;}
    }

    PAY_LAB_INV();
    if(verbose){cout << "Exiting function PAY_LAB_INV in period " << t << endl;}

    COMPET2();
    if(verbose){cout << "Exiting function COMPET2 in period " << t << endl;}

    PROFIT();
    if(verbose){cout << "Exiting function PROFIT in period " << t << endl;}

    MACRO();
    if(verbose){cout << "Exiting function MACRO in period " << t << endl;}
    
    ENTRYEXIT();
    if(verbose){cout << "Exiting function ENTRYEXIT in period " << t << endl;}

    BANKING();
    if(verbose){cout << "Exiting function BANKING in period " << t << endl;}

    BAILOUT();
    if(verbose){cout << "Exiting function BAILOUT in period " << t << endl;}

    GOV_BUDGET();
    if(verbose){cout << "Exiting function GOV_BUDGET in period " << t << endl;}

    TAYLOR();
    if(verbose){cout << "Exiting function TAYLOR in period " << t << endl;}

    SETTLEMENT();
    if(verbose){cout << "Exiting function SETTLEMENT in period " << t << endl;}

    TECHANGEND();
    if(verbose){cout << "Exiting function TECHANGEND in period " << t << endl;}

    if(t>t_start_climbox && t%freqclim==0)
    {
      if (flag_cum_emissions==0)
      {
        
        CLIMATEBOX();
        if(iterations==niterclim)
        {
          std::cerr << "\n\n ERROR: Carbon content loop did not converge in period " << t << endl;
          Errors << "\n Carbon content loop did not converge in period " << t << endl;
        }
        
        if(verbose){cout << "Exiting function CLIMATEBOX in period " << t << endl;}
      }
      else
      {
        CLIMATEBOX_CUM_EMISS();
        if(verbose){cout << "Exiting function CLIMATEBOX_CUM_EMISS in period " << t << endl;}
      }
    }

    if(a2_nord>0)
    {
      SHOCKSNORD();
      if(verbose){cout << "Exiting function SHOCKSNORD in period " << t << endl;}
    }
    else
    {
      SHOCKS();
      if(verbose){cout << "Exiting function SHOCKS in period " << t << endl;}
    }

    DEPOSITCHECK();
    if(verbose){cout << "Exiting function DEPOSITCHECK in period " << t << endl;}

    NEGATIVITYCHECK();
    if(verbose){cout << "Exiting function NEGATIVITYCHECK in period " << t << endl;}

    CHECKSUMS();
    if(verbose){cout << "Exiting function CHECKSUMS in period " << t << endl;}

    ADJUSTSTOCKS();
    if(verbose){cout << "Exiting function ADJUSTSTOCKS in period " << t << endl;}

    SFC_CHECK();
    if(verbose){cout << "Exiting function SFC_CHECK in period " << t << endl;}

    SAVE();
    if(verbose){cout << "Exiting function SAVE in period " << t << endl;}

    UPDATE();
    if(verbose){cout << "Exiting function UPDATE in period " << t << endl;}

    UPDATECLIMATE();
    if(verbose){cout << "Exiting function UPDATECLIMATE in period " << t << endl;}

    OVERBOOST();

    if(verbose){cout << "Exiting function OVERBOOST; end of period " << t << endl;}
  }
  Errors.close();
  if(verbose){cout << "End of simulation loop"<< endl;}
  return(EXIT_SUCCESS);
}

///////////MODEL FUNCTIONS/////////////////////////////

void SETPARAMS(const rapidjson::Document& inputs)
{
      //Read the values of parameters, flags, initial values and one-off shocks from the document "inputs" and set the respective variables to those values
      N1=inputs["params"][0]["N1"].GetInt();
      N2=inputs["params"][0]["N2"].GetInt();
      NB=inputs["params"][0]["NB"].GetInt();
      T=inputs["params"][0]["T"].GetInt();
      varphi=inputs["params"][0]["varphi"].GetDouble();
      nu=inputs["params"][0]["nu"].GetDouble();
      xi=inputs["params"][0]["xi"].GetDouble();
      o1=inputs["params"][0]["o1"].GetDouble();
      o2=inputs["params"][0]["o2"].GetDouble();
      uu11=inputs["params"][0]["uu11"].GetDouble();
      uu21=inputs["params"][0]["uu21"].GetDouble();
      uu12=inputs["params"][0]["uu12"].GetDouble();
      uu22=inputs["params"][0]["uu22"].GetDouble();
      uu31=inputs["params"][0]["uu31"].GetDouble();
      uu41=inputs["params"][0]["uu41"].GetDouble();
      uu32=inputs["params"][0]["uu32"].GetDouble();
      uu42=inputs["params"][0]["uu42"].GetDouble();
      uu51=inputs["params"][0]["uu51"].GetDouble();
      uu61=inputs["params"][0]["uu61"].GetDouble();
      uu52=inputs["params"][0]["uu52"].GetDouble();
      uu62=inputs["params"][0]["uu62"].GetDouble();
      b_a11=inputs["params"][0]["b_a11"].GetDouble();
      b_b11=inputs["params"][0]["b_b11"].GetDouble();
      b_a12=inputs["params"][0]["b_a12"].GetDouble();
      b_b12=inputs["params"][0]["b_b12"].GetDouble();
      b_a2=inputs["params"][0]["b_a2"].GetDouble();
      b_b2=inputs["params"][0]["b_b2"].GetDouble();
      b_a3=inputs["params"][0]["b_a3"].GetDouble();
      b_b3=inputs["params"][0]["b_b3"].GetDouble();
      mi1=inputs["params"][0]["mi1"].GetDouble();
      mi2=inputs["params"][0]["mi2"].GetDouble();
      Gamma=inputs["params"][0]["Gamma"].GetDouble();
      chi=inputs["params"][0]["chi"].GetDouble();
      omega1=inputs["params"][0]["omega1"].GetDouble();
      omega2=inputs["params"][0]["omega2"].GetDouble();
      psi1=inputs["params"][0]["psi1"].GetDouble();
      psi2=inputs["params"][0]["psi2"].GetDouble();
      psi3=inputs["params"][0]["psi3"].GetDouble();
      deltami2=inputs["params"][0]["deltami2"].GetDouble();
      w_min=inputs["params"][0]["w_min"].GetDouble();
      pmin=inputs["params"][0]["pmin"].GetDouble();
      u=inputs["params"][0]["u"].GetDouble();
      alfa=inputs["params"][0]["alfa"].GetDouble();
      b=inputs["params"][0]["b"].GetDouble();
      dim_mach=inputs["params"][0]["dim_mach"].GetDouble();
      agemax=inputs["params"][0]["agemax"].GetDouble();
      floor_default_probability=inputs["params"][0]["floor_default_probability"].GetDouble();
      upsilon=inputs["params"][0]["upsilon"].GetDouble();
      lambdaB1=inputs["params"][0]["lambdaB1"].GetDouble();
      lambdaB2=inputs["params"][0]["lambdaB2"].GetDouble();
      riskWeightLoans=inputs["params"][0]["riskWeightLoans"].GetDouble();	
      riskWeightGovBonds=inputs["params"][0]["riskWeightGovBonds"].GetDouble();	
      capitalAdequacyRatioTarget=inputs["params"][0]["capitalAdequacyRatioTarget"].GetDouble();
      bankmarkdown=inputs["params"][0]["bankmarkdown"].GetDouble();  
      centralbankmarkdown=inputs["params"][0]["centralbankmarkdown"].GetDouble();
      d1=inputs["params"][0]["d1"].GetDouble();
      d2=inputs["params"][0]["d2"].GetDouble();
      db=inputs["params"][0]["db"].GetDouble();
      repayment_share=inputs["params"][0]["repayment_share"].GetDouble();
      bonds_share=inputs["params"][0]["bonds_share"].GetDouble();
      pareto_a=inputs["params"][0]["pareto_a"].GetDouble();
      pareto_k=inputs["params"][0]["pareto_k"].GetDouble();
      pareto_p=inputs["params"][0]["pareto_p"].GetDouble();
      d_cpi_target=inputs["params"][0]["d_cpi_target"].GetDouble();
      ustar=inputs["params"][0]["ustar"].GetDouble();
      w1sup=inputs["params"][0]["w1sup"].GetDouble();
      w1inf=inputs["params"][0]["w1inf"].GetDouble();
      w2sup=inputs["params"][0]["w2sup"].GetDouble();
      w2inf=inputs["params"][0]["w2inf"].GetDouble();
      k_const=inputs["params"][0]["k_const"].GetDouble();
      aliqw=inputs["params"][0]["aliqw"].GetDouble(); 
      taylor1=inputs["params"][0]["taylor1"].GetDouble();
      taylor2=inputs["params"][0]["taylor2"].GetDouble();
      bondsmarkdown=inputs["params"][0]["bondsmarkdown"].GetDouble();
      mdw=inputs["params"][0]["mdw"].GetDouble();
      phi2=inputs["params"][0]["phi2"].GetDouble();
      b1sup=inputs["params"][0]["b1sup"].GetDouble();
      b1inf=inputs["params"][0]["b1inf"].GetDouble();
      b2sup=inputs["params"][0]["b2sup"].GetDouble();
      b2inf=inputs["params"][0]["b2inf"].GetDouble();
      aliq=inputs["params"][0]["aliq"].GetDouble();
      aliqb=inputs["params"][0]["aliqb"].GetDouble();
      wu=inputs["params"][0]["wu"].GetDouble();
      de=inputs["params"][0]["de"].GetDouble();
      a1=inputs["params"][0]["a1"].GetDouble();
      a2=inputs["params"][0]["a2"].GetDouble();
      a3=inputs["params"][0]["a3"].GetDouble();
      a4=inputs["params"][0]["a4"].GetDouble();
      f2_entry_min=inputs["params"][0]["f2_entry_min"].GetDouble();
      en_target_capaciyU=inputs["params"][0]["en_target_capaciyU"].GetDouble();
      kappa=inputs["params"][0]["kappa"].GetDouble();
      taylor=inputs["params"][0]["taylor"].GetDouble();
      omicron=inputs["params"][0]["omicron"].GetDouble();
      I_max=inputs["params"][0]["I_max"].GetDouble();
      omega3=inputs["params"][0]["omega3"].GetDouble();
      de_premium["biomass"]=inputs["params"][0]["iota_e"].GetDouble();
      de_premium["coal"]=0;
      de_premium["gas"]=0;
      de_premium["hydro"]=inputs["params"][0]["iota_e"].GetDouble();
      de_premium["nuclear"]=0;
      de_premium["oil"]=0;
      de_premium["solar"]=inputs["params"][0]["iota_e"].GetDouble();
      de_premium["wind"]=inputs["params"][0]["iota_e"].GetDouble();
      profit_share_energy_inv=inputs["params"][0]["profit_share_energy_inv"].GetDouble();
      redistribute_co2TaxRev=inputs["params"][0]["redistribute_co2TaxRev"].GetDouble();
      d_f=inputs["params"][0]["d_f"].GetDouble();
      g_ls=inputs["params"][0]["g_ls"].GetDouble();

      share_RD_en=inputs["climparams"][0]["share_RD_en"].GetDouble();
      share_de_0=inputs["climparams"][0]["share_de_0"].GetDouble();
      payback_en=inputs["climparams"][0]["payback_en"].GetInt();
      life_plant=inputs["climparams"][0]["life_plant"].GetInt();
      exp_quota=inputs["climparams"][0]["exp_quota"].GetDouble();
      o1_en=inputs["climparams"][0]["o1_en"].GetDouble();
      uu1_en=inputs["climparams"][0]["uu1_en"].GetDouble();
      uu2_en=inputs["climparams"][0]["uu2_en"].GetDouble();
      exp_quota_param=inputs["params"][0]["exp_quota_param"].GetDouble();
      t_start_climbox=inputs["climparams"][0]["t_start_climbox"].GetInt();
      intercept_temp=inputs["climparams"][0]["intercept_temp"].GetDouble();
      slope_temp=inputs["climparams"][0]["slope_temp"].GetDouble();
      tc1=inputs["climparams"][0]["tc1"].GetDouble();
      tc2=inputs["climparams"][0]["tc2"].GetDouble();
      ndep=inputs["climparams"][0]["ndep"].GetInt();
      laydep.ReSize(ndep);
      laydep(1)=inputs["climparams"][0]["laydep1"].GetDouble();
      laydep(2)=inputs["climparams"][0]["laydep2"].GetDouble();
      laydep(3)=inputs["climparams"][0]["laydep3"].GetDouble();
      laydep(4)=inputs["climparams"][0]["laydep4"].GetDouble();
      laydep(5)=inputs["climparams"][0]["laydep5"].GetDouble();
      fertil=inputs["climparams"][0]["fertil"].GetDouble();
      heatstress=inputs["climparams"][0]["heatstress"].GetDouble();
      humtime=inputs["climparams"][0]["humtime"].GetDouble();
      biotime=inputs["climparams"][0]["biotime"].GetDouble();
      humfrac=inputs["climparams"][0]["humfrac"].GetDouble();
      eddydif=inputs["climparams"][0]["eddydif"].GetDouble();
      ConrefT=inputs["climparams"][0]["ConrefT"].GetDouble();
      rev0=inputs["climparams"][0]["rev0"].GetDouble();
      revC=inputs["climparams"][0]["revC"].GetDouble();
      niterclim=inputs["climparams"][0]["niterclim"].GetInt();
      forCO2=inputs["climparams"][0]["forCO2"].GetDouble();
      otherforcefac=inputs["climparams"][0]["otherforcefac"].GetDouble();
      outrad=inputs["climparams"][0]["outrad"].GetDouble();
      secyr=inputs["climparams"][0]["secyr"].GetDouble();
      seasurf=inputs["climparams"][0]["seasurf"].GetDouble();
      heatcap=inputs["climparams"][0]["heatcap"].GetDouble();
      freqclim=inputs["climparams"][0]["freqclim"].GetInt();
      g_emiss_global=inputs["climparams"][0]["g_emiss_global"].GetDouble();
      emiss_share=inputs["climparams"][0]["emiss_share"].GetDouble();
      
      nshocks=inputs["climshockparams"][0]["nshocks"].GetInt();
      a_0.ReSize(nshocks);
      a_0(1)=inputs["climshockparams"][0]["a1_0"].GetDouble();
      a_0(2)=inputs["climshockparams"][0]["a2_0"].GetDouble();
      a_0(3)=inputs["climshockparams"][0]["a3_0"].GetDouble();
      a_0(4)=inputs["climshockparams"][0]["a4_0"].GetDouble();
      a_0(5)=inputs["climshockparams"][0]["a5_0"].GetDouble();
      a_0(6)=inputs["climshockparams"][0]["a6_0"].GetDouble();
      a_0(7)=inputs["climshockparams"][0]["a7_0"].GetDouble();
      b_0.ReSize(nshocks);
      b_0(1)=inputs["climshockparams"][0]["b1_0"].GetDouble();
      b_0(2)=inputs["climshockparams"][0]["b2_0"].GetDouble();
      b_0(3)=inputs["climshockparams"][0]["b3_0"].GetDouble();
      b_0(4)=inputs["climshockparams"][0]["b4_0"].GetDouble();
      b_0(5)=inputs["climshockparams"][0]["b5_0"].GetDouble();
      b_0(6)=inputs["climshockparams"][0]["b6_0"].GetDouble();
      b_0(7)=inputs["climshockparams"][0]["b7_0"].GetDouble();
      shockexponent1.ReSize(nshocks);
      shockexponent1(1)=inputs["climshockparams"][0]["shockexponent1_1"].GetDouble();
      shockexponent1(2)=inputs["climshockparams"][0]["shockexponent2_1"].GetDouble();
      shockexponent1(3)=inputs["climshockparams"][0]["shockexponent3_1"].GetDouble();
      shockexponent1(4)=inputs["climshockparams"][0]["shockexponent4_1"].GetDouble();
      shockexponent1(5)=inputs["climshockparams"][0]["shockexponent5_1"].GetDouble();
      shockexponent1(6)=inputs["climshockparams"][0]["shockexponent6_1"].GetDouble();
      shockexponent1(7)=inputs["climshockparams"][0]["shockexponent7_1"].GetDouble();
      shockexponent2.ReSize(nshocks);
      shockexponent2(1)=inputs["climshockparams"][0]["shockexponent1_2"].GetDouble();
      shockexponent2(2)=inputs["climshockparams"][0]["shockexponent2_2"].GetDouble();
      shockexponent2(3)=inputs["climshockparams"][0]["shockexponent3_2"].GetDouble();
      shockexponent2(4)=inputs["climshockparams"][0]["shockexponent4_2"].GetDouble();
      shockexponent2(5)=inputs["climshockparams"][0]["shockexponent5_2"].GetDouble();
      shockexponent2(6)=inputs["climshockparams"][0]["shockexponent6_2"].GetDouble();
      shockexponent2(7)=inputs["climshockparams"][0]["shockexponent7_2"].GetDouble();
      a2_nord=inputs["climshockparams"][0]["a2_nord"].GetDouble();
      sd_nord=inputs["climshockparams"][0]["sd_nord"].GetDouble();

      A0=inputs["inits"][0]["A0"].GetDouble();
      LS0=inputs["inits"][0]["LS0"].GetDouble();
      W10=inputs["inits"][0]["W10"].GetDouble();
      W20=inputs["inits"][0]["W20"].GetDouble();
      L0=inputs["inits"][0]["L0"].GetDouble();
      w0=inputs["inits"][0]["w0"].GetDouble();
      K0=inputs["inits"][0]["K0"].GetDouble();
      bankmarkup_init=inputs["inits"][0]["bankmarkup_init"].GetDouble();
      FirmDefaultProbability_init=inputs["inits"][0]["FirmDefaultProbability_init"].GetDouble();
      A0_en=inputs["inits"][0]["A0_en"].GetInt();
      A0_ef=inputs["inits"][0]["A0_ef"].GetInt();
      K_ge0_perc=inputs["inits"][0]["K_ge0_perc"].GetDouble();
      pf0=inputs["inits"][0]["pf0"].GetDouble();
      mi_en0=inputs["inits"][0]["mi_en0"].GetDouble();
      A_de0=inputs["inits"][0]["A_de0"].GetDouble();
      EM0=inputs["inits"][0]["EM0"].GetDouble();
      CF_ge0=inputs["inits"][0]["CF_ge0"].GetDouble();
      t_CO2_0=inputs["inits"][0]["t_CO2_0"].GetDouble();
      t_CO2_en_0=inputs["inits"][0]["t_CO2_en_0"].GetDouble();
      r_base=inputs["inits"][0]["r_base"].GetDouble();
      d_Am_init=inputs["inits"][0]["d_Am_init"].GetDouble();
      D_h0=inputs["inits"][0]["D_h0"].GetDouble();
      NW_b0=inputs["inits"][0]["NW_b0"].GetDouble();
      pm=inputs["inits"][0]["pm"].GetDouble();
      D_e0=inputs["inits"][0]["D_e0"].GetDouble();

      Emiss_yearly_0=inputs["climinits"][0]["Emiss_yearly_0"].GetDouble();
      Cum_emissions_0=inputs["climinits"][0]["Cum_emissions_0"].GetDouble();
      T_0_cumemiss=inputs["climinits"][0]["T_0_cumemiss"].GetDouble();
      Con00=inputs["climinits"][0]["Con00"].GetDouble();
      Conref=Con00*laydep(1);
      NPP0=inputs["climinits"][0]["NPP0"].GetDouble();
      Cat0=inputs["climinits"][0]["Cat0"].GetDouble();
      Catinit0=inputs["climinits"][0]["Catinit0"].GetDouble();
      Honinit0.ReSize(ndep);
      Coninit0.ReSize(ndep);
      Honinit0(1)=inputs["climinits"][0]["Honinit01"].GetDouble();
      Honinit0(2)=inputs["climinits"][0]["Honinit02"].GetDouble();
      Honinit0(3)=inputs["climinits"][0]["Honinit03"].GetDouble();
      Honinit0(4)=inputs["climinits"][0]["Honinit04"].GetDouble();
      Honinit0(5)=inputs["climinits"][0]["Honinit05"].GetDouble();
      Coninit0(1)=inputs["climinits"][0]["Coninit01"].GetDouble();
      Coninit0(2)=inputs["climinits"][0]["Coninit02"].GetDouble();
      Coninit0(3)=inputs["climinits"][0]["Coninit03"].GetDouble();
      Coninit0(4)=inputs["climinits"][0]["Coninit04"].GetDouble();
      Coninit0(5)=inputs["climinits"][0]["Coninit05"].GetDouble();
      Tmixedinit0=inputs["climinits"][0]["Tmixedinit0"].GetDouble();
      biominit0=inputs["climinits"][0]["biominit0"].GetDouble();
      huminit0=inputs["climinits"][0]["huminit0"].GetDouble();
      Catinit1=inputs["climinits"][0]["Catinit1"].GetDouble();
      Honinit1.ReSize(ndep);
      Coninit1.ReSize(ndep);
      Honinit1(1)=inputs["climinits"][0]["Honinit11"].GetDouble();
      Honinit1(2)=inputs["climinits"][0]["Honinit12"].GetDouble();
      Honinit1(3)=inputs["climinits"][0]["Honinit13"].GetDouble();
      Honinit1(4)=inputs["climinits"][0]["Honinit14"].GetDouble();
      Honinit1(5)=inputs["climinits"][0]["Honinit15"].GetDouble();
      Coninit1(1)=inputs["climinits"][0]["Coninit11"].GetDouble();
      Coninit1(2)=inputs["climinits"][0]["Coninit12"].GetDouble();
      Coninit1(3)=inputs["climinits"][0]["Coninit13"].GetDouble();
      Coninit1(4)=inputs["climinits"][0]["Coninit14"].GetDouble();
      Coninit1(5)=inputs["climinits"][0]["Coninit15"].GetDouble();
      Tmixedinit1=inputs["climinits"][0]["Tmixedinit1"].GetDouble();
      biominit1=inputs["climinits"][0]["biominit1"].GetDouble();
      huminit1=inputs["climinits"][0]["huminit1"].GetDouble();

      flag_cum_emissions=inputs["flags"][0]["flag_cum_emissions"].GetInt();
      flag_tax_CO2=inputs["flags"][0]["flag_tax_CO2"].GetInt();
      flag_capshocks=inputs["flags"][0]["flag_capshocks"].GetInt();
      flag_outputshocks=inputs["flags"][0]["flag_outputshocks"].GetInt();
      flag_inventshocks=inputs["flags"][0]["flag_inventshocks"].GetInt();
      flag_encapshocks=inputs["flags"][0]["flag_encapshocks"].GetInt();
      flag_popshocks=inputs["flags"][0]["flag_popshocks"].GetInt();
      flag_prodshocks1=inputs["flags"][0]["flag_prodshocks1"].GetInt();
      flag_prodshocks2=inputs["flags"][0]["flag_prodshocks2"].GetInt();
      flag_share_END=inputs["flags"][0]["flag_share_END"].GetInt();
      flag_energy_exp=inputs["flags"][0]["flag_energy_exp"].GetInt();
      flagbailout=inputs["flags"][0]["flagbailout"].GetInt();
      flag_entry=inputs["flags"][0]["flag_entry"].GetInt();
      flag_nonCO2_force=inputs["flags"][0]["flag_nonCO2_force"].GetInt();
      flag_validation=inputs["flags"][0]["flag_validation"].GetInt();
      flag_inventories=inputs["flags"][0]["flag_inventories"].GetInt();
      flag_rate_setting_loans=inputs["flags"][0]["flag_rate_setting_loans"].GetInt();
      flag_endogenous_exp_quota=inputs["flags"][0]["flag_endogenous_exp_quota"].GetInt();
      flag_WITCH_on=inputs["flags"][0]["flag_WITCH_on"].GetInt();
      flag_constant_WITCH_input=inputs["flags"][0]["flag_constant_WITCH_input"].GetInt();

      WITCH_scenario=inputs["WITCH_param"][0]["WITCH_scenario"].GetInt();
      WITCH_t0_year=inputs["WITCH_param"][0]["WITCH_t0_year"].GetInt();
      burnin=inputs["WITCH_param"][0]["burnin"].GetInt();
      rescale_ene_price=inputs["WITCH_param"][0]["rescale_ene_price"].GetDouble();
      rescale_ene_env_filth=inputs["WITCH_param"][0]["rescale_ene_env_filth"].GetDouble();
      ntechs=inputs["WITCH_param"][0]["ntechs"].GetInt();
      ene_tecs.resize(ntechs);
      ene_tecs[0]="biomass";
      ene_tecs[1]="coal";
      ene_tecs[2]="gas";
      ene_tecs[3]="hydro";
      ene_tecs[4]="nuclear";
      ene_tecs[5]="oil";
      ene_tecs[6]="solar";
      ene_tecs[7]="wind";

      e_plant_lifetime["biomass"]=inputs["WITCH_param"][0]["lifetime_biomass"].GetInt();
      e_plant_lifetime["coal"]=inputs["WITCH_param"][0]["lifetime_coal"].GetInt();
      e_plant_lifetime["gas"]=inputs["WITCH_param"][0]["lifetime_gas"].GetInt();
      e_plant_lifetime["hydro"]=inputs["WITCH_param"][0]["lifetime_hydro"].GetInt();
      e_plant_lifetime["nuclear"]=inputs["WITCH_param"][0]["lifetime_nuclear"].GetInt();
      e_plant_lifetime["oil"]=inputs["WITCH_param"][0]["lifetime_oil"].GetInt();
      e_plant_lifetime["solar"]=inputs["WITCH_param"][0]["lifetime_solar"].GetInt();
      e_plant_lifetime["wind"]=inputs["WITCH_param"][0]["lifetime_wind"].GetInt();
}

void RESIZE(void)
{
  //Resize all vectors, matrices and arrays to the needed dimensions
  X_a.ReSize(nshocks);
  X_b.ReSize(nshocks);
  A.ReSize(T,N1);
  C.ReSize(T,N1);
  C_secondhand.ReSize(T,N1);
  A_en.ReSize(T,N1);
  A_ef.ReSize(T,N1);
  A_de.ReSize(T);
  EM_de.ReSize(T);
  C_de.ReSize(T);
  G_de.ReSize(T);
  G_ge.ReSize(T);
  G_ge_n.ReSize(T);
  CF_ge.ReSize(T);
  IC_en_quota.ReSize(T);
  G_de_temp.ReSize(T);  
  shocks_encapstock_de.ReSize(T);
  shocks_encapstock_ge.ReSize(T);
  D2.ReSize(2,N2);
  De.ReSize(N2);
  f2.ReSize(3,N2);
  E2.ReSize(N2);
  c2.ReSize(N2);
  Q2.ReSize(N2);
  N.ReSize(2,N2);
  Inventories.ReSize(2,N2);
  Deposits_2.ReSize(2,N2);
  NW_2.ReSize(2,N2);
  NW_2_c.ReSize(N2);
  CapitalStock.ReSize(2,N2);
  deltaCapitalStock.ReSize(2,N2);
  DebtServiceToSales2.ReSize(N2);
  scrap_age.ReSize(N2);
  I.ReSize(N2);
  EI.ReSize(2,N2);
  EI_n.ReSize(N2);
  SI.ReSize(N2);
  SI_n.ReSize(N2);
  S2.ReSize(2,N2);
  Sales2.ReSize(N2);
  Loans_2.ReSize(2,N2);
  DebtService_2.ReSize(2,N2);
  CreditDemand.ReSize(N2);
  p2.ReSize(N2);
  DebtServiceToSales2_temp.ReSize(N2);
  DS2_order.resize(N2);
  scrap_p1.resize(N2);
  scrap_wage.resize(N2);
  scrap_energy.resize(N2);
  scrap_emission.resize(N2);
  scrap_supplier_ok.resize(N2);
  weight_labprod2.resize(N2);
  weight_eneff2.resize(N2);
  held_supplier.resize(N1);
  held_vintage.resize(N1);
  marked_count.assign(N2,0);
  mu2.ReSize(2,N2);
  LoanInterest_2.ReSize(N2);
  DebtRemittances2.ReSize(N2);
  baddebt_2.ReSize(N2);
  EId.ReSize(N2);
  SId.ReSize(N2);
  SIp.ReSize(N2);
  EIp.ReSize(N2);
  Ip.ReSize(N2);
  A2.ReSize(N2);
  A2_mprod.ReSize(N2);
  A2e.ReSize(N2);
  c2e.ReSize(N2);
  Ld2.ReSize(N2);
  Ld2_control.ReSize(N2);
  l2.ReSize(N2);
  n_mach.ReSize(N2); 
  Qd.ReSize(N2);
  K.ReSize(N2);
  K_cur.ReSize(N2);
  Kd.ReSize(N2);
  Ktrig.ReSize(N2);
  Pi2.ReSize(N2);
  Q2temp.ReSize(N2);
  Wages_2.ReSize(N2);
  Investment_2.ReSize(N2);
  EnergyPayments_2.ReSize(N2);
  InterestDeposits_2.ReSize(N2);
  Taxes_2.ReSize(N2);
  Taxes_CO2_2.ReSize(N2);
  f_temp2.ReSize(N2);
  D_temp2.ReSize(N2);
  dN.ReSize(N2);
  dNm.ReSize(N2);
  supl.ReSize(N2);
  Cmach.ReSize(N2);
  CmachEI.ReSize(N2);
  CmachSI.ReSize(N2);
  Ne.ReSize(N2);
  mol.ReSize(N2);
  Dividends_2.ReSize(N2);
  Match.ReSize(N2,N1);
  BankingSupplier_2.ReSize(N2);  
  r_deb_h.ReSize(N2);
  FirmDefaultProbability.ReSize(N2);
  k.ReSize(N2);
  D2_en.ReSize(N2);
  A2e_en.ReSize(N2);
  A2e_ef.ReSize(N2);
  A2e2.ReSize(N2);
  A2e_en2.ReSize(N2);
  A2e_ef2.ReSize(N2);
  A2_en.ReSize(N2);
  A2_ef.ReSize(N2);
  Emiss2.ReSize(N2);
  Injection_2.ReSize(N2);
  exiting_2.ReSize(N2);
  exit_payments2.ReSize(N2);
  exit_equity2.ReSize(N2);
  exit_marketshare2.ReSize(N2);
  n_mach_entry.ReSize(N2);
  shocks_capstock.ReSize(N2);
  shocks_invent.ReSize(N2);
  Loss_Capital.ReSize(N2);
  Loss_Inventories.ReSize(N2);
  k_entry.ReSize(N2);
  EntryShare.ReSize(N2);
  CompEntry.ReSize(N2);
  K_temp.ReSize(N2);
  K_loss.ReSize(N2);
  C_loss.ReSize(N2);
  I_loss.ReSize(N1);
  marker_age.ReSize(N2);
  Deposits_1.ReSize(2,N1);
  NW_1.ReSize(2,N1);
  NW_1_c.ReSize(N1);
  Wages_1.ReSize(N1);
  S1.ReSize(N1);
  Sales1.ReSize(N1);
  p1.ReSize(N1);
  RD.ReSize(2,N1);
  f1.ReSize(2,N1);
  Q1.ReSize(N1);
  D1.ReSize(N1);
  A1.ReSize(N1);
  A1inn.ReSize(N1);
  A1pinn.ReSize(N1);
  A1imm.ReSize(N1);
  A1pimm.ReSize(N1);
  Pi1.ReSize(N1);
  Ld1.ReSize(N1);
  Ld1rd.ReSize(N1);
  ee1.ReSize(N1);
  nclient.ReSize(N1);
  Dividends_1.ReSize(N1);
  EnergyPayments_1.ReSize(N1);
  InterestDeposits_1.ReSize(N1);
  Taxes_1.ReSize(N1);
  Taxes_CO2_1.ReSize(N1);
  BankingSupplier_1.ReSize(N1);
  c1.ReSize(N1);
  A1p.ReSize(N1);
  RDin.ReSize(N1);
  RDim.ReSize(N1);
  Inn.ReSize(N1);
  Imm.ReSize(N1);
  Td.ReSize(N1+1);
  D1_en.ReSize(N1);
  A1_en.ReSize(N1);
  A1_ef.ReSize(N1);
  A1p_en.ReSize(N1);
  A1p_ef.ReSize(N1);
  EE_inn.ReSize(N1);
  EEp_inn.ReSize(N1);
  EE_imm.ReSize(N1);
  EEp_imm.ReSize(N1);
  EF_inn.ReSize(N1);
  EFp_inn.ReSize(N1);
  EF_imm.ReSize(N1);
  EFp_imm.ReSize(N1);   
  Emiss1.ReSize(N1);
  Injection_1.ReSize(N1);
  Balances_1.ReSize(N1);
  baddebt_1.ReSize(N1);
  exiting_1.ReSize(N1);
  exiting_1_payments.ReSize(N1);
  shocks_machprod.ReSize(N1);
  shocks_techprod.ReSize(N1);
  shocks_labprod1.ReSize(N1);
  shocks_labprod2.ReSize(N2);
  shocks_eneff1.ReSize(N1);
  shocks_eneff2.ReSize(N2);
  shocks_output1.ReSize(N1);
  shocks_output2.ReSize(N2);
  S1_temp.ReSize(2,N1);
  S2_temp.ReSize(2,N2);

  NbClient_1.ReSize(NB);
  NbClient_2.ReSize(NB);
  bonds_dem.ReSize(NB);
  DebtServiceToSales2_bank.ReSize(N2,NB);
  DS2_rating.ReSize(N2,NB);
  fB.ReSize(2,NB); 
  NW_b.ReSize(2,NB); 
  NW_b_c.ReSize(NB);
  Loans_b.ReSize(2,NB);
  BankMatch_1.ReSize(N1,NB); 
  BankMatch_2.ReSize(N2,NB); 
  BankCredit.ReSize(NB);
  Taxes_b.ReSize(NB);
  BaselBankCredit.ReSize(NB);
  InterestDeposits.ReSize(NB);
  Deposits.ReSize(2,NB);
  BankProfits.ReSize(NB);
  Dividends_b.ReSize(NB); 
  Bailout_b.ReSize(NB);
  LossAbsorbed.ReSize(NB);
  BankEquity_temp.ReSize(NB);
  Bank_active.ReSize(NB);
  GB_b.ReSize(2,NB);   
  bonds_purchased.ReSize(NB);
  BondRepayments_b.ReSize(NB);
  NL_1.ReSize(NB);
  NL_2.ReSize(NB);
  BankProfits_temp.ReSize(NB);
  r_deb.ReSize(NB);
  bankmarkup.ReSize(NB);
  riskWeightedAssets.ReSize(NB);	
  capitalAdequacyRatio.ReSize(NB);
  Bond_share.ReSize(NB);
  Deposits_hb.ReSize(2,NB);
  Deposits_eb.ReSize(2,NB);
  Advances_b.ReSize(2,NB);
  Reserves_b.ReSize(2,NB);
  InterestBonds_b.ReSize(NB);
  LoanInterest.ReSize(NB);
  InterestReserves_b.ReSize(NB);
  InterestAdvances_b.ReSize(NB);
  Outflows.ReSize(NB);
  Inflows.ReSize(NB);
  DepositShare_e.ReSize(NB);
  DepositShare_h.ReSize(NB);
  baddebt_b.ReSize(NB);
  capitalRecovered.ReSize(NB);
  capitalRecovered2.ReSize(NB);
  capitalRecoveredShare.ReSize(NB);
  LossEntry_b.ReSize(NB);
  ReserveBalance.ReSize(NB);
  ShareBonds.ReSize(NB);
  ShareReserves.ReSize(NB);
  ShareAdvances.ReSize(NB);
  Adjustment.ReSize(NB);
  prior.ReSize(NB);

  fluxC.ReSize((ndep-1));
  Con.ReSize(2,ndep);
  Hon.ReSize(2,ndep);
  Ton.ReSize(2,ndep);
  Cax.ReSize(niterclim);
  Caxx.ReSize(niterclim);
  Cay.ReSize(niterclim);
  Cayy.ReSize(niterclim);
  Caa.ReSize(niterclim);
  fluxH.ReSize((ndep-1));
  Emiss_TOT.ReSize((freqclim*2));

  //Initialising all resized variables to zero to avoid memory issues
  X_a=0;
  X_b=0;
  A=0;
  C=0;
  C_secondhand=0;
  A_en=0;
  A_ef=0;
  A_de=0;
  EM_de=0;
  C_de=0;
  G_de=0;
  G_ge=0;
  G_ge_n=0;
  CF_ge=0;
  IC_en_quota=0;
  G_de_temp=0;  
  shocks_encapstock_de=0;
  shocks_encapstock_ge=0;
  D2=0;
  De=0;
  f2=0;
  E2=0;
  c2=0;
  Q2=0;
  N=0;
  Inventories=0;
  Deposits_2=0;
  NW_2=0;
  NW_2_c=0;
  CapitalStock=0;
  deltaCapitalStock=0;
  DebtServiceToSales2=0;
  scrap_age=0;
  I=0;
  EI=0;
  EI_n=0;
  SI=0;
  SI_n=0;
  S2=0;
  Sales2=0;
  Loans_2=0;
  DebtService_2=0;
  CreditDemand=0;
  p2=0;
  DebtServiceToSales2_temp=0;
  mu2=0;
  LoanInterest_2=0;
  DebtRemittances2=0;
  baddebt_2=0;
  EId=0;
  SId=0;
  SIp=0;
  EIp=0;
  Ip=0;
  A2=0;
  A2_mprod=0;
  A2e=0;
  c2e=0;
  Ld2=0;
  Ld2_control=0;
  l2=0;
  n_mach=0; 
  Qd=0;
  K=0;
  K_cur=0;
  Kd=0;
  Ktrig=0;
  Pi2=0;
  Q2temp=0;
  Wages_2=0;
  Investment_2=0;
  EnergyPayments_2=0;
  InterestDeposits_2=0;
  Taxes_2=0;
  Taxes_CO2_2=0;
  f_temp2=0;
  D_temp2=0;
  dN=0;
  dNm=0;
  supl=0;
  Cmach=0;
  CmachEI=0;
  CmachSI=0;
  Ne=0;
  mol=0;
  Dividends_2=0;
  Match=0;
  BankingSupplier_2=0;  
  r_deb_h=0;
  FirmDefaultProbability=0;
  k=0;
  D2_en=0;
  A2e_en=0;
  A2e_ef=0;
  A2e2=0;
  A2e_en2=0;
  A2e_ef2=0;
  A2_en=0;
  A2_ef=0;
  Emiss2=0;
  Injection_2=0;
  exiting_2=0;
  exit_payments2=0;
  exit_equity2=0;
  exit_marketshare2=0;
  n_mach_entry=0;
  shocks_capstock=0;
  shocks_invent=0;
  Loss_Capital=0;
  Loss_Inventories=0;
  k_entry=0;
  EntryShare=0;
  CompEntry=0;
  K_temp=0;
  K_loss=0;
  C_loss=0;
  I_loss=0;
  marker_age=0;
  Deposits_1=0;
  NW_1=0;
  NW_1_c=0;
  Wages_1=0;
  S1=0;
  Sales1=0;
  p1=0;
  RD=0;
  f1=0;
  Q1=0;
  D1=0;
  A1=0;
  A1inn=0;
  A1pinn=0;
  A1imm=0;
  A1pimm=0;
  Pi1=0;
  Ld1=0;
  Ld1rd=0;
  ee1=0;
  nclient=0;
  Dividends_1=0;
  EnergyPayments_1=0;
  InterestDeposits_1=0;
  Taxes_1=0;
  Taxes_CO2_1=0;
  BankingSupplier_1=0;
  c1=0;
  A1p=0;
  RDin=0;
  RDim=0;
  Inn=0;
  Imm=0;
  Td=0;
  D1_en=0;
  A1_en=0;
  A1_ef=0;
  A1p_en=0;
  A1p_ef=0;
  EE_inn=0;
  EEp_inn=0;
  EE_imm=0;
  EEp_imm=0;
  EF_inn=0;
  EFp_inn=0;
  EF_imm=0;
  EFp_imm=0;   
  Emiss1=0;
  Injection_1=0;
  Balances_1=0;
  baddebt_1=0;
  exiting_1=0;
  exiting_1_payments=0;
  shocks_machprod=0;
  shocks_techprod=0;
  shocks_labprod1=0;
  shocks_labprod2=0;
  shocks_eneff1=0;
  shocks_eneff2=0;
  shocks_output1=0;
  shocks_output2=0;
  S1_temp=0;
  S2_temp=0;
  NbClient_1=0;
  NbClient_2=0;
  bonds_dem=0;
  DebtServiceToSales2_bank=0;
  DS2_rating=0;
  fB=0; 
  NW_b=0; 
  NW_b_c=0;
  Loans_b=0;
  BankMatch_1=0; 
  BankMatch_2=0; 
  BankCredit=0;
  Taxes_b=0;
  BaselBankCredit=0;
  InterestDeposits=0;
  Deposits=0;
  BankProfits=0;
  Dividends_b=0; 
  Bailout_b=0;
  LossAbsorbed=0;
  BankEquity_temp=0;
  Bank_active=0;
  GB_b=0;   
  bonds_purchased=0;
  BondRepayments_b=0;
  NL_1=0;
  NL_2=0;
  BankProfits_temp=0;
  r_deb=0;
  bankmarkup=0;
  riskWeightedAssets=0;	
  capitalAdequacyRatio=0;
  Bond_share=0;
  Deposits_hb=0;
  Deposits_eb=0;
  Advances_b=0;
  Reserves_b=0;
  InterestBonds_b=0;
  LoanInterest=0;
  InterestReserves_b=0;
  InterestAdvances_b=0;
  Outflows=0;
  Inflows=0;
  DepositShare_e=0;
  DepositShare_h=0;
  baddebt_b=0;
  capitalRecovered=0;
  capitalRecovered2=0;
  capitalRecoveredShare=0;
  LossEntry_b=0;
  ReserveBalance=0;
  ShareBonds=0;
  ShareReserves=0;
  ShareAdvances=0;
  Adjustment=0;
  prior=0;
  fluxC=0;
  Con=0;
  Hon=0;
  Ton=0;
  Cax=0;
  Caxx=0;
  Cay=0;
  Cayy=0;
  Caa=0;
  fluxH=0;
  Emiss_TOT=0;

  //One block each, zeroed, instead of T times N1 separate row allocations.
  age.resize3(T,N1,N2);
  C_pb.resize3(T,N1,N2);
  g_pb.resize3(T,N1,N2);
  g_c.resize3(T,N1,N2);
  g_c2.resize3(T,N1,N2);
  g_c3.resize3(T,N1,N2);
  g.resize3(T,N1,N2);
  gtemp.resize3(T,N1,N2);
  g_price.resize3(T,N1,N2);
  g_secondhand.resize(T);
  g_secondhand_p.resize(T);
  age_secondhand.resize(T);
  for (std::size_t i = 0; i < T; i++)
  {
    g_secondhand[i].resize(N1);
    g_secondhand_p[i].resize(N1);
    age_secondhand[i].resize(N1);
    for (std::size_t j = 0; j < N1; j++)
    {
      g_secondhand[i][j]=0;
      g_secondhand_p[i][j]=0;
      age_secondhand[i][j]=0;
    }
  }

  //Initialising all other global variables to zero to avoid memory issues
  i=0;	                                         
  ii=0;	                                     
  iii=0;					                     
  j=0;                                                                        
  t=0;                                          
  tt=0;                                         
  rni=0;                                        
  t0=0;                                         
  t00=0;                                        
  n=0;                                          
  iterations=0;                                 
  pareto_rv=0;                                  
  tolerance=0;                                  
  deviation=0;                                  
  parber=0;                                     
  rnd=0;                                        
  N1r=0;						                 
  N2r=0;                                        
  step=0;                                       
  stepb=0;                                    
  counter=0;                                       
  age0=0;                                                                    
  D20=0;                                        
  DS2_min_index=0;                              
  newbroch=0;                                   
  indsupl=0;                                    
  flag=0;                                       
  payback=0;                                                                         
  nmachprod=0;                                  
  nmp_temp=0;                                   
  cmin=0;                                       
  imin=0;                                       
  jmin=0;                                       
  tmin=0;                                       
  MaxFunds=0;                                                                  
  p1test=0;                                    
  rated_firm_2=0;                               
  Qpast=0;                                                                        
  scrapmax=0;                                   
  cmax=0;                                       
  ind_i=0;                                      
  ind_tt=0;                                     
  scrap_n=0;                                    
  sendingBank=0;                                
  receivingBank=0;                              
  c_de_min=0;                                   
  cf_min_ge=0;                                  
  Q_de_temp=0;                                  
  idmin=0;                                      
  parber_en_de=0;                               
  parber_en_ge=0;                               
  l2m=0;                                        
  p2m=0;                                        
  Cres=0;                                       
  Cresb=0;                                    
  cpi_temp=0;                                   
  maxbank=0;                                    
  max_equity=0;                                 
  multip_bailout=0;                                                        
  ns1=0;                                        
  ns2=0;                                        
  mD1=0;                                        
  mD2=0;                                        
  multip_entry=0;                               
  injection=0;                                  
  injection2=0;                                 
  n_mach_exit=0;                                
  n_mach_exit2=0;                               
  n_mach_needed=0;                              
  n_mach_resid=0;                               
  n_mach_resid2=0;                                                     
  cpi_init=0;                                   
  GDP_init=0;                                   
  baddebt_2_temp=0;                             
  markdownCapital=0;                            
  post=0;                                       
  prior_cb=0;                                   
  post_cb=0;                                    
  DepositsCheck_1=0;                            
  DepositsCheck_2=0;                            
  p2_entry=0;                                   
  f2_exit=0;                                    
  CurrentDemand=0;                              
  CompEntry_m=0;                                
  K_gap=0;                                      
  K_top=0;                                      
  loss=0;                                       
  lossj=0;                                      
  rani=0;                                       
  rant=0;                                       
  ranj=0;                                       
  reduction=0;                                  
  K_temp_sum=0;                                                                 
  ptemp=0;                                      
  Ldtemp=0;                                     
  Deposits_h=0;                              
  Deposits_e=0;                              
  GB_cb=0;                                   
  GB=0;                                      
  Deposits_fuel=0;                           
  Deposits_fuel_cb=0;                        
  Advances=0;                                
  Reserves=0;                                
  CapitalStock_e=0;                          
  NW_h=0;                                    
  NW_gov=0;                                  
  NW_cb=0;                                   
  NW_e=0;                                    
  NW_f=0;                                    
  NW_h_c=0;                                     
  NW_gov_c=0;                                   
  NW_cb_c=0;                                    
  NW_e_c=0;                                     
  NW_f_c=0;                                     
  NWSum=0;                                      
  RealAssets=0;                                 
  Wages_en=0;                                   
  Wages=0;                                      
  Dividends=0;                               
  Benefits=0;                                   
  EnergyPayments=0;                             
  Dividends_e=0;                                
  Taxes_h=0;                                    
  Taxes_CO2_e=0;                                
  Taxes_CO2=0;                                  
  InterestDeposits_h=0;                         
  InterestDeposits_e=0;                         
  InterestBonds=0;                              
  InterestBonds_cb=0;                           
  BondRepayments_cb=0;                          
  Taxes=0;                                      
  Consumption=0;                                
  FirmTransfers=0;                              
  FirmTransfers_1=0;                            
  FirmTransfers_2=0;                            
  InterestReserves=0;                           
  InterestAdvances=0;                           
  TransferCB=0;                                 
  FuelCost=0;                                   
  TransferFuel=0;                                                                                                                             
  Balance_h=0;                                  
  Balance_1=0;                                  
  Balance_2=0;                                  
  Balance_e=0;                                  
  Balance_b=0;                                  
  Balance_g=0;                                  
  Balance_cb=0;                                 
  Balance_f=0;                                  
  BalanceSum=0;                              
  w=0;				                         
  LS=0;                                         
  U=0;                                       
  Divtot_1=0;                                   
  Divtot_2=0;                                   
  Cons=0;                                       
  Deposits_recovered_1=0;                       
  Deposits_recovered_2=0;                                                       
  r_depo=0;                                     
  bonds_dem_tot=0;                              
  r_bonds=0;                                    
  Bailout=0;                                    
  G=0;                                          
  Deficit=0;                                    
  PSBR=0;                                       
  NewBonds=0;                                   
  EntryCosts=0;                                 
  BankTransfer=0;                               
  r_cbreserves=0;                               
  r_a=0;                                        
  r=0;                                          
  ProfitCB=0;                                
  Adjustment_cb=0;                              
  d_cpi_target_a=0;                             
  inflation_a=0;                                
  pf=0;                                         
  mi_en=0;                                      
  c_en=0;                                    
  D1_en_TOT=0;                                  
  D2_en_TOT=0;                                  
  D_en_TOT=0;                                
  K_ge=0;                                       
  K_de=0;                                       
  K_gelag=0;                                    
  K_delag=0;                                    
  Q_ge=0;                                       
  Q_de=0;                                       
  EI_en=0;                                      
  EI_en_de=0;                                   
  EI_en_ge=0;                                   
  IC_en=0;                                      
  LDexp_en=0;                                   
  PC_en=0;                                      
  c_infra=0;                                    
  share_de=0;                                   
  Rev_en=0;                                     
  RD_en_de=0;                                   
  RD_en_ge=0;                                   
  LDrd_de=0;                                    
  LDrd_ge=0;                                    
  Inn_en_ge=0;                                  
  Inn_en_de=0;                                  
  A_de_inn=0;                                   
  EM_de_inn=0;                                  
  CF_ge_inn=0;                                  
  ProfitEnergy=0;                               
  G_de_0=0;                                     
  G_ge_0=0;                                     
  G_ge_n_0=0;                                   
  Tmixed=0;                                  
  Emiss_yearly_calib=0;                      
  g_rate_em_y=0;                                
  Emiss_yearly=0;                            
  Emiss1_TOT=0;                                 
  Emiss2_TOT=0;                                 
  Emiss_en=0;                                   
  NPP=0;                                                               
  Cum_emissions=0;	                             
  shock_pop=0;                                                                 
  t_CO2=0;                                      
  t_CO2_en=0;                                   
  Emiss_gauge=0;                                
  Cat=0;                                     
  humrelease=0;                                 
  hum=0;                                     
  biorelease=0;                                 
  biom=0;                                    
  Cat1=0;                                       
  dCat1=0;                                      
  Con1=0;                                       
  Ctot1=0;                                      
  FCO2=0;                                       
  Fin=0;                                        
  Fout=0;                                       
  Emiss_global=0;                               
  Am=0;                 	                 
  Am_a=0;                     	                 
  Am2=0;                                        
  Am1=0;                                        
  ftot=0;                                    
  Em2=0;                                     
  cpi=0;                                     
  kpi=0;                                        
  Am_en=0;                                   
  LD1rdtot=0;                                   
  LDentot=0;                                    
  Tdtot=0;                                      
  LD1tot=0;					                 
  LD2tot=0;                                     
  LSe=0;                                        
  LD=0;                                         
  LD2=0;                                        
  Pitot1=0;				                     
  Pitot2=0;                                     
  dNtot=0;                                      
  dNmtot=0;                                     
  ExpansionInvestment_r=0;                      
  ExpansionInvestment_n=0;                      
  ReplacementInvestment_r=0;                    
  ReplacementInvestment_n=0;                    
  Investment_r=0;                               
  Investment_n=0;                               
  Consumption_r=0;                                                    
  Q2tot=0;                                      
  Q1tot=0;                                                                                                   
  A_mi=0;                                       
  A1_mi=0;                                      
  A2_en_mi=0;                                   
  A2_ef_mi=0;                                   
  A1_en_mi=0;                                   
  A1_ef_mi=0;                                   
  A_sd=0;                                       
  H1=0;                                         
  H2=0;                                         
  HB=0;                                         
  GDP_r=0;                                   
  GDP_n=0;                                   
  d_U=0;                                        
  d_cpi=0;                                      
  d_Am=0;                                       
  dw=0;                                         
  dw2=0;                                        
  A2scr=0;                                      
  A1scr=0;                                      
  Utilisation=0;                                
  counter_bankfailure=0;      
  Loan_interest_e=0;
  DebtService_e=0;
  DeafaultedDebtRecovered_e=0;
  BadDebt_e=0;
  DebtRemittances_e=0;
  DebtWrittenOff_e=0;
  DefaultedDeposits_e=0;
  Loans_preDefault_e=0;                  
}

void INITIALIZE(int Exseed)
{
  //Set seed
  seed=Exseed;			
  //Pointer to seed
  p_seed=&seed;		
  //Tolerance level used to check deviations from stock-flow consistency		
  tolerance=1e-06;
  //Numbers of agents as doubles
  N1r=double(N1);
	N2r=double(N2);
  //Generate empty output files
  INTFILE();	

  //Initialize energy sector
  if(flag_WITCH_on==1)
  {
    Deposits_e=0;
  }
  else
  {
    Deposits_e=D_e0;
  }

  Loans_e=0;
  
  for(i=0; i<ene_tecs.size(); i++){
    if(flag_WITCH_on==1)
    {
      deposits_e_mult_tech[ene_tecs[i]]=10;
    }
    else
    {
      deposits_e_mult_tech[ene_tecs[i]]=0;
    }
      Deposits_e(1)+=deposits_e_mult_tech[ene_tecs[i]];
      Deposits_e(2)+=deposits_e_mult_tech[ene_tecs[i]];
      interestsDeposits_e_mult_tech[ene_tecs[i]]=0.0;
      d_en_mult_tech[ene_tecs[i]]=0.0;
      LDexp_en_mult_tech[ene_tecs[i]]=0.0;
      LDprod_en_mult_tech[ene_tecs[i]]=0.0;
      LDmaint_en_mult_tech[ene_tecs[i]]=0.0;
      Rev_en_mult_tech[ene_tecs[i]]=0.0;
      capitalStock_e_mult_tech[ene_tecs[i]]=0.0;
      profitEnergy_e_mult_tech[ene_tecs[i]]=0.0;
      delta_capitalStock_e_mult_tech[ene_tecs[i]]=0.0;
      emiss_e_mult_tech[ene_tecs[i]]=0.0;
      taxes_CO2_e_mult_tech[ene_tecs[i]]=0.0;
      dividends_e_mult_tech[ene_tecs[i]]=0.0;
      effective_dividendRate_e_mult_tech[ene_tecs[i]]=0.0;
      wage_e_mult_tech[ene_tecs[i]]=0.0;
      loans_e_mult_tech[ene_tecs[i]]=0.0;
      Deposits_e(1)+=loans_e_mult_tech[ene_tecs[i]];
      Deposits_e(2)+=loans_e_mult_tech[ene_tecs[i]];
      creditDemand_e_mult_tech[ene_tecs[i]]=0.0;
      loan_interest_e_mult_tech[ene_tecs[i]]=0.0;
      debtRemittances_e_mult_tech[ene_tecs[i]]=0.0;
      badDebt_e_mult_tech[ene_tecs[i]]=0.0;
      shareDefaultedDebt_e_mult_tech[ene_tecs[i]]=0.0;
      exiting_e_mult_tech[ene_tecs[i]]=0;
      unpaidCtax_e_mult_tech[ene_tecs[i]]=0.0;
      nom_inv_mult_tech[ene_tecs[i]]=0.0;
      prodCost_mult_tech[ene_tecs[i]]=0.0;
      maintCost_mult_tech[ene_tecs[i]]=0.0;
      capacity_mult_tech[ene_tecs[i]]=0.0;
      realInv_e_multi_tech[ene_tecs[i]]=0.0;
      unit_cost_inv_e_multi_tech[ene_tecs[i]]=0.0;
      uC_prod_mult_tech[ene_tecs[i]]=0.0;
      uC_cap_mult_tech[ene_tecs[i]]=0.0;
      depreciatedCapital_e_multi_tech[ene_tecs[i]].ReSize(T);
      depreciatedCapital_e_multi_tech[ene_tecs[i]]=0.0;
      inv_redistr_profit_e_multi_tech[ene_tecs[i]]=0.0;
  }

  //Determines number of firm customers of each bank
  ALLOCATEBANKCUSTOMERS();
  
  //Set all variables describing firm-bank networks to zero
  BankingSupplier_2=0;
  BankMatch_2=0;
  NbClient_2=0;
  BankingSupplier_1=0;
  BankMatch_1=0;
  NbClient_1=0;
  //Initial deposits and for firms
  Deposits_1=W10;
  Deposits_2=W20;
  Loans_2=L0;
  //Set deposits and loans from point of view of banks to zero
  Deposits=0;
  Loans_b=0;
  //First, iterate over all C-Firms
  for(j=1; j<=N2; j++)  
  {
    //Match C-firms to banks and give those banks the initial deposits and loans
    //While firm j does not have a supplier of banking services
    while (BankingSupplier_2(j)==0)    
    {
      //Randomly draw a bank
      rni=int(ran1(p_seed)*N1*N2)%NB+1;  
      //If that bank has not yet reached the number of customers determined above
      if (NbClient_2(rni)< NL_2(rni))     
      {
        //j becomes a customer of that bank
        BankMatch_2(j,rni)=1;      
        BankingSupplier_2(j)= rni;     
        NbClient_2(rni)++;             
        Deposits(1,rni)+=Deposits_2(1,j);
        Deposits(2,rni)+=Deposits_2(2,j);
        Loans_b(1,rni)+=Loans_2(1,j);
        Loans_b(2,rni)+=Loans_2(2,j);
      }
      else 
      {
        BankMatch_2(j,rni)=0;
      }
    }                                
  }  
  
  //Iterate over all K-Firms
  for(i=1; i<=N1; i++)  
  {
    //Match K-firms to banks and give those banks the initial deposits
    //While firm i does not have a supplier of banking services
    while (BankingSupplier_1(i)==0)  
    {
      //Randomly draw a bank
      rni=int(ran1(p_seed)*N1*N2)%NB+1;
      //If that bank has not yet reached the number of customers determined above
      if (NbClient_1(rni)< NL_1(rni))     
      {
        BankMatch_1(i,rni)=1;  
        BankingSupplier_1(i)= rni;      
        NbClient_1(rni)++;             
        Deposits(1,rni)+=Deposits_1(1,i);
        Deposits(2,rni)+=Deposits_1(2,i);
      }
      else {
      BankMatch_1(i,rni)=0;
      }
    }                                
  } 

  //Energy price
  if(flag_WITCH_on==1)
  {
    UPDATE_WITCH_YEAR(0);
    WITCH_price_adj=1;
    c_en(2)=FETCH_ELEPRICE_WITCH(WITCH_scenario);
  }
  else
  {
    c_en(2)=mi_en0;
  }
  //Set initial prices
  p2=(1+mi2)*(w0/A0+c_en(2)/A0_en);
  p1=(1+mi1)*(w0/(A0*pm)+c_en(2)/A0_en);
  //Initialise measure of average labour productivity
  Am(2)=(A0*N2+A0*pm*N1)/(N1+N2);	
  //Set initial household deposits
  Deposits_h=D_h0;
  //Finish initialising bank balance sheets
  for(i=1; i<=NB; i++)
  {
    //Bank's market share is equal to its share of firm customers
    fB(1,i)=(NbClient_2(i)+NbClient_1(i))/(N2+N1);
    fB(2,i)=(NbClient_2(i)+NbClient_1(i))/(N2+N1);
    //Bank's initial stock of government bonds is a % of its loan portfolio
    GB_b(1,i)=varphi*Loans_b(1,i);
    GB_b(2,i)=varphi*Loans_b(2,i);
    //Give the bank a share of household and energy sector deposits
    Deposits(1,i)+=fB(1,i)*(Deposits_h(1)+Deposits_e(1));
    Deposits(2,i)+=fB(2,i)*(Deposits_h(2)+Deposits_e(2));
    Deposits_hb(1,i)+=fB(1,i)*(Deposits_h(1));
    Deposits_hb(2,i)+=fB(2,i)*(Deposits_h(2));
    Deposits_eb(1,i)+=fB(1,i)*(Deposits_e(1));
    Deposits_eb(2,i)+=fB(2,i)*(Deposits_e(2));
    //Initialise share of bank in aggregate household and energy sector deposits
    DepositShare_e(i)=Deposits_eb(1,i)/(Deposits_e(1));
    DepositShare_h(i)=Deposits_hb(1,i)/(Deposits_h(1));
    //Bank reserves are given as a residual when setting an exogenous initial value for bank net worth (NW_b0)
    //and assuming that initial central bank advances are zero
    Reserves_b(1,i)=Deposits(1,i)+NW_b0*fB(1,i)-GB_b(1,i)-Loans_b(1,i);
    Reserves_b(2,i)=Deposits(2,i)+NW_b0*fB(2,i)-GB_b(2,i)-Loans_b(2,i);
    NW_b(1,i)=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits(1,i)-Advances_b(1,i);
    NW_b(2,i)=Loans_b(2,i)+Reserves_b(2,i)+GB_b(2,i)-Deposits(2,i)-Advances_b(2,i);
    riskWeightedAssets(i)=riskWeightLoans*Loans_b(1,i)+riskWeightGovBonds*GB_b(1,i);
    capitalAdequacyRatio(i)=NW_b(1,i)/riskWeightedAssets(i);
  }
  //Governmen bonds held by central bank are given by identity from the above
  GB_cb(1)=NW_b0+W10*N1+W20*N2+D_h0+Deposits_e(1)-GB_b.Row(1).Sum()-Loans_b.Row(1).Sum();
  GB_cb(2)=NW_b0+W10*N1+W20*N2+D_h0+Deposits_e(2)-GB_b.Row(2).Sum()-Loans_b.Row(2).Sum();
  //Initial stock of overall government bonds (banks+CB)
  GB(1)=GB_cb(1)+GB_b.Row(1).Sum();
  GB(2)=GB_cb(2)+GB_b.Row(2).Sum();
  //CB advances assumed to be zero initially
  Advances_b=0;
  Advances=0;
  //CB Stock of gov. bonds implies an equal aggregate stock of reserves
  Reserves=GB_cb;
  //Set nominal value of initial capital stock
  CapitalStock=K0/dim_mach*(1+mi1)*(w0/(A0*pm)+c_en(2)/A0_en);
  deltaCapitalStock=0;
  //Initialise net worth of households, K-Firms, government and central bank
  NW_h=Deposits_h;
  NW_1=Deposits_1;
  NW_gov=-GB;
  NW_cb=GB_cb-Reserves;
  //initial central bank rate is set by converting annual r_base to quarterly rate
  r=pow((1+r_base),0.25)-1;
  //convert quarterly inflation target to annual one
  d_cpi_target_a=pow((1+d_cpi_target),4)-1;

  //Initialise deposit, CB reserve and gov. bond interest rates
  r_depo=0;       // For now, interest rate set to zero, change that to "r_depo=r-bankmarkdown;" 
                  // when banks compete on the deposit market
  r_cbreserves=0; // For now, interest rate set to zero, change that to "r_cbreserves=r-centralbankmarkdown;" 
                  // later
  r_bonds=r;      // For now, interest rate set to r, change that to "r_bonds=r-bondsmarkdown;" 
                  // when the bond market is fixed

  //Initial central bank profits
  ProfitCB(2)=r_cbreserves*Reserves(2);
  //Initial wage rate
	w=w0;		
  //Initial change in mean productivity
  d_Am=d_Am_init;
  //Initial wage gross growth rate
  dw2=1;

  //Set initial value for carbon tax
  if(flag_tax_CO2>0)
  {
    t_CO2_en=t_CO2_en_0;
  }

  //Initial productivities, energy efficiencies, environmental friendliness						
	A1=A0;							
	A1p=A0*pm;
	A2=A0;  
  A2_mprod=A0;    
  A=A0;

  A1_en=A0_en;
  A1p_en=A0_en*pm;
  A2_en=A0_en;
  A_en=A0_en;
  
  A1_ef=A0_ef;
  A1p_ef=A0_ef;
  A2_ef=A0_ef;
  A_ef=A0_ef;
  
  //Set initial values for energy sector
  //Brown energy tech.
  A_de=A_de0;                     
  EM_de=EM0;
  //Fossil fuel price
  pf=pf0;
  //Mark-up
  mi_en=mi_en0;
  //Green energy tech
  CF_ge(1)=CF_ge0;

  //Initialise climate module
  if (flag_nonCO2_force==0)
  {
    //Case of no non-CO2 forcing
    //Initial atmospheric carbon
    Cat(1)=Catinit0;     
    //Initial biosphere carbon
    biom(1)=biominit0;
    //Initial humus carbon
    hum(1)=huminit0; 
    //Initial temperature anomaly
    Tmixed(1)=Tmixedinit0; 
    //Initial carbon and heat content of ocean layers
    for (j=1;  j<=ndep; j++)
    {
      Con(1,j)=Coninit0(j);
      Hon(1,j)=Honinit0(j);
    }
  }
  else
  {
    //Case of non-CO2 forcing
    //Initial atmospheric carbon
    Cat(1)=Catinit1;     
    //Initial biosphere carbon
    biom(1)=biominit1;
    //Initial humus carbon
    hum(1)=huminit1; 
    //Initial temperature anomaly
    Tmixed(1)=Tmixedinit1;
    //Initial carbon and heat content of ocean layers
    for (j=1;  j<=ndep; j++)
    {
      Con(1,j)=Coninit1(j);
      Hon(1,j)=Honinit1(j);
    }
  }

  if(flag_cum_emissions==1)
  {
    //If using cumulative emissions module, initialise temperature anomaly
    Tmixed(1)=T_0_cumemiss;
  }
  //Initialise cumulative emissions, yearly emissions, etc.
  Cum_emissions=Cum_emissions_0;                   
  Emiss_yearly_calib=Emiss_yearly_0;
  g_rate_em_y=0;
  Emiss_TOT=0;
  Emiss_yearly=0;

  //Consumer price index
	cpi(2)=(1+mi2)*(w0/A0+c_en(2)/A0_en);			
  cpi(3)=(1+mi2)*(w0/A0+c_en(2)/A0_en);
  cpi(4)=(1+mi2)*(w0/A0+c_en(2)/A0_en);
  cpi(5)=(1+mi2)*(w0/A0+c_en(2)/A0_en);
  initial_price=cpi(2);

  //Initialise cost vectors
	C=w0/A0+c_en(2)/A0_en;							
	c2=w0/A0+c_en(2)/A0_en;
	c1=w0/(A0*pm)+c_en(2)/(A0_en*pm);
  //Competitiveness and ability to serve demand
	Em2=1;							
	l2=1;				
  E2=1;
	f1=1/N1r;
	f2=1/N2r;					
  //Set initial capital stock in terms of productive capacity					
	K=K0;								
  //Initial labour force
	LS=LS0;
  //t0 needed when iterating over technology arrays
	t0=1;

  for (i=1; i<=N2; i++)
	{
    FirmDefaultProbability(i)=FirmDefaultProbability_init;
  }

  Bank_active=1;   
  //Set initial baseline lending rate; additive mark-up over CB rate
  bankmarkup=bankmarkup_init;
  r_deb=r+bankmarkup;

  //initialise debt service
  DebtService_2=L0*(r+bankmarkup_init)+repayment_share*L0;

	EI=0;				
  //Initial investment is given by amount necessary to replace average amount of machines reaching max. age
	I=ROUND((((K0)/(agemax+1)))/dim_mach)*dim_mach;
	mu2=mi2;
	Td.element(0)=0;
	step=N2/N1;							
	counter=0;
	Match=0;

  //Match C-firms to K-firms & distribute initial stocks of machines
	for (i=1; i<=N1; i++)
	{
		counter+=step;
		for (j=0; j<step; j++)
		{
			Match(counter-j,i)=1;
			supl(counter-j)=i;
		}
	}

	for (j=1; j<=N2; j++)					
	{
		n_mach(j)=K(j)/dim_mach;
		while (n_mach(j) > 0)
		{
			i++;	
			if (i > N1)
      {
				i=1;
      }
      //Initial capital stock is assumed to not come from initial supplier
      if (supl(j) != i)
			{
				//Random age
        age0=int(ran1(p_seed)*(agemax+1))%int((agemax+1))+1;
        g[0][i-1][j-1]++;
        g_price[0][i-1][j-1]=p1(i);	
				gtemp[0][i-1][j-1]++;
				g_c[0][i-1][j-1]++;
        g_c2[0][i-1][j-1]++;
        g_c3[0][i-1][j-1]++;
				age[0][i-1][j-1]=age0;
        n_mach(j)--;
			}
		}
	}

  //Set initial demand for consumption goods to a level consistent with unemployment=ustar
  D20=w0*(LS0*(1-ustar)-I.Sum()/dim_mach/(A0*pm));
	D2=D20/N2r;				
  De=D2.Row(1);		
  N=omicron*D20/N2r;
  //Nominal value of C-firm inventories
  Inventories=(omicron*D20/N2r)*p2(1);
  NW_2=Deposits_2+CapitalStock+Inventories-Loans_2;
  //Assume initial firm sales are homogeneous
	S1=(((I(1)/dim_mach)*N2r)/N1r)*p1;
	S2=D20/N2r*p2(1);
  //Rough measure of initial net revenue
  mol=D20/N2r*p2(1)-D20/N2r/A0*w0-D20/N2r/A0_en*c_en(2);
  //Initial dividends
  Dividends(2)=(mol(1)-L0*(r+bankmarkup_init))*N2r*d2+(S1.Sum()-(((I.element(1)/dim_mach)*N2r))*c1(1))*d1+db*(r_bonds*GB_b.Row(1).Sum()+L0*(r+bankmarkup_init)*N2r)+de*c_en(2)*(De.Sum()/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm));
  U(2)=(LS0-((D20)/A0+(((I.element(1)/dim_mach)*N2r))/(A0*pm)))/LS0;
  
  //Give the energy sector an initial stock of productive capacity
  G_ge(1)=K_ge0_perc*(1.1*(D20/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm)));
  G_ge_n(1)=CF_ge(1)*G_ge(1);
  G_de(1)=(1.1*(D20/A0_en+((I(1)/dim_mach)*N2r)/(A0_en*pm)))-G_ge(1);
  C_de(1)=pf/A_de(1)+t_CO2_en*EM_de(1);   
  K_gelag=G_ge(1);
  K_delag=G_de(1);        
  G_de_0=G_de(1);
  G_ge_0=G_ge(1);
  G_ge_n_0=G_ge_n(1);
  CapitalStock_e=G_ge_n_0;
  NW_e=Deposits_e+CapitalStock_e;

  //Initialise technological change
  RD.Row(1)=nu*S1;
  RD.Row(2)=nu*S1;						
  t=0;
  TECHANGEND();						
}

void SETVARS(void)
{
  for (j=1; j<=N2; ++j)
  { 
    Loans_2(1,j)=Loans_2(2,j);
    Deposits_2(1,j)=Deposits_2(2,j);
    deltaCapitalStock(1,j)=0;
    InterestDeposits_2(j)=0;
    CapitalStock(1,j)=CapitalStock(2,j);
    c2(j)=0;
		A2(j)=0;
    A2_mprod(j)=0;
    A2e(j)=0;
    A2e_en(j)=0;
    A2e_ef(j)=0;
    A2e2(j)=0;
    A2e_en2(j)=0;
    A2e_ef2(j)=0;
	  c2e(j)=0;
    EI(1,j)=0;		
		SI(j)=0;
		I(j)=0;
    Dividends_2(j)=0;
    Taxes_2(j)=0;
    Injection_2(j)=0;
    DebtRemittances2(j)=0;
    A2_en(j)=0;
    A2_ef(j)=0;
    baddebt_2(j)=0;
    S2(1,j)=0;
    S2_temp(1,j)=0;
    D2(1,j)=0;
    Q2(j)=0;
    k(j)=0;
    EId(j)=0;
    SId(j)=0;
    exiting_2(j)=0;
    exit_payments2(j)=0;
    exit_equity2(j)=0;
    exit_marketshare2(j)=0;
    n_mach_entry(j)=0;
    Loss_Capital(j)=0;
    Loss_Inventories(j)=0;
    k_entry(j)=0;
    marker_age(j)=0;
    K_loss(j)=0;
    C_loss(j)=0;
  }

  for (i=1; i<=N1; i++)
  {
    Deposits_1(1,i)=Deposits_1(2,i);
    InterestDeposits_1(i)=0;
    Q1(i)=0;
	  D1(i)=0;
    S1(i)=0;
    S1_temp(1,i)=0;
    Dividends_1(i)=0;
    Taxes_1(i)=0;
    Injection_1(i)=0;
    baddebt_1(i)=0;
    exiting_1(i)=0;
    exiting_1_payments(i)=0;
    I_loss(i)=0;
  }

  for (i=1; i<=NB; i++)
  {
    Loans_b(1,i)=Loans_b(2,i);
    Deposits(1,i)=Deposits(2,i);
    Deposits_hb(1,i)=Deposits_hb(2,i);
    Deposits_eb(1,i)=Deposits_eb(2,i);
    GB_b(1,i)=GB_b(2,i);
    Advances_b(1,i)=Advances_b(2,i);
    Reserves_b(1,i)=Reserves_b(2,i);
    InterestDeposits(i)=0;
    baddebt_b(i)=0;
    Outflows(i)=0;
    Inflows(i)=0;
    NW_b(1,i)=NW_b(2,i);
    Dividends_b(i)=0;
    Taxes_b(i)=0;
    LoanInterest(i)=0;
    bonds_dem(i)=0;
    Bailout_b(i)=0;
    LossAbsorbed(i)=0;
    capitalRecovered(i)=0;
    capitalRecovered2(i)=0;
    capitalRecoveredShare(i)=0;
    LossEntry_b(i)=0;
  }

  K_loss_share=0;
  I_loss_share=0;
  C_loss_share=0;
  Taxes=0;
  Taxes_CO2(1)=0;
  Loan_interest_e=0;
  DebtService_e=0;
  DeafaultedDebtRecovered_e=0;
  BadDebt_e=0;
  DebtRemittances_e=0;
  DebtWrittenOff_e=0;
  DefaultedDeposits_e=0;
  Loans_preDefault_e=0;
  Wages=0;
  Deposits_recovered_1=0;
  Deposits_recovered_2=0;
  EntryCosts=0;
  BankTransfer=0;
  Deposits_h(1)=Deposits_h(2);
  Deposits_e(1)=Deposits_e(2);
  GB_cb(1)=GB_cb(2);
  GB(1)=GB(2);
  Advances(1)=Advances(2);
  Reserves(1)=Reserves(2);
  Dividends(1)=0;
  EnergyPayments=0;
  InterestReserves=0;
  InterestAdvances=0;
  FirmTransfers=0;
  FirmTransfers_1=0;
  FirmTransfers_2=0;
  Divtot_1=0;
  Divtot_2=0;
  Pitot1=0;
	Pitot2=0;
  LD2tot=0;
	LD1tot=0;
  dNtot=0;
	dNmtot=0;
  DebtServiceToSales2_bank=0; 
  Bailout=0;
  NewBonds=0;
  cpi(1)=0;
  kpi=0;
  Em2(1)=0;
  ns1=0;
  ns2=0;
  mD1=0;
  mD2=0;
  n_mach_exit=0;
  n_mach_needed=0;
  Emiss_en=0;
  Am(1)=0;
  Am2=0;
  Am1=0;
  Am_a=0;
  Am_en(1)=0;
  ftot=0;
  Consumption_r=0;
  A_mi=0;
  A1_mi=0;
  A2_en_mi=0;
  A2_ef_mi=0;
  A1_en_mi=0;
  A1_ef_mi=0;
	A_sd=0;
  H1=0; 
	H2=0;
	HB=0; 
  PC_en=0;
  IC_en=0;
  p2_entry=0;
  f2_exit=0;
  EntryShare=0;
  CompEntry=0;
  C_secondhand=std::numeric_limits<double>::infinity();
  counter_bankfailure=0;
  FuelCost=0;
  Dividends_e=0;
  Wages_en=0;
  Loans_e(1)=0;
  for(i=0; i<ene_tecs.size(); i++)
  {
    tech=ene_tecs[i];
    exiting_e_mult_tech[tech]=0;
  }
}

void DEPOSITINTEREST(void)
{
  //Firms, Households and the Energy sector receive deposit interest
  for (i=1; i<=N1; i++)
	{
    sendingBank=BankingSupplier_1(i);
    InterestDeposits_1(i)=r_depo*Deposits_1(2,i);
    Deposits_1(1,i)+=InterestDeposits_1(i);
    Deposits(1,sendingBank)+=InterestDeposits_1(i);
    InterestDeposits(sendingBank)+=InterestDeposits_1(i);
  }

  for (j=1; j<=N2; j++)
	{
    sendingBank=BankingSupplier_2(j);
    InterestDeposits_2(j)=r_depo*Deposits_2(2,j);
    InterestDeposits(sendingBank)+=InterestDeposits_2(j);
    Deposits_2(1,j)+=InterestDeposits_2(j);
    Deposits(1,sendingBank)+=InterestDeposits_2(j);
  }

  InterestDeposits_h=r_depo*Deposits_h(2);
  Deposits_h(1)+=InterestDeposits_h;

  if(flag_WITCH_on==0)
  {
    InterestDeposits_e=r_depo*Deposits_e(2);
    Deposits_e(1)+=InterestDeposits_e;
  }
  else
  {
    InterestDeposits_e=0;
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      interestsDeposits_e_mult_tech[tech]=r_depo*deposits_e_mult_tech[tech];
      InterestDeposits_e+=interestsDeposits_e_mult_tech[tech];
      deposits_e_mult_tech[tech]+=interestsDeposits_e_mult_tech[tech];
      Deposits_e(1)+=interestsDeposits_e_mult_tech[tech];
    }
  }

  for(i=1; i<=NB; i++)
  {
    InterestDeposits(i)+=r_depo*Deposits_hb(2,i);
    Deposits_hb(1,i)+=r_depo*Deposits_hb(2,i);
    Deposits(1,i)+=r_depo*Deposits_hb(2,i);
    InterestDeposits(i)+=r_depo*Deposits_eb(2,i);
    Deposits_eb(1,i)+=r_depo*Deposits_eb(2,i);
    Deposits(1,i)+=r_depo*Deposits_eb(2,i);
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
    
    if(Deposits_eb.Row(1).Sum()>0)
    {
      DepositShare_e(i)=Deposits_eb(1,i)/Deposits_eb.Row(1).Sum();
    }
    else
    {
      DepositShare_e(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }
  }
}

void MACH(void)	 
{  
  if ((int)vintage_energy.size() < N1*(t-t0+1))
  {
    vintage_energy.resize(N1*(t-t0+1));
    vintage_emission.resize(N1*(t-t0+1));
  }

  //Iterate over all K-Firms and all machine vintages still in use (newer than t0)
  for (i=1; i<=N1; i++)
  {
    for (tt=t0; tt<=t; tt++)
    {
      if (A(tt,i) > 0 & A_en(tt,i)>0)
      {
        //Calculate the unit cost of production implied by using a machine of vintage tt
        //produced by i to produce consumption goods in the current period
        //The energy and emission-tax halves of the unit cost depend on the
        //vintage and not on the firm using it, so they are kept for COSTPROD
        //rather than divided out again for each of the 200 firms. The sum
        //below is the same one, in the same order, as it always was.
        vintage_energy[(i-1)*(t-t0+1)+(tt-t0)]=c_en(2)/A_en(tt,i);
        vintage_emission[(i-1)*(t-t0+1)+(tt-t0)]=t_CO2*A_ef(tt,i)/A_en(tt,i);
        C(tt,i)=w(2)/A(tt,i)+vintage_energy[(i-1)*(t-t0+1)+(tt-t0)]+vintage_emission[(i-1)*(t-t0+1)+(tt-t0)];
      }
      else
      {
        //Machine labour productivity (A) and energy efficiency (A_en) should be greater than 0
        std::cerr << "\n\n ERROR: A_en(tt,i) or A(tt,i) <= 0 in period " << t << " for K-firm "<< i << endl;
        Errors << "\n A_en(tt,i) or A(tt,i) <= 0 in period " << t << " for K-firm "<< i << endl;
        exit(EXIT_FAILURE);
      }
    }

    //Calculate unit cost of production for K-Firm i (unit labour cost+unit energy cost+emission tax to be paid per unit of output)
    c1(i)=w(2)/((1-shocks_labprod1(i))*A1p(i))+c_en(2)/((1-shocks_eneff1(i))*A1p_en(i))+t_CO2*A1p_ef(i)/((1-shocks_eneff1(i))*A1p_en(i));

    p1(i)=(1+mi1)*c1(i);

    //Ensure that price does not fall below a lower bound
    if (p1(i) < pmin)
    {
      p1(i)=pmin;
    }
  }


  //C-firms receive capital ordered in the last period
  for (j=1; j<=N2; j++)						
	{													        
    //Productive capacity is increased by amount of expansion investment (EI is expressed in units producible)
    //and decreased by capacity which was scrapped in the previous period
    K(j)+=EI(2,j)-scrap_age(j);		
    //Nominal value of the capital stock is updated by changes from last period
    CapitalStock(1,j)+=deltaCapitalStock(2,j);					
    //Re-initialise the variable recording amount of productive capacity to be scrapped
    scrap_age(j)=0;	
  }

  //Update the machine frequency arrays. gtemp incorporates all changes to
  //stocks of individual machine vintages held by each firm; g is set equal to
  //it to take in everything that happened over the previous period, and g_c,
  //g_c2 and g_c3 are the copies used when calculating C-Firms' production
  //costs. All four take the same values, so all four are copied from gtemp,
  //and a whole period's suppliers and firms lie together, so each is one copy
  //rather than one per supplier.
  //g_c2 and g_c3 are read in one place, ADJUSTEMISSENLAB, which is reached only
  //from the two flag_capshocks branches in PRODMACH. With that flag off nothing
  //reads them, so they are refreshed only when something will.
  const int capital_shocks_run=(flag_capshocks!=0);
  for (tt=t0; tt<=t; tt++)
  {
    const double* source=gtemp[tt-1][0].p;
    const std::size_t bytes=(std::size_t)N1*N2*sizeof(double);
    std::memcpy(g[tt-1][0].p, source, bytes);
    std::memcpy(g_c[tt-1][0].p, source, bytes);
    if (capital_shocks_run)
    {
      std::memcpy(g_c2[tt-1][0].p, source, bytes);
      std::memcpy(g_c3[tt-1][0].p, source, bytes);
    }
  }


  //Number of machines owned by each C-firm
  for (j=1; j<=N2; j++)
  {
    n_mach(j)=K(j)/dim_mach;
    if (!(n_mach(j)>0))
    {
      //Every C-Firm should own at least 1 machine. If j does not, this indicates a bug somewhere
      std::cerr << "\n\n ERROR n_mach = 0 in period " << t << " for C-firm "<< j << endl;
      Errors << "\n n_mach = 0 in period " << t << " for C-firm "<< j << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Each C-firm's unit cost is a weighted average of the unit costs implied by
  //all of the machine vintages of which it owns one or more units, and its
  //labour productivity, energy efficiency and environmental friendliness are
  //weighted averages of those of the same vintages. The vintage indices are
  //outside the firm index here rather than inside it: each firm still
  //accumulates its own terms in the same order, but a vintage's productivity
  //is read once for all firms rather than once per firm, and the firm loop
  //runs down contiguous memory.
  //The two shock factors below are the same for every vintage a firm holds, so
  //they are worked out once for the period rather than once for each of the
  //52 million firm-vintage pairs the loop after this one visits.
  //
  //A firm holds units of about two vintages in a hundred, so nearly all of
  //those pairs contribute a productivity multiplied by a machine count of
  //zero, which leaves each running total exactly where it was and can be
  //skipped. That is only true while the divisors are non-zero: with a zero
  //shock factor the term would be a NaN, which does not leave a total where it
  //was, so when any firm has one the loop adds every term as before.
  int weights_sparse=1;
  for (j=1; j<=N2; j++)
  {
    weight_labprod2[j-1]=1-shocks_labprod2(j);
    weight_eneff2[j-1]=1-shocks_eneff2(j);
    if (weight_labprod2[j-1]==0 || weight_eneff2[j-1]==0)
    {
      weights_sparse=0;
    }
  }

  //While the loop below is reading every firm's count of every vintage anyway,
  //it records which firms hold each one. SCRAPPING wants that same set and was
  //reading all 52 million counts a run to find it. Recording them here rather
  //than in a pass of its own costs less than that pass would: it is one store
  //per firm that holds something, against reading every count again. Nothing
  //writes g in between, and the firms are recorded in increasing order, which
  //is the order SCRAPPING met them in.
  const int n_vintage_held=t-t0+1;
  if ((int)holders.size() < N1*n_vintage_held*N2)
  {
    holders.resize(N1*n_vintage_held*N2);
    holder_count.resize(N1*n_vintage_held);
  }

  for (i=1; i<=N1; i++)
  {
    for (tt=t0; tt<=t; tt++)
    {
      const double A_it=A(tt,i);
      const double A_en_it=A_en(tt,i);
      const double A_ef_it=A_ef(tt,i);
      const double w_2=w(2);
      const double c_en_2=c_en(2);
      const VintageRow<double> g_it=g[tt-1][i-1];

      //The five running totals and the firm's machine count, reached through
      //their own storage rather than through a bounds-checked accessor apiece.
      double* c2_s=c2.Store();
      double* A2_s=A2.Store();
      double* A2_mprod_s=A2_mprod.Store();
      double* A2_en_s=A2_en.Store();
      double* A2_ef_s=A2_ef.Store();
      const double* n_mach_s=n_mach.Store();

      const int vintage_row=(i-1)*n_vintage_held+(tt-t0);
      int* row_holders=&holders[(size_t)vintage_row*N2];
      int n_holders=0;

      for (j=1; j<=N2; j++)
      {
        const double g_itj=g_it[j-1];
        if (weights_sparse && g_itj==0)
        {
          continue;
        }
        if (g_itj != 0)
        {
          row_holders[n_holders++]=j;
        }
        const double labprod2=weight_labprod2[j-1];
        const double eneff2=weight_eneff2[j-1];
        const double n_machj=n_mach_s[j-1];

        c2_s[j-1]+=(w_2/(labprod2*A_it)+c_en_2/(eneff2*A_en_it)+t_CO2*A_ef_it/(eneff2*A_en_it))*g_itj/n_machj;
        A2_s[j-1]+=labprod2*A_it*g_itj/n_machj;
        A2_mprod_s[j-1]+=A_it*g_itj/n_machj;
        A2_en_s[j-1]+=eneff2*A_en_it*g_itj/n_machj;
        A2_ef_s[j-1]+=A_ef_it*g_itj/n_machj;
      }
      holder_count[vintage_row]=n_holders;
    }
  }

  //C-firms revise their mark-up and set their price
  for (j=1; j<=N2; j++)
  {
    //C-Firm j updates its mark-up based on the previous period's change in its market share
    if (f2(3,j)>0)
    {
      mu2(1,j)=max(0.0,mu2(2,j)*(1+deltami2*((f2(2,j)-f2(3,j))/f2(3,j))));
    }
    else
    {
      mu2(1,j)=mu2(2,j);
    }

    if (mu2(1,j) <= 0) 
    {
      //Mark-ups should be positive
      std::cerr << "\n\n ERROR: mark-up out of range in period " << t << " for C-firm "<< j << endl;
      Errors << "\n Mark-up out of range in period " << t << " for C-firm "<< j << endl;
      exit(EXIT_FAILURE);
    }

    p2(j)=(1+mu2(1,j))*c2(j);						

    //Ensure that price does not fall below a lower bound
    if (p2(j) < pmin)
    {
      p2(j)=pmin;
    }

  }
}

void BROCHURE(void)			
{
  //K-Firms send brochures to attract customers
  //Re-initialise number of clients of each K-Firm
  nclient=0;

  //If a C-Firm does not have a valid supplier of capital goods, randomly assign one
  //This happens e.g. if j's supplier exited in the previous period
	for (j=1; j<=N2; j++)  
	{											 
    if (supl(j) < 1 || supl(j) > N1)
		{
      supl(j)=int(ran1(p_seed)*N1*N2)%N1+1;
      Match(j,supl(j))=1;
    }
  }
  
  //K-firms send brochures to potential customers
	for (i=1; i<=N1; i++)  
	{							        
		//Count number of C-Firms matched to K-Firm i
    for (j=1; j<=N2; j++)
    {
			nclient(i)+=Match(j,i);  
    }

    //Number of brochures sent is a function of number of existing clients
		newbroch=int(ROUND(nclient(i)*Gamma));
		
    //Ensure that every firm sends at least 1 brochure
    if (newbroch==0)
    {
      newbroch++;
    }

    //Brochures are sent to randomly drawn C-Firms
    //Note that this does not ensure that a randomly drawn C-Firm is not already a customer of i!
    while (newbroch > 0)
		{
			rni=int(ran1(p_seed)*N1*N2)%N2+1;
      Match(rni,i)=1;
      newbroch--;
		}
	}

  //C-firms choose their preferred supplier of machine tools
  for (j=1; j<=N2; j++)  
	{									     
		//current supplier of j as integer
    indsupl=int(supl(j));
		for (i=1; i<=N1; i++)
		{
			if (A1(i) > 0)
			{
        //If j has received a brochure from i and i's technology is more convenient, j switches
        if (Match(j,i)==1 && p1(i)+(w(2)/A1(i)+c_en(2)/A1_en(i)+t_CO2*A1_ef(i)/A1_en(i))*b < p1(indsupl)+(w(2)/A1(indsupl)+ c_en(2)/A1_en(indsupl)+t_CO2*A1_ef(indsupl)/A1_en(indsupl))*b)
        {
          indsupl=i;
        }
			}
			else 
      {
        //Technology offered by i should imply a positive labour productivity
        std::cerr << "\n\n ERROR: A1(i) = 0 in period " << t << " for K-firm "<< i << endl;
        Errors << "\n A1(i) = 0 in period " << t << " for K-firm "<< i << endl;
        exit(EXIT_FAILURE);
      }
    }
    //Update supplier of j
		supl(j)=indsupl;
		
    //Reset the K-Firm C-Firm network by setting the entries of all i from which j has received
    //brochures but which j did not choose as its supplier to zero
    for (i=1; i<=N1; i++)
		{
			if (i != indsupl)
      {
				Match(j,i)=0;
      }
		}
	}

  //Update the number of clients of each K-Firm
	nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			if (Match(j,i) == 1)
      {
				nclient(i)++;
      }
		}
	}
}

void INVEST(void)
{
  //Every firm's scrapping decision, worked out in one sweep before the loop
  //below rather than once per firm inside it.
  SCRAPPING();

	for (j=1; j<=N2; j++)
	{
    //C-firms determine expected demand, desired production, and demand for investment
    De(j)=alfa*De(j)+(1-alfa)*D2(2,j);
    if (De(j)<=0)
    {
        De(j)=1;
    }

    //Desired inventories                    
    Ne(j)=De(j)*omicron; 
    //Desired output                             
		Qd(j)=De(j)+Ne(j)-N(2,j);        

		if (Qd(j) < 0)
    {
			Qd(j)=0;
    }

    //Desired output implies desired capital stock
    Kd(j)=Qd(j)/u;


    //Determine the capital stock which the firm will have available in t+1 once aged machines are scrapped
    //Machines scrapped due to age can still be used in t but will be removed at the beginning of t+1
		Ktrig(j)=ROUND((K(j)-scrap_age(j))/dim_mach)*dim_mach;	 

    //If the desired capital stock is larger than what the firm will have available after scrapping, it wants to engage in expansion investment
		if (Kd(j) >= Ktrig(j))						          
		{								                      
      //If there is a constraint on expansion investment, determine the maximum capital stock which can be achieved
      if(I_max>0)
      {
        K_top=K(j)*(1+I_max);
        K_top=ROUND(K_top/dim_mach)*dim_mach;
        if(K_top<(K(j)+dim_mach))
        {
          K_top+=dim_mach;
        }
      }
      else
      {
        K_top=Kd(j)+1;
      }

      if(Kd(j)>K_top)
      {
        EId(j)=K_top-Ktrig(j);
      }
      else
      {
        EId(j)=floor((Kd(j)-(K(j)-scrap_age(j)))/dim_mach)*dim_mach;
      }
    }
		else 
    {
      EId(j)=0;
    }
    
    if(SId(j)==0 && EId(j)==0 && marker_age(j)==1)
    {
      EId(j)=dim_mach;
    }

		if (Qd(j) > K(j))							
		{											        
			Qd(j)=K(j);								  
		}							 

    //C-firms determine effective production cost
		if (Qd(j) > 0 && Qd(j) < K(j))
		{
			COSTPROD();								  
		}
		else
		{
			A2e(j)=A2(j);
			c2e(j)=c2(j);
      A2e_en(j)=A2_en(j);
      A2e_ef(j)=A2_ef(j);
		} 
	}

	ORD();                   
}

void SCRAPPING(void) 
{
  //Every C-firm's scrapping decision for this period. The firm is the
  //innermost dimension of the vintage arrays, so sweeping the firms inside the
  //vintages reads along memory rather than across it, and reads each vintage's
  //own quantities once for all firms rather than once per firm.
  //
  //Each firm meets its own vintages in the order i then tt, so its running
  //totals accumulate the same terms in the same sequence they would if it were
  //the only firm here. Nothing INVEST does between firms touches what this
  //reads: K, supl, g, age, C, A and the supplier's own productivities are all
  //left alone there.

  //Per firm, what does not depend on which vintage is being looked at. The
  //three subtracted terms stay separate rather than being summed here, so the
  //payback denominator is still evaluated left to right as it was.
  for (j=1; j<=N2; j++)
  {
    K_temp(j)=K(j)/dim_mach;
    const int supplier=int(supl(j));
    scrap_p1[j-1]=p1(supplier);
    scrap_wage[j-1]=w(2)/A1(supplier);
    scrap_energy[j-1]=c_en(2)/A1_en(supplier);
    scrap_emission[j-1]=t_CO2*A1_ef(supplier)/A1_en(supplier);
    scrap_supplier_ok[j-1]=(A1(supplier) > 0 && A1_en(supplier) > 0) ? 1 : 0;
  }
  //The supplier index the loop over firms ends on.
  indsupl=int(supl(N2));

  const int n_vintage_scrap=t-t0+1;

  //Which machines each firm is left wanting to scrap, and what each of them
  //costs it to run. CANCMACH is the only other place g_pb and C_pb are used,
  //and it works through this list, so entries outside the list are never read
  //and the two arrays no longer have to be cleared for every firm and vintage
  //before they are filled in for the few.
  marked_capacity=N1*(t-t0+1);
  if ((int)marked_supplier.size() < N2*marked_capacity)
  {
    marked_supplier.resize(N2*marked_capacity);
    marked_vintage.resize(N2*marked_capacity);
  }
  std::fill(marked_count.begin(), marked_count.end(), 0);

  for (i=1; i<=N1; i++)
	{
		for (tt=t0; tt<=t; tt++)
		{
      const double A_it=A(tt,i);
      const double A_en_it=A_en(tt,i);
      const double C_it=C(tt,i);
      const double w_2=w(2);

      //The firms holding this vintage, in increasing order, as MACH recorded
      //them. A firm holds units of about two vintages in a hundred, and the
      //rest reach neither branch below.
      const int vintage_row=(i-1)*n_vintage_scrap+(tt-t0);
      const int* row_holders=&holders[(size_t)vintage_row*N2];
      const int n_holders=holder_count[vintage_row];

      for (int holder=0; holder<n_holders; holder++)
      {
        j=row_holders[holder];

        double scrapped_for_age=0;

        //If a machine has reached its maximum age, it is scrapped, unless the firm has only 1 machine remaining
        if (g[tt-1][i-1][j-1] > 0 && age[tt-1][i-1][j-1] > (agemax))
			  {
          scrapped_for_age=min(g[tt-1][i-1][j-1],(K_temp(j)-1));
          g_pb[tt-1][i-1][j-1]=scrapped_for_age;
          C_pb[tt-1][i-1][j-1]=C_it;
          if (scrapped_for_age > 0)
          {
            marked_supplier[(j-1)*marked_capacity+marked_count[j-1]]=i;
            marked_vintage[(j-1)*marked_capacity+marked_count[j-1]]=tt;
            marked_count[j-1]++;
          }
          scrap_age(j)+=dim_mach*scrapped_for_age;	  
          K_temp(j)-=scrapped_for_age;
          //If the firm has only one machine left, set the marker_age flag. This will ensure later that
          //its desired investment is at least 1 machine
          if(K_temp(j)==1)
          {
            marker_age(j)=1;
          }
			  }

        //Machines which have not reached their maximum age are scrapped if a superior technology is available from the current supplier
        //Only enter here if the relevant machines have not already been scrapped due to age
        if (g[tt-1][i-1][j-1] > 0 && scrapped_for_age==0)
			  {
          if (w_2 > 0 && A_it > 0 && scrap_supplier_ok[j-1] && A_en_it > 0)
          {
            //calculate payback variable: (price of machine)/(unit cost difference offered by the machine)
            //Payback can be interpreted as number of units which have to be produced with a new machine sucht that
            //the savings in unit cost are equal to the purchase price of the machine
            //The first three terms are the unit cost of this vintage, which MACH
            //put in C(tt,i) this period from the same w(2), c_en(2), t_CO2, A,
            //A_en and A_ef; nothing between the two touches any of them.
            payback=scrap_p1[j-1]/(C_it-scrap_wage[j-1]-scrap_energy[j-1]-scrap_emission[j-1]);
          }
          else 
          {
            //All of the variables in the denominators of the formula above should never become zero or negative
            std::cerr << "\n\n ERROR: payback division by zero in period " << t << " for C-firm "<< j << endl;
            Errors << "\n Payback division by zero in period " << t << " for C-firm "<< j << endl;
            exit(EXIT_FAILURE);
          }

          //If payback is smaller than an exogenous threshold, firm wants to replace the machine in question with the newer vintage
				  if (payback <= b && payback>0)
				  {
					  g_pb[tt-1][i-1][j-1]=g[tt-1][i-1][j-1];
		 			  C_pb[tt-1][i-1][j-1]=C_it;
					  SId(j)+=dim_mach*g_pb[tt-1][i-1][j-1];
            marked_supplier[(j-1)*marked_capacity+marked_count[j-1]]=i;
            marked_vintage[(j-1)*marked_capacity+marked_count[j-1]]=tt;
            marked_count[j-1]++;
				  }
			  }
      }
		}
	}
}

void COSTPROD(void)
{     
  //C-firms determine effective production cost based on desired production; most efficient machines used first
  nmachprod=ceil(Qd(j)/dim_mach);  
	nmp_temp=nmachprod;

  //The vintages j actually holds a machine of, with the unit cost each implies
  //for it. Only these can win the comparison below, so the search runs over
  //them rather than over every vintage in existence, and only these have their
  //cost worked out. A firm holds units of about two vintages in a hundred.
  //
  //They are collected in the order i and then tt, so vintages of equal cost
  //resolve in favour of the same one the full scan resolved in favour of. The
  //cost does not change while j draws its machines down; only the counts in
  //g_c do.
  const int n_vintage=t-t0+1;
  if ((int)vintage_cost.size() < N1*n_vintage)
  {
    vintage_cost.resize(N1*n_vintage);
    held_supplier.resize(N1*n_vintage);
    held_vintage.resize(N1*n_vintage);
  }
  int n_held=0;
  {
    const double w_2=w(2);
    const double labprod2=1-shocks_labprod2(j);
    //j's own count for each vintage, reached by walking rather than by
    //subscripting. The subscript works the address out from the array's base
    //and strides every time, and it has to reload all three each time round,
    //because the writes below are to other arrays the compiler cannot prove sit
    //elsewhere in memory. One vintage on is a fixed distance, so the walk is
    //one addition.
    const double* const machines=g_c.data.data();
    const size_t vintage_stride=(size_t)g_c.d1*g_c.d2;
    for (i=1; i<=N1; i++)
    {
      const double* machines_held=machines+(size_t)(t0-1)*vintage_stride
                                          +(size_t)(i-1)*g_c.d2+(j-1);
      for (tt=t0; tt<=t; tt++, machines_held+=vintage_stride)
      {
        if (*machines_held > 0)
        {
          vintage_cost[n_held]=w_2/(labprod2*A(tt,i))+vintage_energy[(i-1)*n_vintage+(tt-t0)]
                                                    +vintage_emission[(i-1)*n_vintage+(tt-t0)];
          held_supplier[n_held]=i;
          held_vintage[n_held]=tt;
          n_held++;
        }
      }
    }
  }

	while (nmp_temp > 0)
	{
    cmin=std::numeric_limits<double>::infinity();
    imin=0;
    jmin=0;
    tmin=0;

    for (int held=0; held<n_held; held++)
    {
      if (g_c[held_vintage[held]-1][held_supplier[held]-1][j-1] > 0 && vintage_cost[held] < cmin)
      {
        cmin=vintage_cost[held];
        imin=held_supplier[held];
        jmin=j;
        tmin=held_vintage[held];
      }
    }

    if (nmachprod>0)
    {
      if (g_c[tmin-1][imin-1][jmin-1] >= nmp_temp)  
      {
        A2e(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*nmp_temp/nmachprod;
        A2e_en(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*nmp_temp/nmachprod;
        A2e_ef(j)+=A_ef(tmin,imin)*nmp_temp/nmachprod;
        c2e(j)+=(w(2)/((1-shocks_labprod2(j))*A(tmin,imin))+c_en(2)/((1-shocks_eneff2(j))*A_en(tmin,imin))+t_CO2*A_ef(tmin,imin)/((1-shocks_eneff2(j))*A_en(tmin,imin)))*nmp_temp/nmachprod;
        g_c[tmin-1][imin-1][jmin-1]-= nmp_temp;
        nmp_temp=0;                 
      }
      else                      
      {
        A2e(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_en(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_ef(j)+=A_ef(tmin,imin)*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        c2e(j)+=(w(2)/((1-shocks_labprod2(j))*A(tmin,imin))+c_en(2)/((1-shocks_eneff2(j))*A_en(tmin,imin))+t_CO2*A_ef(tmin,imin)/((1-shocks_eneff2(j))*A_en(tmin,imin)))*g_c[tmin-1][imin-1][jmin-1]/nmachprod;
        nmp_temp-=g_c[tmin-1][imin-1][jmin-1];
        g_c[tmin-1][imin-1][jmin-1]=0;        
      }                       
    }
    else 
    {
      std::cerr << "\n\n ERROR: nmachprod = 0!!!" << endl;
      Errors << "\n nmachprod = 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
	} 
}

void ORD(void)  
{
  //C-firms calculate internal funds & the maximum amount of loans they are willing to take on; based on this, investment is possibly scaled back
  for (j=1; j<=N2; ++j)
  {
		MaxFunds=max(0.0,Deposits_2(1,j)+phi2*mol(j)-Loans_2(1,j)-c2e(j)*Qd(j));

    indsupl=int(supl(j));
		p1test=p1(indsupl);					
    if ((EId(j)/dim_mach)*p1test < MaxFunds)
    {
      EIp(j)=EId(j);   
      MaxFunds-=(EId(j)/dim_mach)*p1(indsupl);
    }
    else
    {
      p1test=p1(indsupl);
      EIp(j)=floor(MaxFunds/p1(indsupl))*dim_mach;
      if(EIp(j)<0)
      {
        EIp(j)=0;
        MaxFunds=0;
      }else
      {
        MaxFunds=0;
      }
    }
			
    if ((SId(j)/dim_mach)*p1test < MaxFunds)
    {
      SIp(j)=SId(j);
    }
    else
    {
      SIp(j)=floor(MaxFunds/p1(indsupl))*dim_mach;
    }

    Ip(j)=EIp(j)+SIp(j);  

    //Determine cost of planned investment; NB: EI and SI are expressed in terms of productive capacity, not number of machines, hence need to be divided by dim_mach
		if (Ip(j) > 0)
		{
      CmachEI(j)=p1(indsupl)*EIp(j)/dim_mach;           
			CmachSI(j)=p1(indsupl)*SIp(j)/dim_mach;   
			Cmach(j)=p1(indsupl)*Ip(j)/dim_mach;			
		}
		else
    {
      CmachEI(j)=0;       
			CmachSI(j)=0;  
			Cmach(j)=0;   
    }

  } 
} 

void ALLOCATECREDIT(void)
{
  
  if(flag_WITCH_on==1)
  {
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      if((loans_e_mult_tech[tech]+nom_inv_mult_tech[tech]+prodCost_mult_tech[tech]+maintCost_mult_tech[tech])<=deposits_e_mult_tech[tech])
      {
          creditDemand_e_mult_tech[tech]=0;
          deposits_e_mult_tech[tech]-=loans_e_mult_tech[tech];
          Deposits_e(1)-=loans_e_mult_tech[tech];
          loans_e_mult_tech[tech]=0;
      }
      else
      {
          creditDemand_e_mult_tech[tech]=loans_e_mult_tech[tech]+nom_inv_mult_tech[tech]+prodCost_mult_tech[tech]+maintCost_mult_tech[tech]-deposits_e_mult_tech[tech];
          deposits_e_mult_tech[tech]+=creditDemand_e_mult_tech[tech]-loans_e_mult_tech[tech];
          Deposits_e(1)+=creditDemand_e_mult_tech[tech]-loans_e_mult_tech[tech];
          loans_e_mult_tech[tech]=creditDemand_e_mult_tech[tech];
      }
      Loans_e(1)+=creditDemand_e_mult_tech[tech];  
    }     

    for (i=1; i<=NB; i++)
    {
        Deposits(1,i)=Deposits(1,i)+DepositShare_e(i)*(Loans_e(1)-Loans_e(2)); 
        Deposits_eb(1,i)=Deposits_eb(1,i)+DepositShare_e(i)*(Loans_e(1)-Loans_e(2));    
        Loans_b(1,i)=Loans_b(1,i)+DepositShare_e(i)*(Loans_e(1)-Loans_e(2));
        BankCredit(i)-=DepositShare_e(i)*Loans_e(1);
        BankCredit(i)=max(0.0, BankCredit(i));
    }
  }
  
  //C-firms calculate credit demand based on outstanding loans (need to be rolled over), planned investment and production
  for (j=1; j<=N2; ++j)
  { 
		if (Loans_2(1,j)+Cmach(j)+(c2e(j)*Qd(j))<=Deposits_2(1,j))
    {
      CreditDemand(j)=0;
    }
		else
    {
      CreditDemand(j)=Loans_2(1,j)+Cmach(j)+(c2e(j)*Qd(j))-Deposits_2(1,j);   
    }
	}

  //Banks allocate credit
  for (i=1; i<=NB; i++)
  {
    for (j=1; j<=NL_2(i); j++)
    {
      ColumnMinimum1(DS2_rating,i,rated_firm_2);
      DS2_rating(rated_firm_2,i)=ColumnMaximum(DS2_rating,i)+1;
			if (BankMatch_2(rated_firm_2,i)==1)
      { 
        //If customer does not need credit, they use deposits to repay their outstanding loans
        if (CreditDemand(rated_firm_2)==0)
        {
			    Q2(rated_firm_2)=Qd(rated_firm_2);
          EI(1,rated_firm_2)=EIp(rated_firm_2);
          SI(rated_firm_2)=SIp(rated_firm_2);
          I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
          Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
          Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
          Deposits(1,i)=Deposits(1,i)-Loans_2(1,rated_firm_2);
          Loans_2(1,rated_firm_2)=0;
        }
        else if(CreditDemand(rated_firm_2)>0)
        {
			    //If remaining credit supply of bank is sufficient, demand is fully satisfied
          if (CreditDemand(rated_firm_2) <= BankCredit(i))	 
          {
            Deposits_2(1,rated_firm_2)=max(0.0,Deposits_2(1,rated_firm_2)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2));
            Deposits(1,i)=Deposits(1,i)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2);
            Loans_b(1,i)=Loans_b(1,i)+CreditDemand(rated_firm_2)-Loans_2(1,rated_firm_2);
            Loans_2(1,rated_firm_2)=CreditDemand(rated_firm_2);
				    BankCredit(i)-=CreditDemand(rated_firm_2);
            Q2(rated_firm_2)=Qd(rated_firm_2);
				    EI(1,rated_firm_2)=EIp(rated_firm_2);
    				SI(rated_firm_2)=SIp(rated_firm_2);
    				I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
          }
          else
          { 
				    //If remaining credit supply is insufficient, first remove replacement investment
            if(Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2) <= BankCredit(i))
            {
              if (Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>=0) 
              {
                Loans_b(1,i)=Loans_b(1,i)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=Loans_2(1,rated_firm_2)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
				        Deposits(1,i)=Deposits(1,i)+CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=CmachEI(rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2));
              }
              else 
					    {
                Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
                Deposits(1,i)=Deposits(1,i)-Loans_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=0;
              }

					    BankCredit(i)-=Loans_2(1,rated_firm_2);
              Q2(rated_firm_2)=Qd(rated_firm_2);
              EI(1,rated_firm_2)=EIp(rated_firm_2);
              SI(rated_firm_2)=0;
              I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
            }
            //If remaining credit supply is still insufficient, remove also expansion investment
            else if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)<=BankCredit(i))
            {
              if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>=0)
              {
                Loans_b(1,i)=Loans_b(1,i)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
				        Deposits(1,i)=Deposits(1,i)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=(c2e(rated_firm_2)*Qd(rated_firm_2));
              }
              else
              {
                Loans_b(1,i)=max(0.0,Loans_b(1,i)-Loans_2(1,rated_firm_2));
                Deposits(1,i)=Deposits(1,i)-Loans_2(1,rated_firm_2);
                Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                Loans_2(1,rated_firm_2)=0;
              }
              BankCredit(i)-=Loans_2(1,rated_firm_2);
              Q2(rated_firm_2)=Qd(rated_firm_2);
              EI(1,rated_firm_2)=0;
              SI(rated_firm_2)=0;
              I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
            }
            else if (Loans_2(1,rated_firm_2)+(c2e(rated_firm_2)*Qd(rated_firm_2))-Deposits_2(1,rated_firm_2)>BankCredit(i))
				    {  
              //If possible, scale back desired production to a level which can be financed
              if (Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2)<=BankCredit(i))
              { 
						    Q2(rated_firm_2)=(BankCredit(i)-Loans_2(1,rated_firm_2)+Deposits_2(1,rated_firm_2))/c2e(rated_firm_2);

                //If production needs to be scaled back too much, firm exits
                if (Q2(rated_firm_2) < 1)
						    {
                  if (Loans_2(1,rated_firm_2)>Deposits_2(1,rated_firm_2)) 
                  {
                    Loans_b(1,i)-=Deposits_2(1,rated_firm_2);
                    Deposits(1,i)-=Deposits_2(1,rated_firm_2);
                    Loans_2(1,rated_firm_2)-=Deposits_2(1,rated_firm_2);
                    Deposits_2(1,rated_firm_2)=0;
								    baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2);
                    baddebt_b(i)+=Loans_2(1,rated_firm_2);
    								BankCredit(i)-=Loans_2(1,rated_firm_2);  
                    Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                    Loans_2(1,rated_firm_2)=0;
                  }
                  else
                  {
                    Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                    baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2);
                    Deposits_recovered_2+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                    Deposits(1,i)-=Deposits_2(1,rated_firm_2);
                    Outflows(i)+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                    Deposits_2(1,rated_firm_2)=0;
                    Loans_2(1,rated_firm_2)=0;
                  }
                  Q2(rated_firm_2)=0;
                  exiting_2(rated_firm_2)=1;
    							f2(1,rated_firm_2)=0;
    							f2(2,rated_firm_2)=0;
    							f2(3,rated_firm_2)=0;
                }
                else 
                {
                  Loans_b(1,i)=Loans_b(1,i)+BankCredit(i)-Loans_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2)+BankCredit(i);
                  Deposits(1,i)=Deposits(1,i)-Loans_2(1,rated_firm_2)+BankCredit(i);
                  Loans_2(1,rated_firm_2)=BankCredit(i);
                  BankCredit(i)=0;
                }
    						EI(1,rated_firm_2)=0;
    						SI(rated_firm_2)=0;
    						I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
              }
              //Otherwise, firm is forced to exit
              else 
              {   
  						  Q2(rated_firm_2)=0;
                BankCredit(i)=0;
                if (Loans_2(1,rated_firm_2)>Deposits_2(1,rated_firm_2)) 
                {
                  Loans_b(1,i)-=Deposits_2(1,rated_firm_2);
                  Deposits(1,i)-=Deposits_2(1,rated_firm_2);
                  Loans_2(1,rated_firm_2)-=Deposits_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=0;
                  baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2);
                  baddebt_b(i)+=Loans_2(1,rated_firm_2); 
                  Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                  Loans_2(1,rated_firm_2)=0;
                }
                else
                {
                  Loans_b(1,i)=max(0.0, Loans_b(1,i)-Loans_2(1,rated_firm_2));
                  baddebt_2(rated_firm_2)=Loans_2(1,rated_firm_2)-Deposits_2(1,rated_firm_2);
                  Deposits_recovered_2+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                  Deposits(1,i)-=Deposits_2(1,rated_firm_2);
                  Outflows(i)+=Deposits_2(1,rated_firm_2)-Loans_2(1,rated_firm_2);
                  Deposits_2(1,rated_firm_2)=0;
                  Loans_2(1,rated_firm_2)=0;
                }
    						EI(1,rated_firm_2)=0;
    						SI(rated_firm_2)=0;
    						I(rated_firm_2)=EI(1,rated_firm_2)+SI(rated_firm_2);
                exiting_2(rated_firm_2)=1;
    						f2(1,rated_firm_2)=0;
    						f2(2,rated_firm_2)=0;
    						f2(3,rated_firm_2)=0;
              }
            }
          }
        }
      } 
    }
  }


  //C-firms calculate labour demand based on production which can be financed
	for (j=1; j<=N2; ++j)
  {
      if (A2e(j) > 0)
      {
        Ld2(j)=Q2(j)/A2e(j);
      }
      else 
      {
        std::cerr << "\n\n ERROR: A2e(j) = 0 in period " << t << " for C-firm " << j << endl;
        Errors << "\n ERROR: A2e(j) = 0 in period " << t << " for C-firm " << j << endl;
        exit(EXIT_FAILURE);
      }
  }
}

void PRODMACH(void)
{
  //K-firms receive demand from C-firms and calculate labour demand
  for (i=1; i<=N1; i++)  
  {
		for (j=1; j<=N2; j++)
		{
      if (Match(j,i) == 1)
      {
          D1(i)+=I(j)/dim_mach;
      }
	  } 

    Q1(i)=D1(i);

    if (A1p(i)>0)
    {
			Ld1(i)=Q1(i)/((1-shocks_labprod1(i))*A1p(i));					
    }
    else
    { 
      std::cerr << "\n\n ERROR: A1p(i)=0 in period " << t << " for K-firm " << i << endl;
      Errors << "\n A1p(i)=0 in period " << t << " for K-firm " << i << endl;
      exit(EXIT_FAILURE);
    }
  }

  //NB: These are labour demands from last period!
  if((LD1rdtot+LDentot) > (LS*g_ls))
  {
    std::cerr << "\n\n ERROR: Remaining labour supply is negative in period " << t << endl;
    Errors << "\n Remaining labour supply is negative in period " << t << endl;
    exit(EXIT_FAILURE);
  } 

  //Determine total labour demand
  LABOR();	

  //Calculate Energy demand based on actual production
  for (i=1; i<=N1; i++)
  {
    if (A1p_en(i)>0)
    {
      D1_en(i)=Q1(i)/((1-shocks_eneff1(i))*A1p_en(i));
    }
    else 
    {
      std::cerr << "\n\n ERROR: A1p_en(i)=0 in period " << t << " for K-firm " << i << endl;
      Errors << "\n A1p_en(i) = 0 in period " << t << " for K-firm " << i << endl;
      exit(EXIT_FAILURE);
    }
  }

  for (j=1; j<=N2; j++)    
  {
    if (A2e_en(j) > 0)
    {
      D2_en(j)=Q2(j)/A2e_en(j);
    }
    else
    {
      std::cerr << "\n\n ERROR: A2e_en(j) = 0 in period " << t << " for C-firm " << j <<endl;
      Errors << "\n A2e_en(j) = 0 in period " << t << " for C-firm " << j << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  if((flag_outputshocks==1))
  {
    for (i=1; i<=N1; i++)
    {
      Qpast=Q1(i);
      if (Qpast > 0)
      {
        Q1(i)=floor(Q1(i)*(1-shocks_output1(i)));
        loss=Qpast-Q1(i);
        I_loss(i)+=loss*dim_mach;
        while(loss>0)
        {
          ranj=int(ran1(p_seed)*N1*N2)%N2+1;
          if (Match(ranj,i) == 1 && I(ranj)>0)
          {
            I(ranj)-=dim_mach;
            if (I(ranj)<EI(1,ranj))
            {
              EI(1,ranj)=I(ranj);
            }
            SI(ranj)=I(ranj)-EI(1,ranj);
            loss-=1;
          }
        }
      }
    }
  }

  if((flag_outputshocks==2))
  {
    loss=ceil(Q1.Sum()*shocks_output1(1));
    while(loss>0)
    {
      rani=int(ran1(p_seed)*N1*N2)%N1+1;
      Qpast=Q1(rani);
      if (Qpast > 0)
      {
        if(Q1(rani)>=loss)
        {
          Q1(rani)-=loss;
          loss=0;
        }
        else
        {
          loss-=Q1(rani);
          Q1(rani)=0;
        }
        lossj=Qpast-Q1(rani);
        I_loss(rani)+=lossj*dim_mach;
        while(lossj>0)
        {
          ranj=int(ran1(p_seed)*N1*N2)%N2+1;
          if (Match(ranj,rani) == 1 && I(ranj)>0)
          {
            I(ranj)-=dim_mach;
            if (I(ranj)<EI(1,ranj))
            {
              EI(1,ranj)=I(ranj);
            }
            SI(ranj)=I(ranj)-EI(1,ranj);
            lossj-=1;
          }
        }
      }
    }
  }
  
  //Calculate nominal value of actual investment
  for (j=1; j<=N2; j++)
  {
    indsupl=int(supl(j));
    SI_n(j)=p1(indsupl)*SI(j)/dim_mach;
    EI_n(j)=p1(indsupl)*EI(1,j)/dim_mach;
    S1(indsupl)+=p1(indsupl)*I(j)/dim_mach;
  }

  //Old machines are scrapped; temporary machine frequency arrays are updated based on expansion & replacement investment
  //
  //The three steps below are taken for every firm in turn rather than one firm
  //at a time. Each firm still has its ages cleared after its own machines have
  //been scrapped and before its new machines are added, because CANCMACH and
  //the last loop touch only the firm's own column of these arrays, and
  //clearing them for all firms at once runs down contiguous memory rather than
  //across it once per firm.
	for (j=1; j<=N2; j++)
	{
    CANCMACH();
  }

  for (j=1; j<=N2; j++)
  {
    for (i=1; i<=N1; i++)
    {
      if(I_loss(i)>0)
      {
        I_loss_share+=1;
      }
    }
  }

  //The age of a machine a firm no longer holds is not cleared here. Every place
  //that reads an age asks for it only where the firm's count of that machine is
  //positive - SCRAPPING and CANCMACH under their own count tests, ENTRYEXIT
  //under gtemp - and the three places a count rises from zero all set the age
  //at the same entry: the purchase below, and the two second-hand transfers in
  //ENTRYEXIT. So a stale age is never read, and clearing them cost a sweep of
  //all 52 million firm-vintage pairs a run.

	for (j=1; j<=N2; j++)
	{
		indsupl=int(supl(j));
		gtemp[t-1][indsupl-1][j-1]+=I(j)/dim_mach;
    g_price[t-1][indsupl-1][j-1]=p1(indsupl);
    deltaCapitalStock(1,j)+=EI_n(j);

		if (I(j) > 0)
		{
			if (gtemp[t-1][indsupl-1][j-1] == I(j)/dim_mach){
        age[t-1][indsupl-1][j-1]=0;
      }
		}
	}

  if (flag_capshocks==1)
  {
    //Shock to C-firms' capital stock
    for (j=1; j<=N2; j++)
    {       
        K_cur(j)=K(j);
        loss=0;
        for (i=1; i<=N1; i++)
        {   
            for (tt=t0; tt<=(t); tt++)
            {
                loss+=gtemp[tt-1][i-1][j-1];
            }
        }
        loss=ROUND(loss*shocks_capstock(j));
        K_loss(j)+=loss*dim_mach;
        if(loss>0)
        {
          K_loss_share+=1;
        }
        while(loss>0)
        {
            rani=int(ran1(p_seed)*N1*N2)%N1+1;
            rant=int(ran1(p_seed)*t0*(t))%((t)-t0+1)+t0;
            if (gtemp[rant-1][rani-1][j-1]>0)
            {
                if(loss<gtemp[rant-1][rani-1][j-1])
                {
                    if(rant==t)
                    {
                      deltaCapitalStock(1,j)-=loss*g_price[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]-=loss;
                      if(SI(j)>0)
                      {
                        K(j)-=min(SI(j),loss*dim_mach);
                        loss=max(0.0,loss*dim_mach-SI(j))/dim_mach;
                      }
                      EI(1,j)-=loss*dim_mach;
                      EI(1,j)=max(EI(1,j),0.0);
                      loss=0;
                    }
                    else
                      {
                        Loss_Capital(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*loss);
                        CapitalStock(1,j)-=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*loss);
                        K(j)-=loss*dim_mach;
                        K_cur(j)-=loss*dim_mach;
                        gtemp[rant-1][rani-1][j-1]-=loss;
                        g_c2[rant-1][rani-1][j-1]-=loss;
                        g_c3[rant-1][rani-1][j-1]-=loss;
                        loss=0;
                      }
                }
                else
                {
                    if(rant==t)
                    {
                      deltaCapitalStock(1,j)-=g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1];
                      if(SI(j)>0)
                      {
                        K(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach-EI(1,j);
                      }
                      EI(1,j)=0;
                      loss-=gtemp[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]=0;
                    }
                    else
                    {
                      Loss_Capital(j)+=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1]);
                      CapitalStock(1,j)-=min(CapitalStock(1,j),g_price[rant-1][rani-1][j-1]*gtemp[rant-1][rani-1][j-1]);
                      K(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach;
                      K_cur(j)-=gtemp[rant-1][rani-1][j-1]*dim_mach;
                      loss-=gtemp[rant-1][rani-1][j-1];
                      g_c2[rant-1][rani-1][j-1]-=gtemp[rant-1][rani-1][j-1];
                      g_c3[rant-1][rani-1][j-1]-=gtemp[rant-1][rani-1][j-1];
                      gtemp[rant-1][rani-1][j-1]=0;
                    }
                }
            }
        }
      if(shocks_capstock(j)>0)
      {
        Q2(j)=min(Q2(j),K_cur(j));
        ADJUSTEMISSENLAB();
      }
    }
  }

  if (flag_capshocks==2)
  {
    //Shock to C-firms' capital stock
    loss=0;
    K_temp_sum=0;
    K_cur=K;
    for (j=1; j<=N2; j++)
    {
      for (i=1; i<=N1; i++)
      {   
          for (tt=t0; tt<=(t); tt++)
          {
              loss+=gtemp[tt-1][i-1][j-1];
              K_temp_sum+=gtemp[tt-1][i-1][j-1];
          }
      }
    }

    if(loss==N2)
    {
      loss=0;
    }
    else
    {
      loss=min(K_temp_sum-N2,ROUND(loss*shocks_capstock(1)));
    }

    while(loss>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      lossj=0;
      for (i=1; i<=N1; i++)
      {   
          for (tt=t0; tt<=(t); tt++)
          {
              lossj+=gtemp[tt-1][i-1][ranj-1];
          }
      }

      lossj=min(loss,(lossj-1));
      K_loss(ranj)+=lossj*dim_mach;
      while(lossj>0)
      {
          rani=int(ran1(p_seed)*N1*N2)%N1+1;
          rant=int(ran1(p_seed)*t0*(t))%((t)-t0+1)+t0;
          if (gtemp[rant-1][rani-1][ranj-1]>0)
          {
              if(lossj<gtemp[rant-1][rani-1][ranj-1])
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=lossj*g_price[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=min(SI(ranj),lossj*dim_mach);
                        lossj=max(0.0,lossj*dim_mach-SI(ranj))/dim_mach;
                      }
                      EI(1,ranj)-=lossj*dim_mach;
                      EI(1,ranj)=max(EI(1,ranj),0.0);
                      lossj=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*lossj);
                      K(ranj)-=lossj*dim_mach;
                      K_cur(ranj)-=lossj*dim_mach;
                      gtemp[rant-1][rani-1][ranj-1]-=lossj;
                      g_c2[rant-1][rani-1][ranj-1]-=lossj;
                      g_c3[rant-1][rani-1][ranj-1]-=lossj;
                      loss-=lossj;
                      lossj=0;
                  }
              }
              else
              {
                  if(rant==t)
                  {
                      deltaCapitalStock(1,ranj)-=g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1];
                      if(SI(ranj)>0)
                      {
                        K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach-EI(1,ranj);
                      }
                      EI(1,ranj)=0;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
                  else
                  {
                      Loss_Capital(ranj)+=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      CapitalStock(1,ranj)-=min(CapitalStock(1,ranj),g_price[rant-1][rani-1][ranj-1]*gtemp[rant-1][rani-1][ranj-1]);
                      K(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      K_cur(ranj)-=gtemp[rant-1][rani-1][ranj-1]*dim_mach;
                      lossj-=gtemp[rant-1][rani-1][ranj-1];
                      loss-=gtemp[rant-1][rani-1][ranj-1];
                      g_c2[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      g_c3[rant-1][rani-1][ranj-1]-=gtemp[rant-1][rani-1][ranj-1];
                      gtemp[rant-1][rani-1][ranj-1]=0;
                  }
              }
          }
      }
    }

    if(K_loss.Sum()>0)
    {
      for (j=1; j<=N2; j++)
      {
        if(K_loss(j)>0)
        {
          K_loss_share+=1;
          Q2(j)=min(Q2(j),K_cur(j));
          ADJUSTEMISSENLAB();
        }
      }
    }
  }

  EN_DEM();
}

void ADJUSTEMISSENLAB(void)
{     
  nmachprod=ceil(Q2(j)/dim_mach);  
  nmp_temp=nmachprod;
  while (nmp_temp > 0)       
  {									         
    cmin=std::numeric_limits<double>::infinity();                                                  
    imin=0;                                                                 
    jmin=0;
    tmin=0;

    for (i=1; i<=N1; i++)     
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (g_c2[tt-1][i-1][j-1] > 0 && (w(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)) < cmin)
        {
          cmin=w(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
          imin=i;
          jmin=j;
          tmin=tt;
        }          
      }
    }

    if (nmachprod>0)
    {
      if (g_c2[tmin-1][imin-1][jmin-1] >= nmp_temp)  
      {
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*nmp_temp/nmachprod;
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*nmp_temp/nmachprod;
        A2e_ef2(j)+=A_ef(tmin,imin)*nmp_temp/nmachprod;
        g_c2[tmin-1][imin-1][jmin-1]-= nmp_temp;
        nmp_temp=0;  
      }
      else
      {
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        A2e_ef2(j)+=A_ef(tmin,imin)*g_c2[tmin-1][imin-1][jmin-1]/nmachprod;
        nmp_temp-=g_c2[tmin-1][imin-1][jmin-1];
        g_c2[tmin-1][imin-1][jmin-1]=0;        
      }
    }          
    else 
    {
      std::cerr << "\n\n ERROR: nmachprod = 0!!!" << endl;
      Errors << "\n nmachprod = 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
  } 

  if(A2e2(j)>0)
  {
    Ld2_control(j)=Q2(j)/A2e2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      Ld2_control(j)=0;
    }
  }

  if(Ld2_control(j)>Ld2(j))
  {
    A2e2(j)=0;
    A2e_en2(j)=0;
    A2e_ef2(j)=0;
    Ldtemp=Ld2(j);
    nmp_temp=0;
    for (i=1; i<=N1; i++)     
    {
      for (tt=t0; tt<=t; tt++)
      {
        if (g_c3[tt-1][i-1][j-1] > 0)
        {
          nmp_temp+=g_c3[tt-1][i-1][j-1];
        }          
      }
    }
    nmachprod=0;
    while (Ldtemp > 0 && nmachprod<nmp_temp)       
    {		
      cmin=std::numeric_limits<double>::infinity();                                                  
      imin=0;                                                                 
      jmin=0;
      tmin=0;

      for (i=1; i<=N1; i++)     
      {
        for (tt=t0; tt<=t; tt++)
        {
          if (g_c3[tt-1][i-1][j-1] > 0 && (w(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i)) < cmin)
          {
            cmin=w(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/A_en(tt,i)+t_CO2*A_ef(tt,i)/A_en(tt,i);
            imin=i;
            jmin=j;
            tmin=tt;
          }          
        }
      }

      if ((g_c3[tmin-1][imin-1][jmin-1]*dim_mach/A(tmin,imin)) > Ldtemp)  
      {
        nmachprod+=ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        A2e_ef2(j)+=A_ef(tmin,imin)*ceil(Ldtemp*A(tmin,imin)/dim_mach);
        g_c3[tmin-1][imin-1][jmin-1]-=ceil(Ldtemp*A(tmin,imin)/dim_mach);
        Ldtemp=0;             
      }
      else                      
      {
        Ldtemp-=(g_c3[tmin-1][imin-1][jmin-1]*dim_mach/A(tmin,imin));
        nmachprod+=g_c3[tmin-1][imin-1][jmin-1];
        A2e2(j)+=(1-shocks_labprod2(j))*A(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        A2e_en2(j)+=(1-shocks_eneff2(j))*A_en(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        A2e_ef2(j)+=A_ef(tmin,imin)*g_c3[tmin-1][imin-1][jmin-1];
        g_c3[tmin-1][imin-1][jmin-1]=0;        
      }             
    }
    A2e2(j)/=nmachprod;
    Q2(j)=A2e2(j)*Ld2(j);
    A2e_en2(j)/=nmachprod;
    A2e_ef2(j)/=nmachprod;
  }

  if(A2e_en2(j)>0)
  {
    D2_en(j)=Q2(j)/A2e_en2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e_en2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e_en2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      D2_en(j)=0;
    }
  }

  if(A2e_ef2(j)>0 && A2e_en2(j)>0)
  {
    Emiss2(j)= A2e_ef2(j)/A2e_en2(j)*Q2(j);
  }
  else
  {
    if(Q2(j)>0)
    {
      std::cerr << "\n\n ERROR: Q2>0 and A2e_ef2 or A2e_en2 == 0!!!" << endl;
      Errors << "\n  Q2>0 and A2e_ef2 or A2e_en2 == 0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      Emiss2(j)=0;
    }
  }
}

void CANCMACH(void)
{
  //Based on actual replacement investment carried out, machines are scrapped
  scrapmax=(scrap_age(j)+SI(j))/dim_mach;
  indsupl=int(supl(j));

  scrap_n=0;
  if(scrapmax>0)
  {
    //First scrap machines which are too old, then ones with high production cost
    //Only machines SCRAPPING marked can appear here, and it recorded them in
    //the order i and then tt.
    for (int marked=0; marked<marked_count[j-1]; marked++)
    {
      {
        i=marked_supplier[(j-1)*marked_capacity+marked];
        tt=marked_vintage[(j-1)*marked_capacity+marked];
        if (scrapmax > 0 && g_pb[tt-1][i-1][j-1] > 0 && age[tt-1][i-1][j-1]>(agemax))
        {
          if (g_pb[tt-1][i-1][j-1] >= scrapmax)	  
          {                                       
            gtemp[tt-1][i-1][j-1]-=scrapmax;
            g_pb[tt-1][i-1][j-1]-=scrapmax;
            scrap_n+=scrapmax*g_price[tt-1][i-1][j-1];
            scrapmax=0;                           
          }
          else        
          {						
            scrapmax-=g_pb[tt-1][i-1][j-1];	        
            gtemp[tt-1][i-1][j-1]-=g_pb[tt-1][i-1][j-1];
            scrap_n+=g_pb[tt-1][i-1][j-1]*g_price[tt-1][i-1][j-1];
            g_pb[tt-1][i-1][j-1]=0;
          }
        }
      }
    }
  }
                            
  while (scrapmax > 0)
  {
    cmax=0;

    for (int marked=0; marked<marked_count[j-1]; marked++)
    {
      i=marked_supplier[(j-1)*marked_capacity+marked];
      tt=marked_vintage[(j-1)*marked_capacity+marked];
      if (g_pb[tt-1][i-1][j-1] > 0 && C_pb[tt-1][i-1][j-1] > cmax)
      {
        ind_i=i;					           
        ind_tt=tt;
        cmax=C_pb[tt-1][i-1][j-1];
      }
    }

    if (g_pb[ind_tt-1][ind_i-1][j-1] >= scrapmax)	
    {                                             
      gtemp[ind_tt-1][ind_i-1][j-1]-=scrapmax;    
		  g_pb[ind_tt-1][ind_i-1][j-1]-=scrapmax;
      scrap_n+=scrapmax*g_price[ind_tt-1][ind_i-1][j-1];
		  scrapmax=0;								                 
    }
	  else                                     
	  {                                        
		  scrapmax-=g_pb[ind_tt-1][ind_i-1][j-1];
      scrap_n+=g_pb[ind_tt-1][ind_i-1][j-1]*g_price[ind_tt-1][ind_i-1][j-1];
 		  gtemp[ind_tt-1][ind_i-1][j-1]-=g_pb[ind_tt-1][ind_i-1][j-1];      
		  g_pb[ind_tt-1][ind_i-1][j-1]=0;
	  }
  }

  deltaCapitalStock(1,j)+=(SI_n(j)-scrap_n);
}

void PAY_LAB_INV(void)
{
  //C-firms pay wages
  for (j=1; j<=N2; j++)
	{
    sendingBank=BankingSupplier_2(j);
    Wages_2(j)=w(2)*(Ld2(j));
    if(Deposits_2(1,j)>=Wages_2(j))
    {
      Deposits_2(1,j)-=Wages_2(j);
      Wages+=Wages_2(j);
      Deposits(1,sendingBank)-=Wages_2(j);
      Outflows(sendingBank)+=Wages_2(j);
    }
    else
    {
      if(Wages_2(j)>0)
      {
        deviation=((Wages_2(j)-Deposits_2(1,j))/Wages_2(j));
        if(deviation<=tolerance)
        {
          Wages_2(j)=Deposits_2(1,j);
          Deposits_2(1,j)-=Wages_2(j);
          Wages+=Wages_2(j);
          Deposits(1,sendingBank)-=Wages_2(j);
          Outflows(sendingBank)+=Wages_2(j);
        }
        else
        {
          std::cerr << "\n\n ERROR: C-Firm " << j << " cannot pay wages in period " << t << endl;
          Errors << "\n C-Firm " << j << " cannot pay wages in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
    }

    //C-firms pay for investment
    indsupl=int(supl(j));
    Investment_2(j)=EI_n(j)+SI_n(j);
    receivingBank=BankingSupplier_1(indsupl);
    if(Deposits_2(1,j)>=Investment_2(j))
    {
      Deposits_2(1,j)-=Investment_2(j);
      Deposits(1,sendingBank)-=Investment_2(j);
      Outflows(sendingBank)+=Investment_2(j);

      Deposits_1(1,indsupl)+=Investment_2(j);
      Deposits(1,receivingBank)+=Investment_2(j);
      Inflows(receivingBank)+=Investment_2(j);
    }
    else
    {
      if(Investment_2(j)>0)
      {
        deviation=((Investment_2(j)-Deposits_2(1,j))/Investment_2(j));
        if(deviation<=tolerance)
        {
          S1(indsupl)-=Investment_2(j);
          S1(indsupl)+=Deposits_2(1,j);
          Investment_2(j)=Deposits_2(1,j);
          Deposits_2(1,j)-=Investment_2(j);
          Deposits(1,sendingBank)-=Investment_2(j);
          Outflows(sendingBank)+=Investment_2(j);
          Deposits_1(1,indsupl)+=Investment_2(j);
          Deposits(1,receivingBank)+=Investment_2(j);
          Inflows(receivingBank)+=Investment_2(j);
        }
        else
        {
          std::cerr << "\n\n ERROR: C-Firm " << j <<  " cannot pay for investment in period " << t << endl;
          Errors << "\n C-Firm " << j <<  " cannot pay for investment in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  
  //K-firms pay wages
  for (i=1; i<=N1; i++)
	{
    sendingBank=BankingSupplier_1(i);
    Wages_1(i)=w(2)*(Ld1(i)+Ld1rd(i));
    if(Deposits_1(1,i)>=Wages_1(i))
    {
      Deposits_1(1,i)-=Wages_1(i);
      Wages+=Wages_1(i);
      Deposits(1,sendingBank)-=Wages_1(i);
      Outflows(sendingBank)+=Wages_1(i);
    }
    else
    {
      if(Wages_1(i)>0)
      {
        deviation=((Wages_1(i)-Deposits_1(1,i))/Wages_1(i));
        if(deviation<=tolerance)
        {
          Wages_1(i)=Deposits_1(1,i);
          Deposits_1(1,i)-=Wages_1(i);
          Wages+=Wages_1(i);
          Deposits(1,sendingBank)-=Wages_1(i);
          Outflows(sendingBank)+=Wages_1(i);
        }
        else
        {
          Wages_1(i)=Deposits_1(1,i);
          Deposits_1(1,i)-=Wages_1(i);
          Wages+=Wages_1(i);
          Deposits(1,sendingBank)-=Wages_1(i);
          Outflows(sendingBank)+=Wages_1(i);
          exiting_1(i)=1;
          exiting_1_payments(i)=1;
        }
      }
    }
  }

  //Energy sector pays wages
  if(flag_WITCH_on==0)
  {
    if(Deposits_e(1)>=Wages_en)
    {
      Deposits_e(1)-=Wages_en;
      Wages+=Wages_en;
    }
    else
    {
      deviation=((Wages_en-Deposits_e(1))/Wages_en);
      if(deviation<=tolerance)
      {
        Wages_en=Deposits_e(1);
        Deposits_e(1)-=Wages_en;
        Wages+=Wages_en;
      }
      else
      {
        std::cerr << "\n\n ERROR: Energy sector cannot pay wages in period " << t << endl;
        Errors << "\n Energy sector cannot pay wages in period " << t << endl;
        //exit(EXIT_FAILURE);
        Wages_en=max(0.0,Deposits_e(1));
        Deposits_e(1)-=Wages_en;
        Wages+=Wages_en;
      }
    }
    
    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)-=Wages_en*DepositShare_e(i);
      Outflows(i)+=Wages_en*DepositShare_e(i);
      Deposits_eb(1,i)-=Wages_en*DepositShare_e(i);
    }
  }
  else
  {
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      if(wage_e_mult_tech[tech]>0)
      {
        if(deposits_e_mult_tech[tech]>=wage_e_mult_tech[tech])
        {
            deposits_e_mult_tech[tech]-=wage_e_mult_tech[tech];
            Wages+=wage_e_mult_tech[tech];
            Wages_en+=wage_e_mult_tech[tech];
        }else
        {
            wage_e_mult_tech[tech]=deposits_e_mult_tech[tech];
            deposits_e_mult_tech[tech]-=wage_e_mult_tech[tech];
            Wages+=wage_e_mult_tech[tech];
            Wages_en+=wage_e_mult_tech[tech];
        }
      }
      Deposits_e(1)-=wage_e_mult_tech[tech];   
    }

    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)-=Wages_en*DepositShare_e(i);
      Deposits_eb(1,i)-=Wages_en*DepositShare_e(i);
      Outflows(i)+=Wages_en*DepositShare_e(i);
    }
  }
  
  Deposits_h(1)+=(Wages+Benefits+govTranfers);

  //Households receive wages & Benefits
  for(i=1; i<=NB; i++)
  {
    Deposits(1,i)+=(Wages+Benefits+govTranfers)*DepositShare_h(i);
    Inflows(i)+=(Wages+Benefits+govTranfers)*DepositShare_h(i);
    Deposits_hb(1,i)+=(Wages+Benefits+govTranfers)*DepositShare_h(i);
  }
}

void COMPET2(void)
{
  //C-firms' market shares are updated using quasi-replicator dynamics. Firms with too low market share exit
  p2m=0;
  l2m=0;
  
  for (j=1; j<=N2; j++)
  {
    if (exiting_2(j) == 0)
    // Exclude exiting firms from the computation of average C-firm competitiveness.
    // They are removed from the average price and the average unsatisfied demand (see below).
    // They are de facto removed from ftot since f2 is set to zero when they exit
    // Their competitiveness is calculated but not included in the average C-firm competitiveness since their market share f2 is set to 0 when they exit
    {
      p2m+=p2(j)/(N2r-exiting_2.Sum());
      l2m+=l2(j)/(N2r-exiting_2.Sum());
    }
  }

	
	ftot=0;

	for (j=1; j<=N2; j++)
	{
    E2(j)=-pow(p2(j)/p2m,omega1)-pow(l2(j)/l2m,omega2);
    ftot(2)+=f2(2,j);
	}

	if (ftot(2)>0)
	{
		for (j=1; j<=N2; j++)
		{
			f2(2,j)/=ftot(2);
			Em2(1)+=E2(j)*f2(2,j);
		}
	}
	else
	{
		std::cerr << "\n\n ERROR: f2 = 0 for all firms in period " << t << endl;
    Errors << "\n f2 = 0 for all firms in period " << t << endl;
    exit(EXIT_FAILURE);
	}

	ftot=0;

	if (Em2(1)==0)
  {
		std::cerr << "\n\n ERROR: Em2(1)=0 in period " << t << endl;
    Errors << "\n Em2(1)=0 in period " << t << endl;
    exit(EXIT_FAILURE);
  }
  
	for (j=1; j<=N2; j++)
	{
    f2(1,j)=f2(2,j)*((2*omega3)/(1+exp((-chi)*((E2(j)-Em2(1))/Em2(1))))+(1-omega3));

    if (f2(1,j) <= (1/(N2r*500)))
    {
      f2(1,j)=0;
      f2(2,j)=0;
      f2(3,j)=0;
      if(exiting_2(j)==0 && exit_payments2(j)==0)
      {
        exit_marketshare2(j)=1;
      }
    }
		ftot(1)+=f2(1,j);
		ftot(2)+=f2(2,j);
		ftot(3)+=f2(3,j);
	}

	if (ftot(1)==0 || ftot(2)==0 || ftot(3)==0)
  {
		std::cerr << "\n\n ERROR: ftot=0 in period " << t << endl;
    Errors << "\n ftot=0 in period " << t << endl;
    exit(EXIT_FAILURE);
  }

	for (j=1; j<=N2; j++)
	{
		f2(1,j)/=ftot(1);
		f2(2,j)/=ftot(2);
		f2(3,j)/=ftot(3);
	}
}

void PROFIT(void)
{
  //Calculate K-firm profit
  for (i=1; i<=N1; i++)      
	{
    kpi+=Q1(i)/Q1.Sum()*p1(i);
    EnergyPayments_1(i)=c_en(1)*D1_en(i);
    Pi1(i)=S1(i)+InterestDeposits_1(i)-Wages_1(i)-EnergyPayments_1(i)-t_CO2*Emiss1(i); 
    sendingBank=BankingSupplier_1(i);

    //K-firms pay for energy; exit if they are unable
    if(Deposits_1(1,i)>=EnergyPayments_1(i))
    {
      Deposits_1(1,i)-=EnergyPayments_1(i);
      EnergyPayments+=EnergyPayments_1(i);
      Deposits(1,sendingBank)-=EnergyPayments_1(i);
      Outflows(sendingBank)+=EnergyPayments_1(i);
    }
    else if (Deposits_1(1,i)>=0)
    {
      deviation=((EnergyPayments_1(i)-Deposits_1(1,i))/EnergyPayments_1(i));
      EnergyPayments_1(i)=Deposits_1(1,i);
      Deposits_1(1,i)-=EnergyPayments_1(i);
      EnergyPayments+=EnergyPayments_1(i);
      Deposits(1,sendingBank)-=EnergyPayments_1(i);
      Outflows(sendingBank)+=EnergyPayments_1(i);
      if(deviation>tolerance)
      {
        exiting_1(i)=1;
        exiting_1_payments(i)=1;
      }
    }

    //K-firms pay CO2 tax
    if(exiting_1(i)==0)
    {
      Taxes_CO2_1(i)=t_CO2*Emiss1(i);
      if(Deposits_1(1,i)>=Taxes_CO2_1(i))
      {
        Taxes_CO2(1)+=Taxes_CO2_1(i);
        Deposits_1(1,i)-=Taxes_CO2_1(i);
        Deposits(1,sendingBank)-=Taxes_CO2_1(i);
        Outflows(sendingBank)+=Taxes_CO2_1(i);
      }
      else if(Deposits_1(1,i)>=0)
      {
        Taxes_CO2_1(i)=Deposits_1(1,i);
        Deposits_1(1,i)-=Taxes_CO2_1(i);
        Taxes_CO2(1)+=Taxes_CO2_1(i);
        Deposits(1,sendingBank)-=Taxes_CO2_1(i);
        Outflows(sendingBank)+=Taxes_CO2_1(i);
      }
    }
    else
    {
      Taxes_CO2_1(i)=0;
    }
    
    if (Pi1(i) > 0 && exiting_1(i)==0)
		{
			//If profit is positive, pay tax
      if(Deposits_1(1,i)>=aliq*Pi1(i))
      {
        Taxes_1(i)=aliq*Pi1(i);
        Deposits_1(1,i)-=Taxes_1(i);
        Deposits(1,sendingBank)-=Taxes_1(i);
        Outflows(sendingBank)+=Taxes_1(i);
        Taxes+=Taxes_1(i);
      }
      else if (Deposits_1(1,i)>=0)
      {
        Taxes_1(i)=Deposits_1(1,i);
        Deposits_1(1,i)-=Taxes_1(i);
        Deposits(1,sendingBank)-=Taxes_1(i);
        Outflows(sendingBank)+=Taxes_1(i);
        Taxes+=Taxes_1(i);
      }
      
      //Pay dividend
      if(Deposits_1(1,i)>=d1*Pi1(i))
      {
        Dividends_1(i)=d1*Pi1(i);
        Deposits_1(1,i)-=Dividends_1(i);
        Deposits(1,sendingBank)-=Dividends_1(i);
        Outflows(sendingBank)+=Dividends_1(i);
      }
      else if (Deposits_1(1,i)>=0)
      {
        Dividends_1(i)=Deposits_1(1,i);
        Deposits_1(1,i)-=Dividends_1(i);
        Deposits(1,sendingBank)-=Dividends_1(i);
        Outflows(sendingBank)+=Dividends_1(i);
      }

			Divtot_1+=Dividends_1(i);	
      Dividends(1)+=Dividends_1(i);
		}
		Pitot1+=Pi1(i);			

    if(Deposits_1(1,i)<0)
    {
      std::cerr << "\n\n ERROR: K-firm " << i << " has negative deposits in period " << t << endl;
      Errors << "\n K-firm " << i << " has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }

    //K-Firms with 0 customers exit
    if(nclient(i) < 1 && exiting_1(i)==0)
    {
      exiting_1(i)=1;
    }

    //Calculate mean deposits of non-failing K-firms for entry below
		if (exiting_1(i)==0)
		{
			mD1+=Deposits_1(1,i);
			ns1++;
		}

    //Prepare failing K-firms for exit
    if(exiting_1(i)==1)
    {
      Deposits_recovered_1+=Deposits_1(1,i);
      sendingBank=BankingSupplier_1(i);
      baddebt_1(i)=-Deposits_1(1,i);
      Deposits(1,sendingBank)-=Deposits_1(1,i);
      Outflows(sendingBank)+=Deposits_1(1,i);
      Deposits_1(1,i)=0;
    }
  }

  //Households receive K-firm dividends and pay tax
  Taxes_h=aliqw*Wages;
  Deposits_h(1)+=(Deposits_recovered_1+Divtot_1-Taxes_h);
  Taxes+=Taxes_h;

  for(i=1; i<=NB; i++)
  {
    Deposits(1,i)+=(Deposits_recovered_1+Divtot_1)*DepositShare_h(i);
    Inflows(i)+=(Deposits_recovered_1+Divtot_1)*DepositShare_h(i);
    Deposits_hb(1,i)+=(Deposits_recovered_1+Divtot_1)*DepositShare_h(i);
    Deposits(1,i)-=Taxes_h*DepositShare_h(i);
    Outflows(i)+=Taxes_h*DepositShare_h(i);
    Deposits_hb(1,i)-=Taxes_h*DepositShare_h(i);
  }

  //Consumption demand is determined; Check if consumption can be financed
  Cons=(a1*(Wages+Benefits-Taxes_h)+a2*(InterestDeposits_h+Dividends(2))+a3*Deposits_h(2))+a4*govTranfers;
  
  if(Deposits_h(1)<Cons)
  {
    Cons=Deposits_h(1);
  }

  for (j=1; j<=N2; j++)	
	{
    cpi(1)+=p2(j)*f2(1,j);					
  }
	
  if (cpi(1) < 0.01)									
	{
		std::cerr << "\n\n ERROR: CPI < 0.01 in period " << t << endl;
    Errors << "\n CPI < 0.01 in period " << t << endl;
    exit(EXIT_FAILURE);
	}

  //Consumption takes place
  ALLOC();

  if(Deposits_h(1)<0)
  {
    std::cerr << "\n\n ERROR: Household deposits negative after consumption in period " << t << endl;
    Errors << "\n  Household deposits negative after consumption in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  //Calculate C-firm profits
	for (j=1; j<=N2; j++)
	{
		dN(j)=N(1,j)-N(2,j);
    Inventories(1,j)=N(1,j)*p2(j);
    dNm(j)=Inventories(1,j)-Inventories(2,j);
    dNtot+=dN(j);
    dNmtot+=dNm(j);
    EnergyPayments_2(j)=c_en(1)*D2_en(j);
    mol(j)=S2(1,j)-Wages_2(j)-EnergyPayments_2(j);
    if(exiting_2(j)==0)
    {
      LoanInterest_2(j)=r_deb_h(j)*Loans_2(1,j);
    }
    else
    {
      LoanInterest_2(j)=0;
    }
    
    //NB: Profit here also includes changes in nominal value of tangible assets
    Pi2(j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)-Investment_2(j)-Wages_2(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);

    //C-firm  pays energy; if unable to do so it exits
    if(Deposits_2(1,j)>=EnergyPayments_2(j))
    {
      sendingBank=BankingSupplier_2(j);
      Deposits_2(1,j)-=EnergyPayments_2(j);
      EnergyPayments+=EnergyPayments_2(j);
      Deposits(1,sendingBank)-=EnergyPayments_2(j);
      Outflows(sendingBank)+=EnergyPayments_2(j);
    }
    else if(Deposits_2(1,j)>=0)
    {
      sendingBank=BankingSupplier_2(j);
      EnergyPayments_2(j)=Deposits_2(1,j);
      Deposits_2(1,j)=0;
      Deposits(1,sendingBank)-=EnergyPayments_2(j);
      Outflows(sendingBank)+=EnergyPayments_2(j);
      EnergyPayments+=EnergyPayments_2(j);
      Taxes_2(j)=0;
      LoanInterest_2(j)=0;
      DebtRemittances2(j)=0;
      exiting_2(j)=1;
      exit_payments2(j)=1;
      Pi2(j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)-Investment_2(j)-Wages_2(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);
    }

    //Non-exiting firms make principal & interest payments on loans; if unable to do so they exit
    if(exiting_2(j)==0)
    {
      DebtRemittances2(j)=repayment_share*Loans_2(1,j);
      if(Deposits_2(1,j)>=(LoanInterest_2(j)+DebtRemittances2(j)))
      {
        receivingBank=BankingSupplier_2(j);
        Deposits_2(1,j)-=(LoanInterest_2(j)+DebtRemittances2(j));
        Deposits(1,receivingBank)-=(LoanInterest_2(j)+DebtRemittances2(j));
        DebtService_2(1,j)=LoanInterest_2(j)+DebtRemittances2(j);
        Loans_2(1,j)-=DebtRemittances2(j);
        Loans_b(1,receivingBank)-=DebtRemittances2(j);
        LoanInterest(receivingBank)+=LoanInterest_2(j);
      }
      else
      {
        LoanInterest_2(j)=0;
        DebtRemittances2(j)=0;
        Taxes_2(j)=0;
        exiting_2(j)=1;
        exit_payments2(j)=1;
        Pi2(j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)-Investment_2(j)-Wages_2(j)-EnergyPayments_2(j)-LoanInterest_2(j)-t_CO2*Emiss2(j);
      }
    }

    //C-firms pay CO2 tax
    if(exiting_2(j)==0)
    {
      Taxes_CO2_2(j)=t_CO2*Emiss2(j);
      sendingBank=BankingSupplier_2(j);
      if(Deposits_2(1,j)>=Taxes_CO2_2(j))
      {
        Taxes_CO2(1)+=Taxes_CO2_2(j);
        Deposits_2(1,j)-=Taxes_CO2_2(j);
        Deposits(1,sendingBank)-=Taxes_CO2_2(j);
        Outflows(sendingBank)+=Taxes_CO2_2(j);
      }
      else if(Deposits_2(1,j)>=0)
      {
        Taxes_CO2_2(j)=Deposits_2(1,j);
        Deposits_2(1,j)-=Taxes_CO2_2(j);
        Taxes_CO2(1)+=Taxes_CO2_2(j);
        Deposits(1,sendingBank)-=Taxes_CO2_2(j);
        Outflows(sendingBank)+=Taxes_CO2_2(j);
        Pi2(j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)-Investment_2(j)-Wages_2(j)-EnergyPayments_2(j)-LoanInterest_2(j)-Taxes_CO2_2(j);
      }
    }
    else
    {
      Taxes_CO2_2(j)=0;
      Pi2(j)=S2(1,j)+InterestDeposits_2(j)+dNm(j)+deltaCapitalStock(1,j)-Investment_2(j)-Wages_2(j)-EnergyPayments_2(j)-LoanInterest_2(j)-Taxes_CO2_2(j);
    }

    //C-firm pays tax
    if(Pi2(j)>0 && exiting_2(j)==0)
    {
      Taxes_2(j)=aliq*Pi2(j);
      sendingBank=BankingSupplier_2(j);
      if(Deposits_2(1,j)>=Taxes_2(j))
      {
        Deposits_2(1,j)-=Taxes_2(j);
        Taxes+=Taxes_2(j);
        Deposits(1,sendingBank)-=Taxes_2(j);
        Outflows(sendingBank)+=Taxes_2(j);
      }
      else if(Deposits_2(1,j)>=0)
      {
        Taxes_2(j)=Deposits_2(1,j);
        Deposits_2(1,j)=0;
        Taxes+=Taxes_2(j);
        Deposits(1,sendingBank)-=Taxes_2(j);
        Outflows(sendingBank)+=Taxes_2(j);
      }
    }
  
    //Non-exiting firms pay dividends
		if (Pi2(j) > 0 && exiting_2(j)==0) 
		{
			sendingBank=BankingSupplier_2(j);
      if(Deposits_2(1,j)>=d2*Pi2(j))
      {
        Dividends_2(j)=d2*Pi2(j);
        Deposits_2(1,j)-=Dividends_2(j);
        Deposits(1,sendingBank)-=Dividends_2(j);
        Outflows(sendingBank)+=Dividends_2(j);
      }
      else if(Deposits_2(1,j)>=0)
      {
        Dividends_2(j)=Deposits_2(1,j);
        Deposits_2(1,j)-=Dividends_2(j);
        Deposits(1,sendingBank)-=Dividends_2(j);
        Outflows(sendingBank)+=Dividends_2(j);
      }
			Divtot_2+=Dividends_2(j);
      Dividends(1)+=Dividends_2(j);
		}

		Pitot2+=Pi2(j);

    if(Deposits_2(1,j)<0)
    {
      std::cerr << "\n\n ERROR: C-firm " << j << " has negative deposits in period " << t << endl;
      Errors << "\n C-firm " << j << " has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }
	} 

    
  if(flag_WITCH_on==0)
  {
    //Energy sector receives revenue
    Deposits_e(1)+=EnergyPayments;

    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)+=EnergyPayments*DepositShare_e(i);
      Deposits_eb(1,i)+=EnergyPayments*DepositShare_e(i);
      Inflows(i)+=EnergyPayments*DepositShare_e(i);
    }

    //Profit of energy sector
    ProfitEnergy=EnergyPayments+InterestDeposits_e-Wages_en+(CapitalStock_e(1)-CapitalStock_e(2))-t_CO2_en*Emiss_en-FuelCost;

    if(Deposits_e(1)>=FuelCost)
    {
      Deposits_e(1)-=FuelCost;
      for(i=1; i<=NB; i++)
      {
        Deposits(1,i)-=FuelCost*DepositShare_e(i);
        Outflows(i)+=FuelCost*DepositShare_e(i);
        Deposits_eb(1,i)-=FuelCost*DepositShare_e(i);
      }
    }
    else
    {
      std::cerr << "\n\n Energy sector cannot pay for fuel in period " << t << endl;
      Errors << "\n Energy sector cannot pay for fuel in period " << t << endl;
      exit(EXIT_FAILURE);
    }
    
    //Energy sector pays C02 tax
    if(Deposits_e(1)>=t_CO2_en*Emiss_en)
    {
      Taxes_CO2_e=t_CO2_en*Emiss_en;
      Taxes_CO2(1)+=Taxes_CO2_e;
    }
    else if(Deposits_e(1)>=0)
    {
      Taxes_CO2_e=Deposits_e(1);
      Taxes_CO2(1)+=Taxes_CO2_e;
    }
    else
    {
      std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
      Errors << "\n Energy sector has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }

    Deposits_e(1)-=Taxes_CO2_e;
    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)-=Taxes_CO2_e*DepositShare_e(i);
      Outflows(i)+=Taxes_CO2_e*DepositShare_e(i);
      Deposits_eb(1,i)-=Taxes_CO2_e*DepositShare_e(i);
    }
    
    //Energy sector dividends
    if(ProfitEnergy>0)
    {
      if(Deposits_e(1)>=de*ProfitEnergy)
      {
        Dividends_e=de*ProfitEnergy;
      }
      else if(Deposits_e(1)>=0)
      {
        Dividends_e=Deposits_e(1);
      }
      else
      {
        std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
        Errors << "\n Energy sector has negative deposits in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      Dividends(1)+=Dividends_e;
      Deposits_e(1)-=Dividends_e;
      for(i=1; i<=NB; i++)
      {
        Deposits(1,i)-=Dividends_e*DepositShare_e(i);
        Outflows(i)+=Dividends_e*DepositShare_e(i);
        Deposits_eb(1,i)-=Dividends_e*DepositShare_e(i);
      }
    }
    else
    {
      Dividends_e=0;
    }

    //Fossil fuel agent receives fuel payment
    Deposits_fuel(1)+=FuelCost;
    Deposits_fuel_cb(1)+=FuelCost;
    TransferFuel=d_f*Deposits_fuel(1);
    Deposits_fuel(1)-=TransferFuel;
    Deposits_fuel_cb(1)-=TransferFuel;

    //Households receive energy sector, C-firm dividend and small transfer from fossil fuel agent
    Deposits_h(1)+=(Divtot_2+Dividends_e+TransferFuel);

    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)+=(Divtot_2+Dividends_e+TransferFuel)*DepositShare_h(i);
      Deposits_hb(1,i)+=(Divtot_2+Dividends_e+TransferFuel)*DepositShare_h(i);
      Inflows(i)+=(Divtot_2+Dividends_e+TransferFuel)*DepositShare_h(i);
    }
  }
  else
  {
    for(i=0; i<ene_tecs.size(); i++)
    {
        tech=ene_tecs[i];
        Rev_en_mult_tech[tech]=(d_en_mult_tech[tech]/D_en_TOT)*EnergyPayments;
        deposits_e_mult_tech[tech]+=(d_en_mult_tech[tech]/D_en_TOT)*EnergyPayments;
        Deposits_e(1)+=(d_en_mult_tech[tech]/D_en_TOT)*EnergyPayments;
        loan_interest_e_mult_tech[tech]=r_deb(1)*loans_e_mult_tech[tech];
        profitEnergy_e_mult_tech[tech]=Rev_en_mult_tech[tech]+interestsDeposits_e_mult_tech[tech]-wage_e_mult_tech[tech]+delta_capitalStock_e_mult_tech[tech]-loan_interest_e_mult_tech[tech];
        debtRemittances_e_mult_tech[tech]=repayment_share*loans_e_mult_tech[tech];
        
        if(deposits_e_mult_tech[tech]>=(loan_interest_e_mult_tech[tech]+debtRemittances_e_mult_tech[tech]))
        {
           deposits_e_mult_tech[tech]-=(loan_interest_e_mult_tech[tech]+debtRemittances_e_mult_tech[tech]);
           DebtService_e+=(loan_interest_e_mult_tech[tech]+debtRemittances_e_mult_tech[tech]);
           loans_e_mult_tech[tech]-=debtRemittances_e_mult_tech[tech];
           Loans_e(1)-=debtRemittances_e_mult_tech[tech];
           DebtRemittances_e+=debtRemittances_e_mult_tech[tech];
           Loan_interest_e+=loan_interest_e_mult_tech[tech];
           shareDefaultedDebt_e_mult_tech[tech]=0;
           badDebt_e_mult_tech[tech]=0;
           Loans_preDefault_e+=loans_e_mult_tech[tech];
        }
        else
        {
          loan_interest_e_mult_tech[tech]=0;
          debtRemittances_e_mult_tech[tech]=0;
          taxes_CO2_e_mult_tech[tech]=0;
          exiting_e_mult_tech[tech]=1;
          badDebt_e_mult_tech[tech]=max(0.0, loans_e_mult_tech[tech]-deposits_e_mult_tech[tech]);
          Deposits_e(1)-=deposits_e_mult_tech[tech];
          DefaultedDeposits_e+=deposits_e_mult_tech[tech];
          deposits_e_mult_tech[tech]=0;
          BadDebt_e+=badDebt_e_mult_tech[tech];
          DeafaultedDebtRecovered_e+=max(0.0, loans_e_mult_tech[tech]-badDebt_e_mult_tech[tech]);
          DebtWrittenOff_e+=loans_e_mult_tech[tech];
          Loans_e(1)-=loans_e_mult_tech[tech];
          shareDefaultedDebt_e_mult_tech[tech]=badDebt_e_mult_tech[tech]/loans_e_mult_tech[tech];
          Loans_preDefault_e+=loans_e_mult_tech[tech];
          loans_e_mult_tech[tech]=0;
          profitEnergy_e_mult_tech[tech]=Rev_en_mult_tech[tech]+interestsDeposits_e_mult_tech[tech]-wage_e_mult_tech[tech]+delta_capitalStock_e_mult_tech[tech]-loan_interest_e_mult_tech[tech];   
        }
    }

    Deposits_e(1)-=DebtService_e;
    ShareDefaultedDebt_e=BadDebt_e/Loans_preDefault_e;
  
    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)+=EnergyPayments*DepositShare_e(i);
      Deposits_eb(1,i)+=EnergyPayments*DepositShare_e(i);
      Inflows(i)+=EnergyPayments*DepositShare_e(i);

      Deposits(1,i)-=DebtService_e*DepositShare_e(i);
      Deposits_eb(1,i)-=DebtService_e*DepositShare_e(i);

      Deposits(1,i)-=DefaultedDeposits_e*DepositShare_e(i);
      Deposits_eb(1,i)-=DefaultedDeposits_e*DepositShare_e(i);

      Loans_b(1,i)-=DepositShare_e(i)*DebtRemittances_e;
      LoanInterest(i)+=DepositShare_e(i)*Loan_interest_e;
      
      baddebt_b(i)+=DepositShare_e(i)*BadDebt_e;
      Loans_b(1,i)-=DepositShare_e(i)*DebtWrittenOff_e;
    }
  
    //Energy sector pays C02 tax
    Taxes_CO2_e=0;
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      if(exiting_e_mult_tech[tech]==0)
      {
        if(deposits_e_mult_tech[tech]>=t_CO2_en*emiss_e_mult_tech[tech])
        {
          taxes_CO2_e_mult_tech[tech]=max(0.0, t_CO2_en*emiss_e_mult_tech[tech]);
          Taxes_CO2(1)+=taxes_CO2_e_mult_tech[tech];
          unpaidCtax_e_mult_tech[tech]=0;
        }
        else if(deposits_e_mult_tech[tech]>=0)
        {
          taxes_CO2_e_mult_tech[tech]=max(0.0, t_CO2_en*emiss_e_mult_tech[tech]);
          unpaidCtax_e_mult_tech[tech]=1-(deposits_e_mult_tech[tech]/taxes_CO2_e_mult_tech[tech]);
          taxes_CO2_e_mult_tech[tech]=deposits_e_mult_tech[tech];
          Taxes_CO2(1)+=taxes_CO2_e_mult_tech[tech];
        }
        else
        {
          std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
          Errors << "\n Energy sector has negative deposits in period " << t << endl;
          exit(EXIT_FAILURE);
        }
      }
      else
      {
        taxes_CO2_e_mult_tech[tech]=0;  
      }
      emiss_e_mult_tech[tech]=0;
      deposits_e_mult_tech[tech]-=taxes_CO2_e_mult_tech[tech];
      Deposits_e(1)-=taxes_CO2_e_mult_tech[tech];
      Taxes_CO2_e+=taxes_CO2_e_mult_tech[tech];
    }
  
    for(i=1; i<=NB; i++)
    {
      for(j=0; j<ene_tecs.size(); j++)
      {
          tech=ene_tecs[j];  
          Deposits(1,i)-=taxes_CO2_e_mult_tech[tech]*DepositShare_e(i);
          Outflows(i)+=taxes_CO2_e_mult_tech[tech]*DepositShare_e(i);
          Deposits_eb(1,i)-=taxes_CO2_e_mult_tech[tech]*DepositShare_e(i);
      }
    }
  
    //Energy sector dividends
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      if(exiting_e_mult_tech[tech]==0)
      {
        effective_dividendRate_e_mult_tech[tech]=min(0.99, de+de_premium[tech]);
        if(deposits_e_mult_tech[tech]>=effective_dividendRate_e_mult_tech[tech]*max(0.0, profitEnergy_e_mult_tech[tech])+inv_redistr_profit_e_multi_tech[tech])
        {
          dividends_e_mult_tech[tech]=effective_dividendRate_e_mult_tech[tech]*max(0.0, profitEnergy_e_mult_tech[tech])+inv_redistr_profit_e_multi_tech[tech];
        }
        else if(deposits_e_mult_tech[tech]>=0)
        {
          dividends_e_mult_tech[tech]=deposits_e_mult_tech[tech];
        }
        else
        {
          std::cerr << "\n\n Energy sector has negative deposits in period " << t << endl;
          Errors << "\n Energy sector has negative deposits in period " << t << endl;
          exit(EXIT_FAILURE);
        }
        Dividends(1)+=dividends_e_mult_tech[tech];
        deposits_e_mult_tech[tech]-=dividends_e_mult_tech[tech];
      }
      else
      {
        dividends_e_mult_tech[tech]=0; 
      }
      Dividends_e+=dividends_e_mult_tech[tech];
    }

    Deposits_e(1)-=Dividends_e;
    for(j=1; j<=NB; j++)
    {
        Deposits(1,j)-=Dividends_e*DepositShare_e(j);
        Outflows(j)+=Dividends_e*DepositShare_e(j);
        Deposits_eb(1,j)-=Dividends_e*DepositShare_e(j);
    }
    //Households receive energy sector & C-firm dividend
    Deposits_h(1)+=(Divtot_2+Dividends_e);

    for(i=1; i<=NB; i++)
    {
      Deposits(1,i)+=(Divtot_2+Dividends_e)*DepositShare_h(i);
      Deposits_hb(1,i)+=(Divtot_2+Dividends_e)*DepositShare_h(i);
      Inflows(i)+=(Divtot_2+Dividends_e)*DepositShare_h(i);
    }
  }

  //C-firms which exit due to negative equity, inability to make payments or low market share are prepared for exit
  for(j=1; j<=N2; j++)
  {
    NW_2(1,j)=CapitalStock(1,j)+deltaCapitalStock(1,j)+Inventories(1,j)+Deposits_2(1,j)-Loans_2(1,j);
    if(NW_2(1,j)<0 && exit_payments2(j)==0 && exiting_2(j)==0 && exit_marketshare2(j)==0)
    {
      exit_equity2(j)=1;
      exiting_2(j)=1;
    }
    
    if(exit_payments2(j)==1 || exit_equity2(j)==1)
    {
      sendingBank=BankingSupplier_2(j); 
      baddebt_2(j)=Loans_2(1,j)-Deposits_2(1,j);
      if(Loans_2(1,j)>Deposits_2(1,j))
      {
        baddebt_b(sendingBank)+=(Loans_2(1,j)-Deposits_2(1,j));
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits(1,sendingBank)-=Deposits_2(1,j);
      }
      else
      {
        Deposits_recovered_2+=Deposits_2(1,j)-Loans_2(1,j);
        Outflows(sendingBank)+=Deposits_2(1,j)-Loans_2(1,j);
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits(1,sendingBank)-=Deposits_2(1,j);
      }
      Deposits_2(1,j)=0;
      Loans_2(1,j)=0;
    }

    if(exit_marketshare2(j)==1 && exit_payments2(j)==0 && exit_equity2(j)==0)
    {
      exiting_2(j)=1;
      sendingBank=BankingSupplier_2(j);
      baddebt_2(j)=Loans_2(1,j)-Deposits_2(1,j);
      if(Loans_2(1,j)>Deposits_2(1,j))
      {
        baddebt_b(sendingBank)+=(Loans_2(1,j)-Deposits_2(1,j));
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits(1,sendingBank)-=Deposits_2(1,j);
      }
      else
      {
        Deposits_recovered_2+=Deposits_2(1,j)-Loans_2(1,j);
        Outflows(sendingBank)+=Deposits_2(1,j)-Loans_2(1,j);
        Loans_b(1,sendingBank)=max(0.0,Loans_b(1,sendingBank)-Loans_2(1,j));
        Deposits(1,sendingBank)-=Deposits_2(1,j);
      }
      Deposits_2(1,j)=0;
      Loans_2(1,j)=0;
    }

    if(exiting_2(j)==1 && exit_payments2(j)==0 && exit_equity2(j)==0 && exit_marketshare2(j)==0)
    {
      exit_payments2(j)=1;
    }

	}

  Deposits_h(1)+=Deposits_recovered_2;

  for(i=1; i<=NB; i++)
  {
    Deposits(1,i)+=Deposits_recovered_2*DepositShare_h(i);
    Deposits_hb(1,i)+=Deposits_recovered_2*DepositShare_h(i);
    Inflows(i)+=Deposits_recovered_2*DepositShare_h(i);
  }
}

void ALLOC(void)
{
  n=1;
	ftot=0;
  Utilisation=Q2.Sum()/K.Sum();
	Cres=Cons;
  cpi_temp=cpi(1);

  if(flag_outputshocks==1 || flag_outputshocks==3)
  {
    for (j=1; j<=N2; j++)
	  {
      C_loss(j)=shocks_output2(j)*Q2(j);
      Q2(j)=(1-shocks_output2(j))*Q2(j);
    }
  }

  if(flag_outputshocks==2 || flag_outputshocks==4)
  {
    lossdouble=Q2.Sum()*shocks_output2(1);
    while(loss>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      if(Q2(ranj)>=lossdouble)
      {
        Q2(ranj)-=lossdouble;
        C_loss(ranj)+=lossdouble;
        lossdouble=0;
      }
      else
      {
        lossdouble-=Q2(ranj);
        C_loss(ranj)+=Q2(ranj);
        Q2(ranj)=0;
      }
    }
  }

  if(flag_inventshocks==1)
  {
    for (j=1; j<=N2; j++)
	  {
      if(N(2,j)>0)
      {
        ptemp=Inventories(2,j)/N(2,j);
        Loss_Inventories(j)+=shocks_invent(j)*N(2,j)*ptemp;
        N(2,j)=(1-shocks_invent(j))*N(2,j);
        Inventories(2,j)-=Loss_Inventories(j);
      }
    }
  }
  
  if(flag_inventshocks==2)
  {
    lossdouble=N.Row(2).Sum()*shocks_invent(1);
    while(lossdouble>0)
    {
      ranj=int(ran1(p_seed)*N1*N2)%N2+1;
      if(N(2,ranj)>0)
      {
        ptemp=Inventories(2,ranj)/N(2,ranj);
        if(N(2,ranj)>=lossdouble)
        {
          ptemp=Inventories(2,ranj)/N(2,ranj);
          Loss_Inventories(ranj)+=lossdouble*ptemp;
          Inventories(2,ranj)-=Loss_Inventories(ranj);
          N(2,ranj)-=lossdouble;
          lossdouble=0;
        }
        else
        {
          Loss_Inventories(ranj)+=Inventories(2,ranj);
          Inventories(2,ranj)=0;
          lossdouble-=N(2,ranj);
          N(2,ranj)=0;
        }
      }
    }
  }

	for (j=1; j<=N2; j++)
	{
		if(C_loss(j)>0)
    {
      C_loss_share+=1;
    }
    f_temp2(j)=f2(1,j);
		ftot(1)+=f_temp2(j);
		Q2temp(j)=Q2(j)+N(2,j);
	}

  //Consumption demand is distributed among C-firms based on market shares
	while (Cres >= 1 && ftot(1) > 0)		
	{
    Cresb=Cres;
		for (j=1; j<=N2; j++)
		{
			if (f_temp2(j) > 0)
			{
				D_temp2(j)=Cres/cpi_temp*f_temp2(j);

				if (n==1)
				{
					D2(1,j)+=D_temp2(j);
				}

				if (D_temp2(j) <= Q2temp(j))	
				{
					if (n > 1)
          {
						D2(1,j)+=D_temp2(j);
          }
          S2(1,j)+=p2(j)*D_temp2(j);
          Cresb-=D_temp2(j)*p2(j);
          if(n==1)
          {
						l2(j)=1;
          }		
          Q2temp(j)-=D_temp2(j);		
				}								
				else						
				{
					if (n > 1)
          {
						D2(1,j)+=Q2temp(j);
          }
          S2(1,j)+=p2(j)*Q2temp(j);
          Cresb-=Q2temp(j)*p2(j);
          f_temp2(j)=0;
          if(n==1)
          {
            l2(j)=1+(D_temp2(j)-Q2temp(j));
          }
          Q2temp(j)=0;
				}
			}
		}
		ftot(1)=f_temp2.Sum();
	  f_temp2/=ftot(1);
		Cres=Cresb;
    cpi_temp=0;
    for (j=1; j<=N2; j++)	
    {
      cpi_temp+=p2(j)*f_temp2(j);						
    }
		n++;
	}

  //Nominal consumption is calculated
  Consumption=S2.Row(1).Sum();
  //This is done to ensure that household deposits do not become negative due to consumption (may happen due to rounding issues when liquidity constraint is binding)
  while(Consumption>Cons)
  {
    for (j=1; j<=N2; j++)
	  {
      if(S2(1,j)>(S2(1,j)/S2.Row(1).Sum()*(Consumption-Cons)))
      {
        S2(1,j)-=(S2(1,j)/S2.Row(1).Sum()*(Consumption-Cons));
      }
    }
    Consumption=S2.Row(1).Sum();
  }

  //Real consumption is calculated
  for(j=1; j<=N2; j++)
  {
    Consumption_r+=S2(1,j)/p2(j);
  }
 
  //Households pay for consumption
  Deposits_h(1)-=Consumption;

  for(i=1; i<=NB; i++)
  {
    Deposits(1,i)-=Consumption*DepositShare_h(i);
    Outflows(i)+=Consumption*DepositShare_h(i);
    Deposits_hb(1,i)-=Consumption*DepositShare_h(i);
  }

  //C-firms receive revenue & we re-compute the CPI
  cpi(1)=S2.Row(1).Sum()/Consumption_r;
  for (j=1; j<=N2; j++)
	{
    receivingBank=BankingSupplier_2(j);
    Deposits_2(1,j)+=S2(1,j);
    Deposits(1,receivingBank)+=S2(1,j);
    Inflows(receivingBank)+=S2(1,j);
    N(1,j)=N(2,j)+flag_inventories*(Q2(j)-S2(1,j)/p2(j));
    if(N(1,j)<0)
    {
      if(fabs(N(1,j))/Q2.Sum()<tolerance)
      {
        N(1,j)=0;
      }
      else
      {
        std::cerr << "\n\n ERROR: Inventories of C-firm " << j << " are negative in period " << t << endl;
        Errors << "\n Inventories of C-firm " << j << " are negative in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  S1_temp.Row(1)=S1;
  S2_temp.Row(1)=S2.Row(1);
}

void ENTRYEXIT(void)			
{
  //Save sales as S1 & S2 are reset during entry
  
  Sales1=S1;
  Sales2=S2.Row(1);
  
  for (i=1; i<=N1; i++)
	{
    if(Deposits_1(1,i)<0)
    {
      std::cerr << "\n\n ERROR: K-firm " << i << " has negative deposits in period " << t << endl;
      Errors << "\n  K-firm " << i << " has negative deposits in period " << t << endl;
      exit(EXIT_FAILURE);
    }
	}

  if (ns1 > 0)
  {
		mD1/=ns1;
  }
	else
  { 
    mD1=Deposits_1.Row(2).Sum()/N1r;
  }

  //Exiting K-firms lose customers; K-firm to be copied is chosen
  for (i=1; i<=N1; i++)
  {
    if (exiting_1(i)==1)
    {
      flag=0;
      for (j=1; j<=N2; j++)		
      {									     
        Match(j,i)=0;		
        if(supl(j)==i)
        {
          supl(j)=0;
        }			
      }
      if(exiting_1.Sum()<N1r)
      {
        while (flag == 0)			
        {
          rni=int(ran1(p_seed)*N1*N2)%N1+1;
          if (exiting_1(rni)==0)
          {
            ee1(i)=rni;
            flag=1;
          }
        }
      }
      else
      {
        ee1(i)=i;
      }
    }
  }

  for (i=1; i<=N1; i++)
  {
    if (exiting_1(i)==1)
    {
      //For entering K-firms, most variables are copied from a surviving K-firm as in the original DSK
      iii=int(ee1(i));
      f1(1,i)=0;
      f1(2,i)=0;
      Ld1rd(i)=0;
      shocks_labprod1(i)=shocks_labprod1(iii);
      shocks_eneff1(i)=shocks_eneff1(iii);
      p1(i)=p1(iii);
      A1(i)=A1(iii);  
      A1p(i)=A1p(iii);
      A1_ef(i)=A1_ef(iii);
      A1p_ef(i)=A1p_ef(iii);
      A1_en(i)=A1_en(iii);
      A1p_en(i)=A1p_en(iii);

      //Entering K-firms' deposits, however, are received as a transfer from households
      receivingBank=BankingSupplier_1(i);
      multip_entry=ran1(p_seed);
      multip_entry=w1inf+multip_entry*(w1sup-w1inf);
      if(Deposits_h(1)>=multip_entry*mD1)
      {
        //If households have sufficient deposits, transfer is equal to mean deposits of surviving firms
        injection=multip_entry*mD1;
        Deposits_h(1)-=injection;
        FirmTransfers+=injection;
        FirmTransfers_1+=injection;
        Injection_1(i)=injection;
        for(j=1; j<=NB; j++)
        {
          Deposits_hb(1,j)-=DepositShare_h(j)*injection;
          Deposits(1,j)-=DepositShare_h(j)*injection;
          Outflows(j)+=DepositShare_h(j)*injection;
        }
        Deposits_1(1,i)+=injection;
        Deposits(1,receivingBank)+=injection;
        Inflows(receivingBank)+=injection;
      }
      else
      {
        //If households cannot finance K-firm entry, this is either financed by gov. (flag_entry=0) or banks simply create the deposits (without corresponding loan) and book this as a loss (flag_entry=1)
        std::cerr << "\n\n Households cannot finance K-firm entry in period " << t << endl;
        Errors << "\n Households cannot finance K-firm entry in period " << t << endl;
        injection2=0;
        if(Deposits_h(1)>0)
        {
          injection2=Deposits_h(1);
          Deposits_h(1)-=injection2;
          FirmTransfers+=injection2;
          FirmTransfers_1+=injection2;
          Injection_1(i)=injection2;
          for(j=1; j<=NB; j++)
          {
            Deposits(1,j)-=Deposits_hb(1,j);
            Outflows(j)+=Deposits_hb(1,j);
            Deposits_hb(1,j)=0;
            DepositShare_h(j)=(NL_1(j)+NL_2(j))/(N1+N2);
          }
          Deposits_1(1,i)+=injection2;
          Deposits(1,receivingBank)+=injection2;
          Inflows(receivingBank)+=injection2;
        }
        
        if(flag_entry==1)
        {
          injection=multip_entry*mD1-injection2;
          LossEntry_b(receivingBank)+=injection;
          Injection_1(i)+=injection;
          Deposits_1(1,i)+=injection;
          Deposits(1,receivingBank)+=injection;
        }
        else
        {
          injection=multip_entry*mD1-injection2;
          EntryCosts+=injection;
          FirmTransfers_1+=injection;
          Injection_1(i)+=injection;
          Deposits_1(1,i)+=injection;
          Deposits(1,receivingBank)+=injection;
          Inflows(receivingBank)+=injection;
        }
      }

      S1(i)=p1(i)*step;
      stepb=step;
      while (stepb > 0)							
      {											
        rni=int(ran1(p_seed)*N1*N2)%N2+1;
        if (Match(rni,i) == 0)
        {
          Match(rni,i)=1;
          stepb--;
        }
      }
    }
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      f2_exit+=f2(1,j);
    }
  }

  if(exiting_2.Sum()>0){
    CurrentDemand=D2.Row(1).Sum();
    if(f2_exit>0){
      n_mach_needed=max(exiting_2.Sum(),ceil(f2_exit*CurrentDemand/dim_mach/u));
    }else{
      f2_exit=f2_entry_min*exiting_2.Sum();
      n_mach_needed=max(exiting_2.Sum(),ceil(f2_exit*CurrentDemand/dim_mach/u));
    }
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          if (gtemp[tt-1][i-1][j-1]>0)
          {
            n_mach_exit+=gtemp[tt-1][i-1][j-1];
          }
        }
      }
    }
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      for (i=1; i<=N1; i++)
      {
        for (tt=t; tt>=t0; tt--)
        {
          if (gtemp[tt-1][i-1][j-1]>0)
          {
            C_secondhand(tt,i)=C(tt,i);
          }
        }
      }
    }
    else
    {
      //Calculate mean deposits of surviving C-firms
      mD2+=Deposits_2(1,j);
			ns2++;
    }
  }

  //The capital stock of exiting C-firms is transferred to their respective banks up to the value of bad debt; the value of second-hand capital is marked down depending on age of the machines
  n_mach_exit2=min(n_mach_needed,n_mach_exit);
  //The cheapest second-hand machine still on offer. It only moves when one is
  //taken below, so it is worked out here and again after each is taken rather
  //than rescanned for every vintage of every exiting firm.
  //
  //Only the vintages still in use carry a price: the rest of the matrix is set
  //to infinity for the period and cannot hold the minimum, so the scan covers
  //434 entries rather than all 12,000.
  double min_secondhand=RowRangeMinimum(C_secondhand,t0,t);
  while(n_mach_exit2>0)
  {
    for(j=1; j<=N2; j++)
    {
      if(exiting_2(j)==1)
      {
        receivingBank=BankingSupplier_2(j);
        baddebt_2_temp=baddebt_2(j);
        for (i=1; i<=N1; i++)
        {
          for (tt=t; tt>=t0; tt--)
          {
            if (gtemp[tt-1][i-1][j-1]>0 && n_mach_exit2>0 && C(tt,i)<=min_secondhand)
            {
              markdownCapital=max(0.0,(1-(double)age[tt-1][i-1][j-1]/(agemax)));
              g_secondhand[tt-1][i-1]+=min(n_mach_exit2,gtemp[tt-1][i-1][j-1]);
              age_secondhand[tt-1][i-1]=age[tt-1][i-1][j-1];
              g_secondhand_p[tt-1][i-1]=markdownCapital*g_price[tt-1][i-1][j-1];
              if(baddebt_2_temp>0)
              {
                capitalRecovered(receivingBank)+=min(n_mach_exit2,gtemp[tt-1][i-1][j-1])*g_secondhand_p[tt-1][i-1];
                baddebt_2_temp-=min(n_mach_exit2,gtemp[tt-1][i-1][j-1])*g_secondhand_p[tt-1][i-1];
              }
              n_mach_exit2-=min(n_mach_exit2,gtemp[tt-1][i-1][j-1]);
              C_secondhand(tt,i)=1000000;
              min_secondhand=RowRangeMinimum(C_secondhand,t0,t);
            }
          }
        }
      }
    }
  }
  
	if (ns2 > 0)
	{
		mD2/=ns2;
	}
	else
  {
    std::cerr << "\n\n ERROR: All C-firms are exiting in period " << t << endl;
    Errors << "\n All C-firms are exiting in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  //K-firms lose exiting C-firms as customers
  for (j=1; j<=N2; j++)
  {
    if (exiting_2(j)==1)		
    {
      flag=0;
      indsupl=int(supl(j));
      if (indsupl>=1)
      {
        Match(j,supl(j))=0;				
        supl(j)=0;
      }
    }
  }

  //Determine the number of machines which each entering C-firm will have based on number of available second-hand machines
  n_mach_resid=min(n_mach_needed,n_mach_exit);
  if(exiting_2.Sum()>n_mach_resid)
  {
    std::cerr << "\n\n ERROR: Not enough second hand capital in period " << t << endl;
    Errors << "\n Not enough second hand capital in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1)
    {
      n_mach_entry(j)=1;
      n_mach_resid--;
      k_entry(j)=ran1(p_seed);
    }
  }

  n_mach_resid2=n_mach_resid;

  for(j=1; j<=N2; j++)
  {
    if(exiting_2(j)==1 && n_mach_resid>0)
    {
      n_mach_entry(j)+=floor(k_entry(j)/k_entry.Sum()*n_mach_resid2);
      n_mach_resid-=floor(k_entry(j)/k_entry.Sum()*n_mach_resid2);
    }
  }


  if(n_mach_resid<0 && n_mach_exit<n_mach_needed)
  {
    std::cerr << "\n\n ERROR: Remaining second hand machines are negative in period " << t << endl;
    Errors << "\n Remaining second hand machines are negative " << t << endl;
    exit(EXIT_FAILURE);
  }

  while(n_mach_resid>0)
  {
    rni=int(ran1(p_seed)*N1*N2)%N2+1;
    if(exiting_2(rni)==1 && n_mach_resid>0)
    {
      n_mach_entry(rni)++;
      n_mach_resid--;
    }
  }

  //Second-hand capital is purchased by households; below will be transferred to newly entering firms
  if(Deposits_h(1)>=(capitalRecovered.Sum()))
  {
    Deposits_h(1)-=capitalRecovered.Sum();
    FirmTransfers+=capitalRecovered.Sum();
    FirmTransfers_2+=capitalRecovered.Sum();
    for(i=1; i<=NB; i++)
    {
      Deposits_hb(1,i)-=DepositShare_h(i)*capitalRecovered.Sum();
      Deposits(1,i)-=DepositShare_h(i)*capitalRecovered.Sum();
      if(capitalRecovered(i)>=DepositShare_h(i)*capitalRecovered.Sum())
      {
        Inflows(i)+=(capitalRecovered(i)-DepositShare_h(i)*capitalRecovered.Sum());
      }
      else
      {
        Outflows(i)+=(DepositShare_h(i)*capitalRecovered.Sum()-capitalRecovered(i));
      }
    }
  }
  else if (Deposits_h(1)>=0)
  {
    //Household deposits are insufficient to purchase second hand capital; households buy as much as they can, the rest is financed either by government (flag_entry=0) or booked as a loss by banks (flag_entry=1)
    std::cerr << "\n\n Households cannot purchase second-hand capital in period " << t << endl;
    Errors << "\n Households cannot purchase second-hand capital in period " << t << endl;
    FirmTransfers+=Deposits_h(1);
    FirmTransfers_2+=Deposits_h(1);
    capitalRecovered2=capitalRecovered;
    for(i=1; i<=NB; i++)
    {
      if(capitalRecovered2(i)>=Deposits_hb(1,i))
      {
        capitalRecovered2(i)-=Deposits_hb(1,i);
        Deposits(1,i)-=Deposits_hb(1,i);
        Deposits_h(1)-=Deposits_hb(1,i);
        Deposits_hb(1,i)=0;
      }
      else
      {
        Deposits(1,i)-=Deposits_hb(1,i);
        Deposits_hb(1,i)-=capitalRecovered2(i);
        Deposits_h(1)-=capitalRecovered2(i);
        capitalRecovered2(i)=0;
        Outflows(i)+=Deposits_hb(1,i);
        Deposits_hb(1,i)=0;
      }
      DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
    }

    for(i=1; i<=NB; i++)
    {
      capitalRecoveredShare(i)=capitalRecovered2(i)/capitalRecovered2.Sum();
    }
    
    for(i=1; i<=NB; i++)
    {
      capitalRecovered2(i)-=(capitalRecoveredShare(i)*Deposits_h(1));
      Inflows(i)+=(capitalRecoveredShare(i)*Deposits_h(1));
    }
    
    if(flag_entry==1)
    {
      for(i=1; i<=NB; i++)
      {
        LossEntry_b(i)+=capitalRecovered2(i);
      }
    }
    else
    {
      EntryCosts+=capitalRecovered2.Sum();
      BankTransfer+=capitalRecovered2.Sum();
      for(i=1; i<=NB; i++)
      {
        Inflows(i)+=capitalRecovered2(i);
      }
    }
    Deposits_h(1)=0;
  }
  else
  {
    std::cerr << "\n\n ERROR: Household deposits are negative in period " << t << endl;
    Errors << "\n Household deposits are negative in period " << t << endl;
    exit(EXIT_FAILURE);
  }

  for(j=1; j<=N2; j++)
  {
    if (exiting_2(j)==1)
    { 
      N(1,j)=0;
      N(2,j)=0;
      //Any inventories of exiting firms are destroyed and hence booked as a loss
      Injection_2(j)=-Inventories(1,j);
      Inventories(1,j)=0;
      Inventories(2,j)=0;
      //Households give entering C-firms a transfer of deposits
      receivingBank=BankingSupplier_2(j);
      multip_entry=ran1(p_seed);
      multip_entry=w2inf+multip_entry*(w2sup-w2inf);
      if(Deposits_h(1)>=multip_entry*mD2)
      {
        injection=multip_entry*mD2;
        Deposits_h(1)-=injection;
        FirmTransfers+=injection;
        FirmTransfers_2+=injection;
        Injection_2(j)+=injection;
        for(i=1; i<=NB; i++)
        {
          Deposits_hb(1,i)-=DepositShare_h(i)*injection;
          Deposits(1,i)-=DepositShare_h(i)*injection;
          Outflows(i)+=DepositShare_h(i)*injection;
        }
        Deposits_2(1,j)=injection;
        Deposits(1,receivingBank)+=injection;
        Inflows(receivingBank)+=injection;
      }
      else
      {
        //If households cannot finance C-firm entry, this is done by government or banks as in the case of K-firms
        std::cerr << "\n\n Households cannot finance C-firm entry in period " << t << endl;
        Errors << "\n Households cannot finance C-firm entry in period " << t << endl;
        
        injection2=0;
        if(Deposits_h(1)>0)
        {
          injection2=Deposits_h(1);
          Deposits_h(1)-=injection2;
          FirmTransfers+=injection2;
          FirmTransfers_2+=injection2;
          Injection_2(j)+=injection2;
          for(i=1; i<=NB; i++)
          {
            Deposits(1,i)-=Deposits_hb(1,i);
            Outflows(i)+=Deposits_hb(1,i);
            Deposits_hb(1,i)=0;
            DepositShare_h(i)=(NL_1(i)+NL_2(i))/(N1+N2);
          }
          Deposits_2(1,j)=injection2;
          Deposits(1,receivingBank)+=injection2;
          Inflows(receivingBank)+=injection2;
        }
        
        if(flag_entry==1)
        {
          injection=multip_entry*mD2-injection2;
          LossEntry_b(receivingBank)+=injection;
          Injection_2(j)+=injection;
          Deposits_2(1,j)+=injection;
          Deposits(1,receivingBank)+=injection;
        }
        else
        {
          injection=multip_entry*mD2-injection2;
          EntryCosts+=injection;
          FirmTransfers_2+=injection;
          Injection_2(j)+=injection;
          Deposits_2(1,j)+=injection;
          Deposits(1,receivingBank)+=injection;
          Inflows(receivingBank)+=injection;
        }
      }
      	
      n_mach(j)=0;
      K(j)=0;
      //First subtract the capital stock previously held by the exiting firm
      Injection_2(j)-=(CapitalStock(1,j)+deltaCapitalStock(1,j));
      CapitalStock(1,j)=0;
      //Clear the exiting firms' entries in the frequency arrays
      n_mach_resid=n_mach_entry(j);
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {							
          g[tt-1][i-1][j-1]=0;
          gtemp[tt-1][i-1][j-1]=0;
          g_c[tt-1][i-1][j-1]=0;
          g_c2[tt-1][i-1][j-1]=0;
          g_c3[tt-1][i-1][j-1]=0;
          age[tt-1][i-1][j-1]=0;
        }
      }

      //Give the entering firm an initial capital stock from the pool of second-hand capital
      while(n_mach_resid>0)
      {
        rni=int(ran1(p_seed)*N1*N2)%N1+1;	
        for (tt=t0; tt<=t; tt++)
        {	
          if(g_secondhand[tt-1][rni-1]>0)
          {
            if(g_secondhand[tt-1][rni-1]>=n_mach_resid)
            {
              g[tt-1][rni-1][j-1]+=n_mach_resid;
              gtemp[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c2[tt-1][rni-1][j-1]+=n_mach_resid;
              g_c3[tt-1][rni-1][j-1]+=n_mach_resid;
              age[tt-1][rni-1][j-1]=age_secondhand[tt-1][rni-1];
              n_mach(j)+=n_mach_resid;
              CapitalStock(1,j)+=n_mach_resid*g_secondhand_p[tt-1][rni-1];
              g_price[tt-1][rni-1][j-1]=g_secondhand_p[tt-1][rni-1];
              K(j)+=n_mach_resid*dim_mach;
              g_secondhand[tt-1][rni-1]-=n_mach_resid;
              n_mach_resid=0;
            }
            else
            {
              g[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              gtemp[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c2[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              g_c3[tt-1][rni-1][j-1]+=g_secondhand[tt-1][rni-1];
              age[tt-1][rni-1][j-1]=age_secondhand[tt-1][rni-1];
              n_mach(j)+=g_secondhand[tt-1][rni-1];
              CapitalStock(1,j)+=g_secondhand[tt-1][rni-1]*g_secondhand_p[tt-1][rni-1];
              g_price[tt-1][rni-1][j-1]=g_secondhand_p[tt-1][rni-1];
              K(j)+=g_secondhand[tt-1][rni-1]*dim_mach;
              n_mach_resid-=g_secondhand[tt-1][rni-1];
              g_secondhand[tt-1][rni-1]=0;
            }
          }
        }
      }

      //New capital stock is added to the net worth injection
      Injection_2(j)+=CapitalStock(1,j);
      rni=int(ran1(p_seed)*N1*N2)%N1+1;	
      supl(j)=rni;
      Match(j,rni)=1;
      EI(1,j)=0;
      scrap_age(j)=0;
      deltaCapitalStock(1,j)=0;
      //Set the newly entering firm's cost, mark-up and price
      c2(j)=0;						
      for (i=1; i<=N1; i++)
      {
        for (tt=t0; tt<=t; tt++)
        {
          c2(j)+=(w(2)/((1-shocks_labprod2(j))*A(tt,i))+c_en(2)/((1-shocks_eneff2(j))*A_en(tt,i))+t_CO2*A_ef(tt,i)/((1-shocks_eneff2(j))*A_en(tt,i)))*g[tt-1][i-1][j-1]/n_mach(j);
        }
      }
      mu2(1,j)=mi2;
      p2(j)=(1+mu2(1,j))*c2(j);
      p2_entry+=p2(j);
      DebtService_2(1,j)=0;
    }
  }

  if(exiting_2.Sum()>0)
  {
    CurrentDemand=D2.Row(1).Sum();
    if(f2_exit>0){
      p2_entry/=exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          CompEntry(j)=-p2(j)/p2_entry;
        }
      }
      CompEntry_m=CompEntry.Sum()/exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          EntryShare(j)=(1/exiting_2.Sum())*((2*omega3)/(1+exp((-chi)*((CompEntry(j)-CompEntry_m)/CompEntry_m)))+(1-omega3));
        }
      }
      EntryShare=EntryShare/EntryShare.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          f2(1,j)=max(f2_entry_min,EntryShare(j)*f2_exit);
          f2(2,j)=f2(1,j);
          f2(3,j)=f2(1,j);
        }
      }
      // Market shares are normalised after entering C-firms get a market share
      ftot(1)=f2.Row(1).Sum();
      ftot(2)=f2.Row(2).Sum();
      ftot(3)=f2.Row(3).Sum();

      for(j=1; j<=N2; j++)
      {
        f2(1,j)/=ftot(1);
        f2(2,j)/=ftot(2);
        f2(3,j)/=ftot(3);
        if (exiting_2(j)==1)
        { 
          D2(1,j)=min(K(j),f2(1,j)*CurrentDemand);
          l2(j)=1+(f2(1,j)*CurrentDemand-D2(1,j));
          De(j)=D2(1,j);
          S2(1,j)=p2(j)*D2(1,j);
          mol(j)=S2(1,j)-D2(1,j)*c2(j);
          if(mol(j)<0)
          {
            mol(j)=0;
          }
        }
      }
    }else{
      f2_exit=f2_entry_min*exiting_2.Sum();
      p2_entry/=exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          CompEntry(j)=-p2(j)/p2_entry;
        }
      }
      CompEntry_m=CompEntry.Sum()/exiting_2.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          EntryShare(j)=(1/exiting_2.Sum())*((2*omega3)/(1+exp((-chi)*((CompEntry(j)-CompEntry_m)/CompEntry_m)))+(1-omega3));
        }
      }
      EntryShare=EntryShare/EntryShare.Sum();
      for(j=1; j<=N2; j++)
      {
        if (exiting_2(j)==1)
        { 
          f2(1,j)=max(f2_entry_min,EntryShare(j)*f2_exit);
          f2(2,j)=f2(1,j);
          f2(3,j)=f2(1,j);
        }
      }
      // Market shares are normalised after entering C-firms get a market share
      ftot(1)=f2.Row(1).Sum();
      ftot(2)=f2.Row(2).Sum();
      ftot(3)=f2.Row(3).Sum();

      for(j=1; j<=N2; j++)
      {
        f2(1,j)/=ftot(1);
        f2(2,j)/=ftot(2);
        f2(3,j)/=ftot(3);
        if (exiting_2(j)==1)
        { 
          D2(1,j)=min(K(j),f2(1,j)*CurrentDemand);
          l2(j)=1+(f2(1,j)*CurrentDemand-D2(1,j));
          De(j)=D2(1,j);
          S2(1,j)=p2(j)*D2(1,j);
          mol(j)=S2(1,j)-D2(1,j)*c2(j);
          if(mol(j)<0)
          {
            mol(j)=0;
          }
        }
      }
    }
  }

  //Update C-firm K-firm network
	nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			nclient(i)+=Match(j,i);
		}
	}
	for (i=1; i<=N1; i++)
	{
		if (nclient(i) == 0)
		{
			stepb=step;
			while (stepb > 0)	
			{
				rni=int(ran1(p_seed)*N1*N2)%N2+1;
				if (Match(rni,i) == 0)
				{
					Match(rni,i)=1;
					stepb--;
				}
			}
		}
	}

  nclient=0;
	for (i=1; i<=N1; i++)
	{
		for (j=1; j<=N2; j++)
		{
			nclient(i)+=Match(j,i);
		}
	}
}

void TECHANGEND(void)
{    
  //Endogenous technological change
  Inn=0;
  Imm=0;
  
  A1inn=0.00001;               
  A1pinn=0.00001;                
  A1imm=0.00001;               
  A1pimm=0.00001;
  
  EE_inn=0.00001;                            
  EEp_inn=0.00001;
  EE_imm=0.00001;
  EEp_imm=0.00001;
  
  EF_inn=std::numeric_limits<double>::infinity();
  EFp_inn=std::numeric_limits<double>::infinity();
  EF_imm=std::numeric_limits<double>::infinity();
  EFp_imm=std::numeric_limits<double>::infinity();

  for (i=1; i<=N1; i++)
  {
    //K-firms determine R&D spending and associated labour demand
    RD(1,i)=nu*S1(i);              
    if (S1(i)==0)                  
    {                                              
      RD(1,i)=RD(2,i);                           
    
      if (nclient(i) < 1)            
      {
        std::cerr << "\n\n ERROR: nclient < 1 for K-firm " << i << " in period " << t << endl;
        Errors << "\n nclient < 1 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    if (w(1)>0)
    {
      Ld1rd(i)=RD(1,i)/w(1);               
    }
    else
    {
      std::cerr << "\n\n ERROR: w=0 in period " << t << endl;
      Errors << "\n w=0 in period " << t << endl;
      exit(EXIT_FAILURE);
    }

    //Divide between innovation and imitation
    RDin(i)=Ld1rd(i)*xi;
    RDim(i)=Ld1rd(i)*(1-xi);
    
    //Determine whether firm innovates and/or imitates
    parber=1-exp(-o1*RDin(i));             
    Inn(i)=bnldev(parber,1,p_seed);          

    parber=1-exp(-o2*RDim(i));             
    Imm(i)=bnldev(parber,1,p_seed);          

    if (Inn(i) == 1)                 
    {  
      //If firm innovates, determine characteristics of new technology
      //Labour productivity
      rnd=betadev(b_a11,b_b11,p_seed);         
      rnd=uu11+rnd*(uu21-uu11);                         
      A1inn(i)=A1(i)*(1+rnd);
      
      rnd=betadev(b_a12,b_b12,p_seed);         
      rnd=uu12+rnd*(uu22-uu12);                        
      A1pinn(i)=A1p(i)*(1+rnd);

      //Energy efficiency
      rnd=betadev(b_a2,b_b2,p_seed);         
      rnd=uu31+rnd*(uu41-uu31);                         
      EE_inn(i)=A1_en(i)*(1+rnd);

      rnd=betadev(b_a2,b_b2,p_seed);         
      rnd=uu32+rnd*(uu42-uu32);                         
      EEp_inn(i)=A1p_en(i)*(1+rnd);

      //Environmental friendliness
      rnd=betadev(b_a3,b_b3,p_seed);         
      rnd=uu51+rnd*(uu61-uu51);                         
      EF_inn(i)=A1_ef(i)*(1-rnd);

      rnd=betadev(b_a3,b_b3,p_seed);        
      rnd=uu52+rnd*(uu62-uu52);                       
      EFp_inn(i)=A1p_ef(i)*(1-rnd);
      
      if (A1pinn(i)==0 || A1inn(i)==0 || A1p(i)==0 || A1(i)==0)
      {
        std::cerr << "\n\n ERROR: A1=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EEp_inn(i)==0 || EE_inn(i)==0 || A1p_en(i)==0 || A1_en(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_en=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_en=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }
      
    if (Imm(i) == 1)
    {
      //If K-firm imitates, determine which other firm's technology it will imitate
      Tdtot=0;
      for (ii=1; ii<=N1; ii++)
      {
          Td.element(ii)=sqrt(((A1(ii)-A1(i))*(A1(ii)-A1(i))) + ((A1p(ii)-A1p(i))*(A1p(ii)-A1p(i))) + ((A1_en(ii)-A1_en(i))*(A1_en(ii)-A1_en(i))) + ((A1_ef(ii)-A1_ef(i))*(A1_ef(ii)-A1_ef(i))) + ((A1p_en(ii)-A1p_en(i))*(A1p_en(ii)-A1p_en(i))) + ((A1p_ef(ii)-A1p_ef(i))*(A1p_ef(ii)-A1p_ef(i))));
          if (Td.element(ii)>0)
          {
            Td.element(ii)=1/Td.element(ii);
          }
          else 
          {
            Td.element(ii)=0;
          }
          Tdtot+=Td.element(ii);
      }
      for (ii=1; ii<=N1; ii++)          
      {
        Td.element(ii)/=Tdtot;
        Td.element(ii)+=Td.element(ii-1);
      }
      rnd=ran1(p_seed);                
      for (ii=1; ii<=N1; ii++)           
      {
        if (rnd <= Td.element(ii) && rnd > Td.element(ii-1))
        {
          A1imm(i)=A1(ii);
          A1pimm(i)=A1p(ii);
          EE_imm(i)=A1_en(ii);
          EEp_imm(i)=A1p_en(ii);
          EF_imm(i)=A1_ef(ii);
          EFp_imm(i)=A1p_ef(ii);
        }
      }

      if (A1pimm(i)==0 || A1imm(i)==0 || A1p(i)==0 || A1(i)==0)
      {
        std::cerr << "\n\n ERROR: A1=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EEp_imm(i)==0 || EE_imm(i)==0 || A1p_en(i)==0 || A1_en(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_en=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_en=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
      if (EFp_imm(i)==0 || EF_imm(i)==0 || A1p_ef(i)==0 || A1_ef(i)==0)
      {
        std::cerr << "\n\n ERROR: A1_ef=0 for K-firm " << i << " in period " << t << endl;
        Errors << "\n A1_ef=0 for K-firm " << i << " in period " << t << endl;
        exit(EXIT_FAILURE);
      }
    }  
    
    //If the imitated technology is superior, adopt it
    if ( ((1+mi1)*(w(1)/(A1pimm(i))+c_en(1)/EEp_imm(i)+t_CO2*EFp_imm(i)/EEp_imm(i)))+(w(1)/A1imm(i)+c_en(1)/EE_imm(i)+t_CO2*EF_imm(i)/EE_imm(i))*b < ((1+mi1)*(w(1)/(A1p(i))+c_en(1)/A1p_en(i)+t_CO2*A1p_ef(i)/A1p_en(i))+(w(1)/A1(i)+c_en(1)/A1_en(i)+t_CO2*A1_ef(i)/A1_en(i))*b ))
    {
      A1(i)=A1imm(i);                            
      A1p(i)=A1pimm(i);                        
      A1_en(i)=EE_imm(i);                       
      A1p_en(i)=EEp_imm(i);                      
      A1_ef(i)=EF_imm(i);                        
      A1p_ef(i)=EFp_imm(i);                      
    }

    //If the innovated technology is superior, adopt it
    if ( ((1+mi1)*(w(1)/(A1pinn(i))+c_en(1)/EEp_inn(i)+t_CO2*EFp_inn(i)/EEp_inn(i)))+(w(1)/A1inn(i)+c_en(1)/EE_inn(i)+t_CO2*EF_inn(i)/EE_inn(i))*b < ((1+mi1)*(w(1)/(A1p(i))+c_en(1)/A1p_en(i)+t_CO2*A1p_ef(i)/A1p_en(i))+(w(1)/A1(i)+c_en(1)/A1_en(i)+t_CO2*A1_ef(i)/A1_en(i))*b )) 
    {
      A1(i)=A1inn(i);                            
      A1p(i)=A1pinn(i);                        
      A1_en(i)=EE_inn(i);                        
      A1p_en(i)=EEp_inn(i);                      
      A1_ef(i)=EF_inn(i);                        
      A1p_ef(i)=EFp_inn(i);                      
    }
    
    //Update the productivity matrices
    if (t < T)
    {
      A(t+1,i)=A1(i);                         
      A_en(t+1,i)=A1_en(i);
      A_ef(t+1,i)=A1_ef(i);
    }
  }
  
  LD1rdtot=Ld1rd.Sum();
}

void DEPOSITCHECK(void)
{
  for(i=1; i<=NB; i++)
  {
    deviation=fabs(DepositShare_e(i)-Deposits_eb(1,i)/Deposits_eb.Row(1).Sum());
    if(deviation>tolerance && Deposits_eb(1,i)>tolerance)
    {
      std::cerr<<"Share error Deposits_eb for bank " << i << " in period " << t << endl;
      Errors << "\n Share error Deposits_eb for bank " << i << " in period " << t << endl;
    }
    deviation=fabs(DepositShare_h(i)-Deposits_hb(1,i)/Deposits_hb.Row(1).Sum());
    if(deviation>tolerance)
    {
      std::cerr<<"Share error Deposits_hb for bank " << i << " in period " << t << endl;
      Errors << "\n Share error Deposits_hb for bank " << i << " in period " << t << endl;
    }
    DepositsCheck_1=Deposits(1,i)-Deposits_hb(1,i)-Deposits_eb(1,i);
    DepositsCheck_2=0;
    for(j=1; j<=N1; j++)
    {
      if(BankMatch_1(j,i)==1)
      {
        DepositsCheck_2+=Deposits_1(1,j);
      }
    }
    for(j=1; j<=N2; j++)
    {
      if(BankMatch_2(j,i)==1)
      {
        DepositsCheck_2+=Deposits_2(1,j);
      }
    }
    deviation=fabs((DepositsCheck_1-DepositsCheck_2)/Deposits(1,i));
    if(deviation>tolerance)
    {
      std::cerr<<"Share error firm deposits for bank " << i << " in period " << t << endl;
      Errors << "\n Share error firm deposits for bank " << i << " in period " << t << endl;
    }
  }
}

void NEGATIVITYCHECK(void)
{
  for (j=1; j<=N2; ++j)
  {
    if(Loans_2(1,j)<0)
    {
      std::cerr<<"Error loans for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error loans for C-firm " << j << " in period " << t << endl;
    }
    if(Deposits_2(1,j)<0)
    {
      std::cerr<<"Error deposits for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error deposits for C-firm " << j << " in period " << t << endl;
    }
    if(CapitalStock(1,j)<0)
    {
      std::cerr<<"Error Capital for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Capital for C-firm " << j << " in period " << t << endl;
    }
    if(Inventories(1,j)<0)
    {
      std::cerr<<"Error Inventories for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Inventories for C-firm " << j << " in period " << t << endl;
    }
    if(Investment_2(j)<0)
    {
      std::cerr<<"Error Investment for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Investment for C-firm " << j << " in period " << t << endl;
    }
    if(Taxes_2(j)<0)
    {
      std::cerr<<"Error Taxes for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Taxes for C-firm " << j << " in period " << t << endl;
    }
    if(Wages_2(j)<0)
    {
      std::cerr<<"Error Wages for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Wages for C-firm " << j << " in period " << t << endl;
    }
    if(EnergyPayments_2(j)<0)
    {
      std::cerr<<"Error EnergyPayments for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error EnergyPayments for C-firm " << j << " in period " << t << endl;
    }
    if(Dividends_2(j)<0)
    {
      std::cerr<<"Error Dividends for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Dividends for C-firm " << j << " in period " << t << endl;
    }
    if(LoanInterest_2(j)<0)
    {
      std::cerr<<"Error LoanInterest for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error LoanInterest for C-firm " << j << " in period " << t << endl;
    }
    if(InterestDeposits_2(j)<0)
    {
      std::cerr<<"Error InterestDeposits for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for C-firm " << j << " in period " << t << endl;
    }
    if(S2(1,j)<0)
    {
      std::cerr<<"Error S2 for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error S2 for C-firm " << j << " in period " << t << endl;
    }
    if(DebtRemittances2(j)<0)
    {
      std::cerr<<"Error DebtRemittances for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error DebtRemittances for C-firm " << j << " in period " << t << endl;
    }
    if(Taxes_CO2_2(j)<0)
    {
      std::cerr<<"Error Taxes_CO2 for C-firm " << j << " in period " << t << endl;
      Errors << "\n Error Taxes_CO2 for C-firm " << j << " in period " << t << endl;
    }
  }

  for (i=1; i<=N1; i++)
  {
    if(Deposits_1(1,i)<0)
    {
      std::cerr<<"Error Deposits for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Deposits for K-firm " << i << " in period " << t << endl;
    }
    if(S1(i)<0)
    {
      std::cerr<<"Error S1 for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error S1 for K-firm " << i << " in period " << t << endl;
    }
    if(Taxes_1(i)<0)
    {
      std::cerr<<"Error Taxes for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Taxes for K-firm " << i << " in period " << t << endl;
    }
    if(Wages_1(i)<0)
    {
      std::cerr<<"Error Wages for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Wages for K-firm " << i << " in period " << t << endl;
    }
    if(EnergyPayments_1(i)<0)
    {
      std::cerr<<"Error EnergyPayments for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error EnergyPayments for K-firm " << i << " in period " << t << endl;
    }
    if(Dividends_1(i)<0)
    {
      std::cerr<<"Error Dividends for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Dividends for K-firm " << i << " in period " << t << endl;
    }
    if(InterestDeposits_1(i)<0)
    {
      std::cerr<<"Error InterestDeposits for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for K-firm " << i << " in period " << t << endl;
    }
    if(Taxes_CO2_1(i)<0)
    {
      std::cerr<<"Error Taxes_CO2 for K-firm " << i << " in period " << t << endl;
      Errors << "\n Error Taxes_CO2 for K-firm " << i << " in period " << t << endl;
    }
  }

  for (i=1; i<=NB; i++)
  {
    if(Loans_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Loans for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Loans for Bank " << i << " in period " << t << endl;
    }
    if(Deposits(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits for Bank " << i << " in period " << t << endl;
    }
    if(Deposits_hb(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits_hb for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits_hb  for Bank " << i << " in period " << t << endl;
    }
    if(Deposits_eb(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Deposits_eb for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Deposits_eb for Bank " << i << " in period " << t << endl;
    }
    if(GB_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error GB for Bank " << i << " in period " << t << endl;
      Errors << "\n Error GB for Bank " << i << " in period " << t << endl;
    }
    if(Reserves_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Reserves for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Reserves for Bank " << i << " in period " << t << endl;
    }
    if(Advances_b(1,i)<(-tolerance*cpi(1)))
    {
      std::cerr<<"Error Advances for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Advances for Bank " << i << " in period " << t << endl;
    }
    if(Taxes_b(i)<0)
    {
      std::cerr<<"Error Taxes for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Taxes for Bank " << i << " in period " << t << endl;
    }
    if(Dividends_b(i)<0)
    {
      std::cerr<<"Error Dividends for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Dividends for Bank " << i << " in period " << t << endl;
    }
    if(LoanInterest(i)<0)
    {
      std::cerr<<"Error LoanInterest for Bank " << i << " in period " << t << endl;
      Errors << "\n Error LoanInterest for Bank " << i << " in period " << t << endl;
    }
    if(InterestDeposits(i)<0)
    {
      std::cerr<<"Error InterestDeposits for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestDeposits for Bank " << i << " in period " << t << endl;
    }
    if(InterestBonds_b(i)<0)
    {
      std::cerr<<"Error InterestBonds for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestBonds for Bank " << i << " in period " << t << endl;
    }
    if(BondRepayments_b(i)<0)
    {
      std::cerr<<"Error BondRepayments for Bank " << i << " in period " << t << endl;
      Errors << "\n Error BondRepayments for Bank " << i << " in period " << t << endl;
    }
    if(InterestReserves_b(i)<0)
    {
      std::cerr<<"Error InterestReserves for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestReserves for Bank " << i << " in period " << t << endl;
    }
    if(InterestAdvances_b(i)<0)
    {
      std::cerr<<"Error InterestAdvances for Bank " << i << " in period " << t << endl;
      Errors << "\n Error InterestAdvances for Bank " << i << " in period " << t << endl;
    }
    if(Bailout_b(i)<0)
    {
      std::cerr<<"Error Bailout for Bank " << i << " in period " << t << endl;
      Errors << "\n Error Bailout for Bank " << i << " in period " << t << endl;
    }
  }

  if(Deposits_h(1)<0)
  {
    std::cerr<<"Error Deposits_h in period " << t << endl;
    Errors << "\n Error Deposits_h in period " << t << endl;
  }
  if(Consumption<0)
  {
    std::cerr<<"Error Consumption in period " << t << endl;
    Errors << "\n Error Consumption in period " << t << endl;
  }
  if(Benefits<0)
  {
    std::cerr<<"Error Benefits in period " << t << endl;
    Errors << "\n Error Benefits in period " << t << endl;
  }
  if(Taxes_h<0)
  {
    std::cerr<<"Error Taxes_h in period " << t << endl;
    Errors << "\n Error Taxes_h in period " << t << endl;
  }
  if(Wages<0)
  {
    std::cerr<<"Error Wages in period " << t << endl;
    Errors << "\n Error Wages in period " << t << endl;
  }
  if(Dividends(1)<0)
  {
    std::cerr<<"Error Dividends in period " << t << endl;
    Errors << "\n Error Dividends in period " << t << endl;
  }
  if(InterestDeposits_h<0)
  {
    std::cerr<<"Error InterestDeposits_h in period " << t << endl;
    Errors << "\n Error InterestDeposits_h in period " << t << endl;
  }

  if(G<0)
  {
    std::cerr<<"Error G in period " << t << endl;
    Errors << "\n Error G in period " << t << endl;
  }
  if(Taxes<0)
  {
    std::cerr<<"Error Taxes in period " << t << endl;
    Errors << "\n Error Taxes in period " << t << endl;
  }
  if(Bailout<0)
  {
    std::cerr<<"Error Bailout in period " << t << endl;
    Errors << "\n Error Bailout in period " << t << endl;
  }
  if(InterestReserves<0)
  {
    std::cerr<<"Error InterestReserves in period " << t << endl;
    Errors << "\n Error InterestReserves in period " << t << endl;
  }
  if(InterestAdvances<0)
  {
    std::cerr<<"Error InterestAdvances in period " << t << endl;
    Errors << "\n Error InterestAdvances in period " << t << endl;
  }

  if(EnergyPayments<0)
  {
    std::cerr<<"Error EnergyPayments in period " << t << endl;
    Errors << "\n Error EnergyPayments in period " << t << endl;
  }
  if(Wages_en<0)
  {
    std::cerr<<"Error Wages_en in period " << t << endl;
    Errors << "\n Error Wages_en in period " << t << endl;
  }
  if(Dividends_e<0)
  {
    std::cerr<<"Error Dividends_e in period " << t << endl;
    Errors << "\n Error Dividends_e in period " << t << endl;
  }
  if(InterestDeposits_e<0)
  {
    std::cerr<<"Error InterestDeposits_e in period " << t << endl;
    Errors << "\n Error InterestDeposits_e in period " << t << endl;
  }

  if(flag_WITCH_on==1)
  {
    for(j=0; j<ene_tecs.size(); j++)
    {
      tech=ene_tecs[j];
      if(wage_e_mult_tech[tech]<0)
      {
          std::cerr<<"Error Wages_en for technology " << tech << " in period " << t << endl;
          Errors << "\n Error Wages_en for technology " << tech << " in period " << t << endl;
      }
      if(dividends_e_mult_tech[tech]<0)
      {
          std::cerr<<"Error Dividends_e for technology " << tech << " in period " << t << endl;
          Errors << "\n Error Dividends_e for technology " << tech << " in period " << t << endl;
      }
      if(interestsDeposits_e_mult_tech[tech]<0)
      {
          std::cerr<<"Error InterestDeposits_e for technology " << tech << " in period " << t << endl;
          Errors << "\n Error InterestDeposits_e for technology " << tech << " in period " << t << endl;
      }
    }
  }
}

void CHECKSUMS(void)
{
  deviation=fabs((Deposits_h(1)-Deposits_hb.Row(1).Sum())/Deposits_hb.Row(1).Sum());
  if(deviation>tolerance && Deposits_hb.Row(1).Sum()>tolerance && Deposits_h(1)>tolerance)
  {
    std::cerr<<"Sum error Deposits_h in period " << t << endl;
    Errors << "\n Sum error Deposits_h in period " << t << endl;
  }
  deviation=fabs((Deposits_e(1)-Deposits_eb.Row(1).Sum())/Deposits_eb.Row(1).Sum());
  if(deviation>tolerance && Deposits_eb.Row(1).Sum()>tolerance && Deposits_e(1)>tolerance)
  {
    std::cerr<<"Sum error Deposits_e in period " << t << endl;
    Errors << "\n Sum error Deposits_e in period " << t << endl;
  }
  deviation=fabs((GB_cb(1)+GB_b.Row(1).Sum()-GB(1))/GB(1));
  if(deviation>tolerance)
  {
    std::cerr<<"Sum error GB in period " << t << endl;
    Errors << "\n Sum error GB in period " << t << endl;
  }
  deviation=fabs((Deposits_1.Row(1).Sum()+Deposits_2.Row(1).Sum()+Deposits_hb.Row(1).Sum()+Deposits_eb.Row(1).Sum()-Deposits.Row(1).Sum())/Deposits.Row(1).Sum());
  if(deviation>tolerance && Deposits.Row(1).Sum()>tolerance)
  {
    std::cerr<<"Sum error Deposits in period " << t << endl;
    Errors << "\n Sum error Deposits in period " << t << endl;
  }
  deviation=fabs((Reserves(1)-Reserves_b.Row(1).Sum())/Reserves_b.Row(1).Sum());
  if(deviation>tolerance && Reserves_b.Row(1).Sum()>tolerance && Reserves(1)>tolerance)
  {
    std::cerr<<"Sum error Reserves in period " << t << endl;
    Errors << "\n Sum error Reserves in period " << t << endl;
  }
  deviation=fabs((Advances(1)-Advances_b.Row(1).Sum())/Advances_b.Row(1).Sum());
  if(deviation>tolerance && Advances_b.Row(1).Sum()>tolerance && Advances(1)>tolerance)
  {
    std::cerr<<"Sum error Advances in period " << t << endl;
    Errors << "\n Sum error Advances in period " << t << endl;
  }
  deviation=fabs((Loans_2.Row(1).Sum()+Loans_e(1)-Loans_b.Row(1).Sum())/Loans_b.Row(1).Sum());
  if(deviation>tolerance && Loans_b.Row(1).Sum()>tolerance && Loans_2.Row(1).Sum()>tolerance)
  {
    std::cerr<<"Sum error Loans in period " << t << endl;
    Errors << "\n Sum error Loans in period " << t << endl;
  }
  
  sumDepEn=0;
  for(i=0; i<ene_tecs.size(); i++)
  {
    tech=ene_tecs[i];
    sumDepEn+=deposits_e_mult_tech[tech];
  }         
  deviation=fabs((sumDepEn-Deposits_e(1))/Deposits_e(1));
  if(deviation>tolerance && Deposits_e(1)>tolerance && sumDepEn>tolerance)
  {
    std::cerr<<"Sum error Deposits_e by technology in period " << t << endl;
    Errors << "\n Sum error Deposits_e by technology in period " << t << endl;
  }

}

void ADJUSTSTOCKS(void)
{
  deviation=0;
  for (i=1; i<=NB; i++)
  {
    prior(i)=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits(1,i)-Advances_b(1,i);
  }
  prior_cb=Advances(1)+GB_cb(1)-Reserves(1);
  
  if(Advances(1)<=0 || Advances_b.Row(1).Sum()<=0)
  {
    Advances(1)=0;
    for(i=1; i<=NB; i++)
    {
      Advances_b(1,i)=0;
    }
  }
  
  if(Reserves(1)<=0 || Reserves_b.Row(1).Sum()<=0)
  {
    Reserves(1)=0;
    for(i=1; i<=NB; i++)
    {
      Reserves_b(1,i)=0;
    }
  }

  if(GB(1)<0)
  {
    GB_cb(1)=GB(1);
    for (i=1; i<=NB; i++)
    {
      GB_b(1,i)=0;
    }
  }

  if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)<tolerance && GB_b.Row(1).Sum()<tolerance)
  {
    GB_cb(1)=GB(1);
    for (i=1; i<=NB; i++)
    {
      GB_b(1,i)=0;
    }
  }

  if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)>tolerance && GB_b.Row(1).Sum()>tolerance)
  {
    if(GB_b.Row(1).Sum()>0)
    {
      for (i=1; i<=NB; i++)
      {
        ShareBonds(i)=GB_b(1,i)/GB_b.Row(1).Sum();
      }
    }
    else
    {
      for (i=1; i<=NB; i++)
      {
        ShareBonds(i)=(NL_1(i)+NL_2(i))/(N1+N2);
      }
    }
  }
  
  for(i=1; i<=NB; i++)
  {
    if(Deposits_h(1)>0 && Deposits_hb.Row(1).Sum()>0)
    {
      DepositShare_h(i)=Deposits_hb(1,i)/Deposits_hb.Row(1).Sum();
    }
    else
    {
      DepositShare_h(i)=(NL_2(i)+NL_1(i))/(N1+N2);
    }

    if(Deposits_e(1)>0 && Deposits_eb.Row(1).Sum()>0)
    {
      DepositShare_e(i)=Deposits_eb(1,i)/Deposits_eb.Row(1).Sum();
    }
    else
    {
      DepositShare_e(i)=(NL_2(i)+NL_1(i))/(N1+N2);
    }

    if(Reserves(1)>0)
    {
      ShareReserves(i)=Reserves_b(1,i)/Reserves_b.Row(1).Sum();
    }

    if(Advances(1)>0)
    {
      ShareAdvances(i)=Advances_b(1,i)/Advances_b.Row(1).Sum();
    }
  }
  
  for (i=1; i<=NB; i++)
  {
    Loans_b(1,i)=0;
    Deposits(1,i)=0;

    for(j=1; j<=N2; j++)
    {
      if(BankMatch_2(j,i)==1)
      {
        Loans_b(1,i)+=Loans_2(1,j);
        Deposits(1,i)+=Deposits_2(1,j);
      }
    }
    for(j=1; j<=N1; j++)
    {
      if(BankMatch_1(j,i)==1)
      {
        Deposits(1,i)+=Deposits_1(1,j);
      }
    }

    if(Deposits_h(1)>0)
    {
      Deposits_hb(1,i)=DepositShare_h(i)*Deposits_h(1);
      Deposits(1,i)+=(DepositShare_h(i)*Deposits_h(1));
    }
    else
    {
      Deposits_hb(1,i)=0;
    }

    if(Deposits_e(1)>0)
    {
      Deposits_eb(1,i)=DepositShare_e(i)*Deposits_e(1);
      Deposits(1,i)+=(DepositShare_e(i)*Deposits_e(1));
    }
    else
    {
      Deposits_eb(1,i)=0;
    }

    if(Loans_e(1)>0)
    {
      Loans_b(1,i)+=(DepositShare_e(i)*Loans_e(1));
    }

    if(Reserves(1)>0)
    {
      Reserves_b(1,i)=ShareReserves(i)*Reserves(1);
    }
    else
    {
      Reserves_b(1,i)=0;
    }

    if(Advances(1)>0)
    {
      Advances_b(1,i)=ShareAdvances(i)*Advances(1);
    }
    else
    {
      Advances_b(1,i)=0;
    }

    if(GB(1)>0 && fabs(GB_cb(1)-GB(1))/GB(1)>tolerance)
    {
      GB_b(1,i)=ShareBonds(i)*(GB(1)-GB_cb(1));
    }

    post=Loans_b(1,i)+Reserves_b(1,i)+GB_b(1,i)-Deposits(1,i)-Advances_b(1,i);
    Adjustment(i)=post-prior(i);
    deviation+=fabs(Adjustment(i));
  }

  post_cb=Advances(1)+GB_cb(1)-Reserves(1);
  Adjustment_cb=post_cb-prior_cb;
  deviation+=fabs(Adjustment_cb);

  if(deviation/GDP_n(1)>tolerance)
  {
    std::cerr << "\n\n ERROR: Adjustment in stocks exceeds tolerance in period " << t << endl;
    Errors << "\n Adjustment in stocks exceeds tolerance in period " << t << endl;
  }  
}

void SFC_CHECK(void)
{
  //Calculate the sectoral balances
  Balance_h=Wages+Benefits+InterestDeposits_h+Dividends(1)+TransferFuel-Taxes_h-Consumption-FirmTransfers+govTranfers;
  Balance_1=Sales1.Sum()+InterestDeposits_1.Sum()+FirmTransfers_1-Wages_1.Sum()-EnergyPayments_1.Sum()-Dividends_1.Sum()-Taxes_1.Sum()-Taxes_CO2_1.Sum();
  Balance_2=Sales2.Sum()+InterestDeposits_2.Sum()+FirmTransfers_2-Wages_2.Sum()-Investment_2.Sum()-LoanInterest_2.Sum()-EnergyPayments_2.Sum()-Dividends_2.Sum()-Taxes_2.Sum()-Taxes_CO2_2.Sum();
  Balance_b=LoanInterest.Sum()+InterestBonds_b.Sum()+InterestReserves_b.Sum()+Bailout_b.Row(1).Sum()+BankTransfer-InterestDeposits.Sum()-Taxes_b.Sum()-InterestAdvances_b.Sum()-Dividends_b.Sum();
  Balance_e=EnergyPayments+InterestDeposits_e-Wages_en-Dividends_e-Loan_interest_e-Taxes_CO2_e-FuelCost;
  Balance_cb=InterestBonds_cb+InterestAdvances-InterestReserves-TransferCB;
  Balance_g=Taxes+TransferCB+Taxes_CO2(1)-InterestBonds-Bailout-EntryCosts-G;
  Balance_f=FuelCost-TransferFuel;

  //Sectoral balances should sum to zero
  BalanceSum=Balance_h+Balance_1+Balance_2+Balance_b+Balance_e+Balance_cb+Balance_g+Balance_f;
  //Deviation needs to be scaled somehow since model variables (and hence possibly deviations due to rounding) will grow over time
  deviation=fabs(BalanceSum)/(fabs(Balance_h)+fabs(Balance_1)+fabs(Balance_2)+fabs(Balance_b)+fabs(Balance_e)+fabs(Balance_cb)+fabs(Balance_g)+fabs(Balance_f));
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Sectoral balances do not sum to zero in period " << t << endl;
    Errors << "\n Sectoral balances do not sum to zero in period " << t << endl;
  }

  //Compare stock and flow measures of bank net worth
  for(i=1; i<=NB; i++)
  {
    NW_b(1,i)+=Adjustment(i);
    if(NW_b(1,i)<=0 && Bank_active(i)==1)
    {
      std::cerr << "\n\n ERROR: NW of active bank " << i << " is negative in period " << t << endl;
      Errors << "\n NW of active bank " << i << " is negative in period " << t << endl;
    }
    NW_b_c(i)=Loans_b(1,i)+GB_b(1,i)+Reserves_b(1,i)-Deposits(1,i)-Advances_b(1,i);
  }
  deviation=fabs((NW_b_c.Sum()-NW_b.Row(1).Sum())/NW_b_c.Sum());
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for BANKS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for BANKS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of K-firm net worth
  for(i=1; i<=N1; i++)
  {
    Balances_1(i)=Sales1(i)+InterestDeposits_1(i)-Wages_1(i)-EnergyPayments_1(i)-Dividends_1(i)-Taxes_1(i)-Taxes_CO2_1(i);
    NW_1(1,i)=Deposits_1(1,i);
    NW_1_c(i)=NW_1(2,i)+Balances_1(i)+baddebt_1(i)+Injection_1(i);
  }
  deviation=fabs((NW_1_c.Sum()-NW_1.Row(1).Sum())/NW_1_c.Sum());
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for K-FIRMS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for K-FIRMS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of C-firm net worth
  for(i=1; i<=N2; i++)
  {
    NW_2(1,i)=CapitalStock(1,i)+deltaCapitalStock(1,i)+Inventories(1,i)+Deposits_2(1,i)-Loans_2(1,i);
    NW_2_c(i)=NW_2(2,i)+Pi2(i)+baddebt_2(i)+Injection_2(i)-Dividends_2(i)-Taxes_2(i)-Loss_Capital(i)-Loss_Inventories(i);
  }

  deviation=fabs((NW_2_c.Sum()-NW_2.Row(1).Sum())/NW_2_c.Sum());
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for C-FIRMS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for C-FIRMS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of Household net worth
  NW_h(1)=Deposits_h(1);
  NW_h_c=NW_h(2)+Balance_h+Deposits_recovered_1+Deposits_recovered_2;
  deviation=fabs((NW_h(1)-NW_h_c)/NW_h_c);
  if(deviation>tolerance && NW_h(1)>0 && NW_h_c>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for HOUSEHOLDS are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for HOUSEHOLDS are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of CB net worth
  NW_cb(1)=GB_cb(1)+Advances(1)-Reserves(1)-Deposits_fuel_cb(1);
  NW_cb_c=NW_cb(2)+Balance_cb+Adjustment_cb;
  deviation=fabs((NW_cb(1)-NW_cb_c)/NW_cb_c);
  if(deviation>tolerance && fabs(NW_cb(1)/GDP_n(1))>tolerance && fabs(NW_cb_c/GDP_n(1))>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the CENTRAL BANK are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the CENTRAL BANK are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of government net worth
  NW_gov(1)=-GB(1);
  NW_gov_c=NW_gov(2)+Balance_g;
  deviation=fabs((NW_gov(1)-NW_gov_c)/NW_gov_c);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the GOVERNMENT are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the GOVERNMENT are not consistent in period " << t << endl;
  }

  //Compare stock and flow measures of Energy sector net worth
  if(flag_WITCH_on==1)
  {
    NW_e(1)=0;
    CapitalStock_e(1)=0;
    for(i=0; i<ene_tecs.size(); i++)
    {
      tech=ene_tecs[i];
      NW_e(1)+=deposits_e_mult_tech[tech]+capitalStock_e_mult_tech[tech]-loans_e_mult_tech[tech];
      CapitalStock_e(1)+=capitalStock_e_mult_tech[tech];
    }      
    NW_e_c=NW_e(2)+Balance_e+(CapitalStock_e(1)-CapitalStock_e(2))+BadDebt_e;
    deviation=fabs((NW_e(1)-NW_e_c)/NW_e_c);
    if(deviation>tolerance)
    {
      std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
      Errors << "\n Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
    }
  }
  else
  {
    NW_e(1)=Deposits_e(1)+CapitalStock_e(1);
    NW_e_c=NW_e(2)+Balance_e+CapitalStock_e(1)-CapitalStock_e(2);
    deviation=fabs((NW_e(1)-NW_e_c)/NW_e_c);
    if(deviation>tolerance)
    {
      std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
      Errors << "\n Stock and flow measures of net worth for the ENERGY SECTOR are not consistent in period " << t << endl;
    }
  }

  NW_f(1)=Deposits_fuel(1);
  NW_f_c=Deposits_fuel(2)+FuelCost-TransferFuel;
  deviation=fabs((NW_f(1)-NW_f_c)/NW_f_c);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Stock and flow measures of net worth for the FOSSIL FUEL SECTOR are not consistent in period " << t << endl;
    Errors << "\n Stock and flow measures of net worth for the FOSSIL FUEL SECTOR are not consistent in period " << t << endl;
  }

  //Sum of all sectoral net worths should be equal to nominal value of tangible assets in the economy
  NWSum=NW_h(1)+NW_1.Row(1).Sum()+NW_2.Row(1).Sum()+NW_b.Row(1).Sum()+NW_e(1)+NW_cb(1)+NW_gov(1)+NW_f(1);
  RealAssets=CapitalStock.Row(1).Sum()+deltaCapitalStock.Row(1).Sum()+Inventories.Row(1).Sum()+CapitalStock_e(1);
  deviation=fabs((NWSum-RealAssets)/RealAssets);
  if(deviation>tolerance)
  {
    std::cerr << "\n\n ERROR: Aggregate net worth not equal to tangible assets in period " << t << endl;
    Errors << "\n Aggregate net worth not equal to tangible assets in period " << t << endl;
  }
}

void OVERBOOST(void)
{
  //Reset t0 to shorten time taken to iterate over technology arrays
  t00=t0;
  flag=0;
  for (tt=t00; tt<=t && flag == 0; tt++)
  {
    for (i=1; i<=N1 && flag == 0; i++)
    {
      for (j=1; j<=N2 && flag == 0; j++)
      {
        if (g[tt-1][i-1][j-1] > 0 || gtemp[tt-1][i-1][j-1] > 0)
          flag=1;
      }
    }
    if (flag == 0)
    {
      t0++;
    }
  }
}

void UPDATE(void)
{
  if (cpi(2)<=0)
	{
		std::cerr << "\n\n ERROR: cpi(t-1)<=0 in period " << t << endl;
    Errors << "\n cpi(t-1)<=0 in period " << t << endl;
    exit(EXIT_FAILURE);
	}
	if (Am(2)<=0)
	{
		std::cerr << "\n\n ERROR: Am(t-1)=0 in period " << t << endl;
    Errors << "\n Am(t-1)=0 in period " << t << endl;
    exit(EXIT_FAILURE);
	}

  //Update mark-up in the energy sector
  if(t>1)
  {
    dw2=kappa*dw2+(1-kappa)*w(1)/w(2);
    mi_en*=dw2;
    CF_ge*=dw2;
    pf*=dw2;
  }

  //Update lagged variables needed in next period
  Taxes_CO2(2)=Taxes_CO2(1);
  Deposits_h(2)=Deposits_h(1);
  Deposits_e(2)=Deposits_e(1);
  Deposits_fuel(2)=Deposits_fuel(1);
  Deposits_fuel_cb(2)=Deposits_fuel_cb(1);
  NW_f(2)=NW_f(1);
  GB_cb(2)=GB_cb(1);
  GB(2)=GB(1);
  Advances(2)=Advances(1);
  Reserves(2)=Reserves(1);
  CapitalStock_e(2)=CapitalStock_e(1);
  Loans_e(2)=Loans_e(1);
  NW_h(2)=NW_h(1);
  NW_e(2)=NW_e(1);
  NW_gov(2)=NW_gov(1);
  NW_cb(2)=NW_cb(1);
  Dividends(2)=Dividends(1);
  U(2)=U(1);
  w(2)=w(1);
  Em2(2)=Em2(1);
  ProfitCB(2)=ProfitCB(1);
  c_en(2)=c_en(1);
  GDP_r(2)=GDP_r(1);
  GDP_n(2)=GDP_n(1);
  cpi(5)=cpi(4);
  cpi(4)=cpi(3);
  cpi(3)=cpi(2);
  cpi(2)=cpi(1);
  if(t==2)
  {
    cpi_init=cpi(1);
    GDP_init=GDP_n(1);
  }
  Am(2)=Am(1);
  Am_en(2)=Am_en(1);

  for (i=1; i<=N1; i++)
  {
    Deposits_1(2,i)=Deposits_1(1,i);
    RD(2,i)=RD(1,i);
    f1(2,i)=f1(1,i);
    NW_1(2,i)=NW_1(1,i);
    S1_temp(2,i)=S1_temp(1,i);
  }

  for (j=1; j<=N2; j++)
  {
    Deposits_2(2,j)=Deposits_2(1,j);
    Loans_2(2,j)=Loans_2(1,j); 
    DebtService_2(2,j)=DebtService_2(1,j);
    f2(3,j)=f2(2,j);
    f2(2,j)=f2(1,j);
    D2(2,j)=D2(1,j);
    N(2,j)=N(1,j);
    Inventories(2,j)=Inventories(1,j);
    EI(2,j)=EI(1,j);
    deltaCapitalStock(2,j)=deltaCapitalStock(1,j);
    S2(2,j)=S2(1,j);    
    S2_temp(2,j)=S2_temp(1,j);
    mu2(2,j)=mu2(1,j);
    CapitalStock(2,j)=CapitalStock(1,j);
    NW_2(2,j)=NW_2(1,j);
  }

  //Machines held get a period older. Doing this for every firm at once rather
  //than inside the loop above runs down contiguous memory instead of across it
  //once per firm; each entry is still incremented exactly once.
  //Every entry is aged, not only the ones a firm still holds. An age is read
  //in three places and each asks for it only where the firm's count of that
  //machine is positive, and the three places a count rises from zero all set
  //the age at the same entry, so ageing an entry nobody holds changes no value
  //that is ever read. What it saves is the machine counts: this used to read
  //all 52 million of them a run to decide, and now reads none.
  for (tt=t0; tt<=t; tt++)
  {
    int* ages=age[tt-1][0].p;
    const int n_entries=N1*N2;
    for (int entry=0; entry<n_entries; entry++)
    {
      ages[entry]=ages[entry]+1;
    }
  }

	for (i=1; i<=NB; i++)
  {
    fB(2,i)=fB(1,i);
    Deposits(2,i)=Deposits(1,i);
    Deposits_hb(2,i)=Deposits_hb(1,i);
    Deposits_eb(2,i)=Deposits_eb(1,i);
    GB_b(2,i)=GB_b(1,i);
    Loans_b(2,i)=Loans_b(1,i);
    Advances_b(2,i)=Advances_b(1,i);
    Reserves_b(2,i)=Reserves_b(1,i);
    NW_b(2,i)=NW_b(1,i);
  } 
}

///////////WRITE OUTPUT/////////////////////////////

void WRITEPROD(void)
{
  //When fulloutput==1, save individual productivity values
  if(fulloutput ==1)
	{
		ofstream inv_prodall1(filename4,ios::app);
		inv_prodall1.setf(ios::fixed);
		inv_prodall1.precision(4);
		inv_prodall1.setf(ios::right);
		
    if (t>1)
    {
			inv_prodall1 << "\n";
    }

    for (i=1; i<=N1; i++)
		{
			inv_prodall1.width(20);
			inv_prodall1 << A1p(i);
		}
		inv_prodall1.close();

    ofstream inv_prodall1_en(filename6,ios::app);
    inv_prodall1_en.setf(ios::fixed);
    inv_prodall1_en.precision(4);
    inv_prodall1_en.setf(ios::right);
    if (t>1)
    {
      inv_prodall1_en << "\n";
    }

    for (i=1; i<=N1; i++)
    {
      inv_prodall1_en.width(20);
      inv_prodall1_en << A1p_en(i);
    }
    inv_prodall1_en.close();
          
    ofstream inv_prodall1_ef(filename8,ios::app);
    inv_prodall1_ef.setf(ios::fixed);
    inv_prodall1_ef.precision(4);
    inv_prodall1_ef.setf(ios::right);
    if (t>1)
    {
      inv_prodall1_ef << "\n";
    }

    for (i=1; i<=N1; i++)
    {
      inv_prodall1_ef.width(20);
      inv_prodall1_ef << A1p_ef(i);
    }
    inv_prodall1_ef.close();

		ofstream inv_prodall2(filename5,ios::app);
		inv_prodall2.setf(ios::fixed);
		inv_prodall2.precision(4);
		inv_prodall2.setf(ios::right);

		if (t>1)
    {
			inv_prodall2 << "\n";
    }

    for (j=1; j<=N2; j++)
		{
			inv_prodall2.width(20);
			inv_prodall2 << A2(j);
    }
		inv_prodall2.close();

    ofstream inv_prodall2_en(filename7,ios::app);
    inv_prodall2_en.setf(ios::fixed);
    inv_prodall2_en.precision(4);
    inv_prodall2_en.setf(ios::right);
    
    if (t>1)
    {
      inv_prodall2_en << "\n";
    }

    for (j=1; j<=N2; j++)
    {
      inv_prodall2_en.width(20);
      inv_prodall2_en << A2_en(j);
    }
    inv_prodall2_en.close();
    
    ofstream inv_prodall2_ef(filename9,ios::app);
    inv_prodall2_ef.setf(ios::fixed);
    inv_prodall2_ef.precision(4);
    inv_prodall2_ef.setf(ios::right);

    if (t>1)
    {
      inv_prodall2_ef << "\n";
    }

    for (j=1; j<=N2; j++)
    {
      inv_prodall2_ef.width(20);
      inv_prodall2_ef << A2_ef(j);
    }
    inv_prodall2_ef.close();
	}

	if(fulloutput ==1 && t>=201) 
	{
		ofstream inv_prod2(filename3,ios::app);
		inv_prod2.setf(ios::fixed);
		inv_prod2.precision(4);
		inv_prod2.setf(ios::right);
		
    if (t>201)
    {
			inv_prod2 << "\n";
    }

    for (j=1; j<=N2; j++)
		{
			inv_prod2.width(20);
			A2scr=log(A2(j))-A_mi;
			inv_prod2 << A2scr;
    }
		inv_prod2.close();

		ofstream inv_prod1(filename2,ios::app);
		inv_prod1.setf(ios::fixed);
		inv_prod1.precision(4);
		inv_prod1.setf(ios::right);
		
    if (t>201)
			inv_prod1 << "\n";
		
    for (i=1; i<=N1; i++)
		{
			inv_prod1.width(20);
			if (nclient(i) >= 1)
			{
				A1scr=log(A1p(i))-A1_mi;
				inv_prod1 << A1scr;
			}
			else
			{
				inv_prod1 << "NaN";
			}
    }
		inv_prod1.close();
  }
}

void WRITEDEB(void)									
{
  if(fulloutput ==1)
	{
		// When fulloutput==1, save stock of loans of all individual C-firms
    ofstream inv_deball2(filename13,ios::app);
		inv_deball2.setf(ios::fixed);
		inv_deball2.precision(4);
		inv_deball2.setf(ios::right);
		if (t>1)
    {
			inv_deball2 << "\n";
    }
    for (j=1; j<=N2; j++)
		{
			inv_deball2.width(60);
      if(exiting_2(j)==1)
      {
        inv_deball2 << 0;
      }
      else
      {
        inv_deball2 << Loans_2(1,j);
      }
		}
		inv_deball2.close();
	}
}

void SAVE(void)	
{
  if(flag_validation==0)
  {
    if(fulloutput==1)
    {
      WRITENW();
    }
    ofstream inv_res(filename1,ios::app);
    inv_res.setf(ios::fixed);
    inv_res.precision(10);
    inv_res.setf(ios::right);
    inv_res.width(60);
    inv_res << t;                                                                               // 1
    inv_res.width(60);
    inv_res << GDP_r(1);                                                                        // 2
    inv_res.width(60);
    inv_res << Consumption_r;                                                                   // 3
    inv_res.width(60);
    inv_res << Investment_r;                                                                    // 4
    inv_res.width(60);
    inv_res << 1-U(1);                                                                          // 5
    inv_res.width(60);
    inv_res << cpi(1)/cpi(5);                                                                   // 6
    inv_res.width(60);
    inv_res << Emiss1_TOT+Emiss2_TOT+Emiss_en;                                                  // 7  
    inv_res.width(60);
    inv_res << D_en_TOT;                                                                        // 8
    inv_res.width(60);
    inv_res << LS;                                                                              // 9
    inv_res.width(60);
    inv_res << K.Sum();                                                                         // 10
    inv_res.width(60);
    inv_res << RDin.Sum()+RDim.Sum();                                                           // 11
    inv_res.width(60);
    inv_res << Inn.Sum()+Imm.Sum();                                                             // 12
    inv_res.width(60);
    inv_res << Am_en(1);                                                                        // 13
    inv_res.width(60);
    inv_res << Am_a;                                                                            // 14
    inv_res.width(60);
    inv_res << exit_marketshare2.Sum();                                                         // 15
    inv_res.width(60);
    inv_res << exit_payments2.Sum();                                                            // 16
    inv_res.width(60);
    inv_res << exit_equity2.Sum();                                                              // 17
    inv_res.width(60);
    inv_res << exiting_1.Sum();                                                                 // 18
    inv_res.width(60);
    inv_res << Bailout;                                                                         // 19
    inv_res.width(60);
    inv_res << baddebt_b.Sum()/(GDP_n(1)*4);                                                    // 20
    inv_res.width(60);
    inv_res << counter_bankfailure;                                                             // 21
    inv_res.width(60);
    inv_res << CapitalStock.Row(1).Sum()/(GDP_n(1)*4);                                          // 22
    inv_res.width(60);
    inv_res << NW_cb(1)/(GDP_n(1)*4);                                                           // 23
    inv_res.width(60);
    inv_res << NW_h(1)/(GDP_n(1)*4);                                                            // 24
    inv_res.width(60);
    inv_res << NW_2.Row(1).Sum()/(GDP_n(1)*4);                                                  // 25
    inv_res.width(60);
    inv_res << NW_b.Row(1).Sum()/(GDP_n(1)*4);                                                  // 26
    inv_res.width(60);
    inv_res << BankProfits.Sum();                                                               // 27
    inv_res.width(60);
    inv_res << Loans_2.Row(1).Sum()/(GDP_n(1)*4);                                               // 28
    inv_res.width(60);
    inv_res << CreditDemand.Sum()/BaselBankCredit.Sum();                                        // 29
    inv_res.width(60);
    inv_res << NW_gov(1)/(GDP_n(1)*4);                                                          // 30
    inv_res.width(60);
    inv_res << NW_e(1)/(GDP_n(1)*4);                                                            // 31
    inv_res.width(60);
    inv_res << NW_1.Row(1).Sum()/(GDP_n(1)*4);                                                  // 32
    inv_res.width(60);
    inv_res << exiting_1_payments.Sum();                                                        // 33
    inv_res.width(60);
    inv_res << cpi(1);                                                                          // 34
    inv_res.width(60);
    inv_res << kpi;                                                                             // 35
    inv_res.width(60);
    inv_res << w(1)/cpi(1);                                                                     // 36
    inv_res.width(60);
    inv_res << Am2;                                                                             // 37
    inv_res.width(60);
    inv_res << Am1;                                                                             // 38
    inv_res.width(60);
    inv_res << GDP_n(1);                                                                        // 39
    inv_res.width(60);
    inv_res << H2;                                                                              // 40
    inv_res.width(60);
    inv_res << Consumption;                                                                     // 41
    inv_res.width(60);
    inv_res << t_CO2_en;                                                                        // 42
    inv_res.width(60);
    inv_res << H1;                                                                              // 43
    inv_res.width(60);
    inv_res << Emiss1_TOT;                                                                      // 44
    inv_res.width(60);
    inv_res << Emiss2_TOT;                                                                      // 45
    inv_res.width(60);
    inv_res << Emiss_en;                                                                        // 46
    inv_res.width(60);
    inv_res << Tmixed(1);                                                                       // 47
    inv_res.width(60);
    inv_res << EnergyPayments;                                                                  // 48
    inv_res.width(60);
    inv_res << FuelCost/GDP_n(1);                                                               // 49
    inv_res.width(60);
    inv_res << r;                                                                               // 50
    inv_res.width(60);
    inv_res << Deposits_e(1)/(GDP_n(1)*4);                                                      // 51
    inv_res.width(60);
    inv_res << CapitalStock_e(1)/(GDP_n(1)*4);                                                  // 52
    inv_res.width(60);
    inv_res << Pitot1;                                                                          // 53
    inv_res.width(60);
    inv_res << Pitot2;                                                                          // 54
    inv_res.width(60);
    inv_res << ProfitEnergy;                                                                    // 55
    inv_res.width(60);
    inv_res << Wages/cpi(1);                                                                    // 56
    inv_res.width(60);
    inv_res << Dividends_1.Sum();                                                               // 57
    inv_res.width(60);
    inv_res << Dividends_2.Sum();                                                               // 58
    inv_res.width(60);
    inv_res << I_loss.Sum();                                                                    // 59
    inv_res.width(60);
    inv_res << c_en(1);                                                                         // 60
    inv_res.width(60);
    inv_res << Deposits_2.Row(1).Sum()/(GDP_n(1)*4);                                            // 61
    inv_res.width(60);
    inv_res << Deposits_1.Row(1).Sum()/(GDP_n(1)*4);                                            // 62
    inv_res.width(60);
    inv_res << (Deposits_2.Row(1).Sum()+Deposits_1.Row(1).Sum())/(GDP_n(1)*4);                  // 63
    inv_res.width(60);
    inv_res << K_ge/(K_ge+K_de);                                                                // 64
    inv_res.width(60);
    inv_res << C_loss_share/N2r;                                                                // 65
    inv_res.width(60);
    inv_res << (Dividends_2.Sum()+Dividends_1.Sum())/cpi(1);                                    // 66
    inv_res.width(60);
    inv_res << Dividends_e/cpi(1);                                                              // 67
    inv_res.width(60);
    inv_res << Dividends_b.Sum()/cpi(1);                                                        // 68
    inv_res.width(60);
    inv_res << (NW_1.Row(1).Sum()+NW_2.Row(1).Sum())/(GDP_n(1)*4);                              // 69
    inv_res.width(60);
    inv_res << GB(1)/(GDP_n(1)*4);                                                              // 70
    inv_res.width(60);
    inv_res << K_loss.Sum();                                                                    // 71
    inv_res.width(60);
    inv_res << C_loss.Sum();                                                                    // 72
    inv_res.width(60);
    inv_res << NW_f(1)/(GDP_n(1)*4);                                                            // 73
    inv_res.width(60);
    inv_res << LDexp_en*w(2);                                                                   // 74
    inv_res.width(60);
    inv_res << Pitot1/GDP_n(1);                                                                 // 75
    inv_res.width(60);
    inv_res << Pitot2/GDP_n(1);                                                                 // 76
    inv_res.width(60);
    inv_res << (BankProfits.Sum())/GDP_n(1);                                                    // 77
    inv_res.width(60);
    inv_res << ProfitEnergy/GDP_n(1);                                                           // 78
    inv_res.width(60);
    inv_res << FuelCost/GDP_n(1);                                                               // 79
    inv_res.width(60);
    inv_res << Wages/GDP_n(1);                                                                  // 80
    inv_res.width(60);
    inv_res << exit_payments2.Sum() + exit_equity2.Sum();                                       // 81
    inv_res.width(60);
    inv_res << r_deb.Sum()/NB;                                                                  // 82
    inv_res.width(60);
    inv_res << t_CO2<<endl;                                                                     // 83
    inv_res.close();
  }
  else
  {
    ofstream inv_val1(filename14,ios::app);
    inv_val1.setf(ios::fixed);
    inv_val1.precision(10);
    inv_val1.setf(ios::right);
    inv_val1.width(60);
    inv_val1 << t;                                         // 1
    inv_val1.width(60); 
    inv_val1 << GDP_r(1);                                  // 2
    inv_val1.width(60);
    inv_val1 << Consumption_r;                             // 3
    inv_val1.width(60);
    inv_val1 << Investment_r;                              // 4
    inv_val1.width(60);
    inv_val1 << 1-U(1);                                    // 5
    inv_val1.width(60);
    inv_val1 << cpi(1)/cpi(5);                             // 6
    inv_val1.width(60);
    inv_val1 << Emiss1_TOT+Emiss2_TOT+Emiss_en;            // 7  
    inv_val1.width(60);
    inv_val1 << D_en_TOT;                                  // 8
    inv_val1.width(60);
    inv_val1 << RD.Row(1).Sum()+RD_en_de+RD_en_ge;         // 9
    inv_val1.width(60);
    inv_val1 << Loans_2.Row(1).Sum();                      // 10
    inv_val1.width(60);
    inv_val1 << baddebt_b.Sum();                           // 11
    inv_val1.width(60);
    inv_val1 << Deposits_h(1)/(GDP_n(1)*4);                // 12
    inv_val1.width(60);
    inv_val1 << Deposits_e(1)/(GDP_n(1)*4);                // 13
    inv_val1.width(60);
    inv_val1 << Deposits_1.Row(1).Sum()/(GDP_n(1)*4);      // 14
    inv_val1.width(60);
    inv_val1 << Deposits_2.Row(1).Sum()/(GDP_n(1)*4);      // 15
    inv_val1.width(60);
    inv_val1 << GB(1)/(GDP_n(1)*4);                        // 16
    inv_val1.width(60);
    inv_val1 << GB_cb(1)/(GDP_n(1)*4);                     // 17
    inv_val1.width(60);
    inv_val1 << Loans_2.Row(1).Sum()/(GDP_n(1)*4);         // 18
    inv_val1.width(60);
    inv_val1 << Advances(1)/(GDP_n(1)*4);                  // 19
    inv_val1.width(60);
    inv_val1 << Reserves(1)/(GDP_n(1)*4);                  // 20
    inv_val1.width(60);
    inv_val1 << NW_2.Row(1).Sum()/(GDP_n(1)*4);            // 21
    inv_val1.width(60);
    inv_val1 << NW_b.Row(1).Sum()/(GDP_n(1)*4);            // 22
    inv_val1.width(60);
    inv_val1 << NW_1.Row(1).Sum()/(GDP_n(1)*4);            // 23
    inv_val1.width(60);
    inv_val1 << CapitalStock.Row(1).Sum()/(GDP_n(1)*4);    // 24
    inv_val1.width(60);
    inv_val1 << EnergyPayments/(GDP_n(1));                 // 25
    inv_val1.width(60);
    inv_val1 << GDP_n(1);                                  // 26
    inv_val1.width(60);
    inv_val1 << Am(1);                                     // 27
    inv_val1.width(60);
    inv_val1 << Am_en(1);                                  // 28
    inv_val1.width(60);
    inv_val1 << Am1;                                       // 29
    inv_val1.width(60);
    inv_val1 << Am2;                                       // 30
    inv_val1.width(60);
    inv_val1 << cpi(1);                                    // 31
    inv_val1.width(60);
    inv_val1 << kpi;                                       // 32
    inv_val1.width(60);
    inv_val1 << exit_marketshare2.Sum();                   // 33
    inv_val1.width(60);
    inv_val1 << exit_payments2.Sum();                      // 34
    inv_val1.width(60);
    inv_val1 << exit_equity2.Sum();                        // 35
    inv_val1.width(60);
    inv_val1 << exiting_1.Sum();                           // 36
    inv_val1.width(60);
    inv_val1 << Bailout/GDP_n(1);                          // 37
    inv_val1.width(60);
    inv_val1 << baddebt_b.Sum()/GDP_n(1);                  // 38
    inv_val1.width(60);
    inv_val1 << counter_bankfailure;                       // 39
    inv_val1.width(60);
    inv_val1 << Emiss_en;                                  // 40
    inv_val1.width(60);
    inv_val1 << Tmixed(1);                                 // 41
    inv_val1.width(60);
    inv_val1 << H1;                                        // 42
    inv_val1.width(60);
    inv_val1 << H2;                                        // 43
    inv_val1.width(60);
    inv_val1 << exit_payments2.Sum()+exit_equity2.Sum();   // 44 
    inv_val1.width(60); 
    inv_val1 << NW_e(1)/(GDP_n(1)*4);                      // 45
    inv_val1.width(60);
    inv_val1 << NW_gov(1)/(GDP_n(1)*4);                    // 46
    inv_val1.width(60);
    inv_val1 << NW_f(1)/(GDP_n(1)*4);                      // 47
    inv_val1.width(60);
    inv_val1 << Balance_h/(GDP_n(1));                      // 48
    inv_val1.width(60);
    inv_val1 << Balance_1/(GDP_n(1));                      // 49
    inv_val1.width(60);
    inv_val1 << Balance_2/(GDP_n(1));                      // 50
    inv_val1.width(60);
    inv_val1 << Balance_b/(GDP_n(1));                      // 51
    inv_val1.width(60);
    inv_val1 << Balance_e/(GDP_n(1));                      // 52
    inv_val1.width(60);
    inv_val1 << Balance_cb/(GDP_n(1));                     // 53
    inv_val1.width(60);
    inv_val1 << Balance_g/(GDP_n(1));                      // 54
    inv_val1.width(60);
    inv_val1 << Balance_f/(GDP_n(1))<<endl;                // 55
    inv_val1.close();

    if(t%50==0)
    {
      ofstream inv_val2(filename15,ios::app);
      inv_val2.setf(ios::fixed);
      inv_val2.precision(10);
      inv_val2.setf(ios::right);
      if (t>50)
      {
        inv_val2 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val2.width(60);
        if(exiting_2(j)==1)
        {
          inv_val2 << "NA";
        }
        else
        {
          inv_val2 << S2_temp(1,j);
        }
      }
      inv_val2.close();

      ofstream inv_val3(filename16,ios::app);
      inv_val3.setf(ios::fixed);
      inv_val3.precision(10);
      inv_val3.setf(ios::right);
      if (t>50)
      {
        inv_val3 << "\n";
      }
      for (i=1; i<=N1; i++)
      {
        inv_val3.width(60);
        if(exiting_1(i)==1)
        {
          inv_val3 << "NA";
        }
        else
        {
          inv_val3 << S1_temp(1,i);
        }
      }
      inv_val3.close();

      ofstream inv_val4(filename17,ios::app);
      inv_val4.setf(ios::fixed);
      inv_val4.precision(10);
      inv_val4.setf(ios::right);
      if (t>50)
      {
        inv_val4 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val4.width(60);
        A2scr=log(A2(j))-A_mi;
        inv_val4 << A2scr;
      }
      inv_val4.close();
      
      ofstream inv_val5(filename18,ios::app);
      inv_val5.setf(ios::fixed);
      inv_val5.precision(10);
      inv_val5.setf(ios::right);
      if (t>50)
      {
        inv_val5 << "\n";
      }
      for (i=1; i<=N1; i++)
      {
        inv_val5.width(60);
        A1scr=log(A1p(i))-A1_mi;
        inv_val5 << A1scr;
      }
      inv_val5.close();

      ofstream inv_val6(filename19,ios::app);
      inv_val6.setf(ios::fixed);
      inv_val6.precision(10);
      inv_val6.setf(ios::right);
      if (t>50)
      {
        inv_val6 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val6.width(60);
        inv_val6 << log(A2_en(j))-A2_en_mi;
      }
      inv_val6.close();

      ofstream inv_val7(filename20,ios::app);
      inv_val7.setf(ios::fixed);
      inv_val7.precision(10);
      inv_val7.setf(ios::right);
      if (t>50)
      {
        inv_val7 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val7.width(60);
        inv_val7 << log(A2_ef(j))-A2_ef_mi;
      }
      inv_val7.close();

      ofstream inv_val8(filename21,ios::app);
      inv_val8.setf(ios::fixed);
      inv_val8.precision(10);
      inv_val8.setf(ios::right);
      if (t>50)
      {
        inv_val8 << "\n";
      }
      for (i=1; i<=N1; i++)
      {
        inv_val8.width(60);
        inv_val8 << log(A1p_en(i))-A1_en_mi;
      }
      inv_val8.close();

      ofstream inv_val9(filename22,ios::app);
      inv_val9.setf(ios::fixed);
      inv_val9.precision(10);
      inv_val9.setf(ios::right);
      if (t>50)
      {
        inv_val9 << "\n";
      }
      for (i=1; i<=N1; i++)
      {
        inv_val9.width(60);
        inv_val9 << log(A1p_ef(i))-A1_ef_mi;
      }
      inv_val9.close();

      ofstream inv_val10(filename23,ios::app);
      inv_val10.setf(ios::fixed);
      inv_val10.precision(10);
      inv_val10.setf(ios::right);
      if (t>50)
      {
        inv_val10 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val10.width(60);
        if(exiting_2(j)==1)
        {
          inv_val10 << "NA";
        }
        else
        {
          inv_val10 << I(j);
        }
      }
      inv_val10.close();

      ofstream inv_val11(filename24,ios::app);
      inv_val11.setf(ios::fixed);
      inv_val11.precision(10);
      inv_val11.setf(ios::right);
      if (t>50)
      {
        inv_val11 << "\n";
      }
      for (j=1; j<=N2; j++)
      {
        inv_val11.width(60);
        inv_val11 << (S2_temp(1,j)-S2_temp(2,j))/S2_temp(2,j);
      }
      inv_val11.close();

      ofstream inv_val12(filename25,ios::app);
      inv_val12.setf(ios::fixed);
      inv_val12.precision(10);
      inv_val12.setf(ios::right);
      if (t>50)
      {
        inv_val12 << "\n";
      }
      for (i=1; i<=N1; i++)
      {
        inv_val12.width(60);
        inv_val12 << (S1_temp(1,i)-S1_temp(2,i))/S1_temp(2,i);
      }
      inv_val12.close();
    }
  }
}

void WRITENW(void)									
{
	if(fulloutput ==1)
	{
		//When fulloutput==1, save the individual net worths of all K-firms, C-firms and Banks
    ofstream inv_nwall1(filename10,ios::app);
		inv_nwall1.setf(ios::fixed);
		inv_nwall1.precision(4);
		inv_nwall1.setf(ios::right);
		if (t>1)
    {
			inv_nwall1 << "\n";
    }
    for (i=1; i<=N1; i++)
		{
			inv_nwall1.width(60);
			inv_nwall1 << NW_1(1,i);
		}
		inv_nwall1.close();

		ofstream inv_nwall2(filename11,ios::app);
		inv_nwall2.setf(ios::fixed);
		inv_nwall2.precision(4);
		inv_nwall2.setf(ios::right);
		if (t>1)
    {
			inv_nwall2 << "\n";
    }
    for (j=1; j<=N2; j++)
		{
			inv_nwall2.width(60);
			inv_nwall2 << NW_2(1,j);
		}
		inv_nwall2.close();

		ofstream inv_nwall3(filename12,ios::app);
		inv_nwall3.setf(ios::fixed);
		inv_nwall3.precision(4);
		inv_nwall3.setf(ios::right);
		if (t>1)
    {
			inv_nwall3 << "\n";
    }
		for (i=1; i<=NB; i++)
		{
			inv_nwall3.width(60);
			inv_nwall3 << NW_b(1,i);
		}
		inv_nwall3.close();
	}
}

///////////GENERATE OUTPUT FOLDERS, FILES & NAMES/////////////////////
//These functions generate the directories for saving output and the names of the .txt files in which model output is saved
int make_directory(const char* name)
{
  #ifdef __linux__
    return mkdir(name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  #elif __APPLE__
    return mkdir(name,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  #else
    return mkdir(name);
  #endif
}


void FOLDERS(char *path)
{
  //Create a folder called "output" in the same directory as the executable
  //Also create a subdirectory of "output" called "errors"
  std::string outstr(path);
  for(j=outstr.length(); j>0; j--)
  {
    if(outstr[j-1]=='/')
    {
      break;
    }
    else
    {
      outstr.pop_back();
    }
  }
  outstr+="output";
  std::string errstr=outstr + "/errors";
  //One past each length, for the terminator strcpy writes.
  char out_dir[outstr.length()+1];
  char err_dir[errstr.length()+1];
  strcpy(out_dir,outstr.c_str());
  strcpy(err_dir,errstr.c_str());
	const int out_fol=make_directory(out_dir);
  const int err_fol=make_directory(err_dir);
}

void GENFILEERRORS(char *path, const char *se, char const* desc)
{
  //File to save error messages
  strcpy(errorfilename,path);							
  char* name_error=strcat(errorfilename,se);			
  name_error=strcat(errorfilename,desc);			
  strcat(errorfilename,".txt");		
}

void GENFILEYMC(char *path, const char *s1, char const* desc)				    
{
	//Standard file with selected aggregate variables
  strcpy(filename1,path);								
	char* name1=strcat(filename1,s1);
	name1=strcat(filename1,desc);
	strcat(filename1,".txt");
}

void GENFILEPROD1(char *path, const char *s2, char const* desc)				    
{
	//File to save log deviation of K-Firm productivity from mean
  strcpy(filename2,path);								
	char* name2=strcat(filename2,s2);
	name2=strcat(filename2,desc);
	strcat(filename2,".txt");
}

void GENFILEPROD2(char *path, const char *s3, char const* desc)				    
{
	//File to save log deviation of C-Firm productivity from mean
  strcpy(filename3,path);								
	char* name3=strcat(filename3,s3);
	name3=strcat(filename3,desc);
	strcat(filename3,".txt");
}

void GENFILEPRODALL1(char *path, const char *s4, char const* desc)				    
{
	//File to save untransformed K-Firm productivity
  strcpy(filename4,path);								
	char* name4=strcat(filename4,s4);
	name4=strcat(filename4,desc);
	strcat(filename4,".txt");
}

void GENFILEPRODALL2(char *path, const char *s5, char const* desc)				    
{
	//File to save untransformed C-Firm productivity
  strcpy(filename5,path);								
	char* name5=strcat(filename5,s5);
	name5=strcat(filename5,desc);					
	strcat(filename5,".txt");
}

void GENFILEPRODALL1_en(char *path, const char *s6, char const* desc)           
{
  //File to save energy efficiency of K-Firms
  strcpy(filename6,path);               
  char* name6=strcat(filename6,s6);
  name6=strcat(filename6,desc);         
  strcat(filename6,".txt");
}

void GENFILEPRODALL2_en(char *path, const char *s7, char const* desc)           
{
  //File to save energy efficiency of C-Firms
  strcpy(filename7,path);               
  char* name7=strcat(filename7,s7);
  name7=strcat(filename7,desc);          
  strcat(filename7,".txt");
}


void GENFILEPRODALL1_ef(char *path, const char *s8, char const* desc)           
{
  //File to save environmental friendliness of K-Firms
  strcpy(filename8,path);               
  char* name8=strcat(filename8,s8); 
  name8=strcat(filename8,desc);          
  strcat(filename8,".txt");
}

void GENFILEPRODALL2_ef(char *path, const char *s9, char const* desc)           
{
  //File to save environmental friendliness of C-Firms
  strcpy(filename9,path);               
  char* name9=strcat(filename9,s9);
  name9=strcat(filename9,desc);          
  strcat(filename9,".txt");
}

void GENFILENWALL1(char *path, const char *s10, char const* desc)
{
	//File to save net worth of K-firms
  strcpy(filename10,path);								
	char* name10=strcat(filename10,s10);
  name10=strcat(filename10,desc);
  strcat(filename10,".txt");
}

void GENFILENWALL2(char *path, const char *s11, char const* desc)
{
	//File to save net worth of C-firms
  strcpy(filename11,path);								
	char* name11=strcat(filename11,s11);
  name11=strcat(filename11,desc);
  strcat(filename11,".txt");
}

void GENFILENWALL3(char *path, const char *s12, char const* desc)
{
	//File to save net worth of banks
  strcpy(filename12,path);								
	char* name12=strcat(filename12,s12);
  name12=strcat(filename12,desc);
  strcat(filename12,".txt");
}

void GENFILEDEBALL2(char *path, const char *s13, char const* desc)
{
	//File to save debt of C-Firms
  strcpy(filename13,path);								
	char* name13=strcat(filename13,s13);
  name13=strcat(filename13,desc);
  strcat(filename13,".txt");
}

void GENFILEVALIDATION1(char *path, const char *s14, char const* seednumber)
{
	//File to save aggregate variables needed for validation
  strcpy(filename14,path);								
	char* name14=strcat(filename14,s14);
  name14=strcat(filename14,"_");
  name14=strcat(filename14,seednumber);
  strcat(filename14,".txt");
}

void GENFILEVALIDATION2(char *path, const char *s15, char const* seednumber)
{
	//File to save C-Firms' sales for validation
  strcpy(filename15,path);								
	char* name15=strcat(filename15,s15);
  name15=strcat(filename15,"_");
  name15=strcat(filename15,seednumber);
  strcat(filename15,".txt");
}

void GENFILEVALIDATION3(char *path, const char *s16, char const* seednumber)
{
	//File to save K-Firms' sales for validation
  strcpy(filename16,path);								
	char* name16=strcat(filename16,s16);
  name16=strcat(filename16,"_");
  name16=strcat(filename16,seednumber);
  strcat(filename16,".txt");
}

void GENFILEVALIDATION4(char *path, const char *s17, char const* seednumber)
{
	//File to save log deviation of C-Firms' productivity from mean for validation
  strcpy(filename17,path);								
	char* name17=strcat(filename17,s17);
  name17=strcat(filename17,"_");
  name17=strcat(filename17,seednumber);
  strcat(filename17,".txt");
}

void GENFILEVALIDATION5(char *path, const char *s18, char const* seednumber)
{
	//File to save log deviation of K-Firms' productivity from mean for validation
  strcpy(filename18,path);								
	char* name18=strcat(filename18,s18);
  name18=strcat(filename18,"_");
  name18=strcat(filename18,seednumber);
  strcat(filename18,".txt");
}

void GENFILEVALIDATION6(char *path, const char *s19, char const* seednumber)
{
	//File to save log deviation of C-Firms' energy efficiency from mean for validation
  strcpy(filename19,path);								
	char* name19=strcat(filename19,s19);
  name19=strcat(filename19,"_");
  name19=strcat(filename19,seednumber);
  strcat(filename19,".txt");
}

void GENFILEVALIDATION7(char *path, const char *s20, char const* seednumber)
{
	//File to save log deviation of C-Firms' environmental friendliness from mean for validation
  strcpy(filename20,path);								
	char* name20=strcat(filename20,s20);
  name20=strcat(filename20,"_");
  name20=strcat(filename20,seednumber);
  strcat(filename20,".txt");
}

void GENFILEVALIDATION8(char *path, const char *s21, char const* seednumber)
{
	//File to save log deviation of K-Firms' energy efficiency from mean for validation
  strcpy(filename21,path);								
	char* name21=strcat(filename21,s21);
  name21=strcat(filename21,"_");
  name21=strcat(filename21,seednumber);
  strcat(filename21,".txt");
}

void GENFILEVALIDATION9(char *path, const char *s22, char const* seednumber)
{
	//File to save log deviation of K-Firms' environmental friendliness from mean for validation
  strcpy(filename22,path);								
	char* name22=strcat(filename22,s22);
  name22=strcat(filename22,"_");
  name22=strcat(filename22,seednumber);
  strcat(filename22,".txt");
}

void GENFILEVALIDATION10(char *path, const char *s23, char const* seednumber)
{
	//File to save C-Firms' investment for validation
  strcpy(filename23,path);								
	char* name23=strcat(filename23,s23);
  name23=strcat(filename23,"_");
  name23=strcat(filename23,seednumber);
  strcat(filename23,".txt");
}

void GENFILEVALIDATION11(char *path, const char *s24, char const* seednumber)
{
	//File to save growth rate of individual C-Firm sales for validation
  strcpy(filename24,path);								
	char* name24=strcat(filename24,s24);
  name24=strcat(filename24,"_");
  name24=strcat(filename24,seednumber);
  strcat(filename24,".txt");
}

void GENFILEVALIDATION12(char *path, const char *s25, char const* seednumber)
{
  //File to save growth rate of individual K-Firm sales for validation
	strcpy(filename25,path);								
	char* name25=strcat(filename25,s25);
  name25=strcat(filename25,"_");
  name25=strcat(filename25,seednumber);
  strcat(filename25,".txt");
}

//This function generates the actual output files
void INTFILE(void) 
{
  ofstream Errors(errorfilename);

  if(flag_validation==0)
  {
    ofstream inv_res(filename1);
  }
  else
  {
    //Aggregate output for validation
    ofstream inv_val1(filename14);
    //C-Firm sales for validation
    ofstream inv_val2(filename15);
    //K-Firm sales for validation
    ofstream inv_val3(filename16);
    //C-Firm productivity distribution for validation
    ofstream inv_val4(filename17);
    //K-Firm productivity distribution for validation
    ofstream inv_val5(filename18);
    //C-Firm energy efficiency distribution for validation
    ofstream inv_val6(filename19);
    //C-Firm env. friendliness distribution for validation
    ofstream inv_val7(filename20);
    //K-Firm energy efficiency distribution for validation
    ofstream inv_val8(filename21);
    //K-Firm env. friendliness distribution for validation
    ofstream inv_val9(filename22);
    //C-Firm investment for validation
    ofstream inv_val10(filename23);
    //C-Firm sales growth for validation
    ofstream inv_val11(filename24);
    //K-Firm sales growth for validation
    ofstream inv_val12(filename25);
  }
  
  if(fulloutput==1)
  {
    // When fulloutput==1, a larger set of model variables will be saved
    //Log K-Firm productivity deviation from mean
    ofstream inv_prod1(filename2);
    //Log C-Firm productivity deviation from mean
    ofstream inv_prod2(filename3);
    //K-Firm productivity
    ofstream inv_prodall1(filename4);
    //C-Firm productivity
    ofstream inv_prodall2(filename5);
    //K-Firm energy efficiency
    ofstream inv_prodall1_en(filename6);
    //C-Firm energy efficiency
    ofstream inv_prodall2_en(filename7);
    //K-Firm environmental friendliness
    ofstream inv_prodall1_ef(filename8);
    //C-Firm environmental friendliness
    ofstream inv_prodall2_ef(filename9);
    //K-Firm net worth
    ofstream inv_nwall1(filename10);
    //C-Firm net worth
    ofstream inv_nwall2(filename11);
    //Bank net worth
    ofstream inv_nwall3(filename12);
    //C-firm debt
    ofstream inv_deball2(filename13);
  }
  else
  {
    //Standard aggregate output file
    ofstream inv_res(filename1);
  }
}

///////////AUXILIARY/////////////////////

void catchAlarm(int sig) {
    //If the time taken to perform the run exceeds the threshold set above, abort the simulation
    
    std::cerr << "\n\n Run timed out!" << endl;
    Errors << "\n Run timed out! " << endl;
    
    exit(EXIT_FAILURE);
}

double ROUND(double x)
{
  //Rounds a double to the closest integer
  double x_floor=floor(x);
  double resto=x-x_floor;
  if (resto > 0.5) x=x_floor+1;
  else x=x_floor;
  return x;
}

void ALLOCATEBANKCUSTOMERS(void)
{
  //Used during initialisation to assign C-Firms and K-Firms to banks
  //Initialise number of C-firm customers of each bank to 0
  NL_2=0;
  double sum_NL_2;
  sum_NL_2=0;

  //re-perform the random drawing of customer numbers until the sum of the random values drawn is equal to the number of C-Firms
  while (sum_NL_2!=N2){
       sum_NL_2=0;
       for (i=1; i<=NB; i++)
        {
          //Draw from truncated pareto. pareto_a=shape parameter, pareto_k=lower bound, pareto_p=upper bound
          pareto_rv = bpareto(pareto_a, pareto_k, pareto_p);
          NL_2(i)=pareto_rv;
          sum_NL_2+=pareto_rv;
        }
  }   

  //Initialise number of C-firm customers of each bank to 0
  NL_1=0;
  double sum_NL_1;
  sum_NL_1=0;
  //re-perform the random drawing of customer numbers until the sum of the random values drawn is equal to the number of K-Firms
  while (sum_NL_1!=N1){
       sum_NL_1=0;
       for (i=1; i<=NB; i++)
        {
          //Draw from truncated pareto. pareto_a=shape parameter, pareto_k=lower bound, pareto_p=upper bound
          //Adjust upper and lower bounds to reflect the number of K-Firms relative to C-Firms
          pareto_rv = bpareto(pareto_a, min(pareto_k*N1r/N2r,0.9), ceil(pareto_p*N1r/N2r));
          NL_1(i)=pareto_rv;
          sum_NL_1+=pareto_rv;
        }
  } 
}

double bpareto(double par_a, double par_k, double par_p)
{

  double z;     // Uniform random number from 0 to 1
  double rv;    // RV to be returned

  // Pull a uniform RV (0 < z < 1)
  do
  {
    z=double(ran1(p_seed));
  }
  while ((z == 0) || (z == 1));

  // Generate the bounded Pareto rv using the inversion method
  rv = pow((pow(par_k, par_a) / (z*pow((par_k/par_p), par_a) - z + 1)), (1.0/par_a));
  // make the variable an integer
  rv=ceil(rv);

  return(rv);
}