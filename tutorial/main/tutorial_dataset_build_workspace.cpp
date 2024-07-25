#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooFormula.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooArgList.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#include "GammaComboEngine.h"




int main(int argc, char* argv[])
{

  /////////////////////////////////////////////////////////
  //
  // First, we define all observables and variables and construct the PDF.
  // It will be included in the workspace. 
  //
  /////////////////////////////////////////////////////////

  // First, define the signal peak of the mass model

  // // It consists of a gaussian to describe the signal...
  // RooRealVar mass("mass","mass", 4360., 6360., "MeV");
  // RooRealVar mean("mean","mean",5370);
  // RooRealVar sigma("sigma","sigma", 20.9);
  // RooGaussian signal_model("g","g", mass, mean, sigma);

  // // ... and of an exponential function to describe the background.
  // RooRealVar exponent("exponent","exponent", -1e-3, -1., 1.);
  // RooRealVar n_bkg("Nbkg","Nbkg", 4900, 0, 10000);
  // RooExponential background_model("background_model", "background_model", mass, exponent);
  // RooExponential bkg_only_model("bkg_only_model", "bkg_only_model", mass, exponent);
  // RooExtendPdf extended_bkg_model("extended_bkg_model", "extended_bkg_model", bkg_only_model, n_bkg);

  // // The number of signal events is related to the branching ratio via the formula <branching ratio> = <n_sig> * <normalization factor>
  // // The normalization factor is not exactly known. Instead, it has to be estimated. The estimator for the normalization factor is a global observable
  // // that constrains its value via a Gaussian Constraint. 
  // RooRealVar  norm_constant_obs( "norm_constant_glob_obs", "global observable of normalization constant", 
  //                             1e-8, 1e-20, 1e-6);     // this is the observed value, the global observable
  // RooRealVar  norm_constant("norm_constant","norm_constant", 1e-8, 1e-20,  1e-6);                   // the normalization constant is a nuisance parameter.
  // RooRealVar  norm_constant_sigma("norm_constant_sigma","norm_constant_sigma", 5e-10);
  // RooGaussian norm_constant_constraint("norm_constant_constraint","norm_constant_constraint", norm_constant_obs, norm_constant, norm_constant_sigma);

  // // Now we can build the mass model by adding the signal and background probability density functions
  // RooRealVar branchingRatio("branchingRatio", "branchingRatio", 1e-7, -1e-6,  0.0001);  // this is the branching ratio, the parameter of interest
  // RooFormulaVar n_sig("Nsig", "branchingRatio/norm_constant", RooArgList(branchingRatio, norm_constant));
  // RooExtendPdf extended_sig_model("extended_sig_model", "extended_sig_model", signal_model, n_sig);

  // RooAddPdf mass_model("mass_model","mass_model", RooArgList(signal_model, background_model), RooArgList(n_sig, n_bkg));

  // /////////////////////////////////////////////////////////
  // //
  // // Open/create the dataset with the data that will be analyzed 
  // // and set the global observables to the observed values
  // //
  // /////////////////////////////////////////////////////////

  // /////////////////////////////////////////////////////////
  // //
  // // In this tutorial, we generate our data instead of using a real dataset.
  // // We generate one mass dataset with 100 events.
  // // We also generate a single measurement of the global observable.
  // double observedValueGlobalObservable = norm_constant_constraint.generate(RooArgSet(norm_constant_obs),1)
  //                                                               ->get(0)
  //                                                               ->getRealValue("norm_constant_glob_obs");
  // //
  // RooDataSet& data = *mass_model.generate(RooArgSet(mass),5000); // we use a reference here to avoid copying the object
  // data.SetName("data"); // the name of the dataset MUST be "data" in order for the framework to identify it.
  // //
  // /////////////////////////////////////////////////////////

  // // Set the global observables to the observed values and make them constant for the fit.
  // norm_constant_obs.setVal(observedValueGlobalObservable);
  // norm_constant_obs.setConstant();

  // // The workspace must also include the fit result of an initial fit of the model to data.
  // RooFitResult& rooFitResult = *mass_model.fitTo(data,RooFit::Save(), RooFit::ExternalConstraints(RooArgSet(norm_constant_constraint)));

  // /////////////////////////////////////////////////////////
  // //
  // // Plot it all so we can see what we did
  // //
  // /////////////////////////////////////////////////////////
  
  // RooPlot* plot = mass.frame();
  // data.plotOn(plot);	
  // extended_bkg_model.plotOn(plot, RooFit::LineColor(kRed));
  // mass_model.plotOn(plot);
  // TCanvas c("c","c",1024, 768);
  // plot->Draw();
  // c.SaveAs("plots/pdf/data_and_fit_in_workspace.pdf");
  
  // rooFitResult.Print();

  // /////////////////////////////////////////////////////////
  // //
  // // We also must define some RooArgSets: 
  // //
  // /////////////////////////////////////////////////////////

  // //One contaninting the constraint PDFs,
  // RooArgSet constraint_set(norm_constant_constraint, "constraint_set");

  // //one containing the global Observables,
  // RooArgSet global_observables_set(norm_constant_obs, "global_observables_set");

  // //one containing the normal Observables (the bin variables in the datasets usually) and
  // RooArgSet dataset_observables_set(mass, "datasetObservables");


  // //one containing the parameters
  // RooArgSet parameters_set(branchingRatio, norm_constant, exponent, n_bkg, "parameters");


  // Start the Gammacombo Engine
  GammaComboEngine gc("", argc, argv);

  int bin_number = gc.getArg()->bin;
  TString decay = gc.getArg()->decay;

  // int bin_number = 3;
  // TString decay = "PIPIEE";

  // gc.run();

  std::cout<<"Decay "<<decay<<" in bin "<<bin_number<<std::endl;
  std::cout<<"Decay 2 "<<decay<<" in bin "<<bin_number<<std::endl;


  

  
  
  


  

  // int bin_number = 2;



  TString dir_plots = "/afs/cern.ch/user/a/ascarabo/work/rare-charm-d0-hhee/Optimization/CLs_plots/";
  int min_mass, max_mass;
  TString year, signal_decay, misid_name;
  int yr_idx, ch_idx, nbins, number_of_bins;
  float ratio_brem1_brem0_value;

  float BF_begin = 0, BF_end = 1e-6;
  TString bin0, bin1, name_bin;

    
  TString dilepton_bins_names[5] = {"211-525","525-565","565-950","950-1100","1100-1500"};
  TString dimuon = "211.316749";
  TString dilepton_bins[5][2] = {{dimuon,"525"},{"525","565"},{"565","950"},{"950","1100"},{"1100","1500"}};

  float dilepton_bins_eff[5][2] = {{211.32,525},{525,565},{565,950},{950,1100},{1100,1500}};

  name_bin = dilepton_bins_names[bin_number];

  
  float n_max_1516, n_max_1718, n_start;
  int nbrems = 2, nyears = 2;
  TString years[2] = {"1516","1718"};
  TString brems[2] = {"brem0","brem1"};

  std::string word;

  float BF_signal_per_bin[5];


  if(decay.Contains("PIPIEE")){
    signal_decay = "PIPIEE";
    misid_name = "PIPIPIPI";
    ch_idx = 0;
    ratio_brem1_brem0_value = 0.052;
    nbins = 5;
    n_max_1516 = 300;
    n_max_1718 = 500;
    number_of_bins = 100;
    BF_signal_per_bin[0] = 7.8 * 1e-8;
    BF_signal_per_bin[1] = 2.4 * 1e-8;
     BF_signal_per_bin[2] = 40.6 * 1e-8 ;
     BF_signal_per_bin[3] = 45.4 * 1e-8;
     BF_signal_per_bin[4] = 2.8 * 1e-8; 

    // BF_signal_per_bin[0] = 0.; 
    // BF_signal_per_bin[1] = 0.; 
    // BF_signal_per_bin[2] = 0.;
    // BF_signal_per_bin[3] = 0.;
    // BF_signal_per_bin[4] = 0.;
  }
  else if(decay.Contains("KKEE")){
    signal_decay = "KKEE";
    misid_name = "KKPIPI";
    ratio_brem1_brem0_value = 0.074;
    ch_idx = 1;
    nbins = 3;
    number_of_bins = 100;
    BF_signal_per_bin[0] = 2.6 * 1e-8;
    BF_signal_per_bin[1] = 0.7 * 1e-8;
    BF_signal_per_bin[2] = 12.0 * 1e-8 ;
    
    
  }

float option_diff = 0;

RooRealVar ratio_brem1_brem0 ("ratio_brem1_brem0","ratio_brem1_brem0",ratio_brem1_brem0_value + option_diff);
ratio_brem1_brem0.setConstant(kTRUE);
float ratio_leak_brem1_bin4_value = 0.043;

RooRealVar ratio_leak_brem1_bin4 ("ratio_leak_brem1_bin4","ratio_leak_brem1_bin4", ratio_leak_brem1_bin4_value+ option_diff);
ratio_leak_brem1_bin4.setConstant(kTRUE);

float ratio_leak_brem0_bin2_value = 0.296;
RooRealVar ratio_leak_brem0_bin2 ("ratio_leak_brem0_bin2","ratio_leak_brem0_bin2", ratio_leak_brem0_bin2_value+ option_diff);
ratio_leak_brem0_bin2.setConstant(kTRUE);

float ratio_leak_brem1_bin2_value = 0.144;
RooRealVar ratio_leak_brem1_bin2 ("ratio_leak_brem1_bin2","ratio_leak_brem1_bin2", ratio_leak_brem1_bin2_value+ option_diff);
ratio_leak_brem1_bin2.setConstant(kTRUE);

float ratio_leak_brem1_bin3_value = 0.034;
RooRealVar ratio_leak_brem1_bin3 ("ratio_leak_brem1_bin3","ratio_leak_brem1_bin3", ratio_leak_brem1_bin3_value+ option_diff);
ratio_leak_brem1_bin3.setConstant(kTRUE);


float Y_Norm_values[nbrems][nyears], Y_Norm_errors[nbrems][nyears];
std::ifstream finput_norm;
finput_norm.open("/afs/cern.ch/user/a/ascarabo/work/rare-charm-d0-hhee/Optimization/Optimization_FoM/Norm_yields.txt");

if(finput_norm.is_open()){
finput_norm >> word >> word >> word >> word
>> word >> Y_Norm_values[0][0] >> word >> Y_Norm_errors[0][0] >> word>> Y_Norm_values[1][0] >> word >> Y_Norm_errors[1][0] 
>> word >> word >> word >> word
>> word >> Y_Norm_values[0][1] >> word >> Y_Norm_errors[0][1] >> word>> Y_Norm_values[1][1] >> word >> Y_Norm_errors[1][1];
if(signal_decay == "KKEE"){
finput_norm >> word >> word >> word >> word
>> word >> Y_Norm_values[0][0] >> word >> Y_Norm_errors[0][0] >> word>> Y_Norm_values[1][0] >> word >> Y_Norm_errors[1][0] 
>> word >> word >> word >> word
>> word >> Y_Norm_values[0][1] >> word >> Y_Norm_errors[0][1] >> word>> Y_Norm_values[1][1] >> word >> Y_Norm_errors[1][1];
}
}

finput_norm.close();

std::ifstream finput_syst;
finput_syst.open ("/afs/cern.ch/user/a/ascarabo/work/rare-charm-d0-hhee/Optimization/systematics_"+signal_decay+".txt");

float syst_yield_percentage[nbins], syst_eff_percentage[nbins];
for(int i=0;i<nbins;i++){
finput_syst >> word >> syst_yield_percentage[i] >> word >> syst_eff_percentage[i];
}
finput_syst.close(); 



TH1F* histo_pull_BF_signal[nbins] ;
TH1F* histo_pull_BF_signal_SYST[nbins] ;
TH1F* histo_pull_comb[nbins][nbrems][nyears] ;
TH1F* histo_pull_eta[nbrems][nyears] ;
TH1F* histo_pull_mis_id[nbins][nyears] ;
TString jj;
RooRealVar *Y_Norm[nbins][nbrems][nyears];
RooRealVar *Y_Norm_obs[nbins][nbrems][nyears];
RooRealVar *Y_Norm_sigma[nbins][nbrems][nyears];
RooGaussian *Y_Norm_constraint[nbins][nbrems][nyears];




for(int ibin = 0; ibin < nbins; ibin++){
  jj = dilepton_bins_names[ibin];
  histo_pull_BF_signal[ibin] = new TH1F("histo_pull_BF_signal_"+jj,"histo_pull_BF_signal "+jj+";Pull (BF_signal_{fit}-BF_signal_{gen})/#sigma(BF_signal); number of toys",100,-5,5);
  histo_pull_BF_signal_SYST[ibin] = new TH1F("histo_pull_BF_signal_SYST_"+jj,"histo_pull_BF_signal_SYST "+jj+";BF_signal_{fit}; number of toys",100,BF_signal_per_bin[ibin]/10,BF_signal_per_bin[ibin]*10);
  // histo_pull_BF_signal[ibin] = new TH1F("histo_pull_BF_signal_"+jj,"histo_pull_BF_signal "+jj+";Pull (n_signal_{fit}-n_signal_{gen})/#sigma(n_signal);",20,-5,5);
  
  for(int iyear = 0; iyear < nyears; iyear++){
    histo_pull_mis_id[ibin][iyear] = new TH1F("histo_pull_mis_id_"+years[iyear]+"_"+jj,"histo_pull_mis_id_"+years[iyear]+" "+jj+";Pull (Y_mis-id_{fit}-Y_mis-id_{gen})/#sigma(Y_mis-id); number of toys",100,-5,5);
    
    for(int ibrem = 0; ibrem < nbrems; ibrem++){
      histo_pull_comb[ibin][ibrem][iyear] = new TH1F("histo_pull_comb_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"histo_pull_comb_"+brems[ibrem]+"_"+years[iyear]+" "+jj+";Pull (Y_Comb_brem0_{fit}-Y_Comb_brem0_{gen})/#sigma(Y_Comb_brem0); number of toys",100,-5,5);
      if(ibin == 0) histo_pull_eta[ibrem][iyear] = new TH1F("histo_pull_eta_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"histo_pull_eta_"+brems[ibrem]+"_"+years[iyear]+" "+jj+";Pull (Y_Eta_brem0_{fit}-Y_Eta_brem0_{gen})/#sigma(Y_Eta_brem0); number of toys",100,-5,5);


      Y_Norm[ibin][ibrem][iyear] = new RooRealVar("Y_Norm_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"Y_Norm_"+brems[ibrem]+"_"+years[iyear]+"_"+jj, Y_Norm_values[ibrem][iyear],0,2000);
      Y_Norm_obs[ibin][ibrem][iyear] = new RooRealVar("Y_Norm_obs_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"Y_Norm_obs_"+brems[ibrem]+"_"+years[iyear]+"_"+jj, Y_Norm_values[ibrem][iyear],0,2000);
      // Y_Norm_sigma[ibin][ibrem][iyear] = new RooRealVar("Y_Norm_sigma_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"Y_Norm_sigma_"+brems[ibrem]+"_"+years[iyear]+"_"+jj, Y_Norm_errors[ibrem][iyear]);
      Y_Norm_sigma[ibin][ibrem][iyear] = new RooRealVar("Y_Norm_sigma_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"Y_Norm_sigma_"+brems[ibrem]+"_"+years[iyear]+"_"+jj, Y_Norm_values[ibrem][iyear] * syst_yield_percentage[ibin]/100.);
      Y_Norm_constraint[ibin][ibrem][iyear] = new RooGaussian("Y_Norm_constraint_"+brems[ibrem]+"_"+years[iyear]+"_"+jj,"Y_Norm_constraint_"+brems[ibrem]+"_"+years[iyear]+"_"+jj, *Y_Norm_obs[ibin][ibrem][iyear],*Y_Norm[ibin][ibrem][iyear],*Y_Norm_sigma[ibin][ibrem][iyear]);
      Y_Norm[ibin][ibrem][iyear]->setConstant(kTRUE);
      Y_Norm_obs[ibin][ibrem][iyear]->setConstant(kTRUE);
    }
    
  }
}


std::ifstream finput;
finput.open ("/afs/cern.ch/work/a/ascarabo/rare-charm-d0-hhee/Efficiencies/Final_efficiencies_evaluation/efficiencies_for_FIT/EffRatio_"+signal_decay+"_for_FIT.txt");

float bin0_i,bin1_i, eff, err;
int y, brem;
float eff_ratio_value[nbins][nbrems][nyears];
while(!finput.eof()){
  finput >> y >> bin0_i >> bin1_i >> brem >> eff >> err ;
  for (int j = 0; j < nbins; j++){
    if((bin0_i == dilepton_bins_eff[j][0]) && (bin1_i == dilepton_bins_eff[j][1])){
      eff_ratio_value[j][brem][y] = eff;
    }
  }
}

finput.close();

std::ifstream finputt;
finputt.open ("/afs/cern.ch/user/a/ascarabo/work/rare-charm-d0-hhee/Optimization/BF_fits/blind_yields_"+signal_decay+".txt");

float n_signal_value[nbins][nbrems][nyears], n_comb_value[nbins][nbrems][nyears], n_misid_value[nbins][nbrems][nyears], n_eta_value[nbins][nbrems][nyears], n_leak_value[nbins][nbrems][nyears];
for(int i=0;i<nbins;i++){
finputt >> word >> word >> word >> word >> word >> word >> word >> n_comb_value[i][0][0] >> word >> word >> word >> n_misid_value[i][0][0] >> word >> word >> word >> n_comb_value[i][1][0] >> word >> word >> word >> n_eta_value[i][0][0] >> word >> word >> word >> n_eta_value[i][1][0] >> word >> word ;
finputt >> word >> word >> word >> word >> word >> word >> word >> n_comb_value[i][0][1] >> word >> word >> word >> n_misid_value[i][0][1] >> word >> word >> word >> n_comb_value[i][1][1] >> word >> word >> word >> n_eta_value[i][0][1] >> word >> word >> word >> n_eta_value[i][1][1] >> word >> word ;
}
finputt.close(); 










TString ws_dir = "/afs/cern.ch/user/a/ascarabo/work/rare-charm-d0-hhee/Optimization/workspaces_signal_pdf/";

TFile *f[nbins];
RooWorkspace *w[nbins];
RooRealVar *x_data[nbins];
RooAbsPdf *model_signal[nbins][nbrems][nyears];
// RooCBShape *model_signal[nbins][nbrems][nyears];
RooAbsPdf *model_misid[nbins][nbrems][nyears];
RooAbsPdf *model_comb[nbins][nbrems][nyears];
// RooPolynomial *model_comb[nbins][nbrems][nyears];
RooAbsPdf *model_eta[nbrems][nyears];
RooAbsPdf *model_leak[nbins][nbrems][nyears];

RooAbsPdf *model_misid_only_bkg[nbins][nbrems][nyears];
RooAbsPdf *model_comb_only_bkg[nbins][nbrems][nyears];
RooAbsPdf *model_eta_only_bkg[nbrems][nyears];
RooAbsPdf *model_leak_only_bkg[nbins][nbrems][nyears];

RooExtendPdf *emodel_signal[nbins][nbrems][nyears];
RooExtendPdf *emodel_misid[nbins][nbrems][nyears];
RooExtendPdf *emodel_comb[nbins][nbrems][nyears];
RooExtendPdf *emodel_eta[nbrems][nyears];
RooExtendPdf *emodel_leak[nbins][nbrems][nyears];

RooExtendPdf *emodel_misid_only_bkg[nbins][nbrems][nyears];
RooExtendPdf *emodel_comb_only_bkg[nbins][nbrems][nyears];
RooExtendPdf *emodel_eta_only_bkg[nbrems][nyears];
RooExtendPdf *emodel_leak_only_bkg[nbins][nbrems][nyears];

RooAddPdf *final_model[nbins][nbrems][nyears];
RooAddPdf *final_model_only_bkg[nbins][nbrems][nyears];


RooRealVar *Eff_Ratio[nbins][nbrems][nyears];
RooRealVar *Eff_Ratio_obs[nbins][nbrems][nyears];
RooRealVar *Eff_Ratio_sigma[nbins][nbrems][nyears];
RooGaussian *Eff_Ratio_constraint[nbins][nbrems][nyears];
RooFormulaVar *n_signal[nbins][nbrems][nyears];
// RooFormulaVar *n_signal_BF[nbins][nbrems][nyears];
// RooRealVar *n_signal[nbins][nbrems][nyears];
RooRealVar *BF_signal[nbins];

RooRealVar *n_comb[nbins][nbrems][nyears];
RooRealVar *n_misid_brem0[nbins][nyears];
RooFormulaVar *n_misid_brem1[nbins][nyears];
RooRealVar *n_eta[nbrems][nyears];
RooRealVar *n_leak[nbins][nbrems][nyears];

float BF_norm_value = 4.0 * 1e-6;
float BF_norm_error = 0.5477 * 1e-6;
// RooRealVar *BF_norm[nbins][nbrems][nyears];
// RooRealVar BF_norm("BF_norm","BF_norm",BF_norm_value);
RooRealVar BF_norm("BF_norm","BF_norm",BF_norm_value,1e-20,1e-5);
RooRealVar BF_norm_obs("BF_norm_obs","BF_norm_obs",BF_norm_value,1e-20,1e-5);
RooRealVar BF_norm_sigma("BF_norm_sigma","BF_norm_sigma",BF_norm_error);
RooGaussian BF_norm_constraint("BF_norm_constraint","BF_norm_constraint", BF_norm_obs, BF_norm, BF_norm_sigma);
BF_norm.setConstant(kTRUE);


// RooCategory sample("sample","sample") ;
// RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;

RooRealVar x_data_toys("D_M", "m(D^{0}) [MeV/c^{2}]",1700,2050);
// RooRealVar *x_data_toys[nbins];





for(int ibin = 0; ibin < nbins; ibin++){
  f[ibin] = new TFile(ws_dir+"ws_"+signal_decay+"_"+dilepton_bins_names[ibin]+".root");
  std::cout<<dilepton_bins_names[ibin]<<std::endl;
  w[ibin] = (RooWorkspace *)f[ibin]->Get("w");
  x_data[ibin] = w[ibin]->var("D_M");
  TString ibin_name(std::to_string(ibin));
  // if(BF_signal_per_bin[ibin] > 1e-7) BF_signal[ibin] = new RooRealVar("BF_signal_bin"+ibin_name,"BF_signal_bin"+ibin_name, BF_signal_per_bin[ibin],-1e-10,BF_signal_per_bin[ibin]*20);
  // else BF_signal[ibin] = new RooRealVar("BF_signal_bin"+ibin_name,"BF_signal_bin"+ibin_name, BF_signal_per_bin[ibin],-5 * 1e-8,BF_signal_per_bin[ibin]*20);

  BF_signal[ibin] = new RooRealVar("BF_signal_bin"+ibin_name,"BF_signal_bin"+ibin_name, BF_signal_per_bin[ibin],0,1e-5);
  


  for(int ibrem = 0; ibrem < nbrems; ibrem++){
    for(int iyear = 0; iyear < nyears; iyear++){

        model_signal[ibin][ibrem][iyear] = w[ibin]->pdf("model_signal_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
        model_misid[ibin][ibrem][iyear] = w[ibin]->pdf("model_misid_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
        model_comb[ibin][ibrem][iyear] = w[ibin]->pdf("model_comb_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);

        model_misid_only_bkg[ibin][ibrem][iyear] = w[ibin]->pdf("model_misid_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
        model_misid_only_bkg[ibin][ibrem][iyear]->SetName("model_misid_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
        model_comb_only_bkg[ibin][ibrem][iyear] = w[ibin]->pdf("model_comb_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
        model_comb_only_bkg[ibin][ibrem][iyear]->SetName("model_comb_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);


        Eff_Ratio[ibin][ibrem][iyear] = new RooRealVar("Eff_Ratio_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"Eff_Ratio_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, eff_ratio_value[ibin][ibrem][iyear],0,5);
        Eff_Ratio_obs[ibin][ibrem][iyear] = new RooRealVar("Eff_Ratio_obs"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"Eff_Ratio_obs"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, eff_ratio_value[ibin][ibrem][iyear],0,5);
        Eff_Ratio_sigma[ibin][ibrem][iyear] = new RooRealVar("Eff_Ratio_sigma"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"Eff_Ratio_sigma"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, eff_ratio_value[ibin][ibrem][iyear] * syst_eff_percentage[ibin]/100.);
        Eff_Ratio_constraint[ibin][ibrem][iyear] = new RooGaussian("Eff_Ratio_constraint"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"Eff_Ratio_constraint"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, *Eff_Ratio_obs[ibin][ibrem][iyear],*Eff_Ratio[ibin][ibrem][iyear],*Eff_Ratio_sigma[ibin][ibrem][iyear]);

        Eff_Ratio[ibin][ibrem][iyear]->setConstant(kTRUE);
        Eff_Ratio_obs[ibin][ibrem][iyear]->setConstant(kTRUE);

        // BF_norm[ibin][ibrem][iyear] = new RooRealVar("BF_norm_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"BF_norm_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,BF_norm_value);
        // BF_norm[ibin][ibrem][iyear]->setConstant(kTRUE);

        // n_signal[ibin][ibrem][iyear] = new RooFormulaVar("n_signal_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"x[0]/x[1]*x[2]*x[3]",RooArgList(*BF_signal[ibin], *BF_norm[ibin][ibrem][iyear], *Y_Norm[ibin][ibrem][iyear], *Eff_Ratio[ibin][ibrem][iyear]));
        n_signal[ibin][ibrem][iyear] = new RooFormulaVar("n_signal_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"x[0]/x[1]*x[2]*x[3]",RooArgList(*BF_signal[ibin], BF_norm, *Y_Norm[ibin][ibrem][iyear], *Eff_Ratio[ibin][ibrem][iyear]));
        std::cout<<" BF_signal[ibin] = "<< BF_signal[ibin]->getValV()  <<" Y_Norm[ibin][ibrem][iyear] = "<< Y_Norm[ibin][ibrem][iyear]->getValV() <<" Eff_Ratio[ibin][ibrem][iyear] = "<< Eff_Ratio[ibin][ibrem][iyear]->getValV() <<std::endl;

      
        n_comb[ibin][ibrem][iyear] = new RooRealVar("n_comb_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"n_comb_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, 10, 0, 300);
        if(ibrem == 0) {
          n_misid_brem0[ibin][iyear] = new RooRealVar("n_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"n_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, 10, 0, 200);
          emodel_misid[ibin][ibrem][iyear]  = new RooExtendPdf("emodel_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_misid[ibin][ibrem][iyear],*n_misid_brem0[ibin][iyear]);
          emodel_misid_only_bkg[ibin][ibrem][iyear]  = new RooExtendPdf("emodel_misid_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_misid_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_misid_only_bkg[ibin][ibrem][iyear],*n_misid_brem0[ibin][iyear]);
          }
        if(ibrem == 1) {
          n_misid_brem1[ibin][iyear] = new RooFormulaVar("n_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"x[0] * x[1]",RooArgList(*n_misid_brem0[ibin][iyear], ratio_brem1_brem0));
          emodel_misid[ibin][ibrem][iyear]  = new RooExtendPdf("emodel_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_misid_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_misid[ibin][ibrem][iyear],*n_misid_brem1[ibin][iyear]);
          emodel_misid_only_bkg[ibin][ibrem][iyear]  = new RooExtendPdf("emodel_misid_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_misid_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_misid_only_bkg[ibin][ibrem][iyear],*n_misid_brem1[ibin][iyear]);
        }
        
        
        emodel_signal[ibin][ibrem][iyear]=  new RooExtendPdf("emodel_signal_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_signal_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_signal[ibin][ibrem][iyear],*n_signal[ibin][ibrem][iyear] );
        emodel_comb[ibin][ibrem][iyear] = new RooExtendPdf("emodel_comb_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_comb_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_comb[ibin][ibrem][iyear],*n_comb[ibin][ibrem][iyear]);
        emodel_comb_only_bkg[ibin][ibrem][iyear] = new RooExtendPdf("emodel_comb_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_comb_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_comb_only_bkg[ibin][ibrem][iyear],*n_comb[ibin][ibrem][iyear]);

        if(ibin == 0){
          model_eta[ibrem][iyear] = w[ibin]->pdf("model_eta_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
          model_eta_only_bkg[ibrem][iyear] = w[ibin]->pdf("model_eta_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
          model_eta_only_bkg[ibrem][iyear]->SetName("model_eta_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_"+dilepton_bins_names[ibin]);
          n_eta[ibrem][iyear] = new RooRealVar("n_eta_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"n_eta_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name, 10, 0, 100);
          emodel_eta[ibrem][iyear]=  new RooExtendPdf("emodel_eta_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_eta_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_eta[ibrem][iyear],*n_eta[ibrem][iyear] );
          emodel_eta_only_bkg[ibrem][iyear]=  new RooExtendPdf("emodel_eta_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"emodel_eta_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,*model_eta_only_bkg[ibrem][iyear],*n_eta[ibrem][iyear] );  

          final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear], *emodel_eta[ibrem][iyear]));
          final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear], *emodel_eta_only_bkg[ibrem][iyear]));

           }
        else if(ibin == 2 && signal_decay == "PIPIEE" && ibrem == 1) {
          //model leak ro -> phi
          model_leak[3][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_phi");
          model_leak_only_bkg[3][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_phi");
          model_leak_only_bkg[3][ibrem][iyear]->SetName("model_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_phi");
          n_leak[3][ibrem][iyear] =  new RooRealVar("n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin3","n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin3",10, 0, 100);

          emodel_leak[3][ibrem][iyear]=  new RooExtendPdf("emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin3","emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin3",*model_leak[3][ibrem][iyear],*n_leak[3][ibrem][iyear] );
          emodel_leak_only_bkg[3][ibrem][iyear]=  new RooExtendPdf("emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin3","emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin3",*model_leak_only_bkg[3][ibrem][iyear],*n_leak[3][ibrem][iyear] );
          }
        else if(ibin == 3 && signal_decay == "PIPIEE") {
          
          if(ibrem == 0) n_leak[2][ibrem][iyear] =  new RooRealVar("n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2","n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2",10, 0, 100);
          if(ibrem == 1) {
            //model leak phi -> ro
            n_leak[2][ibrem][iyear] =  new RooRealVar("n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2","n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2",10, 0, 100);
            //model leak phi -> hm
            model_leak[4][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_hm");
            model_leak_only_bkg[4][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_hm");
            model_leak_only_bkg[4][ibrem][iyear]->SetName("model_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_hm");
            n_leak[4][ibrem][iyear] =  new RooRealVar("n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin4","n_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin4",10, 0, 100);
            emodel_leak[4][ibrem][iyear]=  new RooExtendPdf("emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin4","emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin4",*model_leak[4][ibrem][iyear],*n_leak[4][ibrem][iyear] );
            emodel_leak_only_bkg[4][ibrem][iyear]=  new RooExtendPdf("emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin4","emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin4",*model_leak[4][ibrem][iyear],*n_leak[4][ibrem][iyear] );
            }

          //model leak phi -> ro
          model_leak[2][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_ro");
          model_leak_only_bkg[2][ibrem][iyear] = w[ibin]->pdf("model_leak_"+brems[ibrem]+"_"+years[iyear]+"_ro");
          model_leak_only_bkg[2][ibrem][iyear]->SetName("model_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_ro");
          //model leak phi -> ro
          emodel_leak[2][ibrem][iyear]=  new RooExtendPdf("emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2","emodel_leak_"+brems[ibrem]+"_"+years[iyear]+"_bin2",*model_leak[2][ibrem][iyear],*n_leak[2][ibrem][iyear] );
          emodel_leak_only_bkg[2][ibrem][iyear]=  new RooExtendPdf("emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin2","emodel_leak_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin2",*model_leak_only_bkg[2][ibrem][iyear],*n_leak[2][ibrem][iyear] );
          //final model r/o with leak
          final_model[2][ibrem][iyear] = NULL;
          delete final_model[2][ibrem][iyear];
          final_model_only_bkg[2][ibrem][iyear] = NULL;
          delete final_model_only_bkg[2][ibrem][iyear];
          final_model[2][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin2","final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin2",RooArgList(*emodel_signal[2][ibrem][iyear], *emodel_comb[2][ibrem][iyear], *emodel_misid[2][ibrem][iyear], *emodel_leak[2][ibrem][iyear]));
          final_model_only_bkg[2][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin2","final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin2",RooArgList( *emodel_comb_only_bkg[2][ibrem][iyear], *emodel_misid_only_bkg[2][ibrem][iyear], *emodel_leak_only_bkg[2][ibrem][iyear]));

          //final model phi in PIPIEE
          final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear]));
          final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear]));
          }
        else if(ibin == 3 && signal_decay == "PIPIEE" && ibrem == 1) {
          final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear], *emodel_leak[ibin][ibrem][iyear]));
          final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear], *emodel_leak_only_bkg[ibin][ibrem][iyear]));

          }
        else if(ibin == 4 && signal_decay == "PIPIEE" && ibrem == 1) {
          //remove leak in high mass
          // final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear], *emodel_leak[ibin][ibrem][iyear]));
          // final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear], *emodel_leak_only_bkg[ibin][ibrem][iyear]));
          final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear]));
          final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear]));

          }
        else {

          final_model[ibin][ibrem][iyear] = new RooAddPdf("final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList(*emodel_signal[ibin][ibrem][iyear], *emodel_comb[ibin][ibrem][iyear], *emodel_misid[ibin][ibrem][iyear]));
          final_model_only_bkg[ibin][ibrem][iyear] = new RooAddPdf("final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"final_model_only_bkg_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgList( *emodel_comb_only_bkg[ibin][ibrem][iyear], *emodel_misid_only_bkg[ibin][ibrem][iyear]));

        }

        n_signal_value[ibin][ibrem][iyear] = n_signal[ibin][ibrem][iyear]->getValV();
        n_comb[ibin][ibrem][iyear]->setVal(n_comb_value[ibin][ibrem][iyear]);
        if(ibrem == 0) n_misid_brem0[ibin][iyear]->setVal(n_misid_value[ibin][ibrem][iyear]);
        if(ibin == 0) n_eta[ibrem][iyear]->setVal(n_eta_value[ibin][ibrem][iyear]);
        
        // sample.defineType(brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name);
        // if(ibin == 2 && signal_decay == "PIPIEE") std::cout<<"wait to add leak model"<<std::endl;
        // else {
        // simPdf.addPdf(*final_model[ibin][ibrem][iyear],brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name) ;
        // }
        // //adding back leak model once it is defined
        // if(ibin == 3 && signal_decay == "PIPIEE") simPdf.addPdf(*final_model[2][ibrem][iyear],brems[ibrem]+"_"+years[iyear]+"_bin2") ;
    }//year

  }//brem

  
}//bins



RooCategory sample_CLs("sample_CLs","sample_CLs") ;
sample_CLs.defineType("brem0_1516") ;
sample_CLs.defineType("brem1_1516") ;
sample_CLs.defineType("brem0_1718") ;
sample_CLs.defineType("brem1_1718") ;

// Construct a simultaneous pdf using category sample as index
  RooSimultaneous mass_model("mass_model","mass_model",sample_CLs) ;
  // Associate model with the physics state and model_ctl with the control state
  mass_model.addPdf(*final_model[bin_number][0][0],"brem0_1516") ;
  mass_model.addPdf(*final_model[bin_number][1][0],"brem1_1516") ;
  mass_model.addPdf(*final_model[bin_number][0][1],"brem0_1718") ;
  mass_model.addPdf(*final_model[bin_number][1][1],"brem1_1718") ;

// Construct a simultaneous pdf using category sample as index
  RooSimultaneous extended_bkg_model("extended_bkg_model","extended_bkg_model",sample_CLs) ;
  // Associate model with the physics state and model_ctl with the control state
  extended_bkg_model.addPdf(*final_model_only_bkg[bin_number][0][0],"brem0_1516") ;
  extended_bkg_model.addPdf(*final_model_only_bkg[bin_number][1][0],"brem1_1516") ;
  extended_bkg_model.addPdf(*final_model_only_bkg[bin_number][0][1],"brem0_1718") ;
  extended_bkg_model.addPdf(*final_model_only_bkg[bin_number][1][1],"brem1_1718") ;


RooDataSet *data_full[nbins][nbrems][nyears];

RooDataSet *data_signal[nbins][nbrems][nyears];
RooDataSet *data_misid[nbins][nbrems][nyears];
RooDataSet *data_comb[nbins][nbrems][nyears];
RooDataSet *data_eta[nbrems][nyears];
RooDataSet *data_leak[nbins][nbrems][nyears];

std::map<std::string, RooDataSet*> map_combData;

int events_SIGNAL, events_SIGNAL_TOTAL, events_MISID, events_COMB, events_LEAK, events_ETA;
int events_SIGNAL_TOY[nbins][nbrems][nyears];



for(int ibin = 0; ibin < nbins; ibin++){

  if(ibin == bin_number){
    TString ibin_name(std::to_string(ibin));
    for(int iyear = 0; iyear < nyears; iyear++){
      for(int ibrem = 0; ibrem < nbrems; ibrem++){
        TString name_combined_data = brems[ibrem]+"_"+years[iyear];
          
        data_full[ibin][ibrem][iyear] = new RooDataSet("data_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,"data_"+brems[ibrem]+"_"+years[iyear]+"_bin"+ibin_name,RooArgSet(x_data_toys));
        events_SIGNAL = round(n_signal_value[ibin][ibrem][iyear]);
        std::cout<<"events_SIGNAL "<<events_SIGNAL <<std::endl;
        std::cout<<"n_signal_value[ibin][ibrem][iyear] "<<n_signal_value[ibin][ibrem][iyear] <<std::endl;
        events_COMB = round(n_comb_value[ibin][ibrem][iyear]);
        if(ibrem == 0) events_MISID = round(n_misid_value[ibin][ibrem][iyear]);
        else events_MISID = round(events_MISID * ratio_brem1_brem0_value);
        if(events_SIGNAL !=0 ) { 
          data_signal[ibin][ibrem][iyear] = model_signal[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_SIGNAL) ;
          data_full[ibin][ibrem][iyear]->append(*data_signal[ibin][ibrem][iyear]);
          }
        if(events_COMB !=0) {
          data_comb[ibin][ibrem][iyear] = model_comb[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_COMB) ;
          data_full[ibin][ibrem][iyear]->append(*data_comb[ibin][ibrem][iyear]);
          }
        if(events_MISID !=0) { 
          data_misid[ibin][ibrem][iyear] = model_misid[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_MISID) ;
          data_full[ibin][ibrem][iyear]->append(*data_misid[ibin][ibrem][iyear]);
          
          }
        if(ibin == 0){ 
          events_ETA = round(n_eta_value[ibin][ibrem][iyear]);
            if(events_ETA !=0 ){
            data_eta[ibrem][iyear] = model_eta[ibrem][iyear]->generate(RooArgSet(x_data_toys),events_ETA) ;
            data_full[ibin][ibrem][iyear]->append(*data_eta[ibrem][iyear]);
            }
          }
        if(ibin == 2 && signal_decay == "PIPIEE") {
          // events_LEAK = r_leak[ibin][ibrem][iyear]->Poisson(n_leak_value[ibin][ibrem][iyear]);
          if(ibrem == 0) {
            events_LEAK = round(n_signal_value[3][ibrem][iyear] * ratio_leak_brem0_bin2_value);
            n_leak[ibin][ibrem][iyear]->setVal(n_signal_value[3][ibrem][iyear] * ratio_leak_brem0_bin2_value);
            n_leak[ibin][ibrem][iyear]->setConstant(kTRUE);
            }
          if(ibrem == 1) {
            events_LEAK = round(n_signal_value[3][ibrem][iyear] * ratio_leak_brem1_bin2_value);
            n_leak[ibin][ibrem][iyear]->setVal(n_signal_value[3][ibrem][iyear] * ratio_leak_brem1_bin2_value);
            n_leak[ibin][ibrem][iyear]->setConstant(kTRUE);
            }

          if(events_LEAK !=0 ) { 
            data_leak[ibin][ibrem][iyear] = model_leak[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_LEAK) ;
            data_full[ibin][ibrem][iyear]->append(*data_leak[ibin][ibrem][iyear]);
            }
          }
        // if(ibin == 4 && ibrem ==1 && signal_decay == "PIPIEE") {
        //   // events_LEAK = r_leak[ibin][ibrem][iyear]->Poisson(n_leak_value[ibin][ibrem][iyear]);
        //   events_LEAK = round(n_signal_value[3][ibrem][iyear] * ratio_leak_brem1_bin4_value);
        //   n_leak[ibin][ibrem][iyear]->setVal(n_signal_value[3][ibrem][iyear] * ratio_leak_brem1_bin4_value);
        //   n_leak[ibin][ibrem][iyear]->setConstant(kTRUE);
        //   if(events_LEAK !=0 ) { 
        //     data_leak[ibin][ibrem][iyear] = model_leak[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_LEAK) ;
        //     data_full[ibin][ibrem][iyear]->append(*data_leak[ibin][ibrem][iyear]);
        //     }
        // }
        if(ibin == 3 && ibrem ==1 && signal_decay == "PIPIEE") {
          // events_LEAK = r_leak[ibin][ibrem][iyear]->Poisson(n_leak_value[ibin][ibrem][iyear]);
          events_LEAK = round(n_signal_value[2][ibrem][iyear] * ratio_leak_brem1_bin3_value);
          n_leak[ibin][ibrem][iyear]->setVal(n_signal_value[2][ibrem][iyear] * ratio_leak_brem1_bin3_value);
          n_leak[ibin][ibrem][iyear]->setConstant(kTRUE);
          if(events_LEAK !=0 ) { 
            data_leak[ibin][ibrem][iyear] = model_leak[ibin][ibrem][iyear]->generate(RooArgSet(x_data_toys),events_LEAK) ;
            data_full[ibin][ibrem][iyear]->append(*data_leak[ibin][ibrem][iyear]);
            }
        }

        
        
        map_combData.insert(std::pair<std::string, RooDataSet*>(name_combined_data.Data(),data_full[ibin][ibrem][iyear]));


        
        

        // }
        }//year
      }//brem
  }//specific bin bumber
}//nbins

//change name
BF_signal[bin_number]->SetName("BFsignal");

RooDataSet data("data","combined data",x_data_toys,RooFit::Index(sample_CLs),RooFit::Import(map_combData)) ;
data.SetName("data"); // the name of the dataset MUST be "data" in order for the framework to identify it.

// //
//   RooDataSet& data = *mass_model.generate(RooArgSet(*x_data[bin_number]),5000); // we use a reference here to avoid copying the object
//   data.SetName("data"); // the name of the dataset MUST be "data" in order for the framework to identify it.
//   //
//   /////////////////////////////////////////////////////////


  //
  // We also must define some RooArgSets: 
  //
  /////////////////////////////////////////////////////////

  //One contaninting the constraint PDFs,
  RooArgSet constraint_set(BF_norm_constraint,
  *Eff_Ratio_constraint[bin_number][0][0], *Eff_Ratio_constraint[bin_number][0][1], *Eff_Ratio_constraint[bin_number][1][0], *Eff_Ratio_constraint[bin_number][1][1], 
  *Y_Norm_constraint[bin_number][0][0], *Y_Norm_constraint[bin_number][0][1], *Y_Norm_constraint[bin_number][1][0], *Y_Norm_constraint[bin_number][1][1],
  "constraint_set");


  // The workspace must also include the fit result of an initial fit of the model to data.
  RooFitResult& rooFitResult = *mass_model.fitTo(data,RooFit::Save(), RooFit::ExternalConstraints(constraint_set));
  std::cout<<" MY RESULTS ON REAL DATA: "<<std::endl;
  rooFitResult.Print();

  /////////////////////////////////////////////////////////

  //one containing the global Observables,
  RooArgSet global_observables_set(BF_norm_obs,
  *Eff_Ratio_obs[bin_number][0][0], *Eff_Ratio_obs[bin_number][0][1], *Eff_Ratio_obs[bin_number][1][0], *Eff_Ratio_obs[bin_number][1][1],
  *Y_Norm_obs[bin_number][0][0], *Y_Norm_obs[bin_number][0][1], *Y_Norm_obs[bin_number][1][0], *Y_Norm_obs[bin_number][1][1],
   "global_observables_set");


  //one containing the normal Observables (the bin variables in the datasets usually) and
  RooArgSet dataset_observables_set(*x_data[bin_number], sample_CLs, "datasetObservables");


  //one containing the parameters

  RooArgSet parameters_set;

  if(bin_number == 0){
  RooArgSet parameters_set2(*BF_signal[bin_number],
  *n_comb[bin_number][0][0], *n_comb[bin_number][1][0], *n_comb[bin_number][0][1], *n_comb[bin_number][1][1],
  *n_misid_brem0[bin_number][0], *n_misid_brem0[bin_number][1], 
  *n_eta[0][0], *n_eta[1][0], *n_eta[0][1], *n_eta[1][1],
  "parameters");
  parameters_set.add(parameters_set2);
  }
  // RooArgSet parameters_set(*BF_signal[bin_number], norm_constant, exponent, n_bkg, "parameters");

  else if(decay == "PIPIEE" && bin_number > 1){
  if(bin_number == 2){
  RooArgSet parameters_set2(*BF_signal[bin_number],
  *n_comb[bin_number][0][0], *n_comb[bin_number][1][0], *n_comb[bin_number][0][1], *n_comb[bin_number][1][1],
  *n_misid_brem0[bin_number][0], *n_misid_brem0[bin_number][1], 
  *n_leak[bin_number][0][0], *n_leak[bin_number][1][0], *n_leak[bin_number][0][1], *n_leak[bin_number][1][1],
  "parameters");
  parameters_set.add(parameters_set2);
  }

  if(bin_number == 3 || bin_number == 4){
  RooArgSet parameters_set2(*BF_signal[bin_number],
  *n_comb[bin_number][0][0], *n_comb[bin_number][1][0], *n_comb[bin_number][0][1], *n_comb[bin_number][1][1],
  *n_misid_brem0[bin_number][0], *n_misid_brem0[bin_number][1], 
  *n_leak[bin_number][1][0], *n_leak[bin_number][1][1],
  "parameters");
  parameters_set.add(parameters_set2);
  }
  }
  else {
  RooArgSet parameters_set2(*BF_signal[bin_number],
  *n_comb[bin_number][0][0], *n_comb[bin_number][1][0], *n_comb[bin_number][0][1], *n_comb[bin_number][1][1],
  *n_misid_brem0[bin_number][0], *n_misid_brem0[bin_number][1], 
  "parameters");
  parameters_set.add(parameters_set2);
  }


/////////////////////////////////////////////////////////
//
// Import everything into a workspace and save it
//
/////////////////////////////////////////////////////////

RooWorkspace workspace("dataset_workspace");
workspace.import(mass_model);
workspace.import(extended_bkg_model);
workspace.import(data);
workspace.import(rooFitResult, "data_fit_result"); // this MUST be called data_fit_result
workspace.defineSet("constraint_set", constraint_set, true);
workspace.defineSet("global_observables_set", global_observables_set, true);
workspace.defineSet("datasetObservables", dataset_observables_set, true);
workspace.defineSet("parameters", parameters_set, true);

TString bn(std::to_string(bin_number));
// Save the workspace to a file
workspace.SaveAs("workspace_"+decay+"_"+bn+".root");

return 0;
// gc.run();




// RooWorkspace *ww = new RooWorkspace("ww", "workspace");
// ww->import(simPdf_CLs);
// RooStats::ModelConfig modelConfig(ww);
// modelConfig.SetPdf(simPdf_CLs);
// modelConfig.SetParametersOfInterest(*BF_signal[bin_number]);

// // modelConfig.SetNuisanceParameters(RooArgSet( *(RooArgSet*) w->allVars().selectByName("BF_norm,Y_Norm_brem0,Y_Norm_brem1,eff_norm_brem0,eff_norm_brem1,eff_sig_brem0,eff_sig_brem1,nbkg_data_expo,nbkg_data_hadr,nbkg_data1")));
// modelConfig.SetObservables(RooArgSet(x_data_toys));
// modelConfig.SetSnapshot(*BF_signal[bin_number]);
// //modelConfig.SetGlobalObservables(RooArgSet(*m_ws.var("Brmumu"),*m_ws.var("effBrem1516"),*m_ws.var("effBrem1718"),*m_ws.var("effmumu1516"),*m_ws.var("effmumu1718"),*m_ws.var("effnoBrem1516"),*m_ws.var("effnoBrem1718")));
// modelConfig.SetName("ModelConfig");
// ww->import(modelConfig);

// RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*) ww->obj("ModelConfig") ;

// RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) sbModel->Clone("BonlyModel") ;
// RooRealVar* poi = (RooRealVar*) bModel->GetParametersOfInterest()->first();
// poi->setVal(0) ;
// bModel->SetSnapshot( *poi  );

// // BF_signal[bin_number]->setVal(1e-7);
// BF_signal[bin_number]->setMin(1e-9);
// BF_signal[bin_number]->setMax(1e-6);



// //running toys




// //Limit
// RooStats::AsymptoticCalculator  asympCalc(combData0, *bModel, *sbModel);
// asympCalc.SetOneSided(true);


// //MC limit
// // RooStats::FrequentistCalculator  freqCalc(*ds, *bModel, *sbModel);

// // RooStats::ProfileLikelihoodTestStat* plr = new RooStats::ProfileLikelihoodTestStat(*sbModel->GetPdf());

// // plr->SetOneSided(true);

// // RooStats::ToyMCSampler* toymcs = (RooStats::ToyMCSampler*) freqCalc.GetTestStatSampler();
// // toymcs->SetTestStatistic(plr);

// // if (!sbModel->GetPdf()->canBeExtended()) {
// //   toymcs->SetNEventsPerToy(1);
// // }

// // freqCalc.SetToys(1000,1000) ;

// RooStats::HypoTestInverter inverter(asympCalc); 
// // RooStats::HypoTestInverter inverter(freqCalc);

// inverter.SetConfidenceLevel(0.95);
// inverter.UseCLs(true);


// inverter.SetVerbose(false);
// //for setting limit
// inverter.SetFixedScan(20,1e-9,1e-5); // set number of points , xmin and xmax
// //for significance study
// // inverter.SetFixedScan(20,1e-8,12e-7); // set number of points , xmin and xmax



// // Calculation is done here 
// RooStats::HypoTestInverterResult* result =  inverter.GetInterval();



// cout << 100*inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << endl;

// std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
// std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
// std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
// std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
// std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
// std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;

// TCanvas* c1 = new TCanvas();
// RooStats::HypoTestInverterPlot* plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot","CLs BF hypothesis scan: "+signal_decay+" "+name_bin+";Signal BF;p-value",result);
// plot->Draw("exp"); 
// // plot->SetXTitle("Signal BF"); 
// // plot->SetYTitle("p-value"); 

// // plot->Draw("CLb 2CL");  // plot also CLb and CLs+b
// c1->Draw();
// // c1->Write();
// c1 -> SaveAs(dir_plots + "CLs_"+signal_decay+"_"+name_bin+".png");
}
