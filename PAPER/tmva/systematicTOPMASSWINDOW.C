/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void systematicTOPMASSWINDOW( TString myMethodList = "" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

  Float_t ptphoton,etaphoton,ptmuon,etamuon,ptjet,etajet,masstop,mtw,deltaRphotonjet,deltaRphotonmuon,ht,costopphoton,deltaphiphotonmet,cvsdiscriminant;
Float_t jetmultiplicity,bjetmultiplicity,leptoncharge;
   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  reader->AddVariable ("ptphoton", &ptphoton);
//  reader->AddVariable ("etaphoton", &etaphoton);
  reader->AddVariable ("ptmuon", &ptmuon);
//  reader->AddVariable ("etamuon", &etamuon);
  reader->AddVariable ("ptjet", &ptjet);
//  reader->AddVariable ("etajet", &etajet);
//  reader->AddVariable ("masstop", &masstop);
//  reader->AddVariable ("mtw", &mtw);
  reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
  reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
//  reader->AddVariable ("ht", &ht);
//  reader->AddVariable ("photonmuonmass", &photonmuonmass);
  reader->AddVariable ("costopphoton", &costopphoton);
//  reader->AddVariable ("topphotonmass", &topphotonmass);
//reader->AddVariable ("pttop", &pttop);
//reader->AddVariable ("etatop", &etatop);
  reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
//  reader->AddVariable ("bjetmultiplicity", &bjetmultiplicity);
  reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
  reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
//  reader->AddVariable ("leptoncharge", &leptoncharge);


/*
   // Spectator variables declared in the training have to be added to the reader, too
   reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   reader->AddSpectator( "spec2 := var1*3",   &spec2 );

   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }
*/
   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVA";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   UInt_t nbin = 16;
  double min=-0.6;
   double max=1;

   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<TH1F*> datahists;
std::vector<TH1F*> revDATAhists;
std::vector<TH1F*> ihists;
std::vector<TH1F*> iSBhists;

float scales[] = {0.0978,1.491,0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024,19.145*0.0169,19.145*0.0169,19.145*0.0169,19.145*0.0169,19.145*0.0169};
//samples_.push_back("WJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("PHJET200400.root");
samples_.push_back("WPHJET.root");
samples_.push_back("T-W-CH.root");
samples_.push_back("TBAR-W-CH.root");
samples_.push_back("T-S-CH.root");
samples_.push_back("TBAR-S-CH.root");
samples_.push_back("T-T-CH.root");
samples_.push_back("TBAR-T-CH.root");
samples_.push_back("TTBAR1.root");
samples_.push_back("TTBAR2.root");
samples_.push_back("TTBAR3.root");
samples_.push_back("TTG.root");
samples_.push_back("WWG.root");
samples_.push_back("WW.root");
samples_.push_back("WZ.root");
samples_.push_back("ZZ.root");
samples_.push_back("ZGAMMA.root");
samples_.push_back("SINGLE-TOP-PH.root");
samples_.push_back("SINGLE-ANTITOP-PH.root");
samples_.push_back("SIGNALtGu.root");
samples_.push_back("SIGNALtGu.root");
samples_.push_back("SIGNALtGu.root");
samples_.push_back("SIGNALtGu.root");
samples_.push_back("SIGNALtGu.root");


datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");

std::vector<string> datasamplesreverse_;
datasamplesreverse_.push_back("etarev/REALDATA1.root");
datasamplesreverse_.push_back("etarev/REALDATA2.root");
datasamplesreverse_.push_back("etarev/REALDATA3.root");

std::vector<string> systematics;
systematics.push_back("__JES__plus");
systematics.push_back("__JES__minus");
systematics.push_back("__JER__plus");
systematics.push_back("__JER__minus");
systematics.push_back("__PU__plus");
systematics.push_back("__PU__minus");
systematics.push_back("__TRIG__plus");
systematics.push_back("__TRIG__minus");
systematics.push_back("__BTAG__plus");
systematics.push_back("__BTAG__minus");
systematics.push_back("__MISSTAG__plus");
systematics.push_back("__MISSTAG__minus");
systematics.push_back("__MUON__plus");
systematics.push_back("__MUON__minus");
systematics.push_back("__PHOTON__plus");
systematics.push_back("__PHOTON__minus");
systematics.push_back("");

map<string, double> eventwight;
TList* hList = new TList();      // list of histograms to store
double insidewphjet;
double outsidewphjet;
double insidewjet;
double outsidewjet;
double mtopup=220;
double mtopdown=130;

for(unsigned int phi=0; phi<systematics.size(); ++phi){
insidewphjet=0;
outsidewphjet=0;
insidewjet=0;
outsidewjet=0;

std::vector<TH1F*> hists;
std::vector<TH1F*> SBhists;

TH1F   *wphjethist(0), *zjethist(0) , *phjethist(0), *wjethist(0), *twchhist(0), *tbarwhist(0),  *tschhist(0), *tbarschhist(0), *ttchhist(0), *tbartchhist(0), *tt1hist(0) ,*tt2hist(0), *tt3hist(0), *ttphhist(0), *wwphhist(0), *wwhist(0), *wzhist(0), *zzhist(0), *zgammahist(0),*singletopphotonhist(0), *singleantitopphotonhist(0), *signalhist5(0) , *signalhist10(0), *signalhist20(0), *signalhist30(0), *signalhist40(0);

TH1F   *wphjethistSB(0), *zjethistSB(0) , *phjethistSB(0), *wjethistSB(0), *twchhistSB(0), *tbarwhistSB(0),  *tschhistSB(0), *tbarschhistSB(0), *ttchhistSB(0), *tbartchhistSB(0), *tt1histSB(0) ,*tt2histSB(0), *tt3histSB(0), *ttphhistSB(0), *wwphhistSB(0), *wwhistSB(0), *wzhistSB(0), *zzhistSB(0), *zgammahistSB(0),*singletopphotonhistSB(0), *singleantitopphotonhistSB(0), *signalhistSB5(0) , *signalhistSB10(0), *signalhistSB20(0), *signalhistSB30(0), *signalhistSB40(0);

TH1F *data1histrev(0), *data2histrev(0) ,*data3histrev(0),*datahistsideband(0);

wphjethist = new TH1F( std::string("BDT__wphjethists").append(systematics[phi]).c_str(), std::string("BDT__wphjethists").append(systematics[phi]).c_str()  , nbin,  min, max    );
zjethist = new TH1F( std::string("BDT__zjethist").append(systematics[phi]).c_str(), std::string("BDT__zjethist").append(systematics[phi]).c_str(),  nbin,  min, max    );
phjethist = new TH1F( std::string("BDT__phjethist").append(systematics[phi]).c_str(),std::string("BDT__phjethist").append(systematics[phi]).c_str() , nbin,  min, max   );
//wjethist = new TH1F( std::string("BDT__wjethist").append(systematics[phi]).c_str(),std::string("BDT__wjethist").append(systematics[phi]).c_str() , nbin, -0.8, 0.8 );
twchhist = new TH1F( std::string("BDT__twchhist").append(systematics[phi]).c_str(),std::string("BDT__twchhist").append(systematics[phi]).c_str() ,nbin,  min, max   );
tbarwhist = new TH1F( std::string("BDT__tbarwhist").append(systematics[phi]).c_str(),std::string("BDT__tbarwhist").append(systematics[phi]).c_str() ,nbin,  min, max    );
tschhist = new TH1F( std::string("BDT__tschhist").append(systematics[phi]).c_str(), std::string("BDT__tschhist").append(systematics[phi]).c_str(), nbin,  min, max  );
tbarschhist = new TH1F( std::string("BDT__tbarschhist").append(systematics[phi]).c_str(),std::string("BDT__tbarschhist").append(systematics[phi]).c_str(), nbin,  min, max  );
ttchhist = new TH1F( std::string("BDT__ttchhist").append(systematics[phi]).c_str(),std::string("BDT__ttchhist").append(systematics[phi]).c_str(), nbin,  min, max   );
tbartchhist = new TH1F( std::string("BDT__tbartchhist").append(systematics[phi]).c_str(), std::string("BDT__tbartchhist").append(systematics[phi]).c_str(), nbin, min, max   );
tt1hist = new TH1F( std::string("BDT__tt1hist").append(systematics[phi]).c_str(), std::string("BDT__tt1hist").append(systematics[phi]).c_str(), nbin, min, max  );
tt2hist = new TH1F( std::string("BDT__tt2hist").append(systematics[phi]).c_str(), std::string("BDT__tt2hist").append(systematics[phi]).c_str(), nbin, min, max  );
tt3hist = new TH1F( std::string("BDT__tt3hist").append(systematics[phi]).c_str(),std::string("BDT__tt3hist").append(systematics[phi]).c_str() , nbin,  min, max  );
ttphhist = new TH1F( std::string("BDT__ttphhist").append(systematics[phi]).c_str(),std::string("BDT__ttphhist").append(systematics[phi]).c_str() ,nbin, min, max   );
wwphhist = new TH1F( std::string("BDT__wwphhist").append(systematics[phi]).c_str(),std::string("BDT__wwphhist").append(systematics[phi]).c_str(), nbin,min, max   );
wwhist = new TH1F( std::string("BDT__wwhist").append(systematics[phi]).c_str(),std::string("BDT__wwhist").append(systematics[phi]).c_str()  ,nbin,  min, max );
wzhist = new TH1F( std::string("BDT__wzhist").append(systematics[phi]).c_str(),std::string("BDT__wzhist").append(systematics[phi]).c_str(), nbin, min, max   );
zzhist = new TH1F( std::string("BDT__zzhist").append(systematics[phi]).c_str(),std::string("BDT__zzhist").append(systematics[phi]).c_str() ,nbin,  min, max  );
zgammahist = new TH1F( std::string("BDT__zgammahist").append(systematics[phi]).c_str(),std::string("BDT__zgammahist").append(systematics[phi]).c_str() ,nbin, min, max  );
singletopphotonhist = new TH1F( std::string("BDT__singletopphotonhist").append(systematics[phi]).c_str(),  std::string("BDT__singletopphotonhist").append(systematics[phi]).c_str(),nbin, min, max );
singleantitopphotonhist = new TH1F( std::string("BDT__singleantitopphotonhist").append(systematics[phi]).c_str(),           std::string("BDT__singleantitopphotonhist").append(systematics[phi]).c_str(),nbin, min, max );
signalhist5 = new TH1F( std::string("BDT__signal5").append(systematics[phi]).c_str(),std::string("BDT__signal5").append(systematics[phi]).c_str()  ,nbin, min, max    );
signalhist10 = new TH1F( std::string("BDT__signal10").append(systematics[phi]).c_str(),std::string("BDT__signal10").append(systematics[phi]).c_str() ,nbin,  min, max   );
signalhist20 = new TH1F( std::string("BDT__signal20").append(systematics[phi]).c_str(),std::string("BDT__signal20").append(systematics[phi]).c_str() ,nbin, min, max    );
signalhist30 = new TH1F( std::string("BDT__signal30").append(systematics[phi]).c_str(),std::string("BDT__signal30").append(systematics[phi]).c_str() ,nbin, min, max  );
signalhist40 = new TH1F( std::string("BDT__signal40").append(systematics[phi]).c_str(),std::string("BDT__signal40").append(systematics[phi]).c_str() ,nbin,min, max  );


//side band region histograms out of top mass
wphjethistSB = new TH1F( std::string("BDT__wphjethistSB").append(systematics[phi]).c_str(), std::string("BDT__wphjethistSB").append(systematics[phi]).c_str()  , nbin,  min, max   );
zjethistSB = new TH1F( std::string("BDT__zjethistSB").append(systematics[phi]).c_str(), std::string("BDT__zjethistSB").append(systematics[phi]).c_str(),  nbin, min, max    );
phjethistSB = new TH1F( std::string("BDT__phjethistSB").append(systematics[phi]).c_str(),std::string("BDT__phjethistSB").append(systematics[phi]).c_str() , nbin,  min, max  );
//wjethist = new TH1F( std::string("BDT__wjethist").append(systematics[phi]).c_str(),std::string("BDT__wjethist").append(systematics[phi]).c_str() , nbin, -0.8, 0.8 );
twchhistSB = new TH1F( std::string("BDT__twchhistSB").append(systematics[phi]).c_str(),std::string("BDT__twchhistSB").append(systematics[phi]).c_str() ,nbin,  min, max  );
tbarwhistSB = new TH1F( std::string("BDT__tbarwhistSB").append(systematics[phi]).c_str(),std::string("BDT__tbarwhistSB").append(systematics[phi]).c_str() ,nbin,  min, max  );
tschhistSB = new TH1F( std::string("BDT__tschhistSB").append(systematics[phi]).c_str(), std::string("BDT__tschhistSB").append(systematics[phi]).c_str(), nbin, min, max   );
tbarschhistSB = new TH1F( std::string("BDT__tbarschhistSB").append(systematics[phi]).c_str(),std::string("BDT__tbarschhistSB").append(systematics[phi]).c_str(), nbin, min, max   );
ttchhistSB = new TH1F( std::string("BDT__ttchhistSB").append(systematics[phi]).c_str(),std::string("BDT__ttchhistSB").append(systematics[phi]).c_str(), nbin, min, max   );
tbartchhistSB = new TH1F( std::string("BDT__tbartchhistSB").append(systematics[phi]).c_str(), std::string("BDT__tbartchhistSB").append(systematics[phi]).c_str(), nbin,  min, max   );
tt1histSB = new TH1F( std::string("BDT__tt1histSB").append(systematics[phi]).c_str(), std::string("BDT__tt1histSB").append(systematics[phi]).c_str(), nbin, min, max   );
tt2histSB = new TH1F( std::string("BDT__tt2histSB").append(systematics[phi]).c_str(), std::string("BDT__tt2histSB").append(systematics[phi]).c_str(), nbin, min, max   );
tt3histSB = new TH1F( std::string("BDT__tt3histSB").append(systematics[phi]).c_str(),std::string("BDT__tt3histSB").append(systematics[phi]).c_str() , nbin, min, max   );
ttphhistSB = new TH1F( std::string("BDT__ttphhistSB").append(systematics[phi]).c_str(),std::string("BDT__ttphhistSB").append(systematics[phi]).c_str() ,nbin, min, max   );
wwphhistSB = new TH1F( std::string("BDT__wwphhistSB").append(systematics[phi]).c_str(),std::string("BDT__wwphhistSB").append(systematics[phi]).c_str(), nbin,  min, max   );
wwhistSB = new TH1F( std::string("BDT__wwhistSB").append(systematics[phi]).c_str(),std::string("BDT__wwhistSB").append(systematics[phi]).c_str()  ,nbin, min, max   );
wzhistSB = new TH1F( std::string("BDT__wzhistSB").append(systematics[phi]).c_str(),std::string("BDT__wzhistSB").append(systematics[phi]).c_str(), nbin,  min, max  );
zzhistSB = new TH1F( std::string("BDT__zzhistSB").append(systematics[phi]).c_str(),std::string("BDT__zzhistSB").append(systematics[phi]).c_str() ,nbin,  min, max   );
zgammahistSB = new TH1F( std::string("BDT__zgammahistSB").append(systematics[phi]).c_str(),std::string("BDT__zgammahistSB").append(systematics[phi]).c_str() ,nbin,  min, max    );
singletopphotonhistSB = new TH1F( std::string("BDT__singletopphotonhistSB").append(systematics[phi]).c_str(),  std::string("BDT__singletopphotonhistSB").append(systematics[phi]).c_str(),nbin, min, max);
singleantitopphotonhistSB = new TH1F( std::string("BDT__singleantitopphotonhistSB").append(systematics[phi]).c_str(),           std::string("BDT__singleantitopphotonhistSB").append(systematics[phi]).c_str(),nbin, min, max  );
signalhistSB5= new TH1F( std::string("BDT__signal5SB").append(systematics[phi]).c_str(),std::string("BDT__signal5SB").append(systematics[phi]).c_str()  ,nbin,  min, max   );
signalhistSB10 = new TH1F( std::string("BDT__signalSB10").append(systematics[phi]).c_str(),std::string("BDT__signalSB10").append(systematics[phi]).c_str() ,nbin, min, max   );
signalhistSB20 = new TH1F( std::string("BDT__signalSB20").append(systematics[phi]).c_str(),std::string("BDT__signalSB20").append(systematics[phi]).c_str() ,nbin, min, max   );
signalhistSB30 = new TH1F( std::string("BDT__signalSB30").append(systematics[phi]).c_str(),std::string("BDT__signalSB30").append(systematics[phi]).c_str() ,nbin,min, max  );
signalhistSB40 = new TH1F( std::string("BDT__signalSB40").append(systematics[phi]).c_str(),std::string("BDT__signalSB40").append(systematics[phi]).c_str() ,nbin,min, max   );


hists.push_back(zjethist);
hists.push_back(phjethist);
hists.push_back(wphjethist);
//hists.push_back(wjethist);
hists.push_back(twchhist);
hists.push_back(tbarwhist);
hists.push_back(tschhist);
hists.push_back(tbarschhist);
hists.push_back(ttchhist);
hists.push_back(tbartchhist);
hists.push_back(tt1hist);
hists.push_back(tt2hist);
hists.push_back(tt3hist);
hists.push_back(ttphhist);
hists.push_back(wwphhist);
hists.push_back(wwhist);
hists.push_back(wzhist);
hists.push_back(zzhist);
hists.push_back(zgammahist);
hists.push_back(singletopphotonhist);
hists.push_back(singleantitopphotonhist);
hists.push_back(signalhist5);
hists.push_back(signalhist10);
hists.push_back(signalhist20);
hists.push_back(signalhist30);
hists.push_back(signalhist40);


SBhists.push_back(zjethistSB);
SBhists.push_back(phjethistSB);
SBhists.push_back(wphjethistSB);
SBhists.push_back(twchhistSB);
SBhists.push_back(tbarwhistSB);
SBhists.push_back(tschhistSB);
SBhists.push_back(tbarschhistSB);
SBhists.push_back(ttchhistSB);
SBhists.push_back(tbartchhistSB);
SBhists.push_back(tt1histSB);
SBhists.push_back(tt2histSB);
SBhists.push_back(tt3histSB);
SBhists.push_back(ttphhistSB);
SBhists.push_back(wwphhistSB);
SBhists.push_back(wwhistSB);
SBhists.push_back(wzhistSB);
SBhists.push_back(zzhistSB);
SBhists.push_back(zgammahistSB);
SBhists.push_back(singletopphotonhistSB);
SBhists.push_back(singleantitopphotonhistSB);
SBhists.push_back(signalhistSB5);
SBhists.push_back(signalhistSB10);
SBhists.push_back(signalhistSB20);
SBhists.push_back(signalhistSB30);
SBhists.push_back(signalhistSB40);

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Sumw2();}

for(unsigned int idx=0; idx<samples_.size(); ++idx){
   TFile *input(0);
TString fname;
if (phi<4) fname =systematics[phi]+"/"+samples_[idx];
    else fname =samples_[idx];
   if (!gSystem->AccessPathName( fname )) input = TFile::Open( fname ); // check if file in local directory exists
   else    
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
  
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corres[1]ponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
  //Double_t  myptphoton,myetaphoton,myptmuon,myetamuon,myptjet,myetajet,mymasstop,mymtw,mydeltaRphotonjet,mydeltaRphotonmuon,myht,mycostopphoton,mydeltaphiphotonmet,mycvsdiscriminant,myjetmultiplicity,mybjetmultiplicity,myleptoncharge;


std::vector<double> *myptphoton=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myetamuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *myetajet=0;
std::vector<double> *mymasstop=0;
//std::vector<double> *mymtw=0;
std::vector<double> *mydeltaRphotonjet=0;
std::vector<double> *mydeltaRphotonmuon=0;
//std::vector<double> *myht=0;
std::vector<double> *mycostopphoton=0;
std::vector<double> *mydeltaphiphotonmet=0;
std::vector<double> *mycvsdiscriminant=0;
std::vector<double> *myjetmultiplicity=0;
//std::vector<double> *mybjetmultiplicity=0;
//std::vector<double> *myleptoncharge=0;
std::vector<double> *myweight=0;
std::vector<double> *mybtagSF=0;
std::vector<double> *mybtagSFup=0;
std::vector<double> *mybtagSFdown=0;
std::vector<double> *mymistagSFup=0;
std::vector<double> *mymistagSFdown=0;
std::vector<double> *mytriggerSF=0;
std::vector<double> *mytriggerSFup=0;
std::vector<double> *mytriggerSFdown=0;
std::vector<double> *myphotonSF=0;
std::vector<double> *myphotonSFup=0;
std::vector<double> *myphotonSFdown=0;
std::vector<double> *mypileupSF=0;
std::vector<double> *mypileupSFup=0;
std::vector<double> *mypileupSFdown=0;
std::vector<double> *mymuonSFup=0;
std::vector<double> *mymuonSFdown=0;
std::vector<double> *mymuonSF=0;

   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
//   Int_t myjetmultiplicity, mybjetmultiplicity , myleptoncharge;
//   Float_t userVar1, userVar2;

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "etamuon", &myetamuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "etajet", &myetajet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
//   theTree->SetBranchAddress( "mtw", &mymtw );
   theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
   theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
//   theTree->SetBranchAddress( "ht", &myht );
   theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
   theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
//   theTree->SetBranchAddress( "bjetmultiplicity", &mybjetmultiplicity );
   theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
   theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );
   theTree->SetBranchAddress( "weight", &myweight);
   theTree->SetBranchAddress( "btagSF", &mybtagSF);
   theTree->SetBranchAddress( "btagSFup", &mybtagSFup);
   theTree->SetBranchAddress( "btagSFdown", &mybtagSFdown);
   theTree->SetBranchAddress( "mistagSFup", &mymistagSFup);
   theTree->SetBranchAddress( "mistagSFdown", &mymistagSFdown);
   theTree->SetBranchAddress( "triggerSF", &mytriggerSF);
   theTree->SetBranchAddress( "triggerSFup", &mytriggerSFup);
   theTree->SetBranchAddress( "triggerSFdown", &mytriggerSFdown);
   theTree->SetBranchAddress( "photonSF", &myphotonSF);
   theTree->SetBranchAddress( "photonSFup", &myphotonSFup);
   theTree->SetBranchAddress( "photonSFdown", &myphotonSFdown);
   theTree->SetBranchAddress( "muonSF", &mymuonSF);
   theTree->SetBranchAddress( "muonSFup", &mymuonSFup);
   theTree->SetBranchAddress( "muonSFdown", &mymuonSFdown);
   theTree->SetBranchAddress( "pileupSF", &mypileupSF);
   theTree->SetBranchAddress( "pileupSFup", &mypileupSFup);
   theTree->SetBranchAddress( "pileupSFdown", &mypileupSFdown);



 // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

//   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
//   std::cout << "--- ... Processing event: " << ievt << std::endl;
double finalweight;

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
//for (int l=0;l<sizeof(myptphoton);l++){
//std::cout << "--- ... reza: " << myptphoton[l] <<std::endl;
//}
//std::cout << "--- ......................."<< (*mycvsdiscriminant)[0]<<std::endl;
      // --- Return the MVA outputs and fill into histograms
ptphoton=(float)(*myptphoton)[0];
etaphoton=(float)(*myetaphoton )[0];
ptmuon=(float)(*myptmuon )[0];
etamuon=(float)(*myetamuon )[0];
ptjet=(float)(*myptjet )[0];
etajet=(float)(*myetajet )[0];
//masstop=(float)(*mymasstop )[0];
//mtw=(float)(*mymtw )[0];
deltaRphotonjet=(float)(*mydeltaRphotonjet )[0];
deltaRphotonmuon=(float)(*mydeltaRphotonmuon )[0];
//ht=(float)(*myht )[0];
costopphoton=(float)(*mycostopphoton )[0];
jetmultiplicity=(float)(*myjetmultiplicity )[0];
//bjetmultiplicity=(float)(*mybjetmultiplicity )[0];
deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];
finalweight=(*myweight)[0];
//cout<<(*myweight)[0]<<endl;
eventwight["__PU__plus"]=(*myweight)[0]*(*mypileupSFup)[0]/(*mypileupSF)[0];
eventwight["__PU__minus"]=(*myweight)[0]*(*mypileupSFdown)[0]/(*mypileupSF)[0];
eventwight["__TRIG__plus"]=(*myweight)[0]*(*mytriggerSFup)[0]/(*mytriggerSF)[0];
eventwight["__TRIG__minus"]=(*myweight)[0]*(*mytriggerSFdown)[0]/(*mytriggerSF)[0];
eventwight["__BTAG__plus"]=(*myweight)[0]*(*mybtagSFup)[0]/(*mybtagSF)[0];
eventwight["__BTAG__minus"]=(*myweight)[0]*(*mybtagSFdown)[0]/(*mybtagSF)[0];
eventwight["__MISSTAG__plus"]=(*myweight)[0]*(*mymistagSFup)[0]/(*mybtagSF)[0];
eventwight["__MISSTAG__minus"]=(*myweight)[0]*(*mymistagSFdown)[0]/(*mybtagSF)[0];
eventwight["__MUON__plus"]=(*myweight)[0]*(*mymuonSFup)[0]/(*mymuonSF)[0];
eventwight["__MUON__minus"]=(*myweight)[0]*(*mymuonSFdown)[0]/(*mymuonSF)[0];
eventwight["__PHOTON__plus"]=(*myweight)[0]*(*myphotonSFup)[0]/(*myphotonSF)[0];
eventwight["__PHOTON__minus"]=(*myweight)[0]*(*myphotonSFdown)[0]/(*myphotonSF)[0];
eventwight[""]=(*myweight)[0];
if (phi<4) finalweight=(*myweight)[0];
else finalweight=eventwight[systematics[phi].c_str()];
//if (samples_[idx]=="SIGNAL.root") finalweight=(*myweight)[0];
//if (samples_[idx]=="WPHJET")  finalweight=(*mypileupSF)[0]*(*mytriggerSF)[0]*(*mybtagSF)[0]*(*mymuonSF)[0]*(*myphotonSF)[0];
//if (finalweight<0) finalweight=30;
//cout<<"negative event weight"<<finalweight<<"            "<<ptphoton<<endl;
      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup){

hists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight* scales[idx] );
if (samples_[idx]=="WPHJET.root") insidewphjet=insidewphjet+finalweight;
}
else {
SBhists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight* scales[idx] );
if (samples_[idx]=="WPHJET.root") outsidewphjet=outsidewphjet+finalweight;}

//cout<<reader->EvaluateMVA( "BDT method")<<"                            "<<finalweight<<endl;   
//cout<<(*myweight)[0]<<endl;
      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
//delete finalweight;
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

delete myptphoton;
delete myetaphoton;
delete myptmuon;
delete myetamuon;
delete myptjet;
delete myetajet;
delete mymasstop;
//delete mymtw;
delete mydeltaRphotonjet;
delete mydeltaRphotonmuon;
//delete myht;
delete mycostopphoton;
delete mydeltaphiphotonmet;
delete mycvsdiscriminant;
delete myjetmultiplicity;
//delete mybjetmultiplicity;
//delete myleptoncharge;
//delete myplot;
delete mybtagSF;
delete mybtagSFup;
delete mybtagSFdown;
delete mymistagSFup;
delete mytriggerSF;
delete mytriggerSFup;
delete mytriggerSFdown;
delete myphotonSF;
delete myphotonSFup;
delete myphotonSFdown;
delete mypileupSF;
delete mypileupSFup;
delete mypileupSFdown;
delete input;
if (idx==samples_.size()-5) hists[idx]->Scale(5/hists[idx]->Integral());
if (idx==samples_.size()-4) hists[idx]->Scale(10/hists[idx]->Integral());
if (idx==samples_.size()-3) hists[idx]->Scale(20/hists[idx]->Integral());
if (idx==samples_.size()-2) hists[idx]->Scale(30/hists[idx]->Integral());
if (idx==samples_.size()-1) hists[idx]->Scale(40/hists[idx]->Integral());

if (samples_[idx]=="WPHJET.root") hists[idx]->Scale(3173/hists[idx]->Integral());
if (samples_[idx]=="TTBAR2.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="TTBAR3.root") hists[idx]->Add(hists[idx-1]);
if (!(samples_[idx]=="TTBAR1.root" || samples_[idx]=="TTBAR2.root"|| samples_[idx]=="WPHJET.root")) hList->Add(hists[idx]);
if (systematics[phi]==""){
ihists.push_back(hists[idx]);
iSBhists.push_back(SBhists[idx]);}
}
}
TH1F *data1hist(0), *data2hist(0) ,*data3hist(0),*datahistsideband(0);
data1hist = new TH1F( "mu_BDT__data1hist",           "mu_BDT__data1hist",           nbin, min, max );
data2hist = new TH1F( "mu_BDT__data2hist",           "mu_BDT__data2hist",           nbin, min, max );
data3hist = new TH1F( "BDT__DATA",           "BDT__DATA",           nbin, min, max );
datahistsideband = new TH1F( "BDT__wphjethist",           "BDT__wphjethist",           nbin,min, max );

datahists.push_back(data1hist);
datahists.push_back(data2hist);
datahists.push_back(data3hist);
datahists.push_back(datahistsideband);


for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
   TFile *input(0);
   TString fname =datasamples_[idx];
   if (!gSystem->AccessPathName( fname )) input = TFile::Open( fname ); // check if file in local directory exists
   else    
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corres[1]ponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //

std::vector<double> *myptphoton=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myetamuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *myetajet=0;
std::vector<double> *mymasstop=0;
//std::vector<double> *mymtw=0;
std::vector<double> *mydeltaRphotonjet=0;
std::vector<double> *mydeltaRphotonmuon=0;
//std::vector<double> *myht=0;
std::vector<double> *mycostopphoton=0;
std::vector<double> *mydeltaphiphotonmet=0;
std::vector<double> *mycvsdiscriminant=0;
std::vector<double> *myjetmultiplicity=0;
//std::vector<double> *mybjetmultiplicity=0;
//std::vector<double> *myleptoncharge=0;

   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
//   Int_t myjetmultiplicity, mybjetmultiplicity , myleptoncharge;
//   Float_t userVar1, userVar2;

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "etamuon", &myetamuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "etajet", &myetajet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
//   theTree->SetBranchAddress( "mtw", &mymtw );
   theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
   theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
//   theTree->SetBranchAddress( "ht", &myht );
   theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
   theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
//   theTree->SetBranchAddress( "bjetmultiplicity", &mybjetmultiplicity );
   theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
   theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );


 // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
//   std::cout << "--- ... Processing event: " << ievt << std::endl;

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
//for (int l=0;l<sizeof(myptphoton);l++){
//std::cout << "--- ... reza: " << myptphoton[l] <<std::endl;
//}
//std::cout << "--- ......................."<< (*mycvsdiscriminant)[0]<<std::endl;
      // --- Return the MVA outputs and fill into histograms
ptphoton=(float)(*myptphoton)[0];
etaphoton=(float)(*myetaphoton )[0];
ptmuon=(float)(*myptmuon )[0];
etamuon=(float)(*myetamuon )[0];
ptjet=(float)(*myptjet )[0];
etajet=(float)(*myetajet )[0];
//masstop=(float)(*mymasstop )[0];
//mtw=(float)(*mymtw )[0];
deltaRphotonjet=(float)(*mydeltaRphotonjet )[0];
deltaRphotonmuon=(float)(*mydeltaRphotonmuon )[0];
//ht=(float)(*myht )[0];
costopphoton=(float)(*mycostopphoton )[0];
jetmultiplicity=(float)(*myjetmultiplicity )[0];
//bjetmultiplicity=(float)(*mybjetmultiplicity )[0];
deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];




      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );
if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) datahists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ) );
else datahists[3]->Fill( reader->EvaluateMVA( "BDT method"           ) );
}
   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

delete myptphoton;
delete myetaphoton;
delete myptmuon;
delete myetamuon;
delete myptjet;
delete myetajet;
delete mymasstop;
//delete mymtw;
delete mydeltaRphotonjet;
delete mydeltaRphotonmuon;
//delete myht;
delete mycostopphoton;
delete mydeltaphiphotonmet;
delete mycvsdiscriminant;
delete myjetmultiplicity;
//delete mybjetmultiplicity;
//delete myleptoncharge;
//delete myplot;
delete input;
}

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
datahists[idx]->Add(datahists[idx-1]);
}
hList->Add(datahists[2]);


TH1F *data1histrev(0), *data2histrev(0) ,*data3histrev(0), *datahistrevsideband(0);

data1histrev = new TH1F( "BDT__data1histrev",           "BDT__data1histrev",           nbin,min, max );
data2histrev = new TH1F( "BDT__data2histrev",           "BDT__data2histrev",           nbin,min, max );
data3histrev = new TH1F( "BDT__wjet",           "BDT__wjet",           nbin, min, max );
datahistrevsideband = new TH1F( "mu_BDT__DATArevsideband",           "mu_BDT__DATArevsideband",           nbin, min, max );

revDATAhists.push_back(data1histrev);
revDATAhists.push_back(data2histrev);
revDATAhists.push_back(data3histrev);
revDATAhists.push_back(datahistrevsideband);

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
   TFile *input(0);

   TString fname =datasamplesreverse_[idx];
   if (!gSystem->AccessPathName( fname )) input = TFile::Open( fname ); // check if file in local directory exists
   else
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server

   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

   // --- Event loop
   //
   //    // Prepare the event tree
   //       // - here the variable names have to corres[1]ponds to your tree
   //          // - you can use the same variables as above which is slightly faster,
   //             //   but of course you can use different ones and copy the values inside the event loop
   //                //
//   Double_t  myptphoton,myetaphoton,myptmuon,myetamuon,myptjet,myetajet,mymasstop,mymtw,mydeltaRphotonjet,mydeltaRphotonmuon,myht,mycostopphoton,mydeltaphiphotonmet,mycvsdiscriminant,myjetmultiplicity,mybjetmultiplicity,myleptoncharge;
   
std::vector<double> *myptphoton=0;
std::vector<double> *myetaphoton=0;
 std::vector<double> *myptmuon=0;
 std::vector<double> *myetamuon=0;
 std::vector<double> *myptjet=0;
 std::vector<double> *myetajet=0;
 std::vector<double> *mymasstop=0;
 //std::vector<double> *mymtw=0;
 std::vector<double> *mydeltaRphotonjet=0;
std::vector<double> *mydeltaRphotonmuon=0;
 //std::vector<double> *myht=0;
 std::vector<double> *mycostopphoton=0;
 std::vector<double> *mydeltaphiphotonmet=0;
 std::vector<double> *mycvsdiscriminant=0;
 std::vector<double> *myjetmultiplicity=0;
 //std::vector<double> *mybjetmultiplicity=0;
 //std::vector<double> *myleptoncharge=0;
   TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "etamuon", &myetamuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "etajet", &myetajet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
//   theTree->SetBranchAddress( "mtw", &mymtw );
   theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
   theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
//   theTree->SetBranchAddress( "ht", &myht );
    theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
    theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
//   theTree->SetBranchAddress( "bjetmultiplicity", &mybjetmultiplicity );
    theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
    theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );
 // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
        Double_t effS       = 0.7;
 
           std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
 
              std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
                 TStopwatch sw;
                    sw.Start();
                       for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
                       //   std::cout << "--- ... Processing event: " << ievt << std::endl;
                             if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
 
                                   theTree->GetEntry(ievt);
                                   //for (int l=0;l<sizeof(myptphoton);l++){
                                   //std::cout << "--- ... reza: " << myptphoton[l] <<std::endl;
                                   //}
      //                             std::cout << "--- ......................."<< (*mycvsdiscriminant)[0]<<std::endl;
                                         // --- Return the MVA outputs and fill into histograms
                                         ptphoton=(float)(*myptphoton)[0];
                                         etaphoton=(float)(*myetaphoton )[0];
                                         ptmuon=(float)(*myptmuon )[0];
                                         etamuon=(float)(*myetamuon )[0];
                                         ptjet=(float)(*myptjet )[0];
                                         etajet=(float)(*myetajet )[0];
                                       //masstop=(float)(*mymasstop )[0];
                                        //mtw=(float)(*mymtw )[0];
                                         deltaRphotonjet=(float)(*mydeltaRphotonjet )[0];
                                        deltaRphotonmuon=(float)(*mydeltaRphotonmuon )[0];
                                         //ht=(float)(*myht )[0];
                                        costopphoton=(float)(*mycostopphoton )[0];
                                         jetmultiplicity=(float)(*myjetmultiplicity )[0];
                                         //bjetmultiplicity=(float)(*mybjetmultiplicity )[0];
                                         deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
                                         cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
                                         //leptoncharge=(float)(*myleptoncharge )[0];
      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) {
insidewjet=insidewjet+1;
revDATAhists[idx]->Fill( reader->EvaluateMVA( "BDT method"           ) );
}
else {outsidewjet=outsidewjet+1;
revDATAhists[3]->Fill( reader->EvaluateMVA( "BDT method"           ) );
}

}
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

delete myptphoton;
delete myetaphoton;
delete myptmuon;
delete myetamuon;
delete myptjet;
delete myetajet;
delete mymasstop;
//delete mymtw;
delete mydeltaRphotonjet;
delete mydeltaRphotonmuon;
//delete myht;
delete mycostopphoton;
delete mydeltaphiphotonmet;
delete mycvsdiscriminant;
delete myjetmultiplicity;
////delete mybjetmultiplicity;
////delete myleptoncharge;
////delete myplot;
//
delete input;
}


double wphjetscale;
wphjetscale=insidewphjet/(insidewphjet+outsidewphjet);
double wjetscale;
wjetscale=insidewjet/(insidewjet+outsidewjet);
//cout<<insidewphjet<<"insidewphjet"<<"       "<<wphjetscale<<"       "<<insidewjet/(insidewjet+outsidewjet)<<endl;
float lumi = 1;
for(unsigned int idx=1; idx<datasamplesreverse_.size(); ++idx){
revDATAhists[idx]->Add(revDATAhists[idx-1]);
 } 
revDATAhists[2]->Scale(212.609/revDATAhists[2]->Integral());
for(unsigned int idx=1; idx<revDATAhists[2]->GetNbinsX()+1; ++idx){
revDATAhists[2]->SetBinError(idx,(revDATAhists[2]->GetBinContent(idx)/revDATAhists[2]->Integral())*75);}

hList->Add(revDATAhists[2]);
revDATAhists[3]->Scale(212.609/revDATAhists[3]->Integral());
revDATAhists[3]->Scale((1-wjetscale)/wjetscale);

//cout<<revDATAhists[2]->Integral()<<endl;

datahists[3]->Add(revDATAhists[3],-1);
datahists[3]->Add(iSBhists[1],-1);
datahists[3]->Add(iSBhists[1],-1);

for(unsigned int idx=3; idx<samples_.size()-5; ++idx){
datahists[3]->Add(iSBhists[idx],-1);}
datahists[3]->Scale(1101.23/datahists[3]->Integral());
for(unsigned int idx=1; idx<datahists[3]->GetNbinsX()+1; ++idx){
datahists[3]->SetBinError(idx,(datahists[3]->GetBinContent(idx)/datahists[3]->Integral())*140);}
hList->Add(datahists[3]);


/////////////////////////////////////////////////////////////////////////
std::vector<string> signalsystematics;
std::vector<string> signalsamples_;


signalsystematics.push_back("__Qscale__plus");
signalsystematics.push_back("__Qscale__minus");
signalsystematics.push_back("__Topmass__plus");
signalsystematics.push_back("__Topmass__minus");

signalsamples_.push_back("SIGNALtGu_scaleup.root");
signalsamples_.push_back("SIGNALtGu_scaledown.root");
signalsamples_.push_back("SIGNALtGu_massup.root");
signalsamples_.push_back("SIGNALtGu_massdown.root");

for(unsigned int phi=0; phi<signalsystematics.size(); ++phi){

TH1F *shist5(0) , *shist10(0), *shist20(0), *shist30(0), *shist40(0);
std::vector<TH1F*> Shists;

shist5 = new TH1F( std::string("BDT__signal5").append(signalsystematics[phi]).c_str(),std::string("BDT__signal5").append(signalsystematics[phi]).c_str()  ,nbin, min, max    );
shist10 = new TH1F( std::string("BDT__signal10").append(signalsystematics[phi]).c_str(),std::string("BDT__signal10").append(signalsystematics[phi]).c_str() ,nbin,  min, max   );
shist20 = new TH1F( std::string("BDT__signal20").append(signalsystematics[phi]).c_str(),std::string("BDT__signal20").append(signalsystematics[phi]).c_str() ,nbin, min, max    );
shist30 = new TH1F( std::string("BDT__signal30").append(signalsystematics[phi]).c_str(),std::string("BDT__signal30").append(signalsystematics[phi]).c_str() ,nbin, min, max  );
shist40 = new TH1F( std::string("BDT__signal40").append(signalsystematics[phi]).c_str(),std::string("BDT__signal40").append(signalsystematics[phi]).c_str() ,nbin,min, max  );

Shists.push_back(shist5);
Shists.push_back(shist10);
Shists.push_back(shist20);
Shists.push_back(shist30);
Shists.push_back(shist40);

TFile *input(0);
TString fname;
fname =signalsamples_[phi];
input = TFile::Open( fname );
std::vector<double> *myptphoton=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myetamuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *myetajet=0;
std::vector<double> *mymasstop=0;
//std::vector<double> *mymtw=0;
std::vector<double> *mydeltaRphotonjet=0;
std::vector<double> *mydeltaRphotonmuon=0;
//std::vector<double> *myht=0;
std::vector<double> *mycostopphoton=0;
std::vector<double> *mydeltaphiphotonmet=0;
std::vector<double> *mycvsdiscriminant=0;
std::vector<double> *myjetmultiplicity=0;
//std::vector<double> *mybjetmultiplicity=0;
//std::vector<double> *myleptoncharge=0;
std::vector<double> *myweight=0;
std::vector<double> *mybtagSF=0;
std::vector<double> *mybtagSFup=0;
std::vector<double> *mybtagSFdown=0;
std::vector<double> *mymistagSFup=0;
std::vector<double> *mymistagSFdown=0;
std::vector<double> *mytriggerSF=0;
std::vector<double> *mytriggerSFup=0;
std::vector<double> *mytriggerSFdown=0;
std::vector<double> *myphotonSF=0;
std::vector<double> *myphotonSFup=0;
std::vector<double> *myphotonSFdown=0;
std::vector<double> *mypileupSF=0;
std::vector<double> *mypileupSFup=0;
std::vector<double> *mypileupSFdown=0;
std::vector<double> *mymuonSFup=0;
std::vector<double> *mymuonSFdown=0;
std::vector<double> *mymuonSF=0;

   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "etamuon", &myetamuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "etajet", &myetajet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
//   theTree->SetBranchAddress( "mtw", &mymtw );
   theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
   theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
//   theTree->SetBranchAddress( "ht", &myht );
   theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
   theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
//   theTree->SetBranchAddress( "bjetmultiplicity", &mybjetmultiplicity );
   theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
   theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );
   theTree->SetBranchAddress( "weight", &myweight);
   theTree->SetBranchAddress( "btagSF", &mybtagSF);
   theTree->SetBranchAddress( "btagSFup", &mybtagSFup);
   theTree->SetBranchAddress( "btagSFdown", &mybtagSFdown);
   theTree->SetBranchAddress( "mistagSFup", &mymistagSFup);
   theTree->SetBranchAddress( "mistagSFdown", &mymistagSFdown);
   theTree->SetBranchAddress( "triggerSF", &mytriggerSF);
   theTree->SetBranchAddress( "triggerSFup", &mytriggerSFup);
   theTree->SetBranchAddress( "triggerSFdown", &mytriggerSFdown);
   theTree->SetBranchAddress( "photonSF", &myphotonSF);
   theTree->SetBranchAddress( "photonSFup", &myphotonSFup);
   theTree->SetBranchAddress( "photonSFdown", &myphotonSFdown);
   theTree->SetBranchAddress( "muonSF", &mymuonSF);
   theTree->SetBranchAddress( "muonSFup", &mymuonSFup);
   theTree->SetBranchAddress( "muonSFdown", &mymuonSFdown);
   theTree->SetBranchAddress( "pileupSF", &mypileupSF);
   theTree->SetBranchAddress( "pileupSFup", &mypileupSFup);
   theTree->SetBranchAddress( "pileupSFdown", &mypileupSFdown);
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
   double finalweight;
   theTree->GetEntry(ievt);

ptphoton=(float)(*myptphoton)[0];
etaphoton=(float)(*myetaphoton )[0];
ptmuon=(float)(*myptmuon )[0];
etamuon=(float)(*myetamuon )[0];
ptjet=(float)(*myptjet )[0];
etajet=(float)(*myetajet )[0];
//masstop=(float)(*mymasstop )[0];
//mtw=(float)(*mymtw )[0];
deltaRphotonjet=(float)(*mydeltaRphotonjet )[0];
deltaRphotonmuon=(float)(*mydeltaRphotonmuon )[0];
//ht=(float)(*myht )[0];
costopphoton=(float)(*mycostopphoton )[0];
jetmultiplicity=(float)(*myjetmultiplicity )[0];
//bjetmultiplicity=(float)(*mybjetmultiplicity )[0];
deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];
finalweight=(*myweight)[0];
//cout<<(*myweight)[0]<<endl;
finalweight=(*myweight)[0];

if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup){
Shists[0] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight);
Shists[1] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight);
Shists[2] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight);
Shists[3] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight);
Shists[4] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight);
}
}

delete myptphoton;
delete myetaphoton;
delete myptmuon;
delete myetamuon;
delete myptjet;
delete myetajet;
delete mymasstop;
//delete mymtw;
delete mydeltaRphotonjet;
delete mydeltaRphotonmuon;
//delete myht;
delete mycostopphoton;
delete mydeltaphiphotonmet;
delete mycvsdiscriminant;
delete myjetmultiplicity;
//delete mybjetmultiplicity;
//delete myleptoncharge;
//delete myplot;
delete mybtagSF;
delete mybtagSFup;
delete mybtagSFdown;
delete mymistagSFup;
delete mytriggerSF;
delete mytriggerSFup;
delete mytriggerSFdown;
delete myphotonSF;
delete myphotonSFup;
delete myphotonSFdown;
delete mypileupSF;
delete mypileupSFup;
delete mypileupSFdown;
delete input;
Shists[0]->Scale(5/Shists[0]->Integral());
Shists[1]->Scale(10/Shists[1]->Integral());
Shists[2]->Scale(20/Shists[2]->Integral());
Shists[3]->Scale(30/Shists[3]->Integral());
Shists[4]->Scale(40/Shists[4]->Integral());

hList->Add(Shists[0]);
hList->Add(Shists[1]);
hList->Add(Shists[2]);
hList->Add(Shists[3]);
hList->Add(Shists[4]);
}


/////////////////////////////////////////////////////////////////////////
TFile *target  = new TFile( "MVApp.root","RECREATE" );
  hList->Write();
   target->Close();
}


