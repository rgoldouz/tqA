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
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "THStack.h"
#include "TMVAGui.C"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include <cmath> 

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

double Bmodification(double &bvalue, double &id){
if (id==0) return bvalue;
else{
const int N=5;
double X[N]={0,0.244,0.679,0.898,1};
double newlight[N]={0,0.2849,0.85505,0.956763,1 };
double newB[N]={0,0.302109,0.770169,0.916225,1};
double newC[N]={0,0.23492,0.62184,0.88628,1};
int ie=0;

for (int i = 0; i < N; i++) {
    if(bvalue < X[i])  {
    if(i == 0) ie = 0; // less than minimal - back extrapolation with the 1st interval
    else  ie = i-1;
    break;
    }
}
  double y1,y2;
  double x1;
  x1 = X[ie];
  double x2;
  x2= X[ie+1];
  double result;

if (abs(id)==5) {
y1 =newB[ie];
y2 = newB[ie+1];}
else if (abs(id)==4) {
y1 =newC[ie];
y2 = newC[ie+1];}
else {
y1 =newlight[ie];
y2 = newlight[ie+1];}
  result= (y1*(x2-bvalue) + y2*(bvalue-x1))/(x2-x1);
return result;}
}
using namespace TMVA;
using namespace std;
void BDToutputtopmasswindow2systematictotaltopmass() 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

  TString sysname ="MVApp.root";
  TFile *sysinput(0);
  sysinput = TFile::Open( sysname ); // if not: download from ROOT server
  std::vector<string> variables_;
  TString name;
  variables_.push_back("BDT__zjethist");
//  variables_.push_back("BDT__phjethist");
  variables_.push_back("BDT__tbartchhist");
  variables_.push_back("BDT__tt3hist");
  variables_.push_back("BDT__ttphhist");
  variables_.push_back("BDT__wwphhist");
  variables_.push_back("BDT__zzhist");
  variables_.push_back("BDT__zgammahist");
  variables_.push_back("BDT__singleantitopphotonhist");

std::vector<TH1F*> addhists;
std::vector<TH1F*> wjetandwphjet;

wjetandwphjet.push_back((TH1F*) sysinput->Get((std::string("BDT__wjet").c_str())));
wjetandwphjet.push_back((TH1F*) sysinput->Get((std::string("BDT__wphjethist").c_str())));

std::vector<TH1F*> jesuphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__JES__plus";
jesuphists.push_back((TH1F*) sysinput->Get(name));
}
//jesuphists.push_back((TH1F*) sysinput->Get((std::string("BDT__wjet").c_str())));
//jesuphists.push_back((TH1F*) sysinput->Get((std::string("BDT__wphjethist").c_str())));
for(unsigned int idx=1; idx<variables_.size(); ++idx){
jesuphists[idx]->Add(jesuphists[idx-1]);
}
addhists.push_back(jesuphists[variables_.size()-1]);

std::vector<TH1F*> jesdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__JES__minus";
jesdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
jesdownhists[idx]->Add(jesdownhists[idx-1]);}
addhists.push_back(jesdownhists[variables_.size()-1]);

std::vector<TH1F*> jeruphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__JER__plus";
jeruphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
jeruphists[idx]->Add(jeruphists[idx-1]);}
addhists.push_back(jeruphists[variables_.size()-1]);

std::vector<TH1F*> jerdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__JER__minus";
jerdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
jerdownhists[idx]->Add(jerdownhists[idx-1]);}
addhists.push_back(jerdownhists[variables_.size()-1]);

std::vector<TH1F*> phesuphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PhES__plus";
phesuphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
phesuphists[idx]->Add(phesuphists[idx-1]);}
addhists.push_back(phesuphists[variables_.size()-1]);

std::vector<TH1F*> phesdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PhES__minus";
phesdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
phesdownhists[idx]->Add(phesdownhists[idx-1]);}
addhists.push_back(phesdownhists[variables_.size()-1]);

std::vector<TH1F*> puuphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PU__plus";
puuphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
puuphists[idx]->Add(puuphists[idx-1]);}
addhists.push_back(puuphists[variables_.size()-1]);

std::vector<TH1F*> pudownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PU__minus";
pudownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
pudownhists[idx]->Add(pudownhists[idx-1]);}
addhists.push_back(pudownhists[variables_.size()-1]);

std::vector<TH1F*> triguphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__TRIG__plus";
triguphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
triguphists[idx]->Add(triguphists[idx-1]);}
addhists.push_back(triguphists[variables_.size()-1]);

std::vector<TH1F*> trigdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__TRIG__minus";
trigdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
trigdownhists[idx]->Add(trigdownhists[idx-1]);}
addhists.push_back(trigdownhists[variables_.size()-1]);

std::vector<TH1F*> btaguphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__BTAG__plus";
btaguphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
btaguphists[idx]->Add(btaguphists[idx-1]);}
addhists.push_back(btaguphists[variables_.size()-1]);

std::vector<TH1F*> btagdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__BTAG__minus";
btagdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
btagdownhists[idx]->Add(btagdownhists[idx-1]);}
addhists.push_back(btagdownhists[variables_.size()-1]);

std::vector<TH1F*> misstaguphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__MISSTAG__plus";
misstaguphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
misstaguphists[idx]->Add(misstaguphists[idx-1]);}
addhists.push_back(misstaguphists[variables_.size()-1]);

std::vector<TH1F*> misstagdownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__MISSTAG__minus";
misstagdownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
misstagdownhists[idx]->Add(misstagdownhists[idx-1]);}
addhists.push_back(misstagdownhists[variables_.size()-1]);

std::vector<TH1F*> muonuphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__MUON__plus";
muonuphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
muonuphists[idx]->Add(muonuphists[idx-1]);}
addhists.push_back(muonuphists[variables_.size()-1]);

std::vector<TH1F*> muondownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__MUON__minus";
muondownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
muondownhists[idx]->Add(muondownhists[idx-1]);}
addhists.push_back(muondownhists[variables_.size()-1]);

std::vector<TH1F*> photonuphists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PHOTON__plus";
photonuphists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
photonuphists[idx]->Add(photonuphists[idx-1]);}
addhists.push_back(photonuphists[variables_.size()-1]);

std::vector<TH1F*> photondownhists;
for(unsigned int i=0; i<variables_.size(); ++i){
name=variables_[i]+"__PHOTON__minus";
photondownhists.push_back((TH1F*) sysinput->Get(name));}
for(unsigned int idx=1; idx<variables_.size(); ++idx){
photondownhists[idx]->Add(photondownhists[idx-1]);}
addhists.push_back(photondownhists[variables_.size()-1]);

std::vector<std::vector<double_t> > vec(photondownhists[0]->GetNbinsX(), vector<double>(18));
  for(int p = 0; p <photondownhists[0]->GetNbinsX(); p++){ //loop over bins
    for(int m = 0; m < 18; m++){ //loop over systematics
vec[p][m]=addhists[m]->GetBinContent(p+1)+wjetandwphjet[0]->GetBinContent(p+1)+wjetandwphjet[1]->GetBinContent(p+1);
cout<<vec[p][m]<<endl;
}}

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
//  reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
  reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
//  reader->AddVariable ("leptoncharge", &leptoncharge);

   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVA";

   // Book method(s)

         TString methodName =  TString("BDT method");
         TString weightfile = dir + prefix + TString("_") + TString("BDT") + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
   
   // Book output histograms
   UInt_t nbin = 20;
   double min=-1;
   double max=1;

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<TH1F*> hists;
std::vector<TH1F*> datahists;
std::vector<TH1F*> revDATAhists;

float scales[] = {0.628,0.0978,34.01,6.133,1.04,0.32,0.02,0.002,0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139,0.049094905/19.145};
samples_.push_back("WJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("G_Pt_50to80.root");
samples_.push_back("G_Pt_80to120.root");
samples_.push_back("G_Pt_120to170.root");
samples_.push_back("G_Pt_170to300.root");
samples_.push_back("G_Pt_300to470.root");
samples_.push_back("G_Pt_470to800.root");
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
datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");

std::vector<string> datasamplesreverse_;
datasamplesreverse_.push_back("etarev/REALDATA1.root");
datasamplesreverse_.push_back("etarev/REALDATA2.root");
datasamplesreverse_.push_back("etarev/REALDATA3.root");

TH1F   *wphjethist(0), *zjethist(0) , *phjethist(0), *wjethist(0), *twchhist(0), *tbarwhist(0),  *tschhist(0), *tbarschhist(0), *ttchhist(0), *tbartchhist(0), *tt1hist(0) ,*tt2hist(0), *tt3hist(0), *ttphhist(0), *wwphhist(0), *wwhist(0), *wzhist(0), *zzhist(0), *zgammahist(0),*singletopphotonhist(0), *singleantitopphotonhist(0), *signalhist(0), *G_Pt_50to80(0),*G_Pt_80to120(0), *G_Pt_120to170(0), *G_Pt_170to300(0) ,*G_Pt_300to470(0),*G_Pt_470to800(0)  ;

TH1F   *wphjethistSB(0), *zjethistSB(0) , *phjethistSB(0), *wjethistSB(0), *twchhistSB(0), *tbarwhistSB(0),  *tschhistSB(0), *tbarschhistSB(0), *ttchhistSB(0), *tbartchhistSB(0), *tt1histSB(0) ,*tt2histSB(0), *tt3histSB(0), *ttphhistSB(0), *wwphhistSB(0), *wwhistSB(0), *wzhistSB(0), *zzhistSB(0), *zgammahistSB(0),*singletopphotonhistSB(0), *singleantitopphotonhistSB(0), *signalhistSB(0), *G_Pt_50to80SB(0),*G_Pt_80to120SB(0), *G_Pt_120to170SB(0), *G_Pt_170to300SB(0) ,*G_Pt_300to470SB(0),*G_Pt_470to800SB(0)  ;

TH1F *data1hist(0), *data2hist(0) ,*data3hist(0) ,*datahistsideband(0);
TH1F *data1histrev(0), *data2histrev(0) ,*data3histrev(0), *datahistrevsideband(0);

wphjethist = new TH1F( "mu_BDT__wphjethist",           "mu_BDT__wphjethist",           nbin, min, max );
zjethist = new TH1F( "mu_BDT__zjethist",           "mu_BDT__zjethist",           nbin, min, max );
G_Pt_50to80= new TH1F( "mu_BDT__G_Pt_50to80",           "mu_BDT__G_Pt_50to80",           nbin, min, max );
G_Pt_80to120= new TH1F( "mu_BDT__G_Pt_80to120",           "mu_BDT__G_Pt_80to120",           nbin, min, max );
G_Pt_120to170= new TH1F( "mu_BDT__G_Pt_120to170",           "mu_BDT__G_Pt_120to170",           nbin, min, max );
G_Pt_170to300= new TH1F( "mu_BDT__G_Pt_170to300",           "mu_BDT__G_Pt_170to300",           nbin, min, max );
G_Pt_300to470= new TH1F( "mu_BDT__G_Pt_300to470",           "mu_BDT__G_Pt_300to470",           nbin, min, max );
G_Pt_470to800= new TH1F( "mu_BDT__G_Pt_470to800",           "mu_BDT__G_Pt_470to800",           nbin, min, max );
wjethist = new TH1F( "mu_BDT__wjethist",           "mu_BDT__wjethist",           nbin, min, max);
twchhist = new TH1F( "mu_BDT__twchhist",           "mu_BDT__twchhist",           nbin, min, max );
tbarwhist = new TH1F( "mu_BDT__tbarwhist",           "mu_BDT__tbarwhist",           nbin, min, max );
tschhist = new TH1F( "mu_BDT__tschhist",           "mu_BDT__tschhist",           nbin, min, max );
tbarschhist = new TH1F( "mu_BDT__tbarschhist",           "mu_BDT__tbarschhist",           nbin, min, max );
ttchhist = new TH1F( "mu_BDT__ttchhist",           "mu_BDT__ttchhist",           nbin, min, max );
tbartchhist = new TH1F( "mu_BDT__tbartchhist",           "mu_BDT__tbartchhist",           nbin, min, max);
tt1hist = new TH1F( "mu_BDT__tt1hist",           "mu_BDT__tt1hist",           nbin,min, max );
tt2hist = new TH1F( "mu_BDT__tt2hist",           "mu_BDT__tt2hist",           nbin, min, max);
tt3hist = new TH1F( "mu_BDT__tt3hist",           "mu_BDT__tt3hist",           nbin, min, max);
ttphhist = new TH1F( "mu_BDT__ttphhist",           "mu_BDT__ttphhist",           nbin, min, max);
wwphhist = new TH1F( "mu_BDT__wwphhist",           "BDT__wwphhist",           nbin,min, max );
wwhist = new TH1F( "mu_BDT__wwhist",           "mu_BDT__wwhist",           nbin,min, max );
wzhist = new TH1F( "mu_BDT__wzhist",           "mu_BDT__wzhist",           nbin, min, max );
zzhist = new TH1F( "mu_BDT__zzhist",           "mu_BDT__zzhist",           nbin, min, max );
zgammahist = new TH1F( "mu_BDT__zgammahist",           "mu_BDT__zgammahist",           nbin,min, max );
singletopphotonhist = new TH1F( "mu_BDT__singletopphotonhist",           "mu_BDT__singletopphotonhist",           nbin, min, max);
singleantitopphotonhist = new TH1F( "mu_BDT__singleantitopphotonhist",           "mu_BDT__singleantitopphotonhist",           nbin,min, max );
signalhist = new TH1F( "mu_BDT__signal100",           "mu_BDT__signal100",           nbin, min, max );


data1hist = new TH1F( "mu_BDT__data1hist",           "mu_BDT__data1hist",           nbin, min, max );
data2hist = new TH1F( "mu_BDT__data2hist",           "mu_BDT__data2hist",           nbin, min, max );
data3hist = new TH1F( "mu_BDT__DATA",           "mu_BDT__DATA",           nbin, min, max );
datahistsideband = new TH1F( "mu_BDT__DATA_sideband",           "mu_BDT__DATA_sideband",           nbin, min, max);


data1histrev = new TH1F( "mu_BDT__data1histrev",           "mu_BDT__data1histrev",           nbin, min, max );
data2histrev = new TH1F( "mu_BDT__data2histrev",           "mu_BDT__data2histrev",           nbin,min, max );
data3histrev = new TH1F( "mu_BDT__DATArev",           "mu_BDT__DATArev",           nbin, min, max );
datahistrevsideband = new TH1F( "mu_BDT__DATArevsideband",           "mu_BDT__DATArevsideband",           nbin, min, max );

wphjethistSB = new TH1F( "mu_BDT__wphjethist__JES__SB",           "mu_BDT__wphjethist__JES__SB",           nbin,min, max );
zjethistSB = new TH1F( "mu_BDT__zjethist__JES__SB",           "mu_BDT__zjethist__JES__SB",           nbin, min, max );
G_Pt_50to80SB= new TH1F( "mu_BDT__G_Pt_50to80SB",           "mu_BDT__G_Pt_50to80SB",           nbin, min, max );
G_Pt_80to120SB= new TH1F( "mu_BDT__G_Pt_80to120SB",           "mu_BDT__G_Pt_80to120SB",           nbin, min, max );
G_Pt_120to170SB= new TH1F( "mu_BDT__G_Pt_120to170SB",           "mu_BDT__G_Pt_120to170SB",           nbin, min, max );
G_Pt_170to300SB= new TH1F( "mu_BDT__G_Pt_170to300SB",           "mu_BDT__G_Pt_170to300SB",           nbin, min, max );
G_Pt_300to470SB= new TH1F( "mu_BDT__G_Pt_300to470SB",           "mu_BDT__G_Pt_300to470SB",           nbin, min, max );
G_Pt_470to800SB= new TH1F( "mu_BDT__G_Pt_470to800SB",           "mu_BDT__G_Pt_470to800SB",           nbin, min, max );
wjethistSB = new TH1F( "mu_BDT__wjethist__JES__SB",           "mu_BDT__wjethist__JES__SB",           nbin, min, max );
twchhistSB = new TH1F( "mu_BDT__twchhist__JES__SB",           "mu_BDT__twchhist__JES__SB",           nbin,min, max);
tbarwhistSB = new TH1F( "mu_BDT__tbarwhist__JES__SB",           "mu_BDT__tbarwhist__JES__SB",           nbin,min, max );
tschhistSB = new TH1F( "mu_BDT__tschhist__JES__SB",           "mu_BDT__tschhist__JES__SB",           nbin, min, max );
tbarschhistSB = new TH1F( "mu_BDT__tbarschhist__JES__SB",           "mu_BDT__tbarschhist__JES__SB",           nbin, min, max );
ttchhistSB = new TH1F( "mu_BDT__ttchhist__JES__SB",           "mu_BDT__ttchhist__JES__SB",           nbin, min, max );
tbartchhistSB = new TH1F( "mu_BDT__tbartchhist__JES__SB",           "mu_BDT__tbartchhist__JES__SB",           nbin, min, max);
tt1histSB = new TH1F( "mu_BDT__tt1hist__JES__SB",           "mu_BDT__tt1hist__JES__SB",           nbin, min, max );
tt2histSB = new TH1F( "mu_BDT__tt2hist__JES__SB",           "mu_BDT__tt2hist__JES__SB",           nbin, min, max );
tt3histSB = new TH1F( "mu_BDT__tt3hist__JES__SB",           "mu_BDT__tt3hist__JES__SB",           nbin, min, max );
ttphhistSB = new TH1F( "mu_BDT__ttphhist__JES__SB",           "mu_BDT__ttphhist__JES__SB",           nbin,min, max );
wwphhistSB = new TH1F( "mu_BDT__wwphhist__JES__SB",           "BDT__wwphhist__JES__SB",           nbin,min, max );
wwhistSB = new TH1F( "mu_BDT__wwhist__JES__SB",           "mu_BDT__wwhist__JES__SB",           nbin, min, max );
wzhistSB = new TH1F( "mu_BDT__wzhist__JES__SB",           "mu_BDT__wzhist__JES__SB",           nbin, min, max );
zzhistSB = new TH1F( "mu_BDT__zzhist__JES__SB",           "mu_BDT__zzhist__JES__SB",           nbin, min, max);
zgammahistSB = new TH1F( "mu_BDT__zgammahist__JES__SB",           "mu_BDT__zgammahist__JES__SB",           nbin, min, max );
singletopphotonhistSB = new TH1F( "mu_BDT__singletopphotonhistSB",           "mu_BDT__singletopphotonhistSB",           nbin,min, max );
singleantitopphotonhistSB = new TH1F( "mu_BDT__singleantitopphotonhistSB",           "mu_BDT__singleantitopphotonhistSB",           nbin, min, max);
signalhistSB = new TH1F( "mu_BDT__signal100__JES__SB",           "mu_BDT__signal100__JES__SB",           nbin,min, max );

std::vector<TH1F*> SBhists;
SBhists.push_back(wjethistSB);
SBhists.push_back(zjethistSB);
SBhists.push_back(G_Pt_50to80SB);
SBhists.push_back(G_Pt_80to120SB);
SBhists.push_back(G_Pt_120to170SB);
SBhists.push_back(G_Pt_170to300SB);
SBhists.push_back(G_Pt_300to470SB);
SBhists.push_back(G_Pt_470to800SB);
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
SBhists.push_back(signalhistSB);

hists.push_back(wjethist);
hists.push_back(zjethist);
hists.push_back(G_Pt_50to80);
hists.push_back(G_Pt_80to120);
hists.push_back(G_Pt_120to170);
hists.push_back(G_Pt_170to300);
hists.push_back(G_Pt_300to470);
hists.push_back(G_Pt_470to800);
hists.push_back(wphjethist);
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
hists.push_back(signalhist);

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Sumw2();}

datahists.push_back(data1hist);
datahists.push_back(data2hist);
datahists.push_back(data3hist);
datahists.push_back(datahistsideband);
for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
datahists[idx]->Sumw2();}

revDATAhists.push_back(data1histrev);
revDATAhists.push_back(data2histrev);
revDATAhists.push_back(data3histrev);
revDATAhists.push_back(datahistrevsideband);

double insidewphjet=0;
double outsidewphjet=0;
double insidewjet=0;
double outsidewjet=0;
double nsignalevent=0;
double mtopup=220;
double mtopdown=130;
//bool SR=false;
//bool SB=true;
bool SR=true;
bool SB=false;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
   TString fname =samples_[idx];
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
std::vector<double> *myjetmatchinginfo=0;

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
	theTree->SetBranchAddress( "jetmatchinginfo", &myjetmatchinginfo );



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

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

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
//deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float) Bmodification((*mycvsdiscriminant)[0],(*myjetmatchinginfo)[0]);
//cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];
finalweight=(*myweight)[0];
//cout<<(*myweight)[0]<<endl;
//if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup){
hists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight );
if (samples_[idx]=="WPHJET.root")insidewphjet=insidewphjet+finalweight;
if (samples_[idx]=="SIGNALtGu.root")nsignalevent=nsignalevent+1;
//cout<<insidewphjet<<endl;
//}
//else {
SBhists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ),finalweight );
if (samples_[idx]=="WPHJET.root")outsidewphjet=outsidewphjet+finalweight;//}


      // Retrieve also per-event error
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
}
for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
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


   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
//   std::cout << "--- ... Processing event: " << ievt << std::endl;
      theTree->GetEntry(ievt);
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
//deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];

//if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) 
datahists[idx] ->Fill( reader->EvaluateMVA( "BDT method"           ) );
//else datahists[3]->Fill( reader->EvaluateMVA( "BDT method"           ) );
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

}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
   TString fname =datasamplesreverse_[idx];
   if (!gSystem->AccessPathName( fname )) input = TFile::Open( fname ); // check if file in local directory exists
   else
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server

   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

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
 for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
                       //   std::cout << "--- ... Processing event: " << ievt << std::endl;
   if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
   theTree->GetEntry(ievt);

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
//deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
//leptoncharge=(float)(*myleptoncharge )[0];

//if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) {
//revDATAhists[idx]->Fill( reader->EvaluateMVA( "BDT method"           ) );
insidewjet=insidewjet+1;
revDATAhists[idx]->Fill( reader->EvaluateMVA( "BDT method"           ) );
//}
//else {
//revDATAhists[3]->Fill( reader->EvaluateMVA( "BDT method"           ) );
outsidewjet=outsidewjet+1;
revDATAhists[3]->Fill( reader->EvaluateMVA( "BDT method"           ) );
//}


//cout<<insidewjet<<endl;
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
////delete mybjetmultiplicity;
////delete myleptoncharge;
////delete myplot;
//
}

double wphjetscale;
wphjetscale=insidewphjet/(insidewphjet+outsidewphjet);
cout<<"wphjetscale=    "<<wphjetscale<<endl;
double wjetscale;
wjetscale=insidewjet/(insidewjet+outsidewjet);
cout<<"wjetscale=    "<<wjetscale<<endl;
cout<<"nsignalevent=    "<<nsignalevent<<endl;
//cout<<insidewphjet<<"insidewphjet"<<"       "<<wphjetscale<<"       "<<insidewjet/(insidewjet+outsidewjet)<<endl;
float lumi = 1;

if (SR==true){
double ff=0;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
if (idx !=0 && idx!=3){
ff=hists[idx]->Integral()+ff;
cout<<samples_[idx]<<"         =        "<<hists[idx]->Integral()<<"     " <<ff<<endl;}
}


THStack *hs1 = new THStack("hs1","BDT output");
for(unsigned int idx=1; idx<datasamplesreverse_.size(); ++idx){
revDATAhists[idx]->Add(revDATAhists[idx-1]);
}
//cout<<"*********************"<< datahists[3]->Integral()<<"       "<<wphjetscale<<endl;
//cout<<"*********************"<< revDATAhists[2]->Integral()<<"       "<<wjetscale<<endl;

revDATAhists[2]->Scale(620.373/revDATAhists[2]->Integral());
for(unsigned int idx=1; idx<revDATAhists[2]->GetNbinsX()+1; ++idx){
//revDATAhists[2]->SetBinError(idx,(revDATAhists[2]->GetBinContent(idx)/revDATAhists[2]->Integral())*74.84);
revDATAhists[2]->SetBinError(idx,0);
//if (revDATAhists[2]->GetBinError(idx)>revDATAhists[2]->GetBinContent(idx)) revDATAhists[2]->SetBinError(idx, revDATAhists[2]->GetBinContent(idx)/2); 
}
//revDATAhists[2]->Scale(wjetscale);



for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
datahists[idx]->Add(datahists[idx-1]);}
cout<<"     " <<datahists[2]->Integral()<<endl;


hists[8]->Scale(3173/hists[8]->Integral());

TH1F *datatoMC(0);

//datahists[3]->Scale(wphjetscale);

//hists[1]->Add(revDATAhists[2]);
//hists[2]->Add(hists[1]);
//datahists[3]->Add(hists[2]);
//hists[4]->Add(datahists[3]);

//for(unsigned int idx=5; idx<samples_.size()-1; ++idx){
//   hists[idx]->Add(hists[idx-1]);}
//cout<<"**********real data***********"<< datahists[2]->Integral()<<"       "<<wphjetscale<<endl;
//cout<<"********** mc ***********"<< hists[18]->Integral()<<"       "<<wjetscale<<endl;
// setup the canvas and draw the histograms




TCanvas *c1 = new TCanvas("c1","signal region",50,50,865,780);
c1->cd();

//W+jet
revDATAhists[2]->SetFillColor(kBlue-2);
revDATAhists[2]->SetLineColor(kBlack);
hs1->Add(revDATAhists[2]);
//Z+jet
hists[1]->SetFillColor(kOrange-4);
hists[1]->SetLineColor(kBlack);
hs1->Add(hists[1]);
//photon+jet
hists[3]->Add(hists[2]);
hists[4]->Add(hists[3]);
hists[5]->Add(hists[4]);
hists[6]->Add(hists[5]);
hists[7]->Add(hists[6]);
hists[7]->SetFillColor(19);
//hs1->Add(hists[7]);
//W+photon+jet
hists[8]->SetFillColor(kGreen-3);
hists[8]->SetLineColor(kBlack);
hs1->Add(hists[8]);

//single top+singletop photon
hists[5+5]->Add(hists[4+5]);
hists[6+5]->Add(hists[5+5]);
hists[7+5]->Add(hists[6+5]);
hists[8+5]->Add(hists[7+5]);
hists[9+5]->Add(hists[8+5]);
hists[19+5]->Add(hists[9+5]);
hists[20+5]->Add(hists[19+5]);
hists[20+5]->SetFillColor(kRed+3);
hists[20+5]->SetLineColor(kBlack);
hs1->Add(hists[20+5]);
//hists[9+5]->SetFillColor(kAzure+10);
//hs1->Add(hists[9+5]);

hists[11+5]->Add(hists[10+5]);
hists[12+5]->Add(hists[11+5]);
hists[13+5]->Add(hists[12+5]);
hists[13+5]->SetFillColor(kPink+1);
hists[13+5]->SetLineColor(kBlack);
hs1->Add(hists[13+5]);
//hists[13+5]->SetFillColor(17);
//hs1->Add(hists[13+5]);
//hists[14+5]->SetFillColor(kSpring-9);
//hs1->Add(hists[14+5]);
hists[15+5]->Add(hists[14+5]);
hists[16+5]->Add(hists[15+5]);
hists[17+5]->Add(hists[16+5]);
hists[17+5]->SetFillColor(kViolet-7);
hists[17+5]->SetLineColor(kBlack);
hs1->Add(hists[17+5]);
hists[18+5]->SetFillColor(kAzure+10);
hists[18+5]->SetLineColor(kBlack);
hs1->Add(hists[18+5]);
//hists[20+5]->Add(hists[19+5]);
//hists[20+5]->SetFillColor(kYellow+3);
//hs1->Add(hists[20+5]);

hs1->Draw("hist");
hs1->SetMaximum(1.2*datahists[2]->GetMaximum());
//hs1->GetXaxis()->SetTitle("BDT output");
hs1->GetYaxis()->SetTitle("Events / 0.1");
hs1->GetYaxis()->SetTitleSize(0.055*1.25);
hs1->GetYaxis()->SetTitleFont(22);
hs1->GetYaxis()->SetTitleOffset(0.6);
hs1->GetYaxis()->SetLabelSize(0.044);

hists[21+5]->SetLineColor(kRed+3);
hists[21+5]->SetLineWidth(3);
hists[21+5]->Draw("histsame");

datahists[2]->SetLineWidth(3.);
datahists[2]->SetLineColor(kBlack);
datahists[2]->SetMarkerColor(kBlack);
datahists[2]->SetMarkerStyle(20.);
datahists[2]->SetMarkerSize(1.35);
datahists[2]->Draw("esame");

TLegend* leg = new TLegend(0.65,0.55,0.85,0.85);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( datahists[2], "Data"               , "PL");
//  leg->AddEntry( hists[25], "Single top+#gamma"               , "F");
  leg->AddEntry( hists[23], "Z#gamma"               , "F");
  leg->AddEntry( hists[22], "WW,WZ,ZZ,WW#gamma "               , "F");
//  leg->AddEntry( hists[19], "WW#gamma"               , "F");
  leg->AddEntry( hists[18], "t#bar{t}, t#bar{t}#gamma"               , "F");
//  leg->AddEntry( hists[17], "t#bar{t}"               , "F");
  leg->AddEntry( hists[25], "Single top, Single top+#gamma"               , "F");
//  leg->AddEntry( hists[14], "Single top"               , "F");
  leg->AddEntry( hists[8], "#gamma+jets"                           , "F");
  leg->AddEntry( hists[1], "Z+jets"                           , "F");
  leg->AddEntry(revDATAhists[2], "W+jets"                           , "F");
  leg->AddEntry( hists[26], "Signal(tu#gamma) 1 pb"               , "L");
//  leg->AddEntry(sum_h, "Stat uncertainty"               , "F");


//  leg->AddEntry( datahists[2], "CMS Data 2012(19.145/fb)"               , "PL");

leg->Draw("same");

}


if (SB==true){
for(unsigned int idx=0; idx<samples_.size(); ++idx){
SBhists[idx]->Scale(lumi*scales[idx]);}
revDATAhists[3]->Scale(620.32/revDATAhists[3]->Integral());
revDATAhists[3]->Scale(1-wjetscale);
SBhists[1]->Add(revDATAhists[3]);
SBhists[2]->Add(SBhists[1]);
SBhists[4]->Add(SBhists[2]);

for(unsigned int idx=5; idx<samples_.size()-1; ++idx){
   SBhists[idx]->Add(SBhists[idx-1]);}
SBhists[20]->SetMaximum(1.5*datahists[3]->GetMaximum());
SBhists[20]->SetFillColor(kMagenta+2);
SBhists[20]->Draw();
SBhists[18]->SetFillColor(kOrange+4);
SBhists[18]->Draw("same");
SBhists[17]->SetFillColor(kOrange-2);
SBhists[17]->Draw("same");

SBhists[16]->SetFillColor(kRed);
SBhists[16]->Draw("same");
SBhists[15]->SetFillColor(kViolet+1);
SBhists[15]->Draw("same");
SBhists[14]->SetFillColor(kSpring-9);
SBhists[14]->Draw("same");
SBhists[13]->SetFillColor(32);
SBhists[13]->Draw("same");
SBhists[12]->SetFillColor(6);
SBhists[12]->Draw("same");
SBhists[9]->SetFillColor(4);
SBhists[9]->Draw("same");
//hists[8]->SetFillColor(4);
//hists[8]->Draw("same");
//hists[7]->SetFillColor(3);
//hists[7]->Draw("same");
//hists[6]->SetFillColor(3);
//hists[6]->Draw("same");
//hists[5]->SetFillColor(2);
//hists[5]->Draw("same");
//hists[4]->SetFillColor(2);
//hists[4]->Draw("same");
//hists[3]->SetFillColor(5);
//hists[3]->Draw("same");
//datahists[3]->SetFillColor(5);
//datahists[3]->Draw("same");
SBhists[2]->SetFillColor(8);
SBhists[2]->Draw("same");
SBhists[1]->SetFillColor(kOrange+7);
SBhists[1]->Draw("same");
revDATAhists[3]->SetFillColor(7);
revDATAhists[3]->Draw("same");
//hists[0]->SetFillColor(7);
//hists[0]->Draw("same");
SBhists[21]->SetFillColor(1);
SBhists[21]->SetFillStyle(3004);
SBhists[21]->Draw("same");


 // plot data points
datahists[3]->SetLineWidth(3.);
datahists[3]->SetLineColor(kBlack);
datahists[3]->SetMarkerColor(kBlack);
datahists[3]->SetMarkerStyle(20.);
datahists[3]->Draw("esame");
//conv->RedrawAxis();


TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( revDATAhists[2], "W JET"                           , "F");
  leg->AddEntry( SBhists[1], "Z JET"                           , "F");
  leg->AddEntry( SBhists[2], "PH JET"                           , "F");
//  leg->AddEntry( datahists[3], "W PH JET"              , "F");
//  leg->AddEntry( hists[5], "TOP-W-CH"               , "F");
//  leg->AddEntry( hists[5], "T-S-CH"               , "F");
//  leg->AddEntry( hists[7], "TOP-S-CH"               , "F");
//  leg->AddEntry( hists[7], "TTBAR-CH"               , "F");
//  leg->AddEntry( hists[8], "TBAR-W-CH"               , "F");
  leg->AddEntry( SBhists[9], "SINGLE TOP  "               , "F");
  leg->AddEntry( SBhists[12], "TTBAR"               , "F");
  leg->AddEntry( SBhists[13], "TTG"               , "F");
  leg->AddEntry( SBhists[14], "WWG"               , "F");
  leg->AddEntry( SBhists[15], "WW"               , "F");
  leg->AddEntry( SBhists[16], "WZ"               , "F");
  leg->AddEntry( SBhists[17], "ZZ"               , "F");
  leg->AddEntry( SBhists[18], "ZGAMMA"               , "F");
  leg->AddEntry( SBhists[20], "SINGLE TOP+PHOTON"               , "F");
  leg->AddEntry( SBhists[21], "SIGNAL"               , "F");
leg->AddEntry( datahists[3], "CMS Data 2012(19.145/fb)"               , "PL");

  leg->Draw("same");
}
} 
