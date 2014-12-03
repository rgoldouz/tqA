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
#include <cmath> 


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


void Gammaeta1( TString myMethodList = "" ) 
{   

 
   // Book output histograms
   UInt_t nbin = 15;
   double min=-3;
   double max=3;
   TFile *input(0);

std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<TH1F*> datahists;
std::vector<TH1F*> revDATAhists;
std::vector<TH1F*> ihists;
std::vector<TH1F*> iSBhists;

float scales[] = {0.0978,1.491,0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139,0.049094905,0.049094905,0.049094905,0.049094905,0.049094905/19.145};


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
systematics.push_back("__PhES__plus");
systematics.push_back("__PhES__minus");
systematics.push_back("__btagSF__plus");
systematics.push_back("__btagSF__minus");
systematics.push_back("__PU__plus");
systematics.push_back("__PU__minus");
systematics.push_back("__TRIG__plus");
systematics.push_back("__TRIG__minus");

//systematics.push_back("__MISSTAG__plus");
//systematics.push_back("__MISSTAG__minus");
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
if (phi<6) fname =systematics[phi]+"/"+samples_[idx];
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
std::vector<double> *myjetmatchinginfo=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myphiphoton=0;

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
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
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
   theTree->SetBranchAddress( "jetmatchinginfo", &myjetmatchinginfo );
theTree->SetBranchAddress( "phiphoton", &myphiphoton );


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

finalweight=(*myweight)[0];
//cout<<(*myweight)[0]<<endl;
eventwight["__PU__plus"]=((*myweight)[0]*(*mypileupSFup)[0])/(*mypileupSF)[0];
eventwight["__PU__minus"]=((*myweight)[0]*(*mypileupSFdown)[0])/(*mypileupSF)[0];
eventwight["__TRIG__plus"]=((*myweight)[0]*(*mytriggerSFup)[0])/(*mytriggerSF)[0];
eventwight["__TRIG__minus"]=((*myweight)[0]*(*mytriggerSFdown)[0])/(*mytriggerSF)[0];
eventwight["__BTAG__plus"]=((*myweight)[0]*(*mybtagSFup)[0])/(*mybtagSF)[0];
eventwight["__BTAG__minus"]=((*myweight)[0]*(*mybtagSFdown)[0])/(*mybtagSF)[0];
eventwight["__MISSTAG__plus"]=((*myweight)[0]*(*mymistagSFup)[0])/(*mybtagSF)[0];
eventwight["__MISSTAG__minus"]=((*myweight)[0]*(*mymistagSFdown)[0])/(*mybtagSF)[0];
eventwight["__MUON__plus"]=((*myweight)[0]*(*mymuonSFup)[0])/(*mymuonSF)[0];
eventwight["__MUON__minus"]=((*myweight)[0]*(*mymuonSFdown)[0])/(*mymuonSF)[0];
eventwight["__PHOTON__plus"]=((*myweight)[0]*(*myphotonSFup)[0])/(*myphotonSF)[0];
eventwight["__PHOTON__minus"]=((*myweight)[0]*(*myphotonSFdown)[0])/(*myphotonSF)[0];
eventwight[""]=(*myweight)[0];
//if (phi<4) finalweight=(*myweight)[0];
if (phi>7) finalweight=eventwight[systematics[phi].c_str()];
//if (samples_[idx]=="SIGNALtGu.root") cout<<"negative event weight"<<finalweight<<"            "<<endl;
//if (samples_[idx]=="WPHJET")  finalweight=(*mypileupSF)[0]*(*mytriggerSF)[0]*(*mybtagSF)[0]*(*mymuonSF)[0]*(*myphotonSF)[0];
//if (finalweight<0) finalweight=30;
//cout<<"negative event weight"<<finalweight<<"            "<<ptphoton<<endl;
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup){

hists[idx] ->Fill( (*myetaphoton)[0] ,finalweight*scales[idx] );
if (samples_[idx]=="WPHJET.root") insidewphjet=insidewphjet+finalweight;
}
else {
SBhists[idx] ->Fill( (*myetaphoton)[0] ,finalweight*scales[idx] );
if (samples_[idx]=="WPHJET.root") outsidewphjet=outsidewphjet+finalweight;}
}}
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
delete myjetmatchinginfo;
//if (idx==samples_.size()-5) hists[idx]->Scale(5/hists[idx]->Integral());
//if (idx==samples_.size()-4) hists[idx]->Scale(10/hists[idx]->Integral());
//if (idx==samples_.size()-3) hists[idx]->Scale(20/hists[idx]->Integral());
//if (idx==samples_.size()-2) hists[idx]->Scale(30/hists[idx]->Integral());
//if (idx==samples_.size()-1) hists[idx]->Scale(40/hists[idx]->Integral());

if (samples_[idx]=="WPHJET.root") hists[idx]->Scale(3173/hists[idx]->Integral());
if (samples_[idx]=="TTBAR2.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="TTBAR3.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="SINGLE-ANTITOP-PH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="TBAR-W-CH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="T-S-CH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="TBAR-S-CH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="T-T-CH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="TBAR-T-CH.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="WZ.root") hists[idx]->Add(hists[idx-1]);
if (samples_[idx]=="ZZ.root") hists[idx]->Add(hists[idx-1]);

if (!(samples_[idx]=="TTBAR1.root" || samples_[idx]=="TTBAR2.root"|| samples_[idx]=="SINGLE-TOP-PH.root"||samples_[idx]=="SIGNALtGu.root" || samples_[idx]=="WPHJET.root"||samples_[idx]=="T-W-CH.root" ||samples_[idx]=="TBAR-W-CH.root" ||samples_[idx]=="T-S-CH.root" ||samples_[idx]=="TBAR-S-CH.root" ||samples_[idx]=="T-T-CH.root" ||samples_[idx]=="WW.root" ||samples_[idx]=="WZ.root")  ) hList->Add(hists[idx]);
if (idx==samples_.size()-1)hList->Add(hists[idx]);
if (systematics[phi]==""){
ihists.push_back(hists[idx]);
iSBhists.push_back(SBhists[idx]);}
}
}
TH1F *data1hist(0), *data2hist(0) ,*data3hist(0),*datahistsideband(0);
data1hist = new TH1F( "mu_BDT__data1hist",           "mu_BDT__data1hist",           nbin, min, max );
data1hist->Sumw2();
data2hist = new TH1F( "mu_BDT__data2hist",           "mu_BDT__data2hist",           nbin, min, max );
data2hist->Sumw2();
data3hist = new TH1F( "BDT__DATA",           "BDT__DATA",           nbin, min, max );
data3hist->Sumw2();
datahistsideband = new TH1F( "BDT__wphjethist",           "BDT__wphjethist",           nbin,min, max );
datahistsideband->Sumw2();

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
std::vector<double> *mycoswphoton=0;
std::vector<double> *myphiphoton=0;

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
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );
theTree->SetBranchAddress( "phiphoton", &myphiphoton );


   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
//   std::cout << "--- ... Processing event: " << ievt << std::endl;

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);

if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) datahists[idx] ->Fill( (*myetaphoton)[0]);
else datahists[3]->Fill( (*myetaphoton)[0] );
}}
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
std::vector<double> *mycoswphoton=0;
 //std::vector<double> *mybjetmultiplicity=0;
 //std::vector<double> *myleptoncharge=0;
std::vector<double> *myphiphoton=0;
	std::vector<double> *mysigmaietaieta=0;

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
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
//   theTree->SetBranchAddress( "leptoncharge", &myleptoncharge );
 // Efficiency calculator for cut method
theTree->SetBranchAddress( "phiphoton", &myphiphoton );
 	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );

              std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
                 TStopwatch sw;
                    sw.Start();
                       for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
                       //   std::cout << "--- ... Processing event: " << ievt << std::endl;
                             if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
 
                                   theTree->GetEntry(ievt);
                                   //for (int l=0;l<sizeof(myptphoton);l++){
 
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup) {
revDATAhists[idx]->Fill( (*myetaphoton)[0]);
insidewjet=insidewjet+1;}
else {revDATAhists[3]->Fill( (*myetaphoton)[0] );
outsidewjet=outsidewjet+1;}
}}}
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
revDATAhists[2]->Scale(302.91/revDATAhists[2]->Integral());
for(unsigned int idx=1; idx<revDATAhists[2]->GetNbinsX()+1; ++idx){
revDATAhists[2]->SetBinError(idx,(revDATAhists[2]->GetBinContent(idx)/revDATAhists[2]->Integral())*74.84);
//if (revDATAhists[2]->GetBinError(idx)>revDATAhists[2]->GetBinContent(idx)) revDATAhists[2]->SetBinError(idx, revDATAhists[2]->GetBinContent(idx)/2); 
}

hList->Add(revDATAhists[2]);
revDATAhists[3]->Scale(302.91/revDATAhists[3]->Integral());
revDATAhists[3]->Scale((1-wjetscale)/wjetscale);

//cout<<revDATAhists[2]->Integral()<<endl;

datahists[3]->Add(revDATAhists[3],-1);
datahists[3]->Add(iSBhists[0],-1);
datahists[3]->Add(iSBhists[1],-1);

for(unsigned int idx=3; idx<samples_.size()-5; ++idx){
datahists[3]->Add(iSBhists[idx],-1);}
for(unsigned int idx=1; idx<nbin; ++idx){
if (datahists[3]->GetBinContent(idx)<0)datahists[3]->SetBinContent(idx,0);
} 
datahists[3]->Scale(985.18/datahists[3]->Integral());
for(unsigned int idx=1; idx<datahists[3]->GetNbinsX()+1; ++idx){
datahists[3]->SetBinError(idx,(datahists[3]->GetBinContent(idx)/datahists[3]->Integral())*139.11);}
hList->Add(datahists[3]);

cout<<datahists[3]->Integral()<<endl;
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

float signalscales[] ={0.1084/19.145,0.0933/19.145,0.0838/19.145,0.0892/19.145};

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

for(unsigned int idx=0; idx<Shists.size(); ++idx){
Shists[idx]->Sumw2();}

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
std::vector<double> *myjetmatchinginfo=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myphiphoton=0;

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
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
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
   theTree->SetBranchAddress( "jetmatchinginfo", &myjetmatchinginfo );
theTree->SetBranchAddress( "phiphoton", &myphiphoton );

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
   double finalweight;
   theTree->GetEntry(ievt);

finalweight=(*myweight)[0];

//cout<<(*myweight)[0]<<endl;

if((*mymasstop )[0]>mtopdown && (*mymasstop )[0]<mtopup){
Shists[0] ->Fill( (*myetaphoton)[0],finalweight*signalscales[phi]);
Shists[1] ->Fill((*myetaphoton)[0],finalweight*signalscales[phi]);
Shists[2] ->Fill( (*myetaphoton)[0],finalweight*signalscales[phi]);
Shists[3] ->Fill( (*myetaphoton)[0],finalweight*signalscales[phi]);
Shists[4] ->Fill( (*myetaphoton)[0],finalweight*signalscales[phi]);
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
delete myjetmatchinginfo;
/*Shists[0]->Scale(5/Shists[0]->Integral());
Shists[1]->Scale(10/Shists[1]->Integral());
Shists[2]->Scale(20/Shists[2]->Integral());
Shists[3]->Scale(30/Shists[3]->Integral());
*/
//Shists[4]->Scale(40/Shists[4]->Integral());

//hList->Add(Shists[0]);
//hList->Add(Shists[1]);
//hList->Add(Shists[2]);
//hList->Add(Shists[3]);
hList->Add(Shists[4]);
}


/////////////////////////////////////////////////////////////////////////
TFile *target  = new TFile( "GAMMAETA.root","RECREATE" );
  hList->Write();
   target->Close();
}


