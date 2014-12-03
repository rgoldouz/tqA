#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLine.h"
#include "THStack.h"
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

void btagredefinitionplot() {
 
// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<string> datasamplesreverse_;
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

std::vector<string> samplesreverse_;
samplesreverse_.push_back("etarev/WJET.root");
samplesreverse_.push_back("etarev/ZJET.root");
samplesreverse_.push_back("etarev/G_Pt_50to80.root");
samplesreverse_.push_back("etarev/G_Pt_80to120.root");
samplesreverse_.push_back("etarev/G_Pt_120to170.root");
samplesreverse_.push_back("etarev/G_Pt_170to300.root");
samplesreverse_.push_back("etarev/G_Pt_300to470.root");
samplesreverse_.push_back("etarev/G_Pt_470to800.root");
samplesreverse_.push_back("etarev/WPHJET.root");
samplesreverse_.push_back("etarev/T-W-CH.root");
samplesreverse_.push_back("etarev/TBAR-W-CH.root");
samplesreverse_.push_back("etarev/T-S-CH.root");
samplesreverse_.push_back("etarev/TBAR-S-CH.root");
samplesreverse_.push_back("etarev/T-T-CH.root");
samplesreverse_.push_back("etarev/TBAR-T-CH.root");
samplesreverse_.push_back("etarev/TTBAR1.root");
samplesreverse_.push_back("etarev/TTBAR2.root");
samplesreverse_.push_back("etarev/TTBAR3.root");
samplesreverse_.push_back("etarev/TTG.root");
samplesreverse_.push_back("etarev/WWG.root");
samplesreverse_.push_back("etarev/WW.root");
samplesreverse_.push_back("etarev/WZ.root");
samplesreverse_.push_back("etarev/ZZ.root");
samplesreverse_.push_back("etarev/ZGAMMA.root");
samplesreverse_.push_back("etarev/SINGLE-TOP-PH.root");
samplesreverse_.push_back("etarev/SINGLE-ANTITOP-PH.root");


datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");
datasamplesreverse_.push_back("etarev/REALDATA1.root");
datasamplesreverse_.push_back("etarev/REALDATA2.root");
datasamplesreverse_.push_back("etarev/REALDATA3.root");


/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<TH1F*> hists_photonpt;
std::vector<TH1F*> hists_photoneta;
std::vector<TH1F*> hists_muonpt;
std::vector<TH1F*> hists_muoneta;
std::vector<TH1F*> hists_jetpt;
std::vector<TH1F*> hists_jeteta;
std::vector<TH1F*> hists_topmass;
std::vector<TH1F*> hists_toppt;
std::vector<TH1F*> hists_topeta;
std::vector<TH1F*> hists_topphotonmass;
std::vector<TH1F*> hists_Nvertex;
std::vector<TH1F*> hists_NvertexnoPU;

std::vector<TH1F*> hists_CSV;
std::vector<TH1F*> hists_deltaphimetphoton;
std::vector<TH1F*> hists_deltarphotonmuon;
std::vector<TH1F*> hists_deltarphotonbjet;
std::vector<TH1F*> hists_jetmulti;
std::vector<TH1F*> hists_costopphoton;

//be carefull____________ to make code shorter, histograms of MC samples are from zero to sample size and after that 3 histo for real data and after that 3 files for data in controlregion for wjet sample.
 
for(unsigned int idx=0; idx<samples_.size(); ++idx){

hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(samples_[idx]).c_str(), std::string("photon_pt").append(samples_[idx]).c_str() , 10,10., 410.));
hists_photonpt[idx]->Sumw2();
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 12, -3., 3));
hists_photoneta[idx]->Sumw2();
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(samples_[idx]).c_str(), std::string("muon_pt").append(samples_[idx]).c_str() , 12, 0, 312));
hists_muonpt[idx]->Sumw2();
hists_muoneta.push_back(new TH1F(  std::string("muon_eta").append(samples_[idx]).c_str(), std::string("muon_eta").append(samples_[idx]).c_str() , 12, -3., 3));
hists_muoneta[idx]->Sumw2();
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(samples_[idx]).c_str(),  std::string("jet_pt").append(samples_[idx]).c_str(),  12, 0, 312));
hists_jetpt[idx]->Sumw2();
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(samples_[idx]).c_str(),  std::string("bjet_eta").append(samples_[idx]).c_str() , 12, -3., 3));
hists_jeteta[idx]->Sumw2();
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(samples_[idx]).c_str(), std::string("top_mass").append(samples_[idx]).c_str() , 40, 80, 350.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(samples_[idx]).c_str(), std::string("toppt").append(samples_[idx]).c_str() ,35,0., 750.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(samples_[idx]).c_str(), std::string("topeta").append(samples_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(samples_[idx]).c_str(), std::string("topphotonmass").append(samples_[idx]).c_str() ,40, 0., 1000.));
hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(samples_[idx]).c_str(), std::string("Nvertex").append(samples_[idx]).c_str() ,60, -0.5, 59.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(samples_[idx]).c_str(), std::string("NvertexnoPU").append(samples_[idx]).c_str() ,60, -0.5, 59.5));

hists_CSV.push_back(new TH1F(  std::string("CSV").append(samples_[idx]).c_str(), std::string("CSV").append(samples_[idx]).c_str() ,20, 0, 1));
hists_CSV[idx]->Sumw2();

hists_deltaphimetphoton.push_back(new TH1F(  std::string("deltaphimetphoton").append(samples_[idx]).c_str(), std::string("deltaphimetphoton").append(samples_[idx]).c_str() ,12, 0, 4));
hists_deltaphimetphoton[idx]->Sumw2();
hists_deltarphotonmuon.push_back(new TH1F(  std::string("deltarphotonmuon").append(samples_[idx]).c_str(), std::string("deltarphotonmuon").append(samples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonmuon[idx]->Sumw2();
hists_deltarphotonbjet.push_back(new TH1F(  std::string("deltarphotonbjet").append(samples_[idx]).c_str(), std::string("deltarphotonbjet").append(samples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonbjet[idx]->Sumw2();
hists_jetmulti.push_back(new TH1F(  std::string("jetmulti").append(samples_[idx]).c_str(), std::string("jetmulti").append(samples_[idx]).c_str() ,8, -0.5, 7.5));
hists_jetmulti[idx]->Sumw2();
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(samples_[idx]).c_str(), std::string("costopphoton").append(samples_[idx]).c_str() ,12, -1, 1));
hists_costopphoton[idx]->Sumw2();

}

for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(datasamples_[idx]).c_str(), std::string("photon_pt").append(datasamples_[idx]).c_str() , 10,10., 410.));
hists_photonpt[idx]->Sumw2();
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(datasamples_[idx]).c_str(),  std::string("photon_eta").append(datasamples_[idx]).c_str() , 12, -3., 3));
hists_photoneta[idx]->Sumw2();
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamples_[idx]).c_str(), std::string("muon_pt").append(datasamples_[idx]).c_str() , 12, 0, 312));
hists_muoneta.push_back(new TH1F(  std::string("muon_eta").append(datasamples_[idx]).c_str(), std::string("muon_eta").append(datasamples_[idx]).c_str() ,  12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamples_[idx]).c_str(),  std::string("jet_pt").append(datasamples_[idx]).c_str(),  12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamples_[idx]).c_str(),  std::string("bjet_eta").append(datasamples_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamples_[idx]).c_str(), std::string("top_mass").append(datasamples_[idx]).c_str() , 40, 80, 350.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(datasamples_[idx]).c_str(), std::string("toppt").append(datasamples_[idx]).c_str() ,35,0., 750.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(datasamples_[idx]).c_str(), std::string("topeta").append(datasamples_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(datasamples_[idx]).c_str(), std::string("topphotonmass").append(datasamples_[idx]).c_str() ,40, 0., 1000.));
hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(datasamples_[idx]).c_str(), std::string("Nvertex").append(datasamples_[idx]).c_str() ,60, -0.5, 59.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(datasamples_[idx]).c_str(), std::string("NvertexnoPU").append(datasamples_[idx]).c_str() ,60, -0.5, 59.5));
hists_CSV.push_back(new TH1F(  std::string("CSV").append(datasamples_[idx]).c_str(), std::string("CSV").append(datasamples_[idx]).c_str() ,20, 0, 1));

hists_deltaphimetphoton.push_back(new TH1F(  std::string("deltaphimetphoton").append(datasamples_[idx]).c_str(), std::string("deltaphimetphoton").append(datasamples_[idx]).c_str() ,12, 0, 4));
hists_deltaphimetphoton[idx]->Sumw2();
hists_deltarphotonmuon.push_back(new TH1F(  std::string("deltarphotonmuon").append(datasamples_[idx]).c_str(), std::string("deltarphotonmuon").append(datasamples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonmuon[idx]->Sumw2();
hists_deltarphotonbjet.push_back(new TH1F(  std::string("deltarphotonbjet").append(datasamples_[idx]).c_str(), std::string("deltarphotonbjet").append(datasamples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonbjet[idx]->Sumw2();
hists_jetmulti.push_back(new TH1F(  std::string("jetmulti").append(datasamples_[idx]).c_str(), std::string("jetmulti").append(datasamples_[idx]).c_str() ,8, -0.5, 7.5));
hists_jetmulti[idx]->Sumw2();
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(datasamples_[idx]).c_str(), std::string("costopphoton").append(datasamples_[idx]).c_str() ,12, -1, 1));
hists_costopphoton[idx]->Sumw2();
}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
hists_photonpt.push_back(new  TH1F(  std::string("photon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("photon_pt").append(datasamplesreverse_[idx]).c_str() , 10,10., 410.));
hists_photoneta.push_back(new  TH1F(  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str() , 12, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("muon_pt").append(datasamplesreverse_[idx]).c_str() ,12, 0, 312));
hists_muoneta.push_back(new  TH1F(  std::string("muon_eta").append(datasamplesreverse_[idx]).c_str(), std::string("muon_eta").append(datasamplesreverse_[idx]).c_str() , 12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamplesreverse_[idx]).c_str(), std::string("top_mass").append(datasamplesreverse_[idx]).c_str() , 40, 80, 350.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(datasamplesreverse_[idx]).c_str(), std::string("toppt").append(datasamplesreverse_[idx]).c_str() ,35,0., 750.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(datasamplesreverse_[idx]).c_str(), std::string("topeta").append(datasamplesreverse_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(datasamplesreverse_[idx]).c_str(), std::string("topphotonmass").append(datasamplesreverse_[idx]).c_str() ,40, 0., 1000.));
hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(datasamplesreverse_[idx]).c_str(), std::string("Nvertex").append(datasamplesreverse_[idx]).c_str() ,60, -0.5, 59.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(datasamplesreverse_[idx]).c_str(), std::string("NvertexnoPU").append(datasamplesreverse_[idx]).c_str() ,60, -0.5, 59.5));
hists_CSV.push_back(new TH1F(  std::string("CSV").append(datasamplesreverse_[idx]).c_str(), std::string("CSV").append(datasamplesreverse_[idx]).c_str() ,20, 0, 1));

hists_deltaphimetphoton.push_back(new TH1F(  std::string("deltaphimetphoton").append(datasamplesreverse_[idx]).c_str(), std::string("deltaphimetphoton").append(datasamplesreverse_[idx]).c_str() ,12, 0, 4));
hists_deltaphimetphoton[idx]->Sumw2();
hists_deltarphotonmuon.push_back(new TH1F(  std::string("deltarphotonmuon").append(datasamplesreverse_[idx]).c_str(), std::string("deltarphotonmuon").append(datasamplesreverse_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonmuon[idx]->Sumw2();
hists_deltarphotonbjet.push_back(new TH1F(  std::string("deltarphotonbjet").append(datasamplesreverse_[idx]).c_str(), std::string("deltarphotonbjet").append(datasamplesreverse_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonbjet[idx]->Sumw2();
hists_jetmulti.push_back(new TH1F(  std::string("jetmulti").append(datasamplesreverse_[idx]).c_str(), std::string("jetmulti").append(datasamplesreverse_[idx]).c_str() ,8, -0.5, 7.5));
hists_jetmulti[idx]->Sumw2();
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(datasamplesreverse_[idx]).c_str(), std::string("costopphoton").append(datasamplesreverse_[idx]).c_str() ,12, -1, 1));
hists_costopphoton[idx]->Sumw2();
}

for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
hists_photonpt.push_back(new  TH1F(  std::string("photon_pt").append(samplesreverse_[idx]).c_str(), std::string("photon_pt").append(samplesreverse_[idx]).c_str() , 10,10., 410.));
hists_photoneta.push_back(new  TH1F(  std::string("photon_eta").append(samplesreverse_[idx]).c_str(),  std::string("photon_eta").append(samplesreverse_[idx]).c_str() , 12, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(samplesreverse_[idx]).c_str(), std::string("muon_pt").append(samplesreverse_[idx]).c_str() , 12, 0, 312));
hists_muoneta.push_back(new  TH1F(  std::string("muon_eta").append(samplesreverse_[idx]).c_str(), std::string("muon_eta").append(samplesreverse_[idx]).c_str() ,  12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(samplesreverse_[idx]).c_str(),  std::string("jet_pt").append(samplesreverse_[idx]).c_str(),12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(samplesreverse_[idx]).c_str(),  std::string("bjet_eta").append(samplesreverse_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(samplesreverse_[idx]).c_str(), std::string("top_mass").append(samplesreverse_[idx]).c_str() , 40, 80, 350.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(samplesreverse_[idx]).c_str(), std::string("toppt").append(samplesreverse_[idx]).c_str() ,35,0., 750.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(samplesreverse_[idx]).c_str(), std::string("topeta").append(samplesreverse_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(samplesreverse_[idx]).c_str(), std::string("topphotonmass").append(samplesreverse_[idx]).c_str() ,40, 0., 1000.));
hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(samplesreverse_[idx]).c_str(), std::string("Nvertex").append(samplesreverse_[idx]).c_str() ,60, -0.5, 59.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(samplesreverse_[idx]).c_str(), std::string("NvertexnoPU").append(samplesreverse_[idx]).c_str() ,60, -0.5, 59.5));

hists_CSV.push_back(new TH1F(  std::string("CSV").append(samplesreverse_[idx]).c_str(), std::string("CSV").append(samplesreverse_[idx]).c_str() ,20, 0, 1));

hists_deltaphimetphoton.push_back(new TH1F(  std::string("deltaphimetphoton").append(samplesreverse_[idx]).c_str(), std::string("deltaphimetphoton").append(samplesreverse_[idx]).c_str() ,12, 0, 4));
hists_deltaphimetphoton[idx]->Sumw2();
hists_deltarphotonmuon.push_back(new TH1F(  std::string("deltarphotonmuon").append(samplesreverse_[idx]).c_str(), std::string("deltarphotonmuon").append(samplesreverse_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonmuon[idx]->Sumw2();
hists_deltarphotonbjet.push_back(new TH1F(  std::string("deltarphotonbjet").append(samplesreverse_[idx]).c_str(), std::string("deltarphotonbjet").append(samplesreverse_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonbjet[idx]->Sumw2();
hists_jetmulti.push_back(new TH1F(  std::string("jetmulti").append(samplesreverse_[idx]).c_str(), std::string("jetmulti").append(samplesreverse_[idx]).c_str() ,8, -0.5, 7.5));
hists_jetmulti[idx]->Sumw2();
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(samplesreverse_[idx]).c_str(), std::string("costopphoton").append(samplesreverse_[idx]).c_str() ,12, -1, 1));
hists_costopphoton[idx]->Sumw2();
}

for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<endl;
	TFile *input(0);
	TString fname =samples_[idx];
	input = TFile::Open( fname ); 
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myetamuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *myetajet=0;
	std::vector<double> *mymasstop=0;
        std::vector<double> *mypttop=0;
        std::vector<double> *myetatop=0;
        std::vector<double> *mytopphotonmass=0;
        std::vector<double> *myNvertex=0;
	std::vector<double> *mypileupSF=0;
	std::vector<double> *myweight=0;
        std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *myjetmultiplicity=0;
	std::vector<double> *myjetmatchinginfo=0;

	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "etamuon", &myetamuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "etajet", &myetajet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
        theTree->SetBranchAddress( "pttop", &mypttop);
        theTree->SetBranchAddress( "etatop", &myetatop);
        theTree->SetBranchAddress( "topphotonmass", &mytopphotonmass);
        theTree->SetBranchAddress( "Nvertex", &myNvertex);
	theTree->SetBranchAddress( "pileupSF", &mypileupSF);
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
	theTree->SetBranchAddress( "jetmatchinginfo", &myjetmatchinginfo );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
        if (finalweight<0) cout<<"negative weight=  " <<finalweight<<"    "<<(*myptphoton)[0]<<samples_[idx]<<endl;

//        if (samples_[idx]=="WPHJET.root") finalweight=finalweight*(1.46604+0.00923692*(*myptphoton)[0]-1.48871e-06*(*myptphoton)[0]*(*myptphoton)[0]);
//        if (samples_[idx]=="ZGAMMA.root") finalweight=finalweight*(1.3292+0.000952237*(*myptphoton)[0]+2.00623e-05*(*myptphoton)[0]*(*myptphoton)[0]-1.41325e-07 *(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]+2.48614e-10*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]);

	if (finalweight<0) cout<<"negative weight=  " <<finalweight<<"    "<<(*myptphoton)[0]<<samples_[idx]<<endl;
//	if (finalweight<100) {
	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
	hists_photonpt[idx] ->Fill( (*myptphoton)[0],finalweight *scales[idx] );
	hists_photoneta[idx] ->Fill( (*myetaphoton )[0],finalweight*scales[idx] );
//	hists_muonpt[idx] ->Fill((*myptmuon )[0] ,(finalweight/(*mypileupSF)[0])*scales[idx] );
	hists_muonpt[idx] ->Fill((*myptmuon )[0] ,finalweight*scales[idx] );
	hists_muoneta[idx] ->Fill( (*myetamuon )[0],finalweight *scales[idx]);
	hists_jetpt[idx] ->Fill((*myptjet )[0] ,finalweight *scales[idx]);
	hists_jeteta[idx] ->Fill((*myetajet )[0] ,finalweight*scales[idx] );
	hists_topmass[idx] ->Fill( (*mymasstop )[0],finalweight *scales[idx]);
	hists_toppt[idx] ->Fill( (*mypttop )[0],finalweight *scales[idx]);
	hists_topeta[idx] ->Fill( (*myetatop)[0],finalweight *scales[idx]);
	hists_topphotonmass[idx] ->Fill( (*mytopphotonmass)[0],finalweight *scales[idx]);
	hists_Nvertex[idx] ->Fill( (*myNvertex)[0],(*mypileupSF)[0]*scales[idx]);
	hists_NvertexnoPU[idx] ->Fill( (*myNvertex)[0],scales[idx]);
	hists_CSV[idx] ->Fill( Bmodification((*mycvsdiscriminant)[0],(*myjetmatchinginfo)[0]) ,finalweight *scales[idx]);
	hists_deltaphimetphoton[idx] ->Fill( (*mydeltaphiphotonmet)[0],finalweight *scales[idx]);
	hists_deltarphotonmuon[idx] ->Fill( (*mydeltaRphotonmuon)[0],finalweight *scales[idx]);
	hists_deltarphotonbjet[idx] ->Fill( (*mydeltaRphotonjet)[0],finalweight *scales[idx]);
	hists_jetmulti[idx] ->Fill( (*myjetmultiplicity)[0],finalweight *scales[idx]);
	hists_costopphoton[idx] ->Fill( (*mycostopphoton)[0],finalweight *scales[idx]);

}
//	}
	}
	delete myptphoton;
	delete myetaphoton;
	delete myptmuon;
	delete myetamuon;
	delete myptjet;
	delete myetajet;
	delete mymasstop;
	delete mypttop;
        delete myetatop;
        delete mytopphotonmass;
        delete myNvertex;
        delete mypileupSF;
	delete input;
	delete mycvsdiscriminant;
	delete mydeltaRphotonjet;
	delete mydeltaRphotonmuon;
	delete mycostopphoton;
	delete mydeltaphiphotonmet;
	delete myjetmultiplicity;

}

for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
	TFile *input(0);
	TString fname =datasamples_[idx];
	input = TFile::Open( fname ); 
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myetamuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *myetajet=0;
	std::vector<double> *mymasstop=0;
        std::vector<double> *mypttop=0;
        std::vector<double> *myetatop=0;
        std::vector<double> *mytopphotonmass=0;
        std::vector<double> *myNvertex=0;
	std::vector<double> *myweight=0;
        std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *myjetmultiplicity=0;


	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	cout<<theTree->GetEntries()<<endl;
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "etamuon", &myetamuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "etajet", &myetajet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
        theTree->SetBranchAddress( "pttop", &mypttop);
        theTree->SetBranchAddress( "etatop", &myetatop);
        theTree->SetBranchAddress( "topphotonmass", &mytopphotonmass);
        theTree->SetBranchAddress( "Nvertex", &myNvertex);
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );


	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
	hists_photonpt[idx+samples_.size()] ->Fill( (*myptphoton)[0]);
	hists_photoneta[idx+samples_.size()] ->Fill( (*myetaphoton )[0] );
	hists_muonpt[idx+samples_.size()] ->Fill((*myptmuon )[0] );
	hists_muoneta[idx+samples_.size()] ->Fill( (*myetamuon )[0] );
	hists_jetpt[idx+samples_.size()] ->Fill((*myptjet )[0]  );
	hists_jeteta[idx+samples_.size()] ->Fill((*myetajet )[0] );
	hists_topmass[idx+samples_.size()] ->Fill( (*mymasstop )[0] );
        hists_toppt[idx+samples_.size()] ->Fill( (*mypttop )[0]);
        hists_topeta[idx+samples_.size()] ->Fill( (*myetatop)[0]);
        hists_topphotonmass[idx+samples_.size()] ->Fill( (*mytopphotonmass)[0]);
        hists_Nvertex[idx+samples_.size()] ->Fill( (*myNvertex)[0]);
	hists_CSV[idx+samples_.size()] ->Fill( (*mycvsdiscriminant)[0]);
	hists_deltaphimetphoton[idx+samples_.size()] ->Fill( (*mydeltaphiphotonmet)[0]);
	hists_deltarphotonmuon[idx+samples_.size()] ->Fill( (*mydeltaRphotonmuon)[0]);
	hists_deltarphotonbjet[idx+samples_.size()] ->Fill( (*mydeltaRphotonjet)[0]);
	hists_jetmulti[idx+samples_.size()] ->Fill( (*myjetmultiplicity)[0]);
	hists_costopphoton[idx+samples_.size()] ->Fill( (*mycostopphoton)[0]);

}
	}
	delete myptphoton;
	delete myetaphoton;
	delete myptmuon;
	delete myetamuon;
	delete myptjet;
	delete myetajet;
	delete mymasstop;
        delete mypttop;
        delete myetatop;
        delete mytopphotonmass;
        delete myNvertex;
	delete input;
	delete mycvsdiscriminant;
	delete mydeltaRphotonjet;
	delete mydeltaRphotonmuon;
	delete mycostopphoton;
	delete mydeltaphiphotonmet;
	delete myjetmultiplicity;
}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
	TFile *input(0);
	TString fname =datasamplesreverse_[idx];
	input = TFile::Open( fname ); 
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myetamuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *myetajet=0;
	std::vector<double> *mymasstop=0;
        std::vector<double> *mypttop=0;
        std::vector<double> *myetatop=0;
        std::vector<double> *mytopphotonmass=0;
	std::vector<double> *myweight=0;
        std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *myjetmultiplicity=0;


	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "etamuon", &myetamuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "etajet", &myetajet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
        theTree->SetBranchAddress( "pttop", &mypttop);
        theTree->SetBranchAddress( "etatop", &myetatop);
        theTree->SetBranchAddress( "topphotonmass", &mytopphotonmass);
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );


	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
//	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
	hists_photonpt[idx+samples_.size()+3] ->Fill( (*myptphoton)[0]);
	hists_photoneta[idx+samples_.size()+3] ->Fill( (*myetaphoton )[0] );
	hists_muonpt[idx+samples_.size()+3] ->Fill((*myptmuon )[0] );
	hists_muoneta[idx+samples_.size()+3] ->Fill( (*myetamuon )[0] );
	hists_jetpt[idx+samples_.size()+3] ->Fill((*myptjet )[0]  );
	hists_jeteta[idx+samples_.size()+3] ->Fill((*myetajet )[0] );
	hists_topmass[idx+samples_.size()+3] ->Fill( (*mymasstop )[0] );
        hists_toppt[idx+samples_.size()+3] ->Fill( (*mypttop )[0]);
        hists_topeta[idx+samples_.size()+3] ->Fill( (*myetatop)[0]);
        hists_topphotonmass[idx+samples_.size()+3] ->Fill( (*mytopphotonmass)[0]);
	hists_CSV[idx+samples_.size()+3] ->Fill( (*mycvsdiscriminant)[0]);
	hists_deltaphimetphoton[idx+samples_.size()+3] ->Fill( (*mydeltaphiphotonmet)[0]);
	hists_deltarphotonmuon[idx+samples_.size()+3] ->Fill( (*mydeltaRphotonmuon)[0]);
	hists_deltarphotonbjet[idx+samples_.size()+3] ->Fill( (*mydeltaRphotonjet)[0]);
	hists_jetmulti[idx+samples_.size()+3] ->Fill( (*myjetmultiplicity)[0]);
	hists_costopphoton[idx+samples_.size()+3] ->Fill( (*mycostopphoton)[0]);

//}
	}
	delete myptphoton;
	delete myetaphoton;
	delete myptmuon;
	delete myetamuon;
	delete myptjet;
	delete myetajet;
	delete mymasstop;
        delete mypttop;
        delete myetatop;
        delete mytopphotonmass;
	delete input;
	delete mycvsdiscriminant;
	delete mydeltaRphotonjet;
	delete mydeltaRphotonmuon;
	delete mycostopphoton;
	delete mydeltaphiphotonmet;
	delete myjetmultiplicity;
}


for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
	TFile *input(0);
	TString fname =samplesreverse_[idx];
	input = TFile::Open( fname ); 
cout<<fname<<endl;
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myetamuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *myetajet=0;
	std::vector<double> *mymasstop=0;
        std::vector<double> *mypttop=0;
        std::vector<double> *myetatop=0;
        std::vector<double> *mytopphotonmass=0;
	std::vector<double> *myweight=0;
        std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *myjetmultiplicity=0;


	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "etamuon", &myetamuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "etajet", &myetajet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
        theTree->SetBranchAddress( "pttop", &mypttop);
        theTree->SetBranchAddress( "etatop", &myetatop);
        theTree->SetBranchAddress( "topphotonmass", &mytopphotonmass);
        theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
	theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
//	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
	hists_photonpt[idx+samples_.size()+6] ->Fill( (*myptphoton)[0],finalweight *scales[idx] );
	hists_photoneta[idx+samples_.size()+6] ->Fill( (*myetaphoton )[0],finalweight *scales[idx]  );
	hists_muonpt[idx+samples_.size()+6] ->Fill((*myptmuon )[0],finalweight *scales[idx]  );
	hists_muoneta[idx+samples_.size()+6] ->Fill( (*myetamuon )[0] ,finalweight *scales[idx] );
	hists_jetpt[idx+samples_.size()+6] ->Fill((*myptjet )[0],finalweight  *scales[idx]  );
	hists_jeteta[idx+samples_.size()+6] ->Fill((*myetajet )[0],finalweight *scales[idx]  );
	hists_topmass[idx+samples_.size()+6] ->Fill( (*mymasstop )[0] ,finalweight *scales[idx] );
        hists_toppt[idx+samples_.size()+6] ->Fill( (*mypttop )[0],finalweight *scales[idx] );
        hists_topeta[idx+samples_.size()+6] ->Fill( (*myetatop)[0],finalweight *scales[idx] );
        hists_topphotonmass[idx+samples_.size()+6] ->Fill( (*mytopphotonmass)[0],finalweight *scales[idx] );
	hists_CSV[idx+samples_.size()+6] ->Fill( (*mycvsdiscriminant)[0],finalweight *scales[idx]);
	hists_deltaphimetphoton[idx+samples_.size()+6] ->Fill( (*mydeltaphiphotonmet)[0],finalweight *scales[idx]);
	hists_deltarphotonmuon[idx+samples_.size()+6] ->Fill( (*mydeltaRphotonmuon)[0],finalweight *scales[idx]);
	hists_deltarphotonbjet[idx+samples_.size()+6] ->Fill( (*mydeltaRphotonjet)[0],finalweight *scales[idx]);
	hists_jetmulti[idx+samples_.size()+6] ->Fill( (*myjetmultiplicity)[0],finalweight *scales[idx]);
	hists_costopphoton[idx+samples_.size()+6] ->Fill( (*mycostopphoton)[0],finalweight *scales[idx]);

//}
	}
	delete myptphoton;
	delete myetaphoton;
	delete myptmuon;
	delete myetamuon;
	delete myptjet;
	delete myetajet;
	delete mymasstop;
        delete mypttop;
        delete myetatop;
        delete mytopphotonmass;
	delete input;
	delete mycvsdiscriminant;
	delete mydeltaRphotonjet;
	delete mydeltaRphotonmuon;
	delete mycostopphoton;
	delete mydeltaphiphotonmet;
	delete myjetmultiplicity;
}

for(unsigned int idx=2; idx<samplesreverse_.size(); ++idx){
	hists_photonpt[idx+samples_.size()+6]->Add(hists_photonpt[idx+samples_.size()+5]); 
	hists_photoneta[idx+samples_.size()+6]->Add(hists_photoneta[idx+samples_.size()+5]); 
	hists_muonpt[idx+samples_.size()+6]->Add(hists_muonpt[idx+samples_.size()+5]);
	hists_muoneta[idx+samples_.size()+6]->Add(hists_muoneta[idx+samples_.size()+5]); 
	hists_jetpt[idx+samples_.size()+6]->Add(hists_jetpt[idx+samples_.size()+5]);
	hists_jeteta[idx+samples_.size()+6]->Add(hists_jeteta[idx+samples_.size()+5]); 
	hists_topmass[idx+samples_.size()+6]->Add(hists_topmass[idx+samples_.size()+5]); 
        hists_toppt[idx+samples_.size()+6]->Add(hists_toppt[idx+samples_.size()+5]);
        hists_topeta[idx+samples_.size()+6]->Add(hists_topeta[idx+samples_.size()+5]);
        hists_topphotonmass[idx+samples_.size()+6]->Add(hists_topphotonmass[idx+samples_.size()+5]);
        hists_Nvertex[idx+samples_.size()+6]->Add(hists_Nvertex[idx+samples_.size()+5]);
        hists_CSV[idx+samples_.size()+6]->Add(hists_CSV[idx+samples_.size()+5]);
	hists_deltaphimetphoton[idx+samples_.size()+6] ->Add(hists_deltaphimetphoton[idx+samples_.size()+5]);
	hists_deltarphotonmuon[idx+samples_.size()+6]->Add(hists_deltarphotonmuon[idx+samples_.size()+5]);
	hists_deltarphotonbjet[idx+samples_.size()+6]->Add(hists_deltarphotonbjet[idx+samples_.size()+5]);
	hists_jetmulti[idx+samples_.size()+6]->Add(hists_jetmulti[idx+samples_.size()+5]);
	hists_costopphoton[idx+samples_.size()+6]->Add(hists_costopphoton[idx+samples_.size()+5]);
}



for(unsigned int idx=1; idx<datasamplesreverse_.size(); ++idx){
	hists_photonpt[idx+samples_.size()+3]->Add(hists_photonpt[idx+samples_.size()+2]); 
	hists_photoneta[idx+samples_.size()+3]->Add(hists_photoneta[idx+samples_.size()+2]); 
	hists_muonpt[idx+samples_.size()+3]->Add(hists_muonpt[idx+samples_.size()+2]);
	hists_muoneta[idx+samples_.size()+3]->Add(hists_muoneta[idx+samples_.size()+2]); 
	hists_jetpt[idx+samples_.size()+3]->Add(hists_jetpt[idx+samples_.size()+2]);
	hists_jeteta[idx+samples_.size()+3]->Add(hists_jeteta[idx+samples_.size()+2]); 
	hists_topmass[idx+samples_.size()+3]->Add(hists_topmass[idx+samples_.size()+2]); 
        hists_toppt[idx+samples_.size()+3]->Add(hists_toppt[idx+samples_.size()+2]);
        hists_topeta[idx+samples_.size()+3]->Add(hists_topeta[idx+samples_.size()+2]);
        hists_topphotonmass[idx+samples_.size()+3]->Add(hists_topphotonmass[idx+samples_.size()+2]);
        hists_Nvertex[idx+samples_.size()+3]->Add(hists_Nvertex[idx+samples_.size()+2]);
        hists_CSV[idx+samples_.size()+3]->Add(hists_CSV[idx+samples_.size()+2]);
	hists_deltaphimetphoton[idx+samples_.size()+3] ->Add(hists_deltaphimetphoton[idx+samples_.size()+2]);
	hists_deltarphotonmuon[idx+samples_.size()+3]->Add(hists_deltarphotonmuon[idx+samples_.size()+2]);
	hists_deltarphotonbjet[idx+samples_.size()+3]->Add(hists_deltarphotonbjet[idx+samples_.size()+2]);
	hists_jetmulti[idx+samples_.size()+3]->Add(hists_jetmulti[idx+samples_.size()+2]);
	hists_costopphoton[idx+samples_.size()+3]->Add(hists_costopphoton[idx+samples_.size()+2]);
}

cout<<"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(1)<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(2)<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(3)<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(4)<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(5)<<endl;
cout<<hists_photonpt[samples_.size()+samplesreverse_.size()+5]->GetBinContent(6)<<endl;




hists_photonpt[samples_.size()+5]->Add(hists_photonpt[samples_.size()+samplesreverse_.size()+5],-1); 
hists_photoneta[samples_.size()+5]->Add(hists_photoneta[samples_.size()+samplesreverse_.size()+5],-1); 
hists_muonpt[samples_.size()+5]->Add(hists_muonpt[samples_.size()+samplesreverse_.size()+5],-1);
hists_muoneta[samples_.size()+5]->Add(hists_muoneta[samples_.size()+samplesreverse_.size()+5],-1); 
hists_jetpt[samples_.size()+5]->Add(hists_jetpt[samples_.size()+samplesreverse_.size()+5],-1);
hists_jeteta[samples_.size()+5]->Add(hists_jeteta[samples_.size()+samplesreverse_.size()+5],-1); 
hists_topmass[samples_.size()+5]->Add(hists_topmass[samples_.size()+samplesreverse_.size()+5],-1); 
hists_toppt[samples_.size()+5]->Add(hists_toppt[samples_.size()+samplesreverse_.size()+5],-1);
hists_topeta[samples_.size()+5]->Add(hists_topeta[samples_.size()+samplesreverse_.size()+5],-1);
hists_topphotonmass[samples_.size()+5]->Add(hists_topphotonmass[samples_.size()+samplesreverse_.size()+5],-1);
hists_CSV[samples_.size()+5]->Add(hists_CSV[samples_.size()+samplesreverse_.size()+5],-1);
	hists_deltaphimetphoton[samples_.size()+5] ->Add(hists_deltaphimetphoton[samples_.size()+samplesreverse_.size()+5],-1);
	hists_deltarphotonmuon[samples_.size()+5]->Add(hists_deltarphotonmuon[samples_.size()+samplesreverse_.size()+5],-1);
	hists_deltarphotonbjet[samples_.size()+5]->Add(hists_deltarphotonbjet[samples_.size()+samplesreverse_.size()+5],-1);
	hists_jetmulti[samples_.size()+5]->Add(hists_jetmulti[samples_.size()+samplesreverse_.size()+5],-1);
	hists_costopphoton[samples_.size()+5]->Add(hists_costopphoton[samples_.size()+samplesreverse_.size()+5],-1);

double nwjet=219.373/hists_muoneta[samples_.size()+5]->Integral();
hists_photonpt[samples_.size()+5]->Scale(nwjet); 
hists_photoneta[samples_.size()+5]->Scale(nwjet); 
hists_muonpt[samples_.size()+5]->Scale(nwjet);
hists_muoneta[samples_.size()+5]->Scale(nwjet); 
hists_jetpt[samples_.size()+5]->Scale(nwjet);
hists_jeteta[samples_.size()+5]->Scale(nwjet); 
hists_topmass[samples_.size()+5]->Scale(nwjet); 
hists_toppt[samples_.size()+5]->Scale(nwjet);
hists_topeta[samples_.size()+5]->Scale(nwjet);
hists_topphotonmass[samples_.size()+5]->Scale(nwjet);
hists_CSV[samples_.size()+5]->Scale(nwjet);
hists_deltaphimetphoton[samples_.size()+5]->Scale(nwjet);
hists_deltarphotonmuon[samples_.size()+5]->Scale(nwjet);
hists_deltarphotonbjet[samples_.size()+5]->Scale(nwjet);
hists_jetmulti[samples_.size()+5]->Scale(nwjet);
hists_costopphoton[samples_.size()+5]->Scale(nwjet);

cout<<"hists_photonpt[samples_.size()+5]"<<hists_photonpt[samples_.size()+5]->Integral()<<endl;
/*
	hists_photonpt[1]->Add(hists_photonpt[samples_.size()+5]); 
	hists_photoneta[1]->Add(hists_photoneta[samples_.size()+5]); 
	hists_muonpt[1]->Add(hists_muonpt[samples_.size()+5]);
	hists_muoneta[1]->Add(hists_muoneta[samples_.size()+5]); 
	hists_jetpt[1]->Add(hists_jetpt[samples_.size()+5]);
	hists_jeteta[1]->Add(hists_jeteta[samples_.size()+5]); 
	hists_topmass[1]->Add(hists_topmass[samples_.size()+5]); 
        hists_toppt[1]->Add(hists_toppt[samples_.size()+5]);
        hists_topeta[1]->Add(hists_topeta[samples_.size()+5]);
        hists_topphotonmass[1]->Add(hists_topphotonmass[samples_.size()+5]);
*/
double nwphjet=1112.2/hists_muoneta[8]->Integral();
hists_photonpt[8]->Scale(1112.2/hists_photonpt[8]->Integral()); 
hists_photoneta[8]->Scale(nwphjet); 
hists_muonpt[8]->Scale(nwphjet);
hists_muoneta[8]->Scale(nwphjet); 
hists_jetpt[8]->Scale(nwphjet);
hists_jeteta[8]->Scale(nwphjet); 
hists_topmass[8]->Scale(nwphjet); 
hists_toppt[8]->Scale(nwphjet);
hists_topeta[8]->Scale(nwphjet);
hists_topphotonmass[8]->Scale(nwphjet);
hists_Nvertex[8]->Scale(nwphjet);
hists_NvertexnoPU[8]->Scale(nwphjet);
hists_CSV[8]->Scale(nwphjet);
hists_deltaphimetphoton[8]->Scale(nwphjet);
hists_deltarphotonmuon[8]->Scale(nwphjet);
hists_deltarphotonbjet[8]->Scale(nwphjet);
hists_jetmulti[8]->Scale(nwphjet);
hists_costopphoton[8]->Scale(nwphjet);

cout<<"hists_photonpt[8]"<<hists_photonpt[8]->Integral()<<endl;
/*
double nsignal=1914.5/hists_muoneta[26]->Integral();
hists_photonpt[26]->Scale(nsignal); 
hists_photoneta[26]->Scale(nsignal); 
hists_muonpt[26]->Scale(nsignal);
hists_muoneta[26]->Scale(nsignal); 
hists_jetpt[26]->Scale(nsignal);
hists_jeteta[26]->Scale(nsignal); 
hists_topmass[26]->Scale(nsignal); 
hists_toppt[26]->Scale(nsignal);
hists_topeta[26]->Scale(nsignal);
hists_topphotonmass[26]->Scale(nsignal);
hists_Nvertex[26]->Scale(nsignal);
hists_NvertexnoPU[26]->Scale(nsignal);
hists_CSV[26]->Scale(nsignal);
hists_deltaphimetphoton[26]->Scale(nsignal);
hists_deltarphotonmuon[26]->Scale(nsignal);
hists_deltarphotonbjet[26]->Scale(nsignal);
hists_jetmulti[26]->Scale(nsignal);
hists_costopphoton[26]->Scale(nsignal);
*/
/*
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
cout<<"sample=  "<<samples_[idx]<<"integral=   "<<hists_photonpt[idx]->Integral()<<endl;
	hists_photonpt[idx]->Add(hists_photonpt[idx-1]); 
	hists_photoneta[idx]->Add(hists_photoneta[idx-1]); 
	hists_muonpt[idx]->Add(hists_muonpt[idx-1]);
	hists_muoneta[idx]->Add(hists_muoneta[idx-1]); 
	hists_jetpt[idx]->Add(hists_jetpt[idx-1]);
	hists_jeteta[idx]->Add(hists_jeteta[idx-1]); 
	hists_topmass[idx]->Add(hists_topmass[idx-1]); 
        hists_toppt[idx]->Add(hists_toppt[idx-1]);
        hists_topeta[idx]->Add(hists_topeta[idx-1]);
        hists_topphotonmass[idx]->Add(hists_topphotonmass[idx-1]);
cout<<"sample=  "<<samples_[idx]<<"integral=   "<<hists_photonpt[idx]->Integral()<<endl;
}
*/
for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
	hists_photonpt[idx+samples_.size()]->Add(hists_photonpt[idx+samples_.size()-1]); 
	hists_photoneta[idx+samples_.size()]->Add(hists_photoneta[idx+samples_.size()-1]); 
	hists_muonpt[idx+samples_.size()]->Add(hists_muonpt[idx+samples_.size()-1]);
	hists_muoneta[idx+samples_.size()]->Add(hists_muoneta[idx+samples_.size()-1]); 
	hists_jetpt[idx+samples_.size()]->Add(hists_jetpt[idx+samples_.size()-1]);
	hists_jeteta[idx+samples_.size()]->Add(hists_jeteta[idx+samples_.size()-1]); 
	hists_topmass[idx+samples_.size()]->Add(hists_topmass[idx+samples_.size()-1]); 
        hists_toppt[idx+samples_.size()]->Add(hists_toppt[idx+samples_.size()-1]);
        hists_topeta[idx+samples_.size()]->Add(hists_topeta[idx+samples_.size()-1]);
        hists_topphotonmass[idx+samples_.size()]->Add(hists_topphotonmass[idx+samples_.size()-1]);
        hists_Nvertex[idx+samples_.size()]->Add(hists_Nvertex[idx+samples_.size()-1]);
        hists_CSV[idx+samples_.size()]->Add(hists_CSV[idx+samples_.size()-1]);
	hists_deltaphimetphoton[idx+samples_.size()] ->Add(hists_deltaphimetphoton[idx+samples_.size()-1]);
	hists_deltarphotonmuon[idx+samples_.size()]->Add(hists_deltarphotonmuon[idx+samples_.size()-1]);
	hists_deltarphotonbjet[idx+samples_.size()]->Add(hists_deltarphotonbjet[idx+samples_.size()-1]);
	hists_jetmulti[idx+samples_.size()]->Add(hists_jetmulti[idx+samples_.size()-1]);
	hists_costopphoton[idx+samples_.size()]->Add(hists_costopphoton[idx+samples_.size()-1]);
}

TH1F *sum_h= new TH1F ( *hists_photonpt[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h->Add(hists_photonpt[idx],1);
}
//sum_h->Add(hists_photonpt[8],-1);
sum_h->Add(hists_photonpt[samples_.size()+5],1);

cout<<" real data = "<< hists_photonpt[samples_.size()+2]->Integral()<<endl;
cout<<" MC  = "<<hists_photonpt[20]->Integral()<<endl;
cout<<" MC  = "<<hists_topphotonmass[20]->Integral()<<endl;
cout<<" real data = "<< hists_topphotonmass[samples_.size()+2]->Integral()<<endl;

cout<<"************************************************************"<<endl;
cout<<" MC OTHER BG = "<<hists_topphotonmass[20]->Integral()-3397-411<<endl;


TLegend* leg = new TLegend(0.4,0.6,0.6,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( hists_photonpt[25], "Single top+#gamma"               , "F");
  leg->AddEntry( hists_photonpt[23], "Z#gamma"               , "F");
  leg->AddEntry( hists_photonpt[22], "Diboson"               , "F");
  leg->AddEntry( hists_photonpt[19], "WW#gamma"               , "F");
  leg->AddEntry( hists_photonpt[18], "t#bar{t}#gamma"               , "F");
  leg->AddEntry( hists_photonpt[17], "t#bar{t}"               , "F");
  leg->AddEntry( hists_photonpt[14], "Single top"               , "F");
  leg->AddEntry( hists_photonpt[8], "W#gamma"              , "F");
  leg->AddEntry( hists_photonpt[7], "#gamma+jets"                           , "F");
  leg->AddEntry( hists_photonpt[1], "Z+jets"                           , "F");
  leg->AddEntry( hists_photonpt[samples_.size()+5], "W+jets"                           , "F");
  leg->AddEntry( hists_photonpt[26], "Signal(t#gammau)"               , "F");
  leg->AddEntry( hists_photonpt[samples_.size()+2], "CMS Data 2012(19.145/fb)"               , "PL");
  leg->AddEntry(sum_h, "MC uncertainty"               , "F");


TCanvas *c1 = new TCanvas("c1","multipads",900,700);
c1->cd(0);

THStack *hs1 = new THStack("hs1","photonpt");
hists_photonpt[samples_.size()+5]->SetFillColor(kBlue-4);
hs1->Add(hists_photonpt[samples_.size()+5]);
hists_photonpt[1]->SetFillColor(kOrange+7);
hs1->Add(hists_photonpt[1]);
hists_photonpt[7]->SetFillColor(19);
hs1->Add(hists_photonpt[7]);
hists_photonpt[8]->SetFillColor(kGreen+1);
hs1->Add(hists_photonpt[8]);
hists_photonpt[10]->Add(hists_photonpt[9]);
hists_photonpt[11]->Add(hists_photonpt[10]);
hists_photonpt[12]->Add(hists_photonpt[11]);
hists_photonpt[13]->Add(hists_photonpt[12]);
hists_photonpt[14]->Add(hists_photonpt[13]);
hists_photonpt[14]->SetFillColor(kAzure+10);
hs1->Add(hists_photonpt[14]);
hists_photonpt[16]->Add(hists_photonpt[15]);
hists_photonpt[17]->Add(hists_photonpt[16]);
hists_photonpt[17]->SetFillColor(6);
hs1->Add(hists_photonpt[17]);
hists_photonpt[18]->SetFillColor(32);
hs1->Add(hists_photonpt[18]);
hists_photonpt[19]->SetFillColor(kSpring-9);
hs1->Add(hists_photonpt[19]);
hists_photonpt[21]->Add(hists_photonpt[20]);
hists_photonpt[22]->Add(hists_photonpt[21]);
hists_photonpt[22]->SetFillColor(kViolet-7);
hs1->Add(hists_photonpt[22]);
hists_photonpt[23]->SetFillColor(kOrange-2);
hs1->Add(hists_photonpt[23]);
hists_photonpt[25]->Add(hists_photonpt[24]);
hists_photonpt[25]->SetFillColor(kYellow+3);
hs1->Add(hists_photonpt[25]);

hs1->Draw("hist");
hs1->SetMaximum(1.2*hists_photonpt[samples_.size()+2]->GetMaximum());
hs1->GetXaxis()->SetTitle("Photon P_{T} (GeV)");
hs1->GetYaxis()->SetTitle("NEVENTS");

hists_photonpt[26]->SetLineColor(kRed+3);
hists_photonpt[26]->SetLineWidth(3);
hists_photonpt[26]->Draw("histsame");

hists_photonpt[samples_.size()+2]->SetLineWidth(3.);
hists_photonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_photonpt[samples_.size()+2]->Draw("esame");

sum_h->SetFillColor(2);
sum_h->SetFillStyle(3013);
sum_h->Draw("e2same");

leg->Draw("same");


TH1F *sum_h2= new TH1F ( *hists_photoneta[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h2->Add(hists_photoneta[idx],1);
}
sum_h2->Add(hists_photoneta[samples_.size()+5],1);

  TCanvas *c2 = new TCanvas("c2","photoneta",900,700);
  c2->cd(0);
THStack *hs2 = new THStack("hs2","photoneta");
hists_photoneta[samples_.size()+5]->SetFillColor(kBlue-4);
hs2->Add(hists_photoneta[samples_.size()+5]);
hists_photoneta[1]->SetFillColor(kOrange+7);
hs2->Add(hists_photoneta[1]);
hists_photoneta[7]->SetFillColor(19);
hs2->Add(hists_photoneta[7]);
hists_photoneta[8]->SetFillColor(kGreen+1);
hs2->Add(hists_photoneta[8]);
hists_photoneta[10]->Add(hists_photoneta[9]);
hists_photoneta[6+5]->Add(hists_photoneta[5+5]);
hists_photoneta[7+5]->Add(hists_photoneta[6+5]);
hists_photoneta[8+5]->Add(hists_photoneta[7+5]);
hists_photoneta[9+5]->Add(hists_photoneta[8+5]);
hists_photoneta[9+5]->SetFillColor(kAzure+10);
hs2->Add(hists_photoneta[9+5]);
hists_photoneta[11+5]->Add(hists_photoneta[10+5]);
hists_photoneta[12+5]->Add(hists_photoneta[11+5]);
hists_photoneta[12+5]->SetFillColor(6);
hs2->Add(hists_photoneta[12+5]);
hists_photoneta[13+5]->SetFillColor(32);
hs2->Add(hists_photoneta[13+5]);
hists_photoneta[14+5]->SetFillColor(kSpring-9);
hs2->Add(hists_photoneta[14+5]);
hists_photoneta[16+5]->Add(hists_photoneta[15+5]);
hists_photoneta[17+5]->Add(hists_photoneta[16+5]);
hists_photoneta[17+5]->SetFillColor(kViolet-7);
hs2->Add(hists_photoneta[17+5]);
hists_photoneta[18+5]->SetFillColor(kOrange-2);
hs2->Add(hists_photoneta[18+5]);
hists_photoneta[20+5]->Add(hists_photoneta[19+5]);
hists_photoneta[20+5]->SetFillColor(kYellow+3);
hs2->Add(hists_photoneta[20+5]);

hs2->Draw("hist");
hs2->SetMaximum(1.2*hists_photoneta[samples_.size()+2]->GetMaximum());
hs2->GetXaxis()->SetTitle("Photon #eta");
hs2->GetYaxis()->SetTitle("NEVENTS");

hists_photoneta[21+5]->SetLineColor(kRed+3);
hists_photoneta[21+5]->SetLineWidth(3);
hists_photoneta[21+5]->Draw("histsame");

hists_photoneta[samples_.size()+2]->SetLineWidth(3.);
hists_photoneta[samples_.size()+2]->SetLineColor(kBlack);
hists_photoneta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photoneta[samples_.size()+2]->SetMarkerStyle(20.);
hists_photoneta[samples_.size()+2]->Draw("esame");
sum_h2->SetFillColor(2);
sum_h2->SetFillStyle(3013);
sum_h2->Draw("e2same");

TH1F *sum_h3= new TH1F ( *hists_muonpt[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h3->Add(hists_muonpt[idx],1);
}
sum_h3->Add(hists_muonpt[samples_.size()+5],1);

  TCanvas *c3 = new TCanvas("c3","multipads",900,700);
  c3->cd(0);
THStack *hs3 = new THStack("hs3","muonpt");
hists_muonpt[samples_.size()+5]->SetFillColor(kBlue-4);
hs3->Add(hists_muonpt[samples_.size()+5]);
hists_muonpt[1]->SetFillColor(kOrange+7);
hs3->Add(hists_muonpt[1]);
hists_muonpt[2+5]->SetFillColor(19);
hs3->Add(hists_muonpt[2+5]);
hists_muonpt[3+5]->SetFillColor(kGreen+1);
hs3->Add(hists_muonpt[3+5]);
hists_muonpt[5+5]->Add(hists_muonpt[4+5]);
hists_muonpt[6+5]->Add(hists_muonpt[5+5]);
hists_muonpt[7+5]->Add(hists_muonpt[6+5]);
hists_muonpt[8+5]->Add(hists_muonpt[7+5]);
hists_muonpt[9+5]->Add(hists_muonpt[8+5]);
hists_muonpt[9+5]->SetFillColor(kAzure+10);
hs3->Add(hists_muonpt[9+5]);
hists_muonpt[11+5]->Add(hists_muonpt[10+5]);
hists_muonpt[12+5]->Add(hists_muonpt[11+5]);
hists_muonpt[12+5]->SetFillColor(6);
hs3->Add(hists_muonpt[12+5]);
hists_muonpt[13+5]->SetFillColor(32);
hs3->Add(hists_muonpt[13+5]);
hists_muonpt[14+5]->SetFillColor(kSpring-9);
hs3->Add(hists_muonpt[14+5]);
hists_muonpt[16+5]->Add(hists_muonpt[15+5]);
hists_muonpt[17+5]->Add(hists_muonpt[16+5]);
hists_muonpt[17+5]->SetFillColor(kViolet-7);
hs3->Add(hists_muonpt[17+5]);
hists_muonpt[18+5]->SetFillColor(kOrange-2);
hs3->Add(hists_muonpt[18+5]);
hists_muonpt[20+5]->Add(hists_muonpt[19+5]);
hists_muonpt[20+5]->SetFillColor(kYellow+3);
hs3->Add(hists_muonpt[20+5]);

hs3->Draw("hist");
hs3->SetMaximum(1.2*hists_muonpt[samples_.size()+2]->GetMaximum());
hs3->GetXaxis()->SetTitle("muon P_{T} (GeV)");
hs3->GetYaxis()->SetTitle("NEVENTS");

hists_muonpt[21+5]->SetLineColor(kRed+3);
hists_muonpt[21+5]->SetLineWidth(3);
hists_muonpt[21+5]->Draw("histsame");

hists_muonpt[samples_.size()+2]->SetLineWidth(3.);
hists_muonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_muonpt[samples_.size()+2]->Draw("esame");

sum_h3->SetFillColor(2);
sum_h3->SetFillStyle(3013);
sum_h3->Draw("e2same");

TH1F *sum_h4= new TH1F ( *hists_muoneta[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h4->Add(hists_muoneta[idx],1);
}
sum_h4->Add(hists_muoneta[samples_.size()+5],1);

  TCanvas *c4 = new TCanvas("c4","multipads",900,700);
  c4->cd(0);
THStack *hs4 = new THStack("hs4","muoneta");
hists_muoneta[samples_.size()+5]->SetFillColor(kBlue-4);
hs4->Add(hists_muoneta[samples_.size()+5]);
hists_muoneta[1]->SetFillColor(kOrange+7);
hs4->Add(hists_muoneta[1]);
hists_muoneta[2+5]->SetFillColor(19);
hs4->Add(hists_muoneta[2+5]);
hists_muoneta[3+5]->SetFillColor(kGreen+1);
hs4->Add(hists_muoneta[3+5]);
hists_muoneta[5+5]->Add(hists_muoneta[4+5]);
hists_muoneta[6+5]->Add(hists_muoneta[5+5]);
hists_muoneta[7+5]->Add(hists_muoneta[6+5]);
hists_muoneta[8+5]->Add(hists_muoneta[7+5]);
hists_muoneta[9+5]->Add(hists_muoneta[8+5]);
hists_muoneta[9+5]->SetFillColor(kAzure+10);
hs4->Add(hists_muoneta[9+5]);
hists_muoneta[11+5]->Add(hists_muoneta[10+5]);
hists_muoneta[12+5]->Add(hists_muoneta[11+5]);
hists_muoneta[12+5]->SetFillColor(6);
hs4->Add(hists_muoneta[12+5]);
hists_muoneta[13+5]->SetFillColor(32);
hs4->Add(hists_muoneta[13+5]);
hists_muoneta[14+5]->SetFillColor(kSpring-9);
hs4->Add(hists_muoneta[14+5]);
hists_muoneta[16+5]->Add(hists_muoneta[15+5]);
hists_muoneta[17+5]->Add(hists_muoneta[16+5]);
hists_muoneta[17+5]->SetFillColor(kViolet-7);
hs4->Add(hists_muoneta[17+5]);
hists_muoneta[18+5]->SetFillColor(kOrange-2);
hs4->Add(hists_muoneta[18+5]);
hists_muoneta[20+5]->Add(hists_muoneta[19+5]);
hists_muoneta[20+5]->SetFillColor(kYellow+3);
hs4->Add(hists_muoneta[20+5]);

hs4->Draw("hist");
hs4->SetMaximum(1.2*hists_muoneta[samples_.size()+2]->GetMaximum());
hs4->GetXaxis()->SetTitle("muon #eta");
hs4->GetYaxis()->SetTitle("NEVENTS");

hists_muoneta[21+5]->SetLineColor(kRed+3);
hists_muoneta[21+5]->SetLineWidth(3);
hists_muoneta[21+5]->Draw("histsame");

hists_muoneta[samples_.size()+2]->SetLineWidth(3.);
hists_muoneta[samples_.size()+2]->SetLineColor(kBlack);
hists_muoneta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muoneta[samples_.size()+2]->SetMarkerStyle(20.);
hists_muoneta[samples_.size()+2]->Draw("esame");
sum_h4->SetFillColor(2);
sum_h4->SetFillStyle(3013);
sum_h4->Draw("e2same");

TH1F *sum_h5= new TH1F ( *hists_jetpt[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h5->Add(hists_jetpt[idx],1);
}
sum_h5->Add(hists_jetpt[samples_.size()+5],1);
  TCanvas *c5 = new TCanvas("c5","multipads",900,700);
  c5->cd(0);
THStack *hs5 = new THStack("hs5","jetpt");
hists_jetpt[samples_.size()+5]->SetFillColor(kBlue-4);
hs5->Add(hists_jetpt[samples_.size()+5]);
hists_jetpt[1]->SetFillColor(kOrange+7);
hs5->Add(hists_jetpt[1]);
hists_jetpt[2+5]->SetFillColor(19);
hs5->Add(hists_jetpt[2+5]);
hists_jetpt[3+5]->SetFillColor(kGreen+1);
hs5->Add(hists_jetpt[3+5]);
hists_jetpt[5+5]->Add(hists_jetpt[4+5]);
hists_jetpt[6+5]->Add(hists_jetpt[5+5]);
hists_jetpt[7+5]->Add(hists_jetpt[6+5]);
hists_jetpt[8+5]->Add(hists_jetpt[7+5]);
hists_jetpt[9+5]->Add(hists_jetpt[8+5]);
hists_jetpt[9+5]->SetFillColor(kAzure+10);
hs5->Add(hists_jetpt[9+5]);
hists_jetpt[11+5]->Add(hists_jetpt[10+5]);
hists_jetpt[12+5]->Add(hists_jetpt[11+5]);
hists_jetpt[12+5]->SetFillColor(6);
hs5->Add(hists_jetpt[12+5]);
hists_jetpt[13+5]->SetFillColor(32);
hs5->Add(hists_jetpt[13+5]);
hists_jetpt[14+5]->SetFillColor(kSpring-9);
hs5->Add(hists_jetpt[14+5]);
hists_jetpt[16+5]->Add(hists_jetpt[15+5]);
hists_jetpt[17+5]->Add(hists_jetpt[16+5]);
hists_jetpt[17+5]->SetFillColor(kViolet-7);
hs5->Add(hists_jetpt[17+5]);
hists_jetpt[18+5]->SetFillColor(kOrange-2);
hs5->Add(hists_jetpt[18+5]);
hists_jetpt[20+5]->Add(hists_jetpt[19+5]);
hists_jetpt[20+5]->SetFillColor(kYellow+3);
hs5->Add(hists_jetpt[20+5]);


hs5->Draw("hist");
hs5->SetMaximum(1.2*hists_jetpt[samples_.size()+2]->GetMaximum());
hs5->GetXaxis()->SetTitle("b-jet P_{T} (GeV)");
hs5->GetYaxis()->SetTitle("NEVENTS");

hists_jetpt[21+5]->SetLineColor(kRed+3);
hists_jetpt[21+5]->SetLineWidth(3);
hists_jetpt[21+5]->Draw("histsame");

hists_jetpt[samples_.size()+2]->SetLineWidth(3.);
hists_jetpt[samples_.size()+2]->SetLineColor(kBlack);
hists_jetpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_jetpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_jetpt[samples_.size()+2]->Draw("esame");
sum_h5->SetFillColor(2);
sum_h5->SetFillStyle(3013);
sum_h5->Draw("e2same");

 
TH1F *sum_h6= new TH1F ( *hists_jeteta[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h6->Add(hists_jeteta[idx],1);
}
sum_h6->Add(hists_jeteta[samples_.size()+5],1);

 TCanvas *c6 = new TCanvas("c6","multipads",900,700);
  c6->cd(0);
THStack *hs6 = new THStack("hs6","jeteta");
hists_jeteta[samples_.size()+5]->SetFillColor(kBlue-4);
hs6->Add(hists_jeteta[samples_.size()+5]);
hists_jeteta[1]->SetFillColor(kOrange+7);
hs6->Add(hists_jeteta[1]);
hists_jeteta[2+5]->SetFillColor(19);
hs6->Add(hists_jeteta[2+5]);
hists_jeteta[3+5]->SetFillColor(kGreen+1);
hs6->Add(hists_jeteta[3+5]);
hists_jeteta[5+5]->Add(hists_jeteta[4+5]);
hists_jeteta[6+5]->Add(hists_jeteta[5+5]);
hists_jeteta[7+5]->Add(hists_jeteta[6+5]);
hists_jeteta[8+5]->Add(hists_jeteta[7+5]);
hists_jeteta[9+5]->Add(hists_jeteta[8+5]);
hists_jeteta[9+5]->SetFillColor(kAzure+10);
hs6->Add(hists_jeteta[9+5]);
hists_jeteta[11+5]->Add(hists_jeteta[10+5]);
hists_jeteta[12+5]->Add(hists_jeteta[11+5]);
hists_jeteta[12+5]->SetFillColor(6);
hs6->Add(hists_jeteta[12+5]);
hists_jeteta[13+5]->SetFillColor(32);
hs6->Add(hists_jeteta[13+5]);
hists_jeteta[14+5]->SetFillColor(kSpring-9);
hs6->Add(hists_jeteta[14+5]);
hists_jeteta[16+5]->Add(hists_jeteta[15+5]);
hists_jeteta[17+5]->Add(hists_jeteta[16+5]);
hists_jeteta[17+5]->SetFillColor(kViolet-7);
hs6->Add(hists_jeteta[17+5]);
hists_jeteta[18+5]->SetFillColor(kOrange-2);
hs6->Add(hists_jeteta[18+5]);
hists_jeteta[20+5]->Add(hists_jeteta[19+5]);
hists_jeteta[20+5]->SetFillColor(kYellow+3);
hs6->Add(hists_jeteta[20+5]);

hs6->Draw("hist");
hs6->SetMaximum(1.05*hists_jeteta[26]->GetMaximum());
hs6->GetXaxis()->SetTitle("b-jet #eta");
hs6->GetYaxis()->SetTitle("NEVENTS");

hists_jeteta[21+5]->SetLineColor(kRed+3);
hists_jeteta[21+5]->SetLineWidth(3);
hists_jeteta[21+5]->Draw("histsame");

hists_jeteta[samples_.size()+2]->SetLineWidth(3.);
hists_jeteta[samples_.size()+2]->SetLineColor(kBlack);
hists_jeteta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_jeteta[samples_.size()+2]->SetMarkerStyle(20.);
hists_jeteta[samples_.size()+2]->Draw("esame");
sum_h6->SetFillColor(2);
sum_h6->SetFillStyle(3013);
sum_h6->Draw("e2same");


  TCanvas *c7 = new TCanvas("c7","multipads",900,700);
  c7->cd(0);
THStack *hs7 = new THStack("hs7","topmass");
hists_topmass[samples_.size()+5]->SetFillColor(kBlue-4);
hs7->Add(hists_topmass[samples_.size()+5]);
hists_topmass[1]->SetFillColor(kOrange+7);
hs7->Add(hists_topmass[1]);
hists_topmass[2]->SetFillColor(19);
hs7->Add(hists_topmass[2]);
hists_topmass[3]->SetFillColor(kGreen+1);
hs7->Add(hists_topmass[3]);
hists_topmass[5]->Add(hists_topmass[4]);
hists_topmass[6]->Add(hists_topmass[5]);
hists_topmass[7]->Add(hists_topmass[6]);
hists_topmass[8]->Add(hists_topmass[7]);
hists_topmass[9]->Add(hists_topmass[8]);
hists_topmass[9]->SetFillColor(kAzure+10);
hs7->Add(hists_topmass[9]);
hists_topmass[11]->Add(hists_topmass[10]);
hists_topmass[12]->Add(hists_topmass[11]);
hists_topmass[12]->SetFillColor(6);
hs7->Add(hists_topmass[12]);
hists_topmass[13]->SetFillColor(32);
hs7->Add(hists_topmass[13]);
hists_topmass[14]->SetFillColor(kSpring-9);
hs7->Add(hists_topmass[14]);
hists_topmass[16]->Add(hists_topmass[15]);
hists_topmass[17]->Add(hists_topmass[16]);
hists_topmass[17]->SetFillColor(kViolet-7);
hs7->Add(hists_topmass[17]);
hists_topmass[18]->SetFillColor(kOrange-2);
hs7->Add(hists_topmass[18]);
hists_topmass[20]->Add(hists_topmass[19]);
hists_topmass[20]->SetFillColor(kYellow+3);
hs7->Add(hists_topmass[20]);


hs7->Draw();
hs7->SetMaximum(1.2*hists_topmass[samples_.size()+2]->GetMaximum());
hs7->GetXaxis()->SetTitle("m_{top}");
hs7->GetYaxis()->SetTitle("NEVENTS");

hists_topmass[21]->SetLineColor(kRed+3);
hists_topmass[21]->SetLineWidth(3);
hists_topmass[21]->Draw("same");

hists_topmass[samples_.size()+2]->SetLineWidth(3.);
hists_topmass[samples_.size()+2]->SetLineColor(kBlack);
hists_topmass[samples_.size()+2]->SetMarkerColor(kBlack);
hists_topmass[samples_.size()+2]->SetMarkerStyle(20.);
hists_topmass[samples_.size()+2]->Draw("esame");
leg->Draw("same");

TLine *line1 = new TLine(130, 0, 130, hists_topmass[samples_.size()+2]->GetMaximum());
    line1->SetLineColor(1);
    line1->SetLineWidth(2);
    line1->Draw("same");
TLine *line2 = new TLine(220, 0, 220, hists_topmass[samples_.size()+2]->GetMaximum());
    line2->SetLineColor(1);
    line2->SetLineWidth(2);
    line2->Draw("same");
TLine *line3 = new TLine(130, hists_topmass[samples_.size()+2]->GetMaximum(), 220, hists_topmass[samples_.size()+2]->GetMaximum());
    line3->SetLineColor(1);
    line3->SetLineWidth(2);
    line3->Draw("same");

  TCanvas *c8 = new TCanvas("c8","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c8->cd(0);
hists_toppt[20]->SetMaximum(1.5*hists_toppt[samples_.size()+2]->GetMaximum());
hists_toppt[20]->SetFillColor(kYellow+3);
hists_toppt[20]->Draw();
hists_toppt[18]->SetFillColor(kOrange-2);
hists_toppt[18]->Draw("same");
hists_toppt[17]->SetFillColor(kViolet-7);
hists_toppt[17]->Draw("same");
//hists_toppt[16]->SetFillColor(kRed);
//hists_toppt[16]->Draw("same");
//hists_toppt[15]->SetFillColor(kViolet+1);
//hists_toppt[15]->Draw("same");
hists_toppt[14]->SetFillColor(kSpring-9);
hists_toppt[14]->Draw("same");
hists_toppt[13]->SetFillColor(32);
hists_toppt[13]->Draw("same");
hists_toppt[12]->SetFillColor(6);
hists_toppt[12]->Draw("same");
hists_toppt[9]->SetFillColor(kAzure+10);
hists_toppt[9]->Draw("same");
hists_toppt[3]->SetFillColor(kGreen+1);
hists_toppt[3]->Draw("same");
hists_toppt[2]->SetFillColor(19);
hists_toppt[2]->Draw("same");
hists_toppt[1]->SetFillColor(kOrange+7);
hists_toppt[1]->Draw("same");
hists_toppt[samples_.size()+5]->SetFillColor(kBlue-4);
hists_toppt[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_toppt[21]->SetLineColor(kRed+3);
hists_toppt[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_toppt[21]->Draw("same");

hists_toppt[samples_.size()+2]->SetLineWidth(3.);
hists_toppt[samples_.size()+2]->SetLineColor(kBlack);
hists_toppt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_toppt[samples_.size()+2]->SetMarkerStyle(20.);
hists_toppt[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c9 = new TCanvas("c9","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c9->cd(0);
hists_topeta[20]->SetMaximum(1.5*hists_topeta[samples_.size()+2]->GetMaximum());
hists_topeta[20]->SetFillColor(kYellow+3);
hists_topeta[20]->Draw();
hists_topeta[18]->SetFillColor(kOrange-2);
hists_topeta[18]->Draw("same");
hists_topeta[17]->SetFillColor(kViolet-7);
hists_topeta[17]->Draw("same");
//hists_topeta[16]->SetFillColor(kRed);
//hists_topeta[16]->Draw("same");
//hists_topeta[15]->SetFillColor(kViolet+1);
//hists_topeta[15]->Draw("same");
hists_topeta[14]->SetFillColor(kSpring-9);
hists_topeta[14]->Draw("same");
hists_topeta[13]->SetFillColor(32);
hists_topeta[13]->Draw("same");
hists_topeta[12]->SetFillColor(6);
hists_topeta[12]->Draw("same");
hists_topeta[9]->SetFillColor(kAzure+10);
hists_topeta[9]->Draw("same");
hists_topeta[3]->SetFillColor(kGreen+1);
hists_topeta[3]->Draw("same");
hists_topeta[2]->SetFillColor(19);
hists_topeta[2]->Draw("same");
hists_topeta[1]->SetFillColor(kOrange+7);
hists_topeta[1]->Draw("same");
hists_topeta[samples_.size()+5]->SetFillColor(kBlue-4);
hists_topeta[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_topeta[21]->SetLineColor(kRed+3);
hists_topeta[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_topeta[21]->Draw("same");

hists_topeta[samples_.size()+2]->SetLineWidth(3.);
hists_topeta[samples_.size()+2]->SetLineColor(kBlack);
hists_topeta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_topeta[samples_.size()+2]->SetMarkerStyle(20.);
hists_topeta[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c10 = new TCanvas("c10","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c10->cd(0);
hists_topphotonmass[20]->SetMaximum(1.5*hists_topphotonmass[samples_.size()+2]->GetMaximum());
hists_topphotonmass[20]->SetFillColor(kYellow+3);
hists_topphotonmass[20]->Draw();
hists_topphotonmass[18]->SetFillColor(kOrange-2);
hists_topphotonmass[18]->Draw("same");
hists_topphotonmass[17]->SetFillColor(kViolet-7);
hists_topphotonmass[17]->Draw("same");
//hists_topphotonmass[16]->SetFillColor(kRed);
//hists_topphotonmass[16]->Draw("same");
//hists_topphotonmass[15]->SetFillColor(kViolet+1);
//hists_topphotonmass[15]->Draw("same");
hists_topphotonmass[14]->SetFillColor(kSpring-9);
hists_topphotonmass[14]->Draw("same");
hists_topphotonmass[13]->SetFillColor(32);
hists_topphotonmass[13]->Draw("same");
hists_topphotonmass[12]->SetFillColor(6);
hists_topphotonmass[12]->Draw("same");
hists_topphotonmass[9]->SetFillColor(kAzure+10);
hists_topphotonmass[9]->Draw("same");
hists_topphotonmass[3]->SetFillColor(kGreen+1);
hists_topphotonmass[3]->Draw("same");
hists_topphotonmass[2]->SetFillColor(19);
hists_topphotonmass[2]->Draw("same");
hists_topphotonmass[1]->SetFillColor(kOrange+7);
hists_topphotonmass[1]->Draw("same");
hists_topphotonmass[samples_.size()+5]->SetFillColor(kBlue-4);
hists_topphotonmass[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_topphotonmass[21]->SetLineColor(kRed+3);
hists_topphotonmass[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_topphotonmass[21]->Draw("same");

hists_topphotonmass[samples_.size()+2]->SetLineWidth(3.);
hists_topphotonmass[samples_.size()+2]->SetLineColor(kBlack);
hists_topphotonmass[samples_.size()+2]->SetMarkerColor(kBlack);
hists_topphotonmass[samples_.size()+2]->SetMarkerStyle(20.);
hists_topphotonmass[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

TCanvas *c11 = new TCanvas("c11","multipads",900,700);
c11->cd(0);
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
        hists_Nvertex[idx]->Add(hists_Nvertex[idx-1]);
        hists_NvertexnoPU[idx]->Add(hists_NvertexnoPU[idx-1]);
}
hists_Nvertex[samples_.size()+2]->Scale(1/hists_Nvertex[samples_.size()+2]->Integral());
hists_Nvertex[20]->Scale(1/hists_Nvertex[20]->Integral());
hists_NvertexnoPU[20]->Scale(1/hists_NvertexnoPU[20]->Integral());

hists_Nvertex[20]->SetLineColor(kRed+3);
hists_NvertexnoPU[20]->SetLineColor(kGreen+3);
hists_Nvertex[samples_.size()+2]->SetLineColor(kBlue+3);
hists_Nvertex[20]->Draw("");
hists_NvertexnoPU[20]->Draw("same");
hists_Nvertex[samples_.size()+2]->Draw("same");

TH1F *sum_h12= new TH1F ( *hists_CSV[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h12->Add(hists_CSV[idx],1);
}
sum_h12->Add(hists_CSV[samples_.size()+5],1);

  TCanvas *c12 = new TCanvas("c12","CSV",900,700);
  c12->cd(0);
THStack *hs12 = new THStack("hs12","CSV");
hists_CSV[samples_.size()+5]->SetFillColor(kBlue-4);
hs12->Add(hists_CSV[samples_.size()+5]);
hists_CSV[1]->SetFillColor(kOrange+7);
hs12->Add(hists_CSV[1]);
hists_CSV[3]->Add(hists_CSV[2]);
hists_CSV[4]->Add(hists_CSV[3]);
hists_CSV[5]->Add(hists_CSV[4]);
hists_CSV[6]->Add(hists_CSV[5]);
hists_CSV[7]->Add(hists_CSV[6]);
hists_CSV[7]->SetFillColor(19);
hs12->Add(hists_CSV[7]);
hists_CSV[8]->SetFillColor(kGreen+1);
hs12->Add(hists_CSV[8]);
hists_CSV[10]->Add(hists_CSV[9]);
hists_CSV[6+5]->Add(hists_CSV[5+5]);
hists_CSV[7+5]->Add(hists_CSV[6+5]);
hists_CSV[8+5]->Add(hists_CSV[7+5]);
hists_CSV[9+5]->Add(hists_CSV[8+5]);
hists_CSV[9+5]->SetFillColor(kAzure+10);
hs12->Add(hists_CSV[9+5]);
hists_CSV[11+5]->Add(hists_CSV[10+5]);
hists_CSV[12+5]->Add(hists_CSV[11+5]);
hists_CSV[12+5]->SetFillColor(6);
hs12->Add(hists_CSV[12+5]);
hists_CSV[13+5]->SetFillColor(32);
hs12->Add(hists_CSV[13+5]);
hists_CSV[14+5]->SetFillColor(kSpring-9);
hs12->Add(hists_CSV[14+5]);
hists_CSV[16+5]->Add(hists_CSV[15+5]);
hists_CSV[17+5]->Add(hists_CSV[16+5]);
hists_CSV[17+5]->SetFillColor(kViolet-7);
hs12->Add(hists_CSV[17+5]);
hists_CSV[18+5]->SetFillColor(kOrange-2);
hs12->Add(hists_CSV[18+5]);
hists_CSV[20+5]->Add(hists_CSV[19+5]);
hists_CSV[20+5]->SetFillColor(kYellow+3);
hs12->Add(hists_CSV[20+5]);

hs12->Draw("hist");
hs12->SetMaximum(1.2*hists_CSV[samples_.size()+2]->GetMaximum());
hs12->GetXaxis()->SetTitle("CSV discriminator");
hs12->GetYaxis()->SetTitle("NEVENTS");

hists_CSV[21+5]->SetLineColor(kRed+3);
hists_CSV[21+5]->SetLineWidth(3);
hists_CSV[21+5]->Draw("histsame");

hists_CSV[samples_.size()+2]->SetLineWidth(3.);
hists_CSV[samples_.size()+2]->SetLineColor(kBlack);
hists_CSV[samples_.size()+2]->SetMarkerColor(kBlack);
hists_CSV[samples_.size()+2]->SetMarkerStyle(20.);
hists_CSV[samples_.size()+2]->Draw("esame");
sum_h12->SetFillColor(2);
sum_h12->SetFillStyle(3013);
sum_h12->Draw("e2same");

TH1F *sum_h13= new TH1F ( *hists_deltaphimetphoton[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h13->Add(hists_deltaphimetphoton[idx],1);
}
sum_h13->Add(hists_deltaphimetphoton[samples_.size()+5],1);

  TCanvas *c13 = new TCanvas("c13","deltaphimetphoton",900,700);
  c13->cd(0);
THStack *hs13 = new THStack("hs13","deltaphimetphoton");
hists_deltaphimetphoton[samples_.size()+5]->SetFillColor(kBlue-4);
hs13->Add(hists_deltaphimetphoton[samples_.size()+5]);
hists_deltaphimetphoton[1]->SetFillColor(kOrange+7);
hs13->Add(hists_deltaphimetphoton[1]);
hists_deltaphimetphoton[3]->Add(hists_deltaphimetphoton[2]);
hists_deltaphimetphoton[4]->Add(hists_deltaphimetphoton[3]);
hists_deltaphimetphoton[5]->Add(hists_deltaphimetphoton[4]);
hists_deltaphimetphoton[6]->Add(hists_deltaphimetphoton[5]);
hists_deltaphimetphoton[7]->Add(hists_deltaphimetphoton[6]);
hists_deltaphimetphoton[7]->SetFillColor(19);
hs13->Add(hists_deltaphimetphoton[7]);
hists_deltaphimetphoton[8]->SetFillColor(kGreen+1);
hs13->Add(hists_deltaphimetphoton[8]);
hists_deltaphimetphoton[10]->Add(hists_deltaphimetphoton[9]);
hists_deltaphimetphoton[6+5]->Add(hists_deltaphimetphoton[5+5]);
hists_deltaphimetphoton[7+5]->Add(hists_deltaphimetphoton[6+5]);
hists_deltaphimetphoton[8+5]->Add(hists_deltaphimetphoton[7+5]);
hists_deltaphimetphoton[9+5]->Add(hists_deltaphimetphoton[8+5]);
hists_deltaphimetphoton[9+5]->SetFillColor(kAzure+10);
hs13->Add(hists_deltaphimetphoton[9+5]);
hists_deltaphimetphoton[11+5]->Add(hists_deltaphimetphoton[10+5]);
hists_deltaphimetphoton[12+5]->Add(hists_deltaphimetphoton[11+5]);
hists_deltaphimetphoton[12+5]->SetFillColor(6);
hs13->Add(hists_deltaphimetphoton[12+5]);
hists_deltaphimetphoton[13+5]->SetFillColor(32);
hs13->Add(hists_deltaphimetphoton[13+5]);
hists_deltaphimetphoton[14+5]->SetFillColor(kSpring-9);
hs13->Add(hists_deltaphimetphoton[14+5]);
hists_deltaphimetphoton[16+5]->Add(hists_deltaphimetphoton[15+5]);
hists_deltaphimetphoton[17+5]->Add(hists_deltaphimetphoton[16+5]);
hists_deltaphimetphoton[17+5]->SetFillColor(kViolet-7);
hs13->Add(hists_deltaphimetphoton[17+5]);
hists_deltaphimetphoton[18+5]->SetFillColor(kOrange-2);
hs13->Add(hists_deltaphimetphoton[18+5]);
hists_deltaphimetphoton[20+5]->Add(hists_deltaphimetphoton[19+5]);
hists_deltaphimetphoton[20+5]->SetFillColor(kYellow+3);
hs13->Add(hists_deltaphimetphoton[20+5]);

hs13->Draw("hist");
hs13->SetMaximum(1.05*hists_deltaphimetphoton[26]->GetMaximum());
hs13->GetXaxis()->SetTitle("#Delta#phi (Photon,MET)");
hs13->GetYaxis()->SetTitle("NEVENTS");

hists_deltaphimetphoton[21+5]->SetLineColor(kRed+3);
hists_deltaphimetphoton[21+5]->SetLineWidth(3);
hists_deltaphimetphoton[21+5]->Draw("histsame");

hists_deltaphimetphoton[samples_.size()+2]->SetLineWidth(3.);
hists_deltaphimetphoton[samples_.size()+2]->SetLineColor(kBlack);
hists_deltaphimetphoton[samples_.size()+2]->SetMarkerColor(kBlack);
hists_deltaphimetphoton[samples_.size()+2]->SetMarkerStyle(20.);
hists_deltaphimetphoton[samples_.size()+2]->Draw("esame");
sum_h13->SetFillColor(2);
sum_h13->SetFillStyle(3013);
sum_h13->Draw("e2same");

TH1F *sum_h14= new TH1F ( *hists_deltarphotonmuon[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h14->Add(hists_deltarphotonmuon[idx],1);
}
sum_h14->Add(hists_deltarphotonmuon[samples_.size()+5],1);

  TCanvas *c14 = new TCanvas("c14","deltarphotonmuon",900,700);
  c14->cd(0);
THStack *hs14 = new THStack("hs14","deltarphotonmuon");
hists_deltarphotonmuon[samples_.size()+5]->SetFillColor(kBlue-4);
hs14->Add(hists_deltarphotonmuon[samples_.size()+5]);
hists_deltarphotonmuon[1]->SetFillColor(kOrange+7);
hs14->Add(hists_deltarphotonmuon[1]);
hists_deltarphotonmuon[3]->Add(hists_deltarphotonmuon[2]);
hists_deltarphotonmuon[4]->Add(hists_deltarphotonmuon[3]);
hists_deltarphotonmuon[5]->Add(hists_deltarphotonmuon[4]);
hists_deltarphotonmuon[6]->Add(hists_deltarphotonmuon[5]);
hists_deltarphotonmuon[7]->Add(hists_deltarphotonmuon[6]);
hists_deltarphotonmuon[7]->SetFillColor(19);
hs14->Add(hists_deltarphotonmuon[7]);
hists_deltarphotonmuon[8]->SetFillColor(kGreen+1);
hs14->Add(hists_deltarphotonmuon[8]);
hists_deltarphotonmuon[10]->Add(hists_deltarphotonmuon[9]);
hists_deltarphotonmuon[6+5]->Add(hists_deltarphotonmuon[5+5]);
hists_deltarphotonmuon[7+5]->Add(hists_deltarphotonmuon[6+5]);
hists_deltarphotonmuon[8+5]->Add(hists_deltarphotonmuon[7+5]);
hists_deltarphotonmuon[9+5]->Add(hists_deltarphotonmuon[8+5]);
hists_deltarphotonmuon[9+5]->SetFillColor(kAzure+10);
hs14->Add(hists_deltarphotonmuon[9+5]);
hists_deltarphotonmuon[11+5]->Add(hists_deltarphotonmuon[10+5]);
hists_deltarphotonmuon[12+5]->Add(hists_deltarphotonmuon[11+5]);
hists_deltarphotonmuon[12+5]->SetFillColor(6);
hs14->Add(hists_deltarphotonmuon[12+5]);
hists_deltarphotonmuon[13+5]->SetFillColor(32);
hs14->Add(hists_deltarphotonmuon[13+5]);
hists_deltarphotonmuon[14+5]->SetFillColor(kSpring-9);
hs14->Add(hists_deltarphotonmuon[14+5]);
hists_deltarphotonmuon[16+5]->Add(hists_deltarphotonmuon[15+5]);
hists_deltarphotonmuon[17+5]->Add(hists_deltarphotonmuon[16+5]);
hists_deltarphotonmuon[17+5]->SetFillColor(kViolet-7);
hs14->Add(hists_deltarphotonmuon[17+5]);
hists_deltarphotonmuon[18+5]->SetFillColor(kOrange-2);
hs14->Add(hists_deltarphotonmuon[18+5]);
hists_deltarphotonmuon[20+5]->Add(hists_deltarphotonmuon[19+5]);
hists_deltarphotonmuon[20+5]->SetFillColor(kYellow+3);
hs14->Add(hists_deltarphotonmuon[20+5]);

hs14->Draw("hist");
hs14->SetMaximum(1.05*hists_deltarphotonmuon[26]->GetMaximum());
hs14->GetXaxis()->SetTitle("#DeltaR (Photon,muon)");
hs14->GetYaxis()->SetTitle("NEVENTS");

hists_deltarphotonmuon[21+5]->SetLineColor(kRed+3);
hists_deltarphotonmuon[21+5]->SetLineWidth(3);
hists_deltarphotonmuon[21+5]->Draw("histsame");

hists_deltarphotonmuon[samples_.size()+2]->SetLineWidth(3.);
hists_deltarphotonmuon[samples_.size()+2]->SetLineColor(kBlack);
hists_deltarphotonmuon[samples_.size()+2]->SetMarkerColor(kBlack);
hists_deltarphotonmuon[samples_.size()+2]->SetMarkerStyle(20.);
hists_deltarphotonmuon[samples_.size()+2]->Draw("esame");
sum_h14->SetFillColor(2);
sum_h14->SetFillStyle(3013);
sum_h14->Draw("e2same");

TH1F *sum_h15= new TH1F ( *hists_deltarphotonbjet[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h15->Add(hists_deltarphotonbjet[idx],1);
}
sum_h15->Add(hists_deltarphotonbjet[samples_.size()+5],1);

  TCanvas *c15 = new TCanvas("c15","deltarphotonbjet",900,700);
  c15->cd(0);
THStack *hs15 = new THStack("hs15","deltarphotonbjet");
hists_deltarphotonbjet[samples_.size()+5]->SetFillColor(kBlue-4);
hs15->Add(hists_deltarphotonbjet[samples_.size()+5]);
hists_deltarphotonbjet[1]->SetFillColor(kOrange+7);
hs15->Add(hists_deltarphotonbjet[1]);
hists_deltarphotonbjet[3]->Add(hists_deltarphotonbjet[2]);
hists_deltarphotonbjet[4]->Add(hists_deltarphotonbjet[3]);
hists_deltarphotonbjet[5]->Add(hists_deltarphotonbjet[4]);
hists_deltarphotonbjet[6]->Add(hists_deltarphotonbjet[5]);
hists_deltarphotonbjet[7]->Add(hists_deltarphotonbjet[6]);
hists_deltarphotonbjet[7]->SetFillColor(19);
hs15->Add(hists_deltarphotonbjet[7]);
hists_deltarphotonbjet[8]->SetFillColor(kGreen+1);
hs15->Add(hists_deltarphotonbjet[8]);
hists_deltarphotonbjet[10]->Add(hists_deltarphotonbjet[9]);
hists_deltarphotonbjet[6+5]->Add(hists_deltarphotonbjet[5+5]);
hists_deltarphotonbjet[7+5]->Add(hists_deltarphotonbjet[6+5]);
hists_deltarphotonbjet[8+5]->Add(hists_deltarphotonbjet[7+5]);
hists_deltarphotonbjet[9+5]->Add(hists_deltarphotonbjet[8+5]);
hists_deltarphotonbjet[9+5]->SetFillColor(kAzure+10);
hs15->Add(hists_deltarphotonbjet[9+5]);
hists_deltarphotonbjet[11+5]->Add(hists_deltarphotonbjet[10+5]);
hists_deltarphotonbjet[12+5]->Add(hists_deltarphotonbjet[11+5]);
hists_deltarphotonbjet[12+5]->SetFillColor(6);
hs15->Add(hists_deltarphotonbjet[12+5]);
hists_deltarphotonbjet[13+5]->SetFillColor(32);
hs15->Add(hists_deltarphotonbjet[13+5]);
hists_deltarphotonbjet[14+5]->SetFillColor(kSpring-9);
hs15->Add(hists_deltarphotonbjet[14+5]);
hists_deltarphotonbjet[16+5]->Add(hists_deltarphotonbjet[15+5]);
hists_deltarphotonbjet[17+5]->Add(hists_deltarphotonbjet[16+5]);
hists_deltarphotonbjet[17+5]->SetFillColor(kViolet-7);
hs15->Add(hists_deltarphotonbjet[17+5]);
hists_deltarphotonbjet[18+5]->SetFillColor(kOrange-2);
hs15->Add(hists_deltarphotonbjet[18+5]);
hists_deltarphotonbjet[20+5]->Add(hists_deltarphotonbjet[19+5]);
hists_deltarphotonbjet[20+5]->SetFillColor(kYellow+3);
hs15->Add(hists_deltarphotonbjet[20+5]);

hs15->Draw("hist");
hs15->SetMaximum(1.05*hists_deltarphotonbjet[26]->GetMaximum());
hs15->GetXaxis()->SetTitle("#DeltaR (Photon,b-jet)");
hs15->GetYaxis()->SetTitle("NEVENTS");

hists_deltarphotonbjet[21+5]->SetLineColor(kRed+3);
hists_deltarphotonbjet[21+5]->SetLineWidth(3);
hists_deltarphotonbjet[21+5]->Draw("histsame");

hists_deltarphotonbjet[samples_.size()+2]->SetLineWidth(3.);
hists_deltarphotonbjet[samples_.size()+2]->SetLineColor(kBlack);
hists_deltarphotonbjet[samples_.size()+2]->SetMarkerColor(kBlack);
hists_deltarphotonbjet[samples_.size()+2]->SetMarkerStyle(20.);
hists_deltarphotonbjet[samples_.size()+2]->Draw("esame");
sum_h15->SetFillColor(2);
sum_h15->SetFillStyle(3013);
sum_h15->Draw("e2same");

TH1F *sum_h16= new TH1F ( *hists_jetmulti[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h16->Add(hists_jetmulti[idx],1);
}
sum_h16->Add(hists_jetmulti[samples_.size()+5],1);

  TCanvas *c16 = new TCanvas("c16","jetmulti",900,700);
  c16->cd(0);
THStack *hs16 = new THStack("hs16","jetmulti");
hists_jetmulti[samples_.size()+5]->SetFillColor(kBlue-4);
hs16->Add(hists_jetmulti[samples_.size()+5]);
hists_jetmulti[1]->SetFillColor(kOrange+7);
hs16->Add(hists_jetmulti[1]);
hists_jetmulti[3]->Add(hists_jetmulti[2]);
hists_jetmulti[4]->Add(hists_jetmulti[3]);
hists_jetmulti[5]->Add(hists_jetmulti[4]);
hists_jetmulti[6]->Add(hists_jetmulti[5]);
hists_jetmulti[7]->Add(hists_jetmulti[6]);
hists_jetmulti[7]->SetFillColor(19);
hs16->Add(hists_jetmulti[7]);
hists_jetmulti[8]->SetFillColor(kGreen+1);
hs16->Add(hists_jetmulti[8]);
hists_jetmulti[10]->Add(hists_jetmulti[9]);
hists_jetmulti[6+5]->Add(hists_jetmulti[5+5]);
hists_jetmulti[7+5]->Add(hists_jetmulti[6+5]);
hists_jetmulti[8+5]->Add(hists_jetmulti[7+5]);
hists_jetmulti[9+5]->Add(hists_jetmulti[8+5]);
hists_jetmulti[9+5]->SetFillColor(kAzure+10);
hs16->Add(hists_jetmulti[9+5]);
hists_jetmulti[11+5]->Add(hists_jetmulti[10+5]);
hists_jetmulti[12+5]->Add(hists_jetmulti[11+5]);
hists_jetmulti[12+5]->SetFillColor(6);
hs16->Add(hists_jetmulti[12+5]);
hists_jetmulti[13+5]->SetFillColor(32);
hs16->Add(hists_jetmulti[13+5]);
hists_jetmulti[14+5]->SetFillColor(kSpring-9);
hs16->Add(hists_jetmulti[14+5]);
hists_jetmulti[16+5]->Add(hists_jetmulti[15+5]);
hists_jetmulti[17+5]->Add(hists_jetmulti[16+5]);
hists_jetmulti[17+5]->SetFillColor(kViolet-7);
hs16->Add(hists_jetmulti[17+5]);
hists_jetmulti[18+5]->SetFillColor(kOrange-2);
hs16->Add(hists_jetmulti[18+5]);
hists_jetmulti[20+5]->Add(hists_jetmulti[19+5]);
hists_jetmulti[20+5]->SetFillColor(kYellow+3);
hs16->Add(hists_jetmulti[20+5]);

hs16->Draw("hist");
hs16->SetMaximum(1.05*hists_jetmulti[26]->GetMaximum());
hs16->GetXaxis()->SetTitle("jet multiplicity");
hs16->GetYaxis()->SetTitle("NEVENTS");

hists_jetmulti[21+5]->SetLineColor(kRed+3);
hists_jetmulti[21+5]->SetLineWidth(3);
hists_jetmulti[21+5]->Draw("histsame");

hists_jetmulti[samples_.size()+2]->SetLineWidth(3.);
hists_jetmulti[samples_.size()+2]->SetLineColor(kBlack);
hists_jetmulti[samples_.size()+2]->SetMarkerColor(kBlack);
hists_jetmulti[samples_.size()+2]->SetMarkerStyle(20.);
hists_jetmulti[samples_.size()+2]->Draw("esame");
sum_h16->SetFillColor(2);
sum_h16->SetFillStyle(3013);
sum_h16->Draw("e2same");
leg->Draw("same");
/*
	hists_deltaphimetphoton[idx+samples_.size()] ->Add(hists_deltaphimetphoton[idx+samples_.size()-1]);
	hists_deltarphotonmuon[idx+samples_.size()]->Add(hists_deltarphotonmuon[idx+samples_.size()-1]);
	hists_deltarphotonbjet[idx+samples_.size()]->Add(hists_deltarphotonbjet[idx+samples_.size()-1]);
	hists_jetmulti[idx+samples_.size()]->Add(hists_jetmulti[idx+samples_.size()-1]);
	hists_costopphoton[idx+samples_.size()]->Add(hists_costopphoton[idx+samples_.size()-1]);
*/

TH1F *sum_h17= new TH1F ( *hists_costopphoton[1] ) ;
for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
sum_h17->Add(hists_costopphoton[idx],1);
}
sum_h17->Add(hists_costopphoton[samples_.size()+5],1);

  TCanvas *c17 = new TCanvas("c17","costopphoton",900,700);
  c17->cd(0);
THStack *hs17 = new THStack("hs17","costopphoton");
hists_costopphoton[samples_.size()+5]->SetFillColor(kBlue-4);
hs17->Add(hists_costopphoton[samples_.size()+5]);
hists_costopphoton[1]->SetFillColor(kOrange+7);
hs17->Add(hists_costopphoton[1]);
hists_costopphoton[3]->Add(hists_costopphoton[2]);
hists_costopphoton[4]->Add(hists_costopphoton[3]);
hists_costopphoton[5]->Add(hists_costopphoton[4]);
hists_costopphoton[6]->Add(hists_costopphoton[5]);
hists_costopphoton[7]->Add(hists_costopphoton[6]);
hists_costopphoton[7]->SetFillColor(19);
hs17->Add(hists_costopphoton[7]);
hists_costopphoton[8]->SetFillColor(kGreen+1);
hs17->Add(hists_costopphoton[8]);
hists_costopphoton[10]->Add(hists_costopphoton[9]);
hists_costopphoton[6+5]->Add(hists_costopphoton[5+5]);
hists_costopphoton[7+5]->Add(hists_costopphoton[6+5]);
hists_costopphoton[8+5]->Add(hists_costopphoton[7+5]);
hists_costopphoton[9+5]->Add(hists_costopphoton[8+5]);
hists_costopphoton[9+5]->SetFillColor(kAzure+10);
hs17->Add(hists_costopphoton[9+5]);
hists_costopphoton[11+5]->Add(hists_costopphoton[10+5]);
hists_costopphoton[12+5]->Add(hists_costopphoton[11+5]);
hists_costopphoton[12+5]->SetFillColor(6);
hs17->Add(hists_costopphoton[12+5]);
hists_costopphoton[13+5]->SetFillColor(32);
hs17->Add(hists_costopphoton[13+5]);
hists_costopphoton[14+5]->SetFillColor(kSpring-9);
hs17->Add(hists_costopphoton[14+5]);
hists_costopphoton[16+5]->Add(hists_costopphoton[15+5]);
hists_costopphoton[17+5]->Add(hists_costopphoton[16+5]);
hists_costopphoton[17+5]->SetFillColor(kViolet-7);
hs17->Add(hists_costopphoton[17+5]);
hists_costopphoton[18+5]->SetFillColor(kOrange-2);
hs17->Add(hists_costopphoton[18+5]);
hists_costopphoton[20+5]->Add(hists_costopphoton[19+5]);
hists_costopphoton[20+5]->SetFillColor(kYellow+3);
hs17->Add(hists_costopphoton[20+5]);

hs17->Draw("hist");
hs17->SetMaximum(1.05*hists_costopphoton[26]->GetMaximum());
hs17->GetXaxis()->SetTitle("Cos(top,photon)");
hs17->GetYaxis()->SetTitle("NEVENTS");

hists_costopphoton[21+5]->SetLineColor(kRed+3);
hists_costopphoton[21+5]->SetLineWidth(3);
hists_costopphoton[21+5]->Draw("histsame");

hists_costopphoton[samples_.size()+2]->SetLineWidth(3.);
hists_costopphoton[samples_.size()+2]->SetLineColor(kBlack);
hists_costopphoton[samples_.size()+2]->SetMarkerColor(kBlack);
hists_costopphoton[samples_.size()+2]->SetMarkerStyle(20.);
hists_costopphoton[samples_.size()+2]->Draw("esame");
sum_h17->SetFillColor(2);
sum_h17->SetFillStyle(3013);
sum_h17->Draw("e2same");
}






