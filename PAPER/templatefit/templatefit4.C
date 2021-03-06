#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
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
#include <cmath> 

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"


using namespace RooFit;

void templatefit4() {

std::vector<string> samples_;
float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139};

double newweight=19.768/19.145;
//single top and single top + gamma is shifted 50%
//float scales[] = {0.0961,0.0978,1.491,1.3*0.0253,1.3*0.0224,1.3*0.0145,1.3*0.0125,1.3*0.0160,1.3*0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,1.3*0.032,1.3*0.024};

//float scales[] = {0.0961,0.0978,1.491,0.7*0.0253,0.7*0.0224,0.7*0.0145,0.7*0.0125,0.7*0.0160,0.7*0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.7*0.032,0.7*0.024};


//Diboson processes are shifted 50%
//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,1.3*0.0017,1.3*0.0055,1.3*0.0032,0.00084,0.02,0.032,0.024};

//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.7*0.0017,0.7*0.0055,0.7*0.0032,0.00084,0.02,0.032,0.024};


//ttbar and ttbar+gamma is shifted 50%
//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,1.3*0.0341,1.3*0.0341,1.3*0.0341,1.3*0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024};

//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.7*0.0341,0.7*0.0341,0.7*0.0341,0.7*0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024};

//DY is shifted 50%
//float scales[] = {0.0961,1.3*0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024};

//float scales[] = {0.0961,0.7*0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024};

//z + gamma is shifted 50% 
//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,1.3*0.02,0.032,0.024};

//float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.7*0.02,0.032,0.024};

samples_.push_back("WPHJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("PHJET200400.root");
//samples_.push_back("WJET.root");
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

//samples_.push_back("SIGNAL.root");
//samples_.push_back("REALDATA.root");

std::vector<string> samplesreverse_;
samplesreverse_.push_back("etarev/WPHJET.root");
samplesreverse_.push_back("etarev/ZJET.root");
samplesreverse_.push_back("etarev/PHJET200400.root");
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

std::vector<string> datasamples_;
std::vector<string> datasamplesreverse_;
//datasamples_.push_back("REALDATA1.root");
//datasamples_.push_back("REALDATA2.root");
//datasamples_.push_back("REALDATA3.root");
//datasamplesreverse_.push_back("etarev/REALDATA1.root");
//datasamplesreverse_.push_back("etarev/REALDATA2.root");
//datasamplesreverse_.push_back("etarev/REALDATA3.root");

datasamples_.push_back("REALDATAab.root");
datasamples_.push_back("REALDATAc.root");
datasamples_.push_back("REALDATAd.root");

datasamplesreverse_.push_back("etarev/REALDATAab.root");
datasamplesreverse_.push_back("etarev/REALDATAc.root");
datasamplesreverse_.push_back("etarev/REALDATAd.root");



//open files

std::vector<TH1F*> BGSR;
std::vector<TH1F*> DATASR;
std::vector<TH1F*> BGCR;
std::vector<TH1F*> DATACR;
int nbin=22;

TFile *input(0);
Double_t ptphoton,ptmuon,ptjet,masstop,mtw,weight, coswph,pttop;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
TH1F   *recophotonPt, *recomuonPt, *recojetPt, *topmass, *WTmass, *Coswphoton ,*top_Pt;
recophotonPt= new TH1F("photon_Pt", "Pt", 55,50., 500.);
recomuonPt= new TH1F("muon_Pt", "Pt", 25, 26., 400.);
recojetPt= new TH1F("Jet_Pt", "Pt", 25, 30., 500.);
topmass= new TH1F("topmass", "topmass", 50, 0., 400.);
top_Pt= new TH1F("top_Pt", "top_Pt" , 30,0,400);
Coswphoton= new TH1F("Coswphoton","Coswphoton",nbin ,-1 ,1);
//Coswphoton->Sumw2();
WTmass = new TH1F("WTmass", "WTmass", 30, 0., 250.);
TString fname =samples_[idx];
input = TFile::Open( fname );
std::vector<double> *myptphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *mymasstop=0;
std::vector<double> *mypttop=0;
std::vector<double> *mymtw=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myweight=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mysigmaietaieta=0;

TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
   theTree->SetBranchAddress( "pttop" , &mypttop);
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
   theTree->SetBranchAddress( "mtw", &mymtw );
    theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      theTree->GetEntry(ievt);

ptphoton= (*myptphoton)[0];
ptmuon=(*myptmuon)[0];
ptjet=(*myptjet)[0];
masstop=(*mymasstop)[0];
pttop=(*mypttop)[0];
mtw=(*mymtw)[0];
coswph=(*mycoswphoton)[0];
weight=(*myweight)[0];
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {

recophotonPt->Fill(ptphoton,weight);
recomuonPt->Fill( ptmuon,weight);
recojetPt->Fill( ptjet,weight);
topmass->Fill( masstop,weight);
top_Pt->Fill( pttop,weight);
Coswphoton->Fill(coswph ,weight*newweight);
WTmass->Fill(mtw,weight);}

else {if(samples_[idx]=="WPHJET.root") Coswphoton->Fill(coswph ,weight*newweight);}
}
}
}
//BGSR.push_back(WTmass);
//BGSR.push_back(recophotonPt);
BGSR.push_back(Coswphoton);
//BGSR.push_back(top_Pt);
}


for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
TH1F   *recophotonPt, *recomuonPt, *recojetPt, *topmass, *WTmass, *Coswphoton,*top_Pt;
recophotonPt= new TH1F("photon_Pt", "Pt", 55,50., 500.);
recophotonPt->Sumw2();
recomuonPt= new TH1F("muon_Pt", "Pt", 25, 26., 400.);
recomuonPt->Sumw2();
recojetPt= new TH1F("Jet_Pt", "Pt", 25, 30., 500.);
recojetPt->Sumw2();
topmass= new TH1F("topmass", "topmass", 50, 0., 400.);
topmass->Sumw2();
top_Pt= new TH1F("top_Pt", "top_Pt" , 30,0,400);
top_Pt->Sumw2();
Coswphoton= new TH1F("Coswphoton","Coswphoton",nbin ,-1 ,1);
//Coswphoton->Sumw2();
WTmass = new TH1F("WTmass", "WTmass", 30, 0., 250.);
WTmass->Sumw2();
TString fname =datasamples_[idx];
input = TFile::Open( fname );
std::vector<double> *myptphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *mymasstop=0;
std::vector<double> *mypttop=0;
std::vector<double> *mymtw=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myweight=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mysigmaietaieta=0;

TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
   theTree->SetBranchAddress( "pttop" , &mypttop);
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
   theTree->SetBranchAddress( "mtw", &mymtw );
    theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      theTree->GetEntry(ievt);

ptphoton= (*myptphoton)[0];
ptmuon=(*myptmuon)[0];
ptjet=(*myptjet)[0];
masstop=(*mymasstop)[0];
pttop=(*mypttop)[0];
mtw=(*mymtw)[0];
coswph=(*mycoswphoton)[0];
weight=(*myweight)[0];
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){

if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {
recophotonPt->Fill(ptphoton,weight);
recomuonPt->Fill( ptmuon,weight);
recojetPt->Fill( ptjet,weight);
topmass->Fill( masstop,weight);
top_Pt->Fill( pttop,weight);
Coswphoton->Fill(coswph );
WTmass->Fill(mtw,weight);}}}}
//DATASR.push_back(WTmass);
//DATASR.push_back(recophotonPt);
DATASR.push_back(Coswphoton);
//DATASR.push_back(top_Pt);

}

for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
TH1F   *recophotonPt, *recomuonPt, *recojetPt, *topmass, *WTmass, *Coswphoton,*top_Pt;
recophotonPt= new TH1F("photon_Pt", "Pt", 55,50., 500.);
recomuonPt= new TH1F("muon_Pt", "Pt", 25, 26., 400.);
recojetPt= new TH1F("Jet_Pt", "Pt", 25, 30., 500.);
topmass= new TH1F("topmass", "topmass", 50, 0., 400.);
top_Pt= new TH1F("top_Pt", "top_Pt" , 30,0,400);
Coswphoton= new TH1F("Coswphoton","Coswphoton",nbin ,-1 ,1);
WTmass = new TH1F("WTmass", "WTmass", 30, 0., 250.);
TString fname =samplesreverse_[idx];
input = TFile::Open( fname );
std::vector<double> *myptphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *mymasstop=0;
std::vector<double> *mypttop=0;
std::vector<double> *mymtw=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myweight=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mysigmaietaieta=0;

TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
   theTree->SetBranchAddress( "pttop" , &mypttop);
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
   theTree->SetBranchAddress( "mtw", &mymtw );
    theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      theTree->GetEntry(ievt);

ptphoton= (*myptphoton)[0];
ptmuon=(*myptmuon)[0];
ptjet=(*myptjet)[0];
masstop=(*mymasstop)[0];
pttop=(*mypttop)[0];
mtw=(*mymtw)[0];
coswph=(*mycoswphoton)[0];
weight=(*myweight)[0];
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {
recophotonPt->Fill(ptphoton,weight);
recomuonPt->Fill( ptmuon,weight);
recojetPt->Fill( ptjet,weight);
topmass->Fill( masstop,weight);
top_Pt->Fill( pttop,weight);
Coswphoton->Fill(coswph ,weight*newweight);
WTmass->Fill(mtw,weight);}}}}}
//BGCR.push_back(WTmass);
//BGCR.push_back(recophotonPt);
BGCR.push_back(Coswphoton);
//BGCR.push_back(top_Pt);

}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
TH1F   *recophotonPt, *recomuonPt, *recojetPt, *topmass, *WTmass, *Coswphoton,*top_Pt;
recophotonPt= new TH1F("photon_Pt", "Pt", 55,50., 500.);
recomuonPt= new TH1F("muon_Pt", "Pt", 25, 26., 400.);
recojetPt= new TH1F("Jet_Pt", "Pt", 25, 30., 500.);
topmass= new TH1F("topmass", "topmass", 50, 0., 400.);
top_Pt= new TH1F("top_Pt", "top_Pt" , 30,0,400);
Coswphoton= new TH1F("Coswphoton","Coswphoton",nbin ,-1 ,1);
WTmass = new TH1F("WTmass", "WTmass", 30, 0., 250.);
TString fname =datasamplesreverse_[idx];
input = TFile::Open( fname );
std::vector<double> *myptphoton=0;
std::vector<double> *myptmuon=0;
std::vector<double> *myptjet=0;
std::vector<double> *mymasstop=0;
std::vector<double> *mypttop=0;
std::vector<double> *mymtw=0;
std::vector<double> *mycoswphoton=0;
std::vector<double> *myweight=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mysigmaietaieta=0;

TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

   theTree->SetBranchAddress("ptphoton", &myptphoton  );
   theTree->SetBranchAddress( "ptmuon", &myptmuon );
   theTree->SetBranchAddress( "ptjet", &myptjet );
   theTree->SetBranchAddress( "masstop", &mymasstop );
   theTree->SetBranchAddress( "pttop" , &mypttop);
   theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
   theTree->SetBranchAddress( "mtw", &mymtw );
    theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );

   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      theTree->GetEntry(ievt);

ptphoton= (*myptphoton)[0];
ptmuon=(*myptmuon)[0];
ptjet=(*myptjet)[0];
masstop=(*mymasstop)[0];
pttop=(*mypttop)[0];
mtw=(*mymtw)[0];
coswph=(*mycoswphoton)[0];
weight=(*myweight)[0];
if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
//if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.03255)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.01155)) {
//if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.0341)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.0121)) {
recophotonPt->Fill(ptphoton,weight);
recomuonPt->Fill( ptmuon,weight);
recojetPt->Fill( ptjet,weight);
topmass->Fill( masstop,weight);
top_Pt->Fill( pttop,weight);
Coswphoton->Fill(coswph );
WTmass->Fill(mtw,weight);}}
}}
}

//DATACR.push_back(WTmass);
//DATACR.push_back(recophotonPt);
DATACR.push_back(Coswphoton);
//TACR.push_back(top_Pt);



}

float lumi = 1;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
BGSR[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
BGCR[idx]->Scale(lumi*scales[idx]);
}

TH1F* wphjethists = BGSR[0];

for(unsigned int idx=2; idx<samples_.size(); ++idx){
   BGSR[idx]->Add(BGSR[idx-1]);
 }
for(unsigned int idx=1; idx<samplesreverse_.size(); ++idx){
BGCR[idx]->Add(BGCR[idx-1]);
 }
for(unsigned int idx=1; idx<datasamplesreverse_.size(); ++idx){
DATACR[idx]->Add(DATACR[idx-1]);
 }
for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
DATASR[idx]->Add(DATASR[idx-1]);
 }
cout<<"real data num  =  "<<DATASR[2]->Integral()<<endl;
TH1F* SIGhists;

//this is W+Jet shape which is driven from data
SIGhists=DATACR[2];
//SIGhists->Add(revBGhists[samples_.size()-1],-1);
//SIGhists->Add(BGCR[samplesreverse_.size()-1],-1);

//TCanvas* conv = new TCanvas("BGBG", "BG" , 600, 600);
//conv->cd(0);
//DATAhists->Draw();
//revBGhists[samples_.size()-1]->Draw("same");

//TCanvas* nv = new TCanvas("data", "data" , 600, 600);
//nv->cd(0);
//SIGhists->Draw();

for(unsigned int idx=1; idx<nbin+1; ++idx){
//SIGhists->SetBinContent(idx,SIGhists->GetBinContent(idx)+SIGhists->GetBinError(idx));
}
for(unsigned int idx=1; idx<nbin+1; ++idx){
//if (wphjethists->GetBinContent(idx)-wphjethists->GetBinError(idx)>0) wphjethists->SetBinContent(idx,wphjethists->GetBinContent(idx)+wphjethists->GetBinError(idx));
}


int norm=1;
double scalesig =norm/SIGhists->Integral();
SIGhists->Scale(scalesig);



//double scaledata =norm/DATAhists->Integral();
//DATAhists->Scale(scaledata);
double NOBG=BGSR[samples_.size()-1]->Integral();
double scalebg =norm/BGSR[samples_.size()-1]->Integral();
BGSR[samples_.size()-1]->Scale(scalebg);

double scalewphjet =norm/wphjethists->Integral();
wphjethists->Scale(scalewphjet);
//TCanvas* conv = new TCanvas("signal", "signal" , 600, 600);
//conv->cd(0);
//SIGhists->Draw();

//TCanvas* nv = new TCanvas("data", "data" , 600, 600);
//nv->cd(0);
//DATAhists->Draw();

//TCanvas* con = new TCanvas("BG", "BG" , 600, 600);
//con->cd(0);
//BGhists[samples_.size()-1]->Draw();


//conv->cd(0);
////Ghists->Draw();

//    RooRealVar x("x", "x", 50.,500.);
//    RooRealVar x("x", "x", -3,3);
//    RooRealVar x("x", "x", 26.,500.);
//      RooRealVar x("x", "x", 0.,250.);
//      RooRealVar x("x", "x", 0.,400.);
    RooRealVar x("x", "x", -1,1);

//wjet template
    RooDataHist* template1 = new RooDataHist("template1", "template1", RooArgList(x), SIGhists);
    RooHistPdf modelTemplate1("modelTemplate1", "modelTemplate1", RooArgSet(x), *template1);
//other BG template
    RooDataHist* template2 = new RooDataHist("template2", "template2", RooArgList(x), BGSR[samples_.size()-1]);
    RooHistPdf modelTemplate2("modelTemplate2", "modelTemplate2", RooArgSet(x), *template2);

//w+ph+jet template
    RooDataHist* template3 = new RooDataHist("template3", "template3", RooArgList(x), wphjethists);
    RooHistPdf modelTemplate3("modelTemplate3", "modelTemplate3", RooArgSet(x), *template3);
//DATA histogram
    RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), DATASR[2]);

    std::cout << " executing fit..." << std::endl;
    std::cout << " Number of data events ..." <<DATASR[2]->Integral()<< std::endl;
    std::cout << " Number of MC events ..." <<1/scalebg<< std::endl;
    std::cout << " Number of other BG ..." <<NOBG<< std::endl;

    RooRealVar coeffTemplate1wjet("coeffTemplate1wjet", "coeffTemplate1wjet", 0, 5000);
    RooRealVar coeffTemplate2BG("coeffTemplate2BG", "coeffTemplate2BG",NOBG, NOBG-(NOBG/3), NOBG+(NOBG/3));
    RooRealVar coeffTemplate3wphjet("coeffTemplate3wphjet", "coeffTemplate3wphjet", 0, 5000);
 //   coeffTemplate2BG.setConstant(kTRUE);

    RooAddPdf sumTemplate("sumTemplate", "coeff*modelTemplate1 + (1 - coeff)*modelTemplate2",RooArgList(modelTemplate1, modelTemplate2, modelTemplate3),RooArgList(coeffTemplate1wjet,coeffTemplate2BG, coeffTemplate3wphjet));
    sumTemplate.fitTo(*sumData);

std::cout << "--> coeffTemplate1wjet = " << coeffTemplate1wjet.getVal() << std::endl;
std::cout << "--> coeffTemplate2BG = " << coeffTemplate2BG.getVal() << std::endl;
std::cout << "--> coeffTemplate3wphjet = " << coeffTemplate3wphjet.getVal() << std::endl;

//SIGhists->Draw();
//wphjethists->SetLineColor(kGreen);
//wphjethists->Draw("same");
RooPlot* frame = x.frame() ;
sumData->plotOn(frame,LineStyle(kDashed),LineWidth(3)) ;
sumTemplate.plotOn(frame,LineWidth(3)) ;

//modelTemplate2.plotOn(frame) ;

//frame->Draw() ;

            TH1D *  h = new TH1D("dummy_fit","dummy_fit",100,0,100);
	    //	    h->SetLineColor(kblue);
	    h->SetLineWidth(3);
	    //sumData->SetLineWidth(3);

BGSR[samples_.size()-1]->Scale(NOBG);
SIGhists->Scale(coeffTemplate1wjet.getVal());
wphjethists->Scale(coeffTemplate3wphjet.getVal());
SIGhists->Add(BGSR[samples_.size()-1]);
wphjethists->Add(SIGhists);
TCanvas *c1 = new TCanvas("c1","signal region",200,10,700,500);
wphjethists->SetMinimum(0);
wphjethists->SetMaximum(1.2*DATASR[2]->GetMaximum());
wphjethists->SetFillColor(kGreen+1);
wphjethists->GetXaxis()->SetTitle("cos(w,#gamma )");
wphjethists->GetYaxis()->SetTitle("Nevents");
wphjethists->SetTitle("");
wphjethists->Draw() ;
SIGhists->SetFillColor(kBlue-4);
SIGhists->Draw("same") ;
BGSR[samples_.size()-1]->SetFillColor(12);
BGSR[samples_.size()-1]->Draw("same") ;
DATASR[2]->SetLineWidth(3.);
DATASR[2]->SetLineColor(kBlack);
DATASR[2]->SetMarkerColor(kBlack);
DATASR[2]->SetMarkerStyle(20.);
DATASR[2]->Draw("esame");
TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry(wphjethists , "W PH JET"               , "F");
  leg->AddEntry(SIGhists , "W JET" , "F");
  leg->AddEntry(BGSR[samples_.size()-1] , "OTHER BG" , "F");
  leg->AddEntry(DATASR[2], "CMS Data 2012(19.145/fb)" , "PL");
  leg->Draw("same");
}
