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

void dataMCinCR() {
 	  	  	
// list of valid histogram names


// list of input samples

std::vector<string> datasamplesreverse_;
float scales[] = {1,0.0978,34.01,6.133,1.04,0.32,0.02,0.002,0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,
0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139,5*(0.049094905/19.145)};


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

std::vector<TH1F*> hists_sigmaietaieta;
std::vector<TH1F*> hists_HOE;




for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
hists_photonpt.push_back(new  TH1F(  std::string("photon_pt").append(samplesreverse_[idx]).c_str(), std::string("photon_pt").append(samplesreverse_[idx]).c_str() , 12,0., 300.));
hists_photoneta.push_back(new  TH1F(  std::string("photon_eta").append(samplesreverse_[idx]).c_str(),  std::string("photon_eta").append(samplesreverse_[idx]).c_str() , 22, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(samplesreverse_[idx]).c_str(), std::string("muon_pt").append(samplesreverse_[idx]).c_str() , 12, 0, 312));
hists_muoneta.push_back(new  TH1F(  std::string("muon_eta").append(samplesreverse_[idx]).c_str(), std::string("muon_eta").append(samplesreverse_[idx]).c_str() ,  12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(samplesreverse_[idx]).c_str(),  std::string("jet_pt").append(samplesreverse_[idx]).c_str(),12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(samplesreverse_[idx]).c_str(),  std::string("bjet_eta").append(samplesreverse_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(samplesreverse_[idx]).c_str(), std::string("top_mass").append(samplesreverse_[idx]).c_str() , 20, 70, 370.));
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
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(samplesreverse_[idx]).c_str(), std::string("costopphoton").append(samplesreverse_[idx]).c_str() ,22, -1, 1));
hists_costopphoton[idx]->Sumw2();

hists_sigmaietaieta.push_back(new TH1F(  std::string("sigmaietaieta").append(samplesreverse_[idx]).c_str(), std::string("sigmaietaieta").append(samplesreverse_[idx]).c_str() ,44, 0.001, 0.111));
hists_sigmaietaieta[idx]->Sumw2();

hists_HOE.push_back(new TH1F(  std::string("HOE").append(samplesreverse_[idx]).c_str(), std::string("HOE").append(samplesreverse_[idx]).c_str() ,8, 0, 0.05));
hists_HOE[idx]->Sumw2();
}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
hists_photonpt.push_back(new  TH1F(  std::string("photon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("photon_pt").append(datasamplesreverse_[idx]).c_str() , 12,0., 300.));
hists_photoneta.push_back(new  TH1F(  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str() , 22, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("muon_pt").append(datasamplesreverse_[idx]).c_str() ,12, 0, 312));
hists_muoneta.push_back(new  TH1F(  std::string("muon_eta").append(datasamplesreverse_[idx]).c_str(), std::string("muon_eta").append(datasamplesreverse_[idx]).c_str() , 12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamplesreverse_[idx]).c_str(), std::string("top_mass").append(datasamplesreverse_[idx]).c_str() , 20, 70, 370.));
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
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(datasamplesreverse_[idx]).c_str(), std::string("costopphoton").append(datasamplesreverse_[idx]).c_str() ,22, -1, 1));
hists_costopphoton[idx]->Sumw2();

hists_sigmaietaieta.push_back(new TH1F(  std::string("sigmaietaieta").append(datasamplesreverse_[idx]).c_str(), std::string("sigmaietaieta").append(datasamplesreverse_[idx]).c_str() ,44, 0.001, 0.111));
hists_sigmaietaieta[idx]->Sumw2();
hists_HOE.push_back(new TH1F(  std::string("HOE").append(datasamplesreverse_[idx]).c_str(), std::string("HOE").append(datasamplesreverse_[idx]).c_str() ,8, 0, 0.05));
hists_HOE[idx]->Sumw2();
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
	std::vector<double> *mysigmaietaieta=0;
	std::vector<double> *myHOE=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mycoswphoton=0;



	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
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
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );
	theTree->SetBranchAddress( "HOE", &myHOE );
        theTree->SetBranchAddress( "coswphoton", &mycoswphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
//	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {

if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
	hists_photonpt[idx] ->Fill( (*myptphoton)[0],finalweight *scales[idx] );
	hists_photoneta[idx] ->Fill( (*myetaphoton )[0],finalweight *scales[idx]  );
	hists_muonpt[idx] ->Fill((*myptmuon )[0],finalweight *scales[idx]  );
	hists_muoneta[idx] ->Fill( (*myetamuon )[0] ,finalweight *scales[idx] );
	hists_jetpt[idx] ->Fill((*myptjet )[0],finalweight  *scales[idx]  );
	hists_jeteta[idx] ->Fill((*myetajet )[0],finalweight *scales[idx]  );
	hists_topmass[idx] ->Fill( (*mymasstop )[0] ,finalweight *scales[idx] );
        hists_toppt[idx] ->Fill( (*mypttop )[0],finalweight *scales[idx] );
        hists_topeta[idx] ->Fill( (*myetatop)[0],finalweight *scales[idx] );
        hists_topphotonmass[idx] ->Fill( (*mytopphotonmass)[0],finalweight *scales[idx] );
	hists_CSV[idx] ->Fill( (*mycvsdiscriminant)[0],finalweight *scales[idx]);
	hists_deltaphimetphoton[idx] ->Fill( (*mydeltaphiphotonmet)[0],finalweight *scales[idx]);
	hists_deltarphotonmuon[idx] ->Fill( (*mydeltaRphotonmuon)[0],finalweight *scales[idx]);
	hists_deltarphotonbjet[idx] ->Fill( (*mydeltaRphotonjet)[0],finalweight *scales[idx]);
	hists_jetmulti[idx] ->Fill( (*myjetmultiplicity)[0],finalweight *scales[idx]);
	hists_costopphoton[idx] ->Fill( (*mycoswphoton)[0],finalweight *scales[idx]);
	hists_sigmaietaieta[idx] ->Fill( (*mysigmaietaieta)[0],finalweight *scales[idx]);
	hists_HOE[idx] ->Fill( (*myHOE)[0],finalweight *scales[idx]);
}
}
}
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
	delete mysigmaietaieta;
	delete myHOE;
}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
	TFile *input(0);
	TString fname =datasamplesreverse_[idx];
	input = TFile::Open( fname ); 
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;
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
	std::vector<double> *mysigmaietaieta=0;
	std::vector<double> *myHOE=0;
	std::vector<double> *mycoswphoton=0;	

	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
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
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );
	theTree->SetBranchAddress( "HOE", &myHOE );
        theTree->SetBranchAddress( "coswphoton", &mycoswphoton );


	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
//	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
	hists_photonpt[samplesreverse_.size()] ->Fill( (*myptphoton)[0]);
	hists_photoneta[samplesreverse_.size()] ->Fill( (*myetaphoton )[0] );
	hists_muonpt[samplesreverse_.size()] ->Fill((*myptmuon )[0] );
	hists_muoneta[samplesreverse_.size()] ->Fill( (*myetamuon )[0] );
	hists_jetpt[samplesreverse_.size()] ->Fill((*myptjet )[0]  );
	hists_jeteta[samplesreverse_.size()] ->Fill((*myetajet )[0] );
	hists_topmass[samplesreverse_.size()] ->Fill( (*mymasstop )[0] );
        hists_toppt[samplesreverse_.size()] ->Fill( (*mypttop )[0]);
        hists_topeta[samplesreverse_.size()] ->Fill( (*myetatop)[0]);
        hists_topphotonmass[samplesreverse_.size()] ->Fill( (*mytopphotonmass)[0]);
	hists_CSV[samplesreverse_.size()] ->Fill( (*mycvsdiscriminant)[0]);
	hists_deltaphimetphoton[samplesreverse_.size()] ->Fill( (*mydeltaphiphotonmet)[0]);
	hists_deltarphotonmuon[samplesreverse_.size()] ->Fill( (*mydeltaRphotonmuon)[0]);
	hists_deltarphotonbjet[samplesreverse_.size()] ->Fill( (*mydeltaRphotonjet)[0]);
	hists_jetmulti[samplesreverse_.size()] ->Fill( (*myjetmultiplicity)[0]);
	hists_costopphoton[samplesreverse_.size()] ->Fill( (*mycoswphoton)[0]);
	hists_sigmaietaieta[samplesreverse_.size()] ->Fill( (*mysigmaietaieta)[0]);
	hists_HOE[samplesreverse_.size()] ->Fill( (*myHOE)[0]);
//}
}
}
}	}
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
	delete mysigmaietaieta;
	delete myHOE;
}


TH1F *wjet2ptphoton = new TH1F("pTphotonwjet","pTphotonwjet",12,0., 300.);
TH1F *wjet2sigmaietaieta = new TH1F("sigmaietaietawjet","sigmaietaietawjet",44, 0.001, 0.111);
TH1F *wjet2coswphoton = new TH1F("sigmaietaietawjet","sigmaietaietawjet",22, -1, 1);

	TFile *input(0);
	TString fname ="WJETREVERSE1.root";
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
	std::vector<double> *mysigmaietaieta=0;
	std::vector<double> *myHOE=0;
	std::vector<double> *myphiphoton=0;
	std::vector<double> *mycoswphoton=0;

	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );
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
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );
	theTree->SetBranchAddress( "HOE", &myHOE );
        theTree->SetBranchAddress( "coswphoton", &mycoswphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
//	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031)||(abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011)) {
//	wjet2->Fill( (*mysigmaietaieta)[0],finalweight  );
	wjet2ptphoton->Fill( (*myptphoton)[0],finalweight  );
	wjet2sigmaietaieta->Fill( (*mysigmaietaieta)[0],finalweight  );
	wjet2coswphoton->Fill( (*mycoswphoton)[0],finalweight  );
}
}
}
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
	delete mysigmaietaieta;
	delete myHOE;


TLegend* leg = new TLegend(0.5,0.5,0.75,0.85);
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
  leg->AddEntry( hists_photonpt[0], "W+jets"                           , "F");
  leg->AddEntry( hists_photonpt[samplesreverse_.size()], "CMS Data 2012(19.145/fb)"               , "PL");


TCanvas *c1 = new TCanvas("c1","multipads",900,700);
c1->cd(0);

hists_photonpt[0]->Add(wjet2ptphoton);
hists_photonpt[0]->Scale(0.436);
//hists_photonpt[0]->Draw();

THStack *hs1 = new THStack("hs1","photonpt");
hists_photonpt[0]->SetFillColor(kBlue-4);
hs1->Add(hists_photonpt[0]);
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
hs1->SetMaximum(1.2*hists_photonpt[samplesreverse_.size()]->GetMaximum());
hs1->GetXaxis()->SetTitle("Photon P_{T} (GeV)");
hs1->GetYaxis()->SetTitle("NEVENTS");

hists_photonpt[samplesreverse_.size()]->SetLineWidth(3.);
hists_photonpt[samplesreverse_.size()]->SetLineColor(kBlack);
hists_photonpt[samplesreverse_.size()]->SetMarkerColor(kBlack);
hists_photonpt[samplesreverse_.size()]->SetMarkerStyle(20.);
hists_photonpt[samplesreverse_.size()]->Draw("esame");

leg->Draw("same");


  TCanvas *c18 = new TCanvas("c18","sigmaietaieta",900,700);
  c18->cd(0);
hists_sigmaietaieta[0]->Add(wjet2sigmaietaieta);
hists_sigmaietaieta[0]->Scale(0.436);
//hists_sigmaietaieta[0]->Scale(0.68);

THStack *hs18 = new THStack("hs18","sigmaietaieta");
hists_sigmaietaieta[0]->SetFillColor(kBlue-4);
hs18->Add(hists_sigmaietaieta[0]);
hists_sigmaietaieta[1]->SetFillColor(kOrange+7);
hs18->Add(hists_sigmaietaieta[1]);
hists_sigmaietaieta[3]->Add(hists_sigmaietaieta[2]);
hists_sigmaietaieta[4]->Add(hists_sigmaietaieta[3]);
hists_sigmaietaieta[5]->Add(hists_sigmaietaieta[4]);
hists_sigmaietaieta[6]->Add(hists_sigmaietaieta[5]);
hists_sigmaietaieta[7]->Add(hists_sigmaietaieta[6]);
hists_sigmaietaieta[7]->SetFillColor(19);
hs18->Add(hists_sigmaietaieta[7]);
hists_sigmaietaieta[8]->SetFillColor(kGreen+1);
hs18->Add(hists_sigmaietaieta[8]);
hists_sigmaietaieta[10]->Add(hists_sigmaietaieta[9]);
hists_sigmaietaieta[6+5]->Add(hists_sigmaietaieta[5+5]);
hists_sigmaietaieta[7+5]->Add(hists_sigmaietaieta[6+5]);
hists_sigmaietaieta[8+5]->Add(hists_sigmaietaieta[7+5]);
hists_sigmaietaieta[9+5]->Add(hists_sigmaietaieta[8+5]);
hists_sigmaietaieta[9+5]->SetFillColor(kAzure+10);
hs18->Add(hists_sigmaietaieta[9+5]);
hists_sigmaietaieta[11+5]->Add(hists_sigmaietaieta[10+5]);
hists_sigmaietaieta[12+5]->Add(hists_sigmaietaieta[11+5]);
hists_sigmaietaieta[12+5]->SetFillColor(6);
hs18->Add(hists_sigmaietaieta[12+5]);
hists_sigmaietaieta[13+5]->SetFillColor(32);
hs18->Add(hists_sigmaietaieta[13+5]);
hists_sigmaietaieta[14+5]->SetFillColor(kSpring-9);
hs18->Add(hists_sigmaietaieta[14+5]);
hists_sigmaietaieta[16+5]->Add(hists_sigmaietaieta[15+5]);
hists_sigmaietaieta[17+5]->Add(hists_sigmaietaieta[16+5]);
hists_sigmaietaieta[17+5]->SetFillColor(kViolet-7);
hs18->Add(hists_sigmaietaieta[17+5]);
hists_sigmaietaieta[18+5]->SetFillColor(kOrange-2);
hs18->Add(hists_sigmaietaieta[18+5]);
hists_sigmaietaieta[20+5]->Add(hists_sigmaietaieta[19+5]);
hists_sigmaietaieta[20+5]->SetFillColor(kYellow+3);
hs18->Add(hists_sigmaietaieta[20+5]);

hs18->Draw("hist");
hs18->SetMaximum(1.2*hists_sigmaietaieta[samplesreverse_.size()]->GetMaximum());
hs18->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
hs18->GetYaxis()->SetTitle("NEVENTS");


hists_sigmaietaieta[samplesreverse_.size()]->SetLineWidth(3.);
hists_sigmaietaieta[samplesreverse_.size()]->SetLineColor(kBlack);
hists_sigmaietaieta[samplesreverse_.size()]->SetMarkerColor(kBlack);
hists_sigmaietaieta[samplesreverse_.size()]->SetMarkerStyle(20.);
hists_sigmaietaieta[samplesreverse_.size()]->Draw("esame");
leg->Draw("same");

  TCanvas *c17 = new TCanvas("c17","costopphoton",900,700);
  c17->cd(0);
hists_costopphoton[0]->Add(wjet2coswphoton);
hists_costopphoton[0]->Scale(0.436);
THStack *hs17 = new THStack("hs17","costopphoton");
hists_costopphoton[0]->SetFillColor(kBlue-4);
hs17->Add(hists_costopphoton[0]);
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
hs17->SetMaximum(1.2*hists_costopphoton[samplesreverse_.size()]->GetMaximum());
hs17->GetXaxis()->SetTitle("Cos(W,#gamma)");
hs17->GetYaxis()->SetTitle("NEVENTS");

hists_costopphoton[samplesreverse_.size()]->SetLineWidth(3.);
hists_costopphoton[samplesreverse_.size()]->SetLineColor(kBlack);
hists_costopphoton[samplesreverse_.size()]->SetMarkerColor(kBlack);
hists_costopphoton[samplesreverse_.size()]->SetMarkerStyle(20.);
hists_costopphoton[samplesreverse_.size()]->Draw("esame");
leg->Draw("same");

}






