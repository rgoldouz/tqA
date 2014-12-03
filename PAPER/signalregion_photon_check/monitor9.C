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
void monitor9() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<string> datasamplesreverse_;
float scales[] = {0.628,0.0978,1.491,0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024,0.025};
samples_.push_back("WJET.root");
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
//be carefull____________ to make code shorter, histograms of MC samples are from zero to sample size and after that 3 histo for real data and after that 3 files for data in controlregion for wjet sample.
 
for(unsigned int idx=0; idx<samples_.size(); ++idx){

hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(samples_[idx]).c_str(), std::string("photon_pt").append(samples_[idx]).c_str() , 20,50., 600.));
hists_photonpt[idx]->Sumw2();
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 30, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(samples_[idx]).c_str(), std::string("muon_pt").append(samples_[idx]).c_str() , 25, 26., 400.));
hists_muoneta.push_back(new TH1F(  std::string("muon_eta").append(samples_[idx]).c_str(), std::string("muon_eta").append(samples_[idx]).c_str() ,  30, -3., 3.));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(samples_[idx]).c_str(),  std::string("jet_pt").append(samples_[idx]).c_str(),  25, 30., 500.));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(samples_[idx]).c_str(),  std::string("bjet_eta").append(samples_[idx]).c_str() , 30, -3., 3.));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(samples_[idx]).c_str(), std::string("top_mass").append(samples_[idx]).c_str() , 40, 0., 400.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(samples_[idx]).c_str(), std::string("toppt").append(samples_[idx]).c_str() ,50,0., 500.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(samples_[idx]).c_str(), std::string("topeta").append(samples_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(samples_[idx]).c_str(), std::string("topphotonmass").append(samples_[idx]).c_str() ,40, 0., 1000.));
}

for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(datasamples_[idx]).c_str(), std::string("photon_pt").append(datasamples_[idx]).c_str() , 20,50., 600.));
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(datasamples_[idx]).c_str(),  std::string("photon_eta").append(datasamples_[idx]).c_str() , 30, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamples_[idx]).c_str(), std::string("muon_pt").append(datasamples_[idx]).c_str() , 25, 26., 400.));
hists_muoneta.push_back(new TH1F(  std::string("muon_eta").append(datasamples_[idx]).c_str(), std::string("muon_eta").append(datasamples_[idx]).c_str() ,  30, -3., 3.));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamples_[idx]).c_str(),  std::string("jet_pt").append(datasamples_[idx]).c_str(),  25, 30., 500.));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamples_[idx]).c_str(),  std::string("bjet_eta").append(datasamples_[idx]).c_str() , 30, -3., 3.));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamples_[idx]).c_str(), std::string("top_mass").append(datasamples_[idx]).c_str() , 40, 0., 400.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(datasamples_[idx]).c_str(), std::string("toppt").append(datasamples_[idx]).c_str() ,50,0., 500.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(datasamples_[idx]).c_str(), std::string("topeta").append(datasamples_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(datasamples_[idx]).c_str(), std::string("topphotonmass").append(datasamples_[idx]).c_str() ,40, 0., 1000.));

}

for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
hists_photonpt.push_back(new  TH1F(  std::string("photon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("photon_pt").append(datasamplesreverse_[idx]).c_str() , 20,50., 600.));
hists_photoneta.push_back(new  TH1F(  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("photon_eta").append(datasamplesreverse_[idx]).c_str() , 30, -3., 3));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamplesreverse_[idx]).c_str(), std::string("muon_pt").append(datasamplesreverse_[idx]).c_str() , 25, 26., 400.));
hists_muoneta.push_back(new  TH1F(  std::string("muon_eta").append(datasamplesreverse_[idx]).c_str(), std::string("muon_eta").append(datasamplesreverse_[idx]).c_str() ,  30, -3., 3.));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  std::string("jet_pt").append(datasamplesreverse_[idx]).c_str(),  25, 30., 500.));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str(),  std::string("bjet_eta").append(datasamplesreverse_[idx]).c_str() , 30, -3., 3.));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamplesreverse_[idx]).c_str(), std::string("top_mass").append(datasamplesreverse_[idx]).c_str() , 40, 0., 400.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(datasamplesreverse_[idx]).c_str(), std::string("toppt").append(datasamplesreverse_[idx]).c_str() ,50,0., 500.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(datasamplesreverse_[idx]).c_str(), std::string("topeta").append(datasamplesreverse_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(datasamplesreverse_[idx]).c_str(), std::string("topphotonmass").append(datasamplesreverse_[idx]).c_str() ,40, 0., 1000.));

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
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
        if (finalweight<0) cout<<"negative weight=  " <<finalweight<<"    "<<(*myptphoton)[0]<<samples_[idx]<<endl;

//        if (samples_[idx]=="WPHJET.root") finalweight=finalweight*(1.46604+0.00923692*(*myptphoton)[0]-1.48871e-06*(*myptphoton)[0]*(*myptphoton)[0]);
//        if (samples_[idx]=="ZGAMMA.root") finalweight=finalweight*(1.3292+0.000952237*(*myptphoton)[0]+2.00623e-05*(*myptphoton)[0]*(*myptphoton)[0]-1.41325e-07 *(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]+2.48614e-10*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]);

	if (finalweight<0) cout<<"negative weight=  " <<finalweight<<"    "<<(*myptphoton)[0]<<samples_[idx]<<endl;
//	if (samples_[idx]=="WPHJET.root" && finalweight>100) 
//{cout<<finalweight*scales[idx]<<"   "<<(*myptphoton)[0]<<endl;
//cout<< "btagSF"<<"   "<<(*mybtagSF)[0]<<endl;
//cout<<"photonSF"<<"   "<<(*myphotonSF)[0]<<endl;
//cout<<"muonSF"<<"   "<<(*mymuonSF)[0]<<endl;
//cout<<"pileupSF"<<"   "<<(*mypileupSF)[0]<<endl;
//}
	if ((*mymasstop )[0]>130 && (*mymasstop )[0]<220 ) {
	hists_photonpt[idx] ->Fill( (*myptphoton)[0],finalweight *scales[idx] );
	hists_photoneta[idx] ->Fill( (*myetaphoton )[0],finalweight*scales[idx] );
	hists_muonpt[idx] ->Fill((*myptmuon )[0] ,finalweight*scales[idx] );
	hists_muoneta[idx] ->Fill( (*myetamuon )[0],finalweight *scales[idx]);
	hists_jetpt[idx] ->Fill((*myptjet )[0] ,finalweight *scales[idx]);
	hists_jeteta[idx] ->Fill((*myetajet )[0] ,finalweight*scales[idx] );
	hists_topmass[idx] ->Fill( (*mymasstop )[0],finalweight *scales[idx]);
	hists_toppt[idx] ->Fill( (*mypttop )[0],finalweight *scales[idx]);
	hists_topeta[idx] ->Fill( (*myetatop)[0],finalweight *scales[idx]);
	hists_topphotonmass[idx] ->Fill( (*mytopphotonmass)[0],finalweight *scales[idx]);
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
	delete input;


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
	std::vector<double> *myweight=0;

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
	theTree->SetBranchAddress( "weight", &myweight);
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
	delete input;
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
	delete input;
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

}
cout<<"WJET       "<<hists_photoneta[0]->Integral()<<endl;
cout<<"W  ph  JET       "<<hists_photoneta[3]->Integral()<<endl;

hists_photonpt[samples_.size()+5]->Scale((212.609)/hists_photonpt[samples_.size()+5]->Integral()); 
hists_photoneta[samples_.size()+5]->Scale((207.609)/hists_photoneta[samples_.size()+5]->Integral()); 
hists_muonpt[samples_.size()+5]->Scale((207.609)/hists_muonpt[samples_.size()+5]->Integral());
hists_muoneta[samples_.size()+5]->Scale((207.609)/hists_muoneta[samples_.size()+5]->Integral()); 
hists_jetpt[samples_.size()+5]->Scale((207.609)/hists_jetpt[samples_.size()+5]->Integral());
hists_jeteta[samples_.size()+5]->Scale((207.609)/hists_jeteta[samples_.size()+5]->Integral()); 
hists_topmass[samples_.size()+5]->Scale((207.609)/hists_topmass[samples_.size()+5]->Integral()); 
hists_toppt[samples_.size()+5]->Scale((207.609)/hists_toppt[samples_.size()+5]->Integral());
hists_topeta[samples_.size()+5]->Scale((207.609)/hists_topeta[samples_.size()+5]->Integral());
hists_topphotonmass[samples_.size()+5]->Scale((207.609)/hists_topphotonmass[samples_.size()+5]->Integral());

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

hists_photonpt[3]->Scale(1101.23/hists_photonpt[3]->Integral()); 
hists_photoneta[3]->Scale(1065.23/hists_photoneta[3]->Integral()); 
hists_muonpt[3]->Scale(1065.23/hists_muonpt[3]->Integral());
hists_muoneta[3]->Scale(1065.23/hists_muoneta[3]->Integral()); 
hists_jetpt[3]->Scale(1065.23/hists_jetpt[3]->Integral());
hists_jeteta[3]->Scale(1065.23/hists_jeteta[3]->Integral()); 
hists_topmass[3]->Scale(1065.23/hists_topmass[3]->Integral()); 
hists_toppt[3]->Scale(1065.23/hists_toppt[3]->Integral());
hists_topeta[3]->Scale(1065.23/hists_topeta[3]->Integral());
hists_topphotonmass[3]->Scale(1065.23/hists_topphotonmass[3]->Integral());


for(unsigned int idx=2; idx<samples_.size()-1; ++idx){
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

}

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

}
TH1F *sum_h= new TH1F ( *hists_photonpt[20] ) ;

cout<<" real data = "<< hists_photonpt[samples_.size()+2]->Integral()<<endl;
cout<<" MC  = "<<hists_photonpt[20]->Integral()<<endl;
cout<<" MC  = "<<hists_topphotonmass[20]->Integral()<<endl;
cout<<" real data = "<< hists_topphotonmass[samples_.size()+2]->Integral()<<endl;
  TCanvas *c1 = new TCanvas("c1","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c1->cd(0);
hists_photonpt[20]->SetMaximum(1.2*hists_photonpt[samples_.size()+2]->GetMaximum());
hists_photonpt[20]->SetFillColor(kYellow+3);
hists_photonpt[20]->Draw();
hists_photonpt[18]->SetFillColor(kOrange-2);
hists_photonpt[18]->Draw("histsame");
hists_photonpt[17]->SetFillColor(kViolet-7);
hists_photonpt[17]->Draw("histsame");
//hists_photonpt[16]->SetFillColor(kRed);
//hists_photonpt[16]->Draw("same");
//hists_photonpt[15]->SetFillColor(kViolet+1);
//hists_photonpt[15]->Draw("same");
hists_photonpt[14]->SetFillColor(kSpring-9);
hists_photonpt[14]->Draw("histsame");
hists_photonpt[13]->SetFillColor(32);
hists_photonpt[13]->Draw("histsame");
hists_photonpt[12]->SetFillColor(6);
hists_photonpt[12]->Draw("histsame");
hists_photonpt[9]->SetFillColor(kAzure+10);
hists_photonpt[9]->Draw("histsame");
hists_photonpt[3]->SetFillColor(kGreen+1);
hists_photonpt[3]->Draw("histsame");
hists_photonpt[2]->SetFillColor(19);
hists_photonpt[2]->Draw("histsame");
hists_photonpt[1]->SetFillColor(kOrange+7);
hists_photonpt[1]->Draw("histsame");
hists_photonpt[samples_.size()+5]->SetFillColor(kBlue-4);
hists_photonpt[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_photonpt[21]->SetLineColor(kRed+3);
hists_photonpt[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_photonpt[21]->Draw("same");

hists_photonpt[samples_.size()+2]->SetLineWidth(3.);
hists_photonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_photonpt[samples_.size()+2]->Draw("esame");

sum_h->SetFillColor(kOrange-7);
sum_h->SetFillStyle(3022);
sum_h->Draw("e2same");

TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( hists_photonpt[samples_.size()+5], "W JET"                           , "F");
  leg->AddEntry( hists_photonpt[1], "Z JET"                           , "F");
  leg->AddEntry( hists_photonpt[2], "PH JET"                           , "F");
  leg->AddEntry( hists_photonpt[3], "W PH JET"              , "F");
  leg->AddEntry( hists_photonpt[9], "SINGLE TOP  "               , "F");
  leg->AddEntry( hists_photonpt[12], "TTBAR"               , "F");
  leg->AddEntry( hists_photonpt[13], "TTG"               , "F");
  leg->AddEntry( hists_photonpt[14], "WWG"               , "F");
//  leg->AddEntry( hists_photonpt[15], "WW"               , "F");
//  leg->AddEntry( hists_photonpt[16], "WZ"               , "F");
  leg->AddEntry( hists_photonpt[17], "DIBOSON"               , "F");
  leg->AddEntry( hists_photonpt[18], "ZGAMMA"               , "F");
  leg->AddEntry( hists_photonpt[20], "SINGLE TOP+PHOTON"               , "F");
  leg->AddEntry( hists_photonpt[21], "SIGNAL"               , "F");
  leg->AddEntry( hists_photonpt[samples_.size()+2], "CMS Data 2012(19.145/fb)"               , "PL");

  leg->Draw("same");

  TCanvas *c2 = new TCanvas("c2","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c2->cd(0);
hists_photoneta[20]->SetMaximum(1.2*hists_photoneta[samples_.size()+2]->GetMaximum());
hists_photoneta[20]->SetFillColor(kYellow+3);
hists_photoneta[20]->Draw();
hists_photoneta[18]->SetFillColor(kOrange-2);
hists_photoneta[18]->Draw("same");
hists_photoneta[17]->SetFillColor(kViolet-7);
hists_photoneta[17]->Draw("same");
//hists_photoneta[16]->SetFillColor(kRed);
//hists_photoneta[16]->Draw("same");
//hists_photoneta[15]->SetFillColor(kViolet+1);
//hists_photoneta[15]->Draw("same");
hists_photoneta[14]->SetFillColor(kSpring-9);
hists_photoneta[14]->Draw("same");
hists_photoneta[13]->SetFillColor(32);
hists_photoneta[13]->Draw("same");
hists_photoneta[12]->SetFillColor(6);
hists_photoneta[12]->Draw("same");
hists_photoneta[9]->SetFillColor(kAzure+10);
hists_photoneta[9]->Draw("same");
hists_photoneta[3]->SetFillColor(kGreen+1);
hists_photoneta[3]->Draw("same");
hists_photoneta[2]->SetFillColor(19);
hists_photoneta[2]->Draw("same");
hists_photoneta[1]->SetFillColor(kOrange+7);
hists_photoneta[1]->Draw("same");
hists_photoneta[samples_.size()+5]->SetFillColor(kBlue-4);
hists_photoneta[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_photoneta[21]->SetLineColor(kRed+3);
hists_photoneta[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_photoneta[21]->Draw("same");

hists_photoneta[samples_.size()+2]->SetLineWidth(3.);
hists_photoneta[samples_.size()+2]->SetLineColor(kBlack);
hists_photoneta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photoneta[samples_.size()+2]->SetMarkerStyle(20.);
hists_photoneta[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c3 = new TCanvas("c3","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c3->cd(0);
hists_muonpt[20]->SetMaximum(1.2*hists_muonpt[samples_.size()+2]->GetMaximum());
hists_muonpt[20]->SetFillColor(kYellow+3);
hists_muonpt[20]->Draw();
hists_muonpt[18]->SetFillColor(kOrange-2);
hists_muonpt[18]->Draw("same");
hists_muonpt[17]->SetFillColor(kViolet-7);
hists_muonpt[17]->Draw("same");
//hists_muonpt[16]->SetFillColor(kRed);
//hists_muonpt[16]->Draw("same");
//hists_muonpt[15]->SetFillColor(kViolet+1);
//hists_muonpt[15]->Draw("same");
hists_muonpt[14]->SetFillColor(kSpring-9);
hists_muonpt[14]->Draw("same");
hists_muonpt[13]->SetFillColor(32);
hists_muonpt[13]->Draw("same");
hists_muonpt[12]->SetFillColor(6);
hists_muonpt[12]->Draw("same");
hists_muonpt[9]->SetFillColor(kAzure+10);
hists_muonpt[9]->Draw("same");
hists_muonpt[3]->SetFillColor(kGreen+1);
hists_muonpt[3]->Draw("same");
hists_muonpt[2]->SetFillColor(19);
hists_muonpt[2]->Draw("same");
hists_muonpt[1]->SetFillColor(kOrange+7);
hists_muonpt[1]->Draw("same");
hists_muonpt[samples_.size()+5]->SetFillColor(kBlue-4);
hists_muonpt[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_muonpt[21]->SetLineColor(kRed+3);
hists_muonpt[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_muonpt[21]->Draw("same");

hists_muonpt[samples_.size()+2]->SetLineWidth(3.);
hists_muonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_muonpt[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c4 = new TCanvas("c4","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c4->cd(0);
hists_muoneta[20]->SetMaximum(1.2*hists_muoneta[samples_.size()+2]->GetMaximum());
hists_muoneta[20]->SetFillColor(kYellow+3);
hists_muoneta[20]->Draw();
hists_muoneta[18]->SetFillColor(kOrange-2);
hists_muoneta[18]->Draw("same");
hists_muoneta[17]->SetFillColor(kViolet-7);
hists_muoneta[17]->Draw("same");
//hists_muoneta[16]->SetFillColor(kRed);
//hists_muoneta[16]->Draw("same");
//hists_muoneta[15]->SetFillColor(kViolet+1);
//hists_muoneta[15]->Draw("same");
hists_muoneta[14]->SetFillColor(kSpring-9);
hists_muoneta[14]->Draw("same");
hists_muoneta[13]->SetFillColor(32);
hists_muoneta[13]->Draw("same");
hists_muoneta[12]->SetFillColor(6);
hists_muoneta[12]->Draw("same");
hists_muoneta[9]->SetFillColor(kAzure+10);
hists_muoneta[9]->Draw("same");
hists_muoneta[3]->SetFillColor(kGreen+1);
hists_muoneta[3]->Draw("same");
hists_muoneta[2]->SetFillColor(19);
hists_muoneta[2]->Draw("same");
hists_muoneta[1]->SetFillColor(kOrange+7);
hists_muoneta[1]->Draw("same");
hists_muoneta[samples_.size()+5]->SetFillColor(kBlue-4);
hists_muoneta[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_muoneta[21]->SetLineColor(kRed+3);
hists_muoneta[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_muoneta[21]->Draw("same");

hists_muoneta[samples_.size()+2]->SetLineWidth(3.);
hists_muoneta[samples_.size()+2]->SetLineColor(kBlack);
hists_muoneta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muoneta[samples_.size()+2]->SetMarkerStyle(20.);
hists_muoneta[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c5 = new TCanvas("c5","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c5->cd(0);
hists_jetpt[20]->SetMaximum(1.2*hists_jetpt[samples_.size()+2]->GetMaximum());
hists_jetpt[20]->SetFillColor(kYellow+3);
hists_jetpt[20]->Draw();
hists_jetpt[18]->SetFillColor(kOrange-2);
hists_jetpt[18]->Draw("same");
hists_jetpt[17]->SetFillColor(kViolet-7);
hists_jetpt[17]->Draw("same");
//hists_jetpt[16]->SetFillColor(kRed);
//hists_jetpt[16]->Draw("same");
//hists_jetpt[15]->SetFillColor(kViolet+1);
//hists_jetpt[15]->Draw("same");
hists_jetpt[14]->SetFillColor(kSpring-9);
hists_jetpt[14]->Draw("same");
hists_jetpt[13]->SetFillColor(32);
hists_jetpt[13]->Draw("same");
hists_jetpt[12]->SetFillColor(6);
hists_jetpt[12]->Draw("same");
hists_jetpt[9]->SetFillColor(kAzure+10);
hists_jetpt[9]->Draw("same");
hists_jetpt[3]->SetFillColor(kGreen+1);
hists_jetpt[3]->Draw("same");
hists_jetpt[2]->SetFillColor(19);
hists_jetpt[2]->Draw("same");
hists_jetpt[1]->SetFillColor(kOrange+7);
hists_jetpt[1]->Draw("same");
hists_jetpt[samples_.size()+5]->SetFillColor(kBlue-4);
hists_jetpt[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_jetpt[21]->SetLineColor(kRed+3);
hists_jetpt[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_jetpt[21]->Draw("same");

hists_jetpt[samples_.size()+2]->SetLineWidth(3.);
hists_jetpt[samples_.size()+2]->SetLineColor(kBlack);
hists_jetpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_jetpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_jetpt[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c6 = new TCanvas("c6","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c6->cd(0);
hists_jeteta[20]->SetMaximum(1.2*hists_jeteta[samples_.size()+2]->GetMaximum());
hists_jeteta[20]->SetFillColor(kYellow+3);
hists_jeteta[20]->Draw();
hists_jeteta[18]->SetFillColor(kOrange-2);
hists_jeteta[18]->Draw("same");
hists_jeteta[17]->SetFillColor(kViolet-7);
hists_jeteta[17]->Draw("same");
//hists_jeteta[16]->SetFillColor(kRed);
//hists_jeteta[16]->Draw("same");
//hists_jeteta[15]->SetFillColor(kViolet+1);
//hists_jeteta[15]->Draw("same");
hists_jeteta[14]->SetFillColor(kSpring-9);
hists_jeteta[14]->Draw("same");
hists_jeteta[13]->SetFillColor(32);
hists_jeteta[13]->Draw("same");
hists_jeteta[12]->SetFillColor(6);
hists_jeteta[12]->Draw("same");
hists_jeteta[9]->SetFillColor(kAzure+10);
hists_jeteta[9]->Draw("same");
hists_jeteta[3]->SetFillColor(kGreen+1);
hists_jeteta[3]->Draw("same");
hists_jeteta[2]->SetFillColor(19);
hists_jeteta[2]->Draw("same");
hists_jeteta[1]->SetFillColor(kOrange+7);
hists_jeteta[1]->Draw("same");
hists_jeteta[samples_.size()+5]->SetFillColor(kBlue-4);
hists_jeteta[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_jeteta[21]->SetLineColor(kRed+3);
hists_jeteta[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
hists_jeteta[21]->Draw("same");

hists_jeteta[samples_.size()+2]->SetLineWidth(3.);
hists_jeteta[samples_.size()+2]->SetLineColor(kBlack);
hists_jeteta[samples_.size()+2]->SetMarkerColor(kBlack);
hists_jeteta[samples_.size()+2]->SetMarkerStyle(20.);
hists_jeteta[samples_.size()+2]->Draw("esame");
  leg->Draw("same");

  TCanvas *c7 = new TCanvas("c7","multipads",900,700);
//  c1->Divide(2,2,0,0);
  c7->cd(0);
hists_topmass[20]->SetMaximum(1.5*hists_topmass[samples_.size()+2]->GetMaximum());
hists_topmass[20]->SetFillColor(kYellow+3);
hists_topmass[20]->Draw();
hists_topmass[18]->SetFillColor(kOrange-2);
hists_topmass[18]->Draw("same");
hists_topmass[17]->SetFillColor(kViolet-7);
hists_topmass[17]->Draw("same");
//hists_topmass[16]->SetFillColor(kRed);
//hists_topmass[16]->Draw("same");
//hists_topmass[15]->SetFillColor(kViolet+1);
//hists_topmass[15]->Draw("same");
hists_topmass[14]->SetFillColor(kSpring-9);
hists_topmass[14]->Draw("same");
hists_topmass[13]->SetFillColor(32);
hists_topmass[13]->Draw("same");
hists_topmass[12]->SetFillColor(6);
hists_topmass[12]->Draw("same");
hists_topmass[9]->SetFillColor(kAzure+10);
hists_topmass[9]->Draw("same");
hists_topmass[3]->SetFillColor(kGreen+1);
hists_topmass[3]->Draw("same");
hists_topmass[2]->SetFillColor(19);
hists_topmass[2]->Draw("same");
hists_topmass[1]->SetFillColor(kOrange+7);
hists_topmass[1]->Draw("same");
hists_topmass[samples_.size()+5]->SetFillColor(kBlue-4);
hists_topmass[samples_.size()+5]->Draw("same");
//hists_photonpt[21]->SetFillColor(kRed+3);
hists_topmass[21]->SetLineColor(kRed+3);
hists_topmass[21]->SetLineWidth(3);
//hists_photonpt[21]->SetFillStyle(3022);
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

}






