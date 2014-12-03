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

void cutflowerrorSIGNAL() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> step_;
double NOE[2];
NOE[0]=312463;
NOE[1]=389959;

samples_.push_back("SIGNALtGu.root");
samples_.push_back("SIGNALtGC.root");

step_.push_back("STEP1");
step_.push_back("STEP2");
step_.push_back("STEP3");
step_.push_back("STEP4");
step_.push_back("STEP5");
step_.push_back("STEP6");
step_.push_back("STEP7");


std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;

	TTree* theTree = (TTree*)input->Get("STEP3/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
	hists_photoneta[idxx] ->Fill( 0.5);
}
}
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   efficiency               "<<hists_photoneta[idx]->Integral()/NOE[idx]<<endl;
}



std::vector<TH1F*> hists_photoneta4;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta4.push_back(new TH1F(  std::string("photon_eta4").append(samples_[idx]).c_str(),std::string("photon_eta4").append(samples_[idx]).c_str() , 1, -3., 3));

}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;

	TTree* theTree = (TTree*)input->Get("STEP4/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
	hists_photoneta4[idxx] ->Fill( 0.5);
}
}
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta4[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   efficiency               "<<hists_photoneta4[idx]->Integral()/NOE[idx]<<endl;
}


std::vector<TH1F*> hists_photoneta7;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta7.push_back(new TH1F(  std::string("photon_eta7").append(samples_[idx]).c_str(),std::string("photon_eta7").append(samples_[idx]).c_str() , 1, -3., 3));

}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;

	TTree* theTree = (TTree*)input->Get("STEP7/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
	hists_photoneta7[idxx] ->Fill( 0.5);
}
}
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP7!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta7[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   efficiency               "<<hists_photoneta7[idx]->Integral()/NOE[idx]<<endl;
}


std::vector<TH1F*> hists_photoneta8;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta8.push_back(new TH1F(  std::string("photon_eta8").append(samples_[idx]).c_str(),std::string("photon_eta8").append(samples_[idx]).c_str() , 1, -3., 3));

}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;
	std::vector<double> *mymasstop=0;
	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "masstop", &mymasstop );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {
//	hists_photoneta8[idxx] ->Fill( 0.5,finalweight*(0.049094905/19.145));}
	hists_photoneta8[idxx] ->Fill( 0.5);}
}
}
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!final!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta8[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   efficiency               "<<hists_photoneta8[idx]->Integral()/NOE[idx]<<endl;
}


}
