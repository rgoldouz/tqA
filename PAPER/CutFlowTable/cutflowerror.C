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
#include <cmath> 

void cutflowerror() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<string> step_;

double newweight=19.768/19.145;


float scales[] = {0.0961,0.0978,1.491,0.628,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139,34.01/19.145,6.133/19.145,1.04/19.145,0.32/19.145,0.02/19.145,0.002/19.145,19.145*0.0169};
samples_.push_back("WPHJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("PHJET200400.root");
samples_.push_back("WJET.root");
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
samples_.push_back("G_Pt_50to80.root");
samples_.push_back("G_Pt_80to120.root");
samples_.push_back("G_Pt_120to170.root");
samples_.push_back("G_Pt_170to300.root");
samples_.push_back("G_Pt_300to470.root");
samples_.push_back("G_Pt_470to800.root");
samples_.push_back("SIGNALtGu.root");

//datasamples_.push_back("REALDATA1.root");
//datasamples_.push_back("REALDATA2.root");
//datasamples_.push_back("REALDATA3.root");

datasamples_.push_back("REALDATAab.root");
datasamples_.push_back("REALDATAc.root");
datasamples_.push_back("REALDATAd.root");


step_.push_back("STEP1");
step_.push_back("STEP2");
step_.push_back("STEP3");
step_.push_back("STEP4");
step_.push_back("STEP5");
step_.push_back("STEP6");
step_.push_back("STEP7");


//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}

std::vector<TFile*> datafiles;
for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
datafiles.push_back(new TFile(datasamples_[idx].c_str()));
}

/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("/Photon_Pt");


for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP1").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP1").append(variables2_[i]).c_str()));
}
float lumi = 1;
//float lumi = 4.429;

std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
hists_photoneta[idx]->Sumw2();
}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;

	TTree* theTree = (TTree*)input->Get("STEP1/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
	hists_photoneta[idxx] ->Fill( 0.5,finalweight*scales[idxx] );
}
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   error               "<<hists_photoneta[idx]->GetBinError(1)<<endl;
}
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}



for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP2").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP2").append(variables2_[i]).c_str()));
}
float lumi = 1;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}


cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int i=0; i<variables2_.size(); ++i){
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP3").append(variables2_[i]).c_str()));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP3").append(variables2_[i]).c_str()));
}
float lumi = 1;
std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
hists_photoneta[idx]->Sumw2();
}

TH1F * histstotal;
histstotal= new TH1F(  "histstotal","histstotal", 1, -3., 3);

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myphiphoton=0;

	TTree* theTree = (TTree*)input->Get("STEP3/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);

	double finalweight;
	finalweight=(*myweight)[0];
 if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
	hists_photoneta[idxx] ->Fill( 0.5,finalweight*scales[idxx]*newweight );
}}
}
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
histstotal->Add(hists_photoneta[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   error               "<<hists_photoneta[idx]->GetBinError(1)<<endl;
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}


for(unsigned int idx=5; idx<10; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    total single top      "<<hists_photoneta[9]->Integral()<<endl;
cout<<"   error   total single top     "<<hists_photoneta[9]->GetBinError(1)<<endl;
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;

for(unsigned int idx=11; idx<13; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   total ttbar     "<<hists_photoneta[12]->Integral()<<endl;
cout<<"   error  total ttbar   "<<hists_photoneta[12]->GetBinError(1)<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

for(unsigned int idx=16; idx<18; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   diboson     "<<hists_photoneta[17]->Integral()<<endl;
cout<<"   error diboson  "<<hists_photoneta[17]->GetBinError(1)<<endl;
cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;

hists_photoneta[20]->Add(hists_photoneta[19]);
cout<<"  integral   singletop+photon     "<<hists_photoneta[20]->Integral()<<endl;
cout<<"   error singletop+photon   "<<hists_photoneta[20]->GetBinError(1)<<endl;
cout<<"singletop+photon "<<"              "<<hists[19]->Integral()+hists[20]->Integral()<<endl;

for(unsigned int idx=22; idx<27; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    gamma jet      "<<hists_photoneta[26]->Integral()<<endl;
cout<<"   error   gamma jet    "<<hists_photoneta[26]->GetBinError(1)<<endl;
cout<<"totalgamma jet"<<"              "<<hists[21]->Integral()+hists[22]->Integral()+hists[23]->Integral()+hists[24]->Integral()+hists[25]->Integral()+hists[26]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"integral totalMC     "<<histstotal->Integral()<<endl;
cout<<" error  totalMC   "<<histstotal->GetBinError(1)<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()  <<"error" << sqrt(datahists[2]->Integral())<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}

cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int i=0; i<variables2_.size(); ++i){
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP4").append(variables2_[i]).c_str()));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP4").append(variables2_[i]).c_str()));
}
float lumi = 1;
std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
hists_photoneta[idx]->Sumw2();
}

TH1F * histstotal;
histstotal= new TH1F(  "histstotal","histstotal", 1, -3., 3);

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myphiphoton=0;


	TTree* theTree = (TTree*)input->Get("STEP4/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);
       theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
 if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
	hists_photoneta[idxx] ->Fill( 0.5,finalweight*scales[idxx]*newweight );
}}
}
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
histstotal->Add(hists_photoneta[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   error               "<<hists_photoneta[idx]->GetBinError(1)<<endl;
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}


for(unsigned int idx=5; idx<10; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    total single top      "<<hists_photoneta[9]->Integral()<<endl;
cout<<"   error   total single top     "<<hists_photoneta[9]->GetBinError(1)<<endl;
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;

for(unsigned int idx=11; idx<13; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   total ttbar     "<<hists_photoneta[12]->Integral()<<endl;
cout<<"   error  total ttbar   "<<hists_photoneta[12]->GetBinError(1)<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

for(unsigned int idx=16; idx<18; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   diboson     "<<hists_photoneta[17]->Integral()<<endl;
cout<<"   error diboson  "<<hists_photoneta[17]->GetBinError(1)<<endl;
cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;

hists_photoneta[20]->Add(hists_photoneta[19]);
cout<<"  integral   singletop+photon     "<<hists_photoneta[20]->Integral()<<endl;
cout<<"   error singletop+photon   "<<hists_photoneta[20]->GetBinError(1)<<endl;
cout<<"singletop+photon "<<"              "<<hists[19]->Integral()+hists[20]->Integral()<<endl;

for(unsigned int idx=22; idx<27; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    gamma jet      "<<hists_photoneta[26]->Integral()<<endl;
cout<<"   error   gamma jet    "<<hists_photoneta[26]->GetBinError(1)<<endl;
cout<<"totalgamma jet"<<"              "<<hists[21]->Integral()+hists[22]->Integral()+hists[23]->Integral()+hists[24]->Integral()+hists[25]->Integral()+hists[26]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"integral totalMC     "<<histstotal->Integral()<<endl;
cout<<" error  totalMC   "<<histstotal->GetBinError(1)<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()  <<"error" << sqrt(datahists[2]->Integral())<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}


for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP5").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP5").append(variables2_[i]).c_str()));
}
float lumi = 1;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}


for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP6").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP6").append(variables2_[i]).c_str()));
}
float lumi = 1;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP6!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}



cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP7!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP7").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP7").append(variables2_[i]).c_str()));
}
float lumi = 1;
std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
hists_photoneta[idx]->Sumw2();
}

TH1F * histstotal;
histstotal= new TH1F(  "histstotal","histstotal", 1, -3., 3);

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;
std::vector<double> *myetaphoton=0;
std::vector<double> *myphiphoton=0;


	TTree* theTree = (TTree*)input->Get("STEP7/FCNC");
	theTree->SetBranchAddress( "weight", &myweight);
       theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
 if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
	hists_photoneta[idxx] ->Fill( 0.5,finalweight*scales[idxx]*newweight );
}}}
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
histstotal->Add(hists_photoneta[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   error               "<<hists_photoneta[idx]->GetBinError(1)<<endl;
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}


for(unsigned int idx=5; idx<10; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    total single top      "<<hists_photoneta[9]->Integral()<<endl;
cout<<"   error   total single top     "<<hists_photoneta[9]->GetBinError(1)<<endl;
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;

for(unsigned int idx=11; idx<13; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   total ttbar     "<<hists_photoneta[12]->Integral()<<endl;
cout<<"   error  total ttbar   "<<hists_photoneta[12]->GetBinError(1)<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

for(unsigned int idx=16; idx<18; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   diboson     "<<hists_photoneta[17]->Integral()<<endl;
cout<<"   error diboson  "<<hists_photoneta[17]->GetBinError(1)<<endl;
cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;

hists_photoneta[20]->Add(hists_photoneta[19]);
cout<<"  integral   singletop+photon     "<<hists_photoneta[20]->Integral()<<endl;
cout<<"   error singletop+photon   "<<hists_photoneta[20]->GetBinError(1)<<endl;
cout<<"singletop+photon "<<"              "<<hists[19]->Integral()+hists[20]->Integral()<<endl;

for(unsigned int idx=22; idx<27; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    gamma jet      "<<hists_photoneta[26]->Integral()<<endl;
cout<<"   error   gamma jet    "<<hists_photoneta[26]->GetBinError(1)<<endl;
cout<<"totalgamma jet"<<"              "<<hists[21]->Integral()+hists[22]->Integral()+hists[23]->Integral()+hists[24]->Integral()+hists[25]->Integral()+hists[26]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"integral totalMC     "<<histstotal->Integral()<<endl;
cout<<" error  totalMC   "<<histstotal->GetBinError(1)<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()  <<"error" << sqrt(datahists[2]->Integral())<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}
/*
for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((step_[mdx].append(variables2_[i])).c_str()));
hists.push_back((TH1F*)files[idx]->Get(std::string("STEP8").append(variables2_[i]).c_str()));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}

std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get(std::string("STEP8").append(variables2_[i]).c_str()));
}
float lumi = 1;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP8!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}
cout<<"total single top"<<"              "<<hists[4]->Integral()+hists[5]->Integral()+hists[6]->Integral()+hists[7]->Integral()+hists[8]->Integral()+hists[9]->Integral()<<endl;
cout<<"total ttbar"<<"              "<<hists[10]->Integral()+hists[11]->Integral()+hists[12]->Integral()<<endl;

cout<<"diboson"<<"              "<<hists[15]->Integral()+hists[16]->Integral()+hists[17]->Integral()<<endl;


for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
cout<<"****************************************************************"<<endl;
cout<<"totalMC"<<"              "<<hists[samples_.size()-2]->Integral()<<endl;

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }
cout<<"DATA"<<"              "<<datahists[2]->Integral()<<endl;
cout<<"--------------------------------------------------------------------------"<<endl;

}*/
TH1F * histstotal;
histstotal= new TH1F(  "histstotal","histstotal", 1, -3., 3);
std::vector<TH1F*> hists_photoneta;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(samples_[idx]).c_str(),std::string("photon_eta").append(samples_[idx]).c_str() , 1, -3., 3));
hists_photoneta[idx]->Sumw2();
}

for(unsigned int idxx=0; idxx<samples_.size(); ++idxx){
	TFile *input(0);
	TString fname =samples_[idxx];
	input = TFile::Open( fname ); 
	std::vector<double> *myweight=0;
	std::vector<double> *mymasstop=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myphiphoton=0;

	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "masstop", &mymasstop );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
	theTree->SetBranchAddress( "phiphoton", &myphiphoton );

	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];
if((*mymasstop )[0]>130 && (*mymasstop )[0]<220) {
if(sqrt(pow((*myetaphoton )[0]+1.76,2)+pow((*myphiphoton )[0]-1.37,2))>0.05 && sqrt(pow((*myetaphoton )[0]-2.37,2)+pow((*myphiphoton )[0]-2.69,2))>0.05){
if(sqrt(pow((*myetaphoton )[0]- 1.61,2)+pow((*myphiphoton )[0]+2.05,2))>0.05 && sqrt(pow((*myetaphoton )[0]-1.75,2)+pow((*myphiphoton )[0]-2.12,2))>0.05){
	hists_photoneta[idxx] ->Fill( 0.5,finalweight*scales[idxx]*newweight );}}}
}
}


////////////////////////////////////////////////////////
///if you want to have MC only turn this part off
hists_photoneta[3]->Scale(newweight*302.91/hists_photoneta[3]->Integral());
hists_photoneta[0]->Scale(newweight*985.18/hists_photoneta[0]->Integral());

hists_photoneta[3]->SetBinError(1,90.86);
hists_photoneta[0]->SetBinError(1,159.16);

for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"  integral                "<<hists_photoneta[idx]->Integral()<<endl;
cout<<samples_[idx]<<"   error               "<<hists_photoneta[idx]->GetBinError(1)<<endl;
}
//////////////////////////////////////////////////////////
for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
histstotal->Add(hists_photoneta[idx]);
}
for(unsigned int idx=5; idx<10; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    total single top      "<<hists_photoneta[9]->Integral()<<endl;
cout<<"   error   total single top     "<<hists_photoneta[9]->GetBinError(1)<<endl;


for(unsigned int idx=11; idx<13; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   total ttbar     "<<hists_photoneta[12]->Integral()<<endl;
cout<<"   error  total ttbar   "<<hists_photoneta[12]->GetBinError(1)<<endl;

for(unsigned int idx=16; idx<18; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral   diboson     "<<hists_photoneta[17]->Integral()<<endl;
cout<<"   error diboson  "<<hists_photoneta[17]->GetBinError(1)<<endl;

hists_photoneta[20]->Add(hists_photoneta[19]);
cout<<"  integral   singletop+photon     "<<hists_photoneta[20]->Integral()<<endl;
cout<<"   error singletop+photon   "<<hists_photoneta[20]->GetBinError(1)<<endl;

for(unsigned int idx=22; idx<27; ++idx){
hists_photoneta[idx]->Add(hists_photoneta[idx-1]);
}
cout<<"  integral    gamma jet      "<<hists_photoneta[26]->Integral()<<endl;
cout<<"   error   gamma jet    "<<hists_photoneta[26]->GetBinError(1)<<endl;

cout<<"****************************************************************"<<endl;
cout<<"integral totalMC     "<<histstotal->Integral()<<endl;
cout<<" error  totalMC   "<<histstotal->GetBinError(1)<<endl;

}
