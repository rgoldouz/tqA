#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void top() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;

float scales[] = {0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341};
samples_.push_back("T-W-CH.root");
samples_.push_back("TBAR-W-CH.root");
samples_.push_back("T-S-CH.root");
samples_.push_back("TBAR-S-CH.root");
samples_.push_back("T-T-CH.root");
samples_.push_back("TBAR-T-CH.root");
samples_.push_back("TTBAR1.root");
samples_.push_back("TTBAR2.root");
samples_.push_back("TTBAR3.root");



//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}


/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("photon_Pt");
//  variables2_.push_back("photon_Eta");
//  variables2_.push_back("muon_Pt");
//  variables2_.push_back("muon_Eta");
//  variables2_.push_back("Jet_Pt");
 // /variables2_.push_back("Jet_Eta");
//  variables2_.push_back("Jet_Multiplicity");
  variables2_.push_back("topmass");
  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
//  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
//  variables2_.push_back("Deltaphiphotonjet");
// variables2_.push_back("Deltaphiphotonmuon");
//  variables2_.push_back("Deltaphimuonjet");
//  variables2_.push_back("DeltaRphotonmuon");
//  variables2_.push_back("DeltaRphotonjet");
//  variables2_.push_back("DeltaRmuonjet");
//  variables2_.push_back("HT");
//  variables2_.push_back("Photonmuonmass");
//  variables2_.push_back("Costopphoton");
//  variables2_.push_back("sigmaetaetaendcap");
//  variables2_.push_back("sigmaetaetabarrel");
//  variables2_.push_back("top_Pt");
//  variables2_.push_back("top_Eta");
//  variables2_.push_back("b_tag_info");
//  variables2_.push_back("top_photon_mass");
 // variables2_.push_back("bJet_Multiplicity");
//  variables2_.push_back("muoncharge");
//  variables2_.push_back("Deltaphiphotonmet");


for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
float lumi = 10;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }

// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);

hists[8]->SetFillColor(kSpring-9);
hists[8]->Draw();
}}
