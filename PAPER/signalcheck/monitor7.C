#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void monitor7() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;

samples_.push_back("TESTSIGNAL005.root");
samples_.push_back("TESTSIGNAL015.root");
//samples_.push_back("TESTSIGNAL02.root");
samples_.push_back("TESTSIGNAL01R.root");



//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}

/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("GENphoton_Pt");
  variables2_.push_back("GENphoton_Eta");
  variables2_.push_back("GENmuon_Pt");
  variables2_.push_back("GENmuon_Eta");
  variables2_.push_back("GENJet_Pt");
  variables2_.push_back("GENJet_Eta");

for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get((std::string("TEST/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}

for(unsigned int idx=0; idx<samples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

hists[idx]->Scale(1/hists[idx]->Integral());
 }

// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);
hists[0]->SetLineColor(kRed);
hists[0]->Draw(); 
hists[1]->SetLineColor(kGreen);
hists[1]->Draw("same");
hists[2]->SetLineColor(kBlue);
hists[2]->Draw("same");
TLegend* legmuo =new TLegend(0.35,0.6,0.85,0.90);

  legmuo->SetFillStyle ( 0);
  legmuo->SetFillColor ( 0);
  legmuo->SetBorderSize( 0);
  legmuo->AddEntry(hists[0] , "005"                  , "F");
  legmuo->AddEntry( hists[1], "015"                           , "F");
  legmuo->AddEntry( hists[2], "02"              , "F");
  legmuo->Draw("same");



}
}
