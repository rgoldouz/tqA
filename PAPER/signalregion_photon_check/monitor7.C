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
std::vector<string> datasamples_;

float scales[] = {0.0961,0.0978,1.491,0.628,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024,4*19.145*0.0169};
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
samples_.push_back("SIGNAL.root");
datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");



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
/*  variables2_.push_back("photon_Pt");
  variables2_.push_back("photon_Eta");
  variables2_.push_back("muon_Pt");
  variables2_.push_back("muon_Eta");
  variables2_.push_back("Jet_Pt");
  variables2_.push_back("Jet_Eta");
  variables2_.push_back("Jet_Multiplicity");
  variables2_.push_back("topmass");
  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
  variables2_.push_back("Deltaphiphotonjet");
 variables2_.push_back("Deltaphiphotonmuon");
  variables2_.push_back("Deltaphimuonjet");
  variables2_.push_back("DeltaRphotonmuon");
  variables2_.push_back("DeltaRphotonjet");
  variables2_.push_back("DeltaRmuonjet");
  variables2_.push_back("HT");
  variables2_.push_back("Photonmuonmass");
  variables2_.push_back("Costopphoton");
  variables2_.push_back("Coswphoton");
  variables2_.push_back("sigmaetaetaendcap");
  variables2_.push_back("sigmaetaetabarrel");
  variables2_.push_back("top_Pt");
  variables2_.push_back("top_Eta");
  variables2_.push_back("b_tag_info");
  variables2_.push_back("top_photon_mass");
 variables2_.push_back("bJet_Multiplicity");
  variables2_.push_back("muoncharge");
  variables2_.push_back("Deltaphiphotonmet");
*/
  variables2_.push_back("Nvertex");
  variables2_.push_back("NvertexNOPU");



for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get((std::string("STEP1/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
std::vector<TH1F*> datahists;
for(unsigned int idx=0; idx<datafiles.size(); ++idx){
datahists.push_back((TH1F*)datafiles[idx]->Get((std::string("STEP1/").append(variables2_[i]).c_str())));
}
float lumi = 19.145;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};

for(unsigned int idx=0; idx<samples_.size(); ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}
//hists[0]->Scale(3173.42/hists[0]->Integral());
//hists[3]->Scale(620.326/hists[3]->Integral());

for(unsigned int idx=1; idx<samples_.size()-1; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

datahists[idx]->Add(datahists[idx-1]);
 }

// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);

hists[20]->SetMinimum(0.01);
hists[20]->SetMaximum(3*datahists[1]->GetMaximum());
hists[20]->SetFillColor(kMagenta+2);
hists[20]->Draw();
//hists[19]->SetFillColor(kOrange-2);
//hists[19]->Draw("same");
hists[18]->SetFillColor(41);
hists[18]->Draw("same");
hists[17]->SetFillColor(kOrange-2);
hists[17]->Draw("same");
hists[16]->SetFillColor(kRed);
hists[16]->Draw("same");
hists[15]->SetFillColor(kViolet+1);
hists[15]->Draw("same");
hists[14]->SetFillColor(kSpring-9);
hists[14]->Draw("same");
hists[13]->SetFillColor(32);
hists[13]->Draw("same");
hists[12]->SetFillColor(6);
hists[12]->Draw("same");
hists[9]->SetFillColor(4);
hists[9]->Draw("same");
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
hists[3]->SetFillColor(7);
hists[3]->Draw("same");
hists[2]->SetFillColor(8);
hists[2]->Draw("same");
hists[1]->SetFillColor(kOrange+7);
hists[1]->Draw("same");
hists[0]->SetFillColor(5);
hists[0]->Draw("same");
hists[21]->SetFillColor(1);
hists[21]->SetFillStyle(3004);
hists[21]->Draw("same");

 // plot data points
datahists[2]->SetLineWidth(3.);
datahists[2]->SetLineColor(kBlack);
datahists[2]->SetMarkerColor(kBlack);
datahists[2]->SetMarkerStyle(20.);
datahists[2]->Draw("esame");
conv->RedrawAxis();

TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( hists[0], "W PH JET"                           , "F");
  leg->AddEntry( hists[1], "Z JET"                           , "F");
  leg->AddEntry( hists[2], "PH JET"                           , "F");
  leg->AddEntry( hists[3], "W JET"              , "F");
//  leg->AddEntry( hists[5], "TOP-W-CH"               , "F");
//  leg->AddEntry( hists[5], "T-S-CH"               , "F");
//  leg->AddEntry( hists[7], "TOP-S-CH"               , "F");
//  leg->AddEntry( hists[7], "TTBAR-CH"               , "F");
//  leg->AddEntry( hists[8], "TBAR-W-CH"               , "F");
  leg->AddEntry( hists[9], "SINGLE TOP	"               , "F");
  leg->AddEntry( hists[12], "TTBAR"               , "F");
  leg->AddEntry( hists[13], "TTG"               , "F");
  leg->AddEntry( hists[14], "WWG"               , "F");
  leg->AddEntry( hists[15], "WW"               , "F");
  leg->AddEntry( hists[16], "WZ"               , "F");
  leg->AddEntry( hists[17], "ZZ"               , "F");
  leg->AddEntry( hists[18], "ZGAMMA"               , "F");
  leg->AddEntry( hists[20], "SINGLE TOP + PHOTON"               , "F");
  leg->AddEntry( hists[21], "SIGNAL"               , "F");
  leg->AddEntry( datahists[2], "CMS Data 2012(19.145 /fb)" , "F"); 

  leg->Draw("same");
}
}
