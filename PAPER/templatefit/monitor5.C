#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void monitor5() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
float scales[] = {1.2*0.0961,0.0978,1.491,0.628,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,1.2*0.02,10*0.0169};
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
samples_.push_back("SIGNAL.root");
samples_.push_back("REALDATA.root");



//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}

/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("photon Pt");
  variables2_.push_back("photon Eta");
  variables2_.push_back("muon Pt");
  variables2_.push_back("muon Eta");
//  variables2_.push_back("Jet Pt");
//  variables2_.push_back("Jet Eta");
//  variables2_.push_back("Jet Multiplicity");
  variables2_.push_back("topmass");
  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
//  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
 // variables2_.push_back("Deltaphiphotonjet");
//  variables2_.push_back("Deltaphiphotonmuon");
//  variables2_.push_back("Deltaphimuonjet");
//  variables2_.push_back("DeltaRphotonmuon");
//  variables2_.push_back("DeltaRphotonjet");
//  variables2_.push_back("DeltaRmuonjet");
//  variables2_.push_back("HT");
//  variables2_.push_back("Photonmuonmass");
  variables2_.push_back("Costopphoton");
  variables2_.push_back("sigmaetaetaendcap");
  variables2_.push_back("sigmaetaetabarrel");



for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
float lumi = 5.319;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};

for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=1; idx<samples_.size()-2; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);

hists[18]->SetMinimum(0.01);
hists[18]->SetMaximum(2*hists[20]->GetMaximum());
hists[18]->SetFillColor(41);
hists[18]->Draw();
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
hists[19]->SetFillColor(1);
hists[19]->SetFillStyle(3004);
hists[19]->Draw("same");

 // plot data points
hists[20]->SetLineWidth(3.);
hists[20]->SetLineColor(kBlack);
hists[20]->SetMarkerColor(kBlack);
hists[20]->SetMarkerStyle(20.);
hists[20]->Draw("esame");
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
  leg->AddEntry( hists[19], "SIGNAL"               , "F");
  leg->AddEntry( hists[20], "CMS Data 2012(5.31/fb)"               , "PL");

  leg->Draw("same");
}
////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables1_;
//  variables1_.push_back("photon Pt");
//  variables1_.push_back("photon Eta");
//  variables2_.push_back("muon Pt");
//  variables2_.push_back("muon Eta");
//  variables2_.push_back("Jet Pt");
//  variables2_.push_back("Jet Eta");
//  variables2_.push_back("Jet Multiplicity");
//  variables1_.push_back("topmass");
//  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
//  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
 // variables2_.push_back("Deltaphiphotonjet");
//  variables2_.push_back("Deltaphiphotonmuon");
//  variables2_.push_back("Deltaphimuonjet");
//  variables2_.push_back("DeltaRphotonmuon");
//  variables2_.push_back("DeltaRphotonjet");
//  variables2_.push_back("DeltaRmuonjet");
//  variables2_.push_back("HT");
//  variables2_.push_back("Photonmuonmass");
//  variables1_.push_back("Costopphoton");


for(unsigned int i=0; i<variables1_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables1_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
float lumi = 5.319;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};


for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=1; idx<samples_.size()-2; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("Signalregion/").append(variables1_[i]).c_str(), std::string("Signalregion/").append(variables1_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);

hists[18]->SetMinimum(0.01);
hists[18]->SetMaximum(2*hists[20]->GetMaximum());
hists[18]->SetFillColor(41);
hists[18]->Draw();
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
hists[19]->SetFillColor(1);
hists[19]->SetFillStyle(3004);
hists[19]->Draw("same");

 // plot data points
hists[20]->SetLineWidth(3.);
hists[20]->SetLineColor(kBlack);
hists[20]->SetMarkerColor(kBlack);
hists[20]->SetMarkerStyle(20.);
hists[20]->Draw("esame");
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
  leg->AddEntry( hists[19], "SIGNAL"               , "F");
  leg->AddEntry( hists[20], "CMS Data 2012(5.31/fb)"               , "PL");

  leg->Draw("same");
}
////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables0_;
//  variables0_.push_back("photon Pt");
//  variables0_.push_back("photon Eta");
//  variables2_.push_back("muon Pt");
//  variables2_.push_back("muon Eta");
//  variables2_.push_back("Jet Pt");
//  variables2_.push_back("Jet Eta");
//  variables2_.push_back("Jet Multiplicity");
//  variables0_.push_back("topmass");
//  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
//  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
//  variables2_.push_back("Deltaphiphotonjet");
//  variables2_.push_back("Deltaphiphotonmuon");
//  variables2_.push_back("Deltaphimuonjet");
//  variables2_.push_back("DeltaRphotonmuon");
//  variables2_.push_back("DeltaRphotonjet");
//  variables2_.push_back("DeltaRmuonjet");
//  variables2_.push_back("HT");
//  variables2_.push_back("Photonmuonmass");
//  variables2_.push_back("Costopphoton");


for(unsigned int i=0; i<variables0_.size(); ++i){

// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
//hists.push_back((TH1F*)files[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables0_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
}
float lumi = 5.319;
//float lumi = 4.429;

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};

for(unsigned int idx=0; idx<samples_.size()-1; ++idx){
hists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=1; idx<samples_.size()-2; ++idx){
//for(unsigned int idx=1; idx<2 ; ++idx){

   hists[idx]->Add(hists[idx-1]);
 }
// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("Sidebandregion/").append(variables0_[i]).c_str(), std::string("Sidebandregion/").append(variables0_[i]).c_str() , 600, 600);
conv->cd(0);
//conv->SetLogy(1);

hists[18]->SetMinimum(0.01);
hists[18]->SetMaximum(2*hists[20]->GetMaximum());
hists[18]->SetFillColor(41);
hists[18]->Draw();
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
hists[19]->SetFillColor(1);
hists[19]->SetFillStyle(3004);
hists[19]->Draw("same");

 // plot data points
hists[20]->SetLineWidth(3.);
hists[20]->SetLineColor(kBlack);
hists[20]->SetMarkerColor(kBlack);
hists[20]->SetMarkerStyle(20.);
hists[20]->Draw("esame");
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
  leg->AddEntry( hists[19], "SIGNAL"               , "F");
  leg->AddEntry( hists[20], "CMS Data 2012(5.31/fb)"               , "PL");

  leg->Draw("same");
}

}

