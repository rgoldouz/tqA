#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void CompareDifferentSigmaIeta() {
 	  	  	
// list of valid histogram names
// list of input samples
std::vector<string> datasamplesmean_;
std::vector<string> datasamplesup10_;
std::vector<string> datasamplesup20_;


datasamplesmean_.push_back("CRmean/REALDATA1.root");
datasamplesmean_.push_back("CRmean/REALDATA2.root");
datasamplesmean_.push_back("CRmean/REALDATA3.root");

datasamplesup10_.push_back("CR10up/REALDATA1.root");
datasamplesup10_.push_back("CR10up/REALDATA2.root");
datasamplesup10_.push_back("CR10up/REALDATA3.root");

datasamplesup20_.push_back("CR20up/REALDATA1.root");
datasamplesup20_.push_back("CR20up/REALDATA2.root");
datasamplesup20_.push_back("CR20up/REALDATA3.root");

//open files

std::vector<TFile*> datafilesmean;
for(unsigned int idx=0; idx<datasamplesmean_.size(); ++idx){
datafilesmean.push_back(new TFile(datasamplesmean_[idx].c_str()));
}

std::vector<TFile*> datafilesup10;
for(unsigned int idx=0; idx<datasamplesup10_.size(); ++idx){
datafilesup10.push_back(new TFile(datasamplesup10_[idx].c_str()));
}

std::vector<TFile*> datafilesup20;
for(unsigned int idx=0; idx<datasamplesup20_.size(); ++idx){
datafilesup20.push_back(new TFile(datasamplesup20_[idx].c_str()));
}

/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("photon_Pt");

for(unsigned int i=0; i<variables2_.size(); ++i){

std::vector<TH1F*> datahistsmean;
std::vector<TH1F*> datahistsup10;
std::vector<TH1F*> datahistsup20;

for(unsigned int idx=0; idx<datafilesmean.size(); ++idx){
datahistsmean.push_back((TH1F*)datafilesmean[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
}
for(unsigned int idx=0; idx<datafilesup10.size(); ++idx){
datahistsup10.push_back((TH1F*)datafilesup10[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
}
for(unsigned int idx=0; idx<datafilesup20.size(); ++idx){
datahistsup20.push_back((TH1F*)datafilesup20[idx]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));
}

for(unsigned int idx=1; idx<datasamplesmean_.size(); ++idx){
datahistsmean[idx]->Add(datahistsmean[idx-1]);
 }
for(unsigned int idx=1; idx<datasamplesup10_.size(); ++idx){
datahistsup10[idx]->Add(datahistsup10[idx-1]);
 }
for(unsigned int idx=1; idx<datasamplesup20_.size(); ++idx){
datahistsup20[idx]->Add(datahistsup20[idx-1]);
 }
cout<<"etarev  "<<datahistsmean[2]->Integral()<<endl;
cout<<"etarev10  "<<datahistsup10[2]->Integral()<<endl;
cout<<"etarev20  "<<datahistsup20[2]->Integral()<<endl;

// setup the canvas and draw the histograms
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
double norm=1;
double scalemean =norm/datahistsmean[2]->Integral();
datahistsmean[2]->Scale(scalemean);
double scaleup10 =norm/datahistsup10[2]->Integral();
datahistsup10[2]->Scale(scaleup10);
double scaleup20 =norm/datahistsup20[2]->Integral();
datahistsup20[2]->Scale(scaleup20);

 // plot data points
datahistsmean[2]->SetLineWidth(2);
datahistsmean[2]->SetLineColor(kGreen);
datahistsmean[2]->SetMarkerColor(kGreen);
datahistsmean[2]->SetMarkerStyle(20.);
datahistsmean[2]->Draw();

datahistsup10[2]->SetLineWidth(2);
datahistsup10[2]->SetLineColor(kRed);
datahistsup10[2]->SetMarkerColor(kRed);
datahistsup10[2]->SetMarkerStyle(20.);
datahistsup10[2]->Draw("same");

datahistsup20[2]->SetLineWidth(2);
datahistsup20[2]->SetLineColor(kBlue);
datahistsup20[2]->SetMarkerColor(kBlue);
datahistsup20[2]->SetMarkerStyle(20.);
datahistsup20[2]->Draw("same");

conv->RedrawAxis();

TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( datahistsmean[2], "Sigmaieta reverced  "                           , "F");
  leg->AddEntry( datahistsup10[2], "Sigmaieta + 10% * Sigmaieta"                           , "F");
  leg->AddEntry( datahistsup20[2], "Sigmaieta + 20% * Sigmaieta"                           , "F");

  leg->Draw("same");
}
}

