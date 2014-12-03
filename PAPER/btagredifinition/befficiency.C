#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

void bjetefficiency() { 	  	  	
std::vector<string> samples_;
samples_.push_back("WPHJET.root");
//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
variables2_.push_back("btageffL");
variables2_.push_back("btageffM");
variables2_.push_back("btageffT");
std::vector<TH1F*> hists;

for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms

//for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[0]->Get((std::string("analyzestep2/").append(variables2_[i]).c_str())));

//hists.push_back((TH1F*)files[idx]->Get((std::string("Signalregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get((std::string("Sidebandregion/").append(variables2_[i]).c_str())));
//hists.push_back((TH1F*)files[idx]->Get("analyzestep1/topmass"));
//}
}
cout<<hists[1]->GetBinContent(5)<<endl;

std::vector<double> effL;
std::vector<double> effM;
std::vector<double> effT;
int bin; 
for(unsigned int i=0; i<variables2_.size(); ++i){
for ( bin=1; bin<21; ++bin){
cout<<hists[i]->GetBinContent(bin)<<endl; 
if (i==0) effL.push_back(hists[i]->GetBinContent(bin));
if (i==1) effM.push_back(hists[i]->GetBinContent(bin));
if (i==2) effT.push_back(hists[i]->GetBinContent(bin));}}
Double_t bdiscriminator[30];
cout<<"  ***********************************************  "<<endl;
for (unsigned int i=0; i<10; ++i){
bdiscriminator[i] =0.244-0.1+i*0.02; }
for (unsigned int i=0; i<10; ++i){
bdiscriminator[10+i] = 0.679-0.1+i*0.02;}
for (unsigned int i=0; i<10; ++i){
bdiscriminator[20+i] = 0.898-0.1+i*0.02;}
cout<<"  ***********************************************  "<<endl;

Double_t finaleff[30];
for (unsigned int i=0; i<10; ++i){
finaleff[i] =effL[2*i+1]/(effL[2*i]+effL[2*i+1]); }
cout<<"  ***********************************************  "<<endl;
for (unsigned int i=0; i<10; ++i){
finaleff[10+i] =effM[2*i+1]/(effM[2*i]+effM[2*i+1]); }
cout<<"  ***********************************************  "<<endl;
for (unsigned int i=0; i<10; ++i){
finaleff[20+i] =effT[2*i+1]/(effT[2*i]+effT[2*i+1]); }
cout<<"  ***********************************************  "<<endl;

//   TGraph *gr = new TGraph(30,bdiscriminator,finaleff);


   TGraph *gr = new TGraph(30,finaleff,bdiscriminator);
   cout<<gr->Eval(0.1071)<<endl;
   cout<<gr->Eval(0.0190)<<endl;
   cout<<gr->Eval(0.00091)<<endl;
   gr->SetLineColor(4);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);

   gr->Draw("AP");

}
 


 

