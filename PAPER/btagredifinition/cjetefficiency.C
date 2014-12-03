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

void cjetefficiency() { 	  	  	
std::vector<string> samples_;
samples_.push_back("WPHJET.root");
//samples_.push_back("WW.root");
//samples_.push_back("WZ.root");
//samples_.push_back("ZZ.root");
//samples_.push_back("SIGNALtGC.root");
//samples_.push_back("SIGNALtGu.root");

//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));}
std::vector<string> variables2_;
/*
variables2_.push_back("btageffL");
variables2_.push_back("btageffM");
variables2_.push_back("btageffT");

variables2_.push_back("ctageffL");
variables2_.push_back("ctageffM");
variables2_.push_back("ctageffT");
*/
variables2_.push_back("ctageffL");
variables2_.push_back("ctageffM");
variables2_.push_back("ctageffT");


std::vector<double> finaleff;

Double_t bdiscriminator[35];
for (unsigned int i=0; i<15; ++i){
bdiscriminator[i] =0.244-0.1+i*0.02; }
for (unsigned int i=0; i<10; ++i){
bdiscriminator[15+i] = 0.679-0.1+i*0.02;}
for (unsigned int i=0; i<10; ++i){
bdiscriminator[25+i] = 0.898-0.1+i*0.02;}

Double_t count[3];
count[0]=15;
count[1]=10;
count[2]=10;

for(unsigned int ii=0; ii<variables2_.size(); ++ii){
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){hists.push_back((TH1F*)files[idx]->Get((std::string("STEP1/").append(variables2_[ii]).c_str())));}

for(unsigned int idx=1; idx<files.size(); ++idx){
hists[idx]->Add(hists[idx-1]);}

//hists[files.size()-1]->Draw();
std::vector<double> eff;
int bin; 

for ( bin=1; bin<count[ii]*2+1; ++bin){eff.push_back(hists[files.size()-1]->GetBinContent(bin));}
for (unsigned int i=0; i<count[ii]; ++i){finaleff.push_back(eff[2*i+1]/(eff[2*i]+eff[2*i+1])); 
//cout<< eff[2*i+1]<<"                "<<eff[2*i+1]/(eff[2*i]+eff[2*i+1])<<endl;
}
}

//   TGraph *gr = new TGraph(30,bdiscriminator,finaleff);

Double_t effi[35];
for (unsigned int i=0; i<35; ++i){
effi[i] =finaleff[i]; }
   TGraph *gr = new TGraph(35,effi,bdiscriminator);
   TGraph *grREV = new TGraph(35,bdiscriminator,effi);



   cout<<"the new working points for nominal SF"<<endl;
   cout<<"loose        "<<gr->Eval(1.008*grREV->Eval(0.244))<<endl;
   cout<<"medium        "<<gr->Eval(0.963*grREV->Eval(0.679))<<endl;
   cout<<"tight        "<<gr->Eval(0.947*grREV->Eval(0.898))<<endl;


   cout<<"the new working points for SF UP"<<endl;
   cout<<"loose        "<<gr->Eval((1.008+2*0.023)*grREV->Eval(0.244))<<endl;
   cout<<"medium        "<<gr->Eval((0.963+2*0.020)*grREV->Eval(0.679))<<endl;
   cout<<"tight        "<<gr->Eval((0.947+2*0.025)*grREV->Eval(0.898))<<endl;


   cout<<"the new working points for SF DOWN"<<endl;
   cout<<"loose        "<<gr->Eval((1.008-2*0.023)*grREV->Eval(0.244))<<endl;
   cout<<"medium        "<<gr->Eval((0.963-2*0.020)*grREV->Eval(0.679))<<endl;
   cout<<"tight        "<<gr->Eval((0.947-2*0.025)*grREV->Eval(0.898))<<endl;

//   cout<<gr->Eval(0.979*grREV->Eval(0.244))<<endl;
//   cout<<gr->Eval(0.966*grREV->Eval(0.679))<<endl;
//   cout<<gr->Eval(0.931*grREV->Eval(0.898))<<endl;

//   cout<<gr->Eval(0.0990)<<endl;
 //  cout<<gr->Eval(0.0142)<<endl;
//   cout<<gr->Eval(0.00155)<<endl;
//   cout<<gr->Eval(0.0142)<<endl;
 //  cout<<gr->Eval(0.0016)<<endl;
   gr->SetLineColor(4);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);

   gr->Draw("AP");

}
 


 

