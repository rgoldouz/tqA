#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void pileupdistribution() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<string> step_;

float scales[] = {0.0961,1.491,0.628,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.032,0.024,19.145*0.0169};
//samples_.push_back("WPHJET.root");
samples_.push_back("ZJET.root");
//samples_.push_back("PHJET200400.root");
//samples_.push_back("WJET.root");
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
//samples_.push_back("SINGLE-TOP-PH.root");
//samples_.push_back("SINGLE-ANTITOP-PH.root");
//samples_.push_back("SIGNAL.root");
datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");
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
  variables2_.push_back("/Nvertex");
  variables2_.push_back("/NvertexNOPU");


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

//float scales[] = {0.0961,0.0978,1.491,1.708,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.1};
//float scales[] = {0.0961,0.0978,1.491,1.708,0.0244,0.0216,0.0108,0.01121,0.0133,0.0129,0.0199,0.1};
//float scales[] = {0.0961,0.116,1.484,2.108,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.006,0.0036,0.00099,1*0.0169};
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!STEP1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

for(unsigned int idx=0; idx<samples_.size(); ++idx){
//hists[idx]->Scale(lumi*scales[idx]);
}
for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<"              "<<hists[idx]->Integral()<<endl;
}

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
TCanvas* conv = new TCanvas(std::string("analyzestep2/").append(variables2_[i]).c_str(), std::string("analyzestep2/").append(variables2_[i]).c_str() , 600, 600);
conv->cd(0);
hists[samples_.size()-1]->Scale(1/hists[samples_.size()-1]->Integral());
datahists[2]->Scale(1/datahists[2]->Integral());

datahists[2]->SetLineColor(kRed+3);
datahists[2]->SetLineWidth(2);
datahists[2]->Draw();
hists[samples_.size()-1]->SetLineColor(kBlue+3);
hists[samples_.size()-1]->SetLineWidth(2);
hists[samples_.size()-1]->Draw("same");
}

}
