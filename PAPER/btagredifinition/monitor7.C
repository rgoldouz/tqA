#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void monitor7() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;

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
//samples_.push_back("SINGLE-TOP-PH.root");
//samples_.push_back("SINGLE-ANTITOP-PH.root");
//samples_.push_back("SIGNAL.root");

//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}


/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;
  variables2_.push_back("btagtotal");
  variables2_.push_back("btagRECO");
  variables2_.push_back("ctagtotal");
  variables2_.push_back("ctagRECO");
  variables2_.push_back("lighttotal");
  variables2_.push_back("lightRECO");

  TFile* outfile = new TFile("testd.root", "RECREATE");

std::vector<TH2F*> finalhists;

for(unsigned int i=0; i<variables2_.size(); ++i){

// load histograms
std::vector<TH2F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH2F*)files[idx]->Get((std::string("STEP4/").append(variables2_[i]).c_str())));
}

for(unsigned int idx=1; idx<samples_.size(); ++idx){
   hists[idx]->Add(hists[idx-1]);
 }
finalhists.push_back(hists[samples_.size()-1]);
 }
finalhists[1]->Divide(finalhists[0]);
finalhists[3]->Divide(finalhists[2]);
finalhists[5]->Divide(finalhists[4]);

/*for(unsigned int i=0; i<finalhists.size(); ++i){
cout<<"1                1"<<endl;
finalhists[i]->Write();
}*/
finalhists[1]->Write();
finalhists[3]->Write();
finalhists[5]->Write();


// setup the canvas and draw the histograms
}
