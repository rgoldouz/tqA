#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

void pileupdata(){
std::vector<string> samples_;
std::vector<string> datasamples_;
//samples_.push_back("MyDataPileupHistogram.root");
samples_.push_back("MyJAN2013DataPileupHistogram.root");
//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}
std::vector<string> variables2_;
variables2_.push_back("pileup");
// load histograms
std::vector<TH1F*> hists;
for(unsigned int idx=0; idx<files.size(); ++idx){
hists.push_back((TH1F*)files[idx]->Get(variables2_[0].c_str()));}
hists.push_back(new TH1F(  "NV", "NV", 60,0, 60));
hists.push_back(new TH1F(  "NVcor", "NVcor", 60,0, 60));
hists[0]->Scale(1/hists[0]->Integral());
float npuSummer12_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
      5.005E-06};
for(unsigned int idx=1; idx<61; ++idx){cout<<idx<<"              "<<hists[0]->GetBinContent(idx)/npuSummer12_S10[idx-1]<<endl;
//cout<<idx<<"              "<<hists[0]->GetBinContent(idx)<<endl;
//hists[1]->SetBinContent(idx,npuSummer12_S10[idx-1]);
hists[1]->SetBinContent(idx,npuSummer12_S10[idx-1]);
}
cout<<"integralMC"<<hists[1]->Integral()<<endl;
cout<<"integralDATA"<<hists[0]->Integral()<<endl;
hists[0]->SetMaximum(0.1);
hists[0]->SetLineColor(kRed);
hists[0]->SetLineWidth(2);
hists[0]->SetMarkerStyle(20);
hists[0]->SetMarkerSize(1);
hists[0]->SetMarkerColor(kRed);
hists[0]->GetXaxis()->SetTitle("number of true PU interaction");
hists[0]->GetYaxis()->SetRangeUser(-0.005, 0.08);
hists[0]->Draw("P");
hists[1]->SetLineColor(kBlue);
hists[1]->SetLineWidth(1);
hists[1]->SetMarkerStyle(22);
hists[1]->SetMarkerSize(1);
hists[1]->SetMarkerColor(kBlue);
hists[1]->Draw("Psame");

    TLegend *leg1 = new TLegend(0.6, 0.6, 0.7, 0.8);
    leg1->SetTextSize(0.03);
    leg1->SetBorderSize(0);
   leg1->SetLineColor(0);
    leg1->SetLineWidth(0);
    leg1->SetFillColor(kWhite);
    leg1->AddEntry(hists[0], "data number of true interaction", "p");
    leg1->AddEntry(hists[1], "summer 12 MC number of interaction", "p");
 
    leg1->Draw();

}



