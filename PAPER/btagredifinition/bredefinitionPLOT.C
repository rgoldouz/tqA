#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TPaveText.h"
void bredefinitionPLOT(){
const int N=5;
double X[N]={0,0.244,0.679,0.898,1};
double newlight[N]={0,0.2305,0.6117,0.8777,1 };
double newlightup[N]={0,0.2185,0.5921,0.8724,1 };
double newlightdown[N]={0,0.2425,0.6641,0.8981,1 };

double newB[N]={0,0.2237,0.7756,0.9125,1};
double newBup[N]={0,0.1739,0.7342,0.9057,1};
double newBdown[N]={0,0.2856,0.8066,0.9193,1};

double newC[N]={0,0.2414,0.7026,0.9022,1};
double newCup[N]={0,0.2267,0.6751,0.8982,1};
double newCdown[N]={0,0.2521,0.7213,0.9062,1};

TGraph *bjets= new TGraph(N,X,newB);
TGraph *bjetsup= new TGraph(N,X,newBup);
TGraph *bjetsdown= new TGraph(N,X,newBdown);

TGraph *cjets= new TGraph(N,X,newC);
TGraph *cjetsup= new TGraph(N,X,newCup);
TGraph *cjetsdown= new TGraph(N,X,newCdown);

TGraph *lightjets= new TGraph(N,X,newlight);
TGraph *lightjetsup= new TGraph(N,X,newlightup);
TGraph *lightjetsdown= new TGraph(N,X,newlightdown);

  bjets->SetLineColor(4);
   bjets->SetLineWidth(3);
   bjets->SetMarkerColor(1);
   bjets->SetMarkerStyle(20);
//   bjets->SetLineStyle(2);
   bjets->GetYaxis()->SetTitle("CSV modified");
   bjets->GetXaxis()->SetTitle("CSV original");
   bjets->SetTitle("b-discriminator ");
   bjets->SetFillColor(10);
   bjets->GetXaxis()->SetRangeUser(0, 1);
   bjets->GetYaxis()->SetRangeUser(0, 1);

   bjetsup->SetLineColor(4);
 //  bjetsup->SetLineWidth(3);
   bjetsup->SetMarkerColor(1);
   bjetsup->SetMarkerStyle(20);
   bjetsup->SetLineStyle(2);

   bjetsdown->SetLineColor(4);
 //  bjetsdown->SetLineWidth(3);
   bjetsdown->SetMarkerColor(1);
   bjetsdown->SetMarkerStyle(20);
   bjetsdown->SetLineStyle(2);

   cjets->SetLineColor(3);
   cjets->SetLineWidth(3);
   cjets->SetMarkerColor(1);
   cjets->SetMarkerStyle(20);
//   cjets->SetLineStyle(2);
   cjets->GetYaxis()->SetTitle("CSV modified");
   cjets->GetXaxis()->SetTitle("CSV original");
   cjets->SetFillColor(10);

   cjetsup->SetLineColor(3);
  // cjetsup->SetLineWidth(3);
   cjetsup->SetMarkerColor(1);
   cjetsup->SetMarkerStyle(20);
   cjetsup->SetLineStyle(2);

   cjetsdown->SetLineColor(3);
  // cjetsdown->SetLineWidth(3);
   cjetsdown->SetMarkerColor(1);
   cjetsdown->SetMarkerStyle(20);
   cjetsdown->SetLineStyle(2);

   lightjets->SetLineColor(2);
   lightjets->SetLineWidth(3);
   lightjets->SetMarkerColor(1);
   lightjets->SetMarkerStyle(20);
 //  lightjets->SetLineStyle(2);
   lightjets->GetYaxis()->SetTitle("Cross Section [pb]");
   lightjets->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   lightjets->SetFillColor(10);

   lightjetsup->SetLineColor(2);
 //  lightjetsup->SetLineWidth(3);
   lightjetsup->SetMarkerColor(1);
   lightjetsup->SetMarkerStyle(20);
   lightjetsup->SetLineStyle(2);

   lightjetsdown->SetLineColor(2);
 //  lightjetsdown->SetLineWidth(3);
   lightjetsdown->SetMarkerColor(1);
   lightjetsdown->SetMarkerStyle(20);
   lightjetsdown->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1","upper limit results ",200,10,700,500);
  c1->SetGrid();

   bjets->Draw("ACP");
   bjetsup->Draw("Csame");
   bjetsdown->Draw("Csame");

   lightjets->Draw("CPsame");
   lightjetsup->Draw("Csame");
   lightjetsdown->Draw("Csame");

   cjets->Draw("CPsame");
   cjetsup->Draw("Csame");
   cjetsdown->Draw("Csame");

    TLegend *leg1 = new TLegend(0.2, 0.6, 0.4, 0.8);
    leg1->SetTextSize(0.05);
    leg1->SetBorderSize(0);
    leg1->SetLineColor(0);
    leg1->SetLineWidth(0);
    leg1->SetFillColor(kWhite);

    leg1->AddEntry(bjets, "b-jet", "L");
    leg1->AddEntry(cjets, "c-jet", "L");
    leg1->AddEntry(lightjets, "light jet", "L");
    leg1->Draw();
}


