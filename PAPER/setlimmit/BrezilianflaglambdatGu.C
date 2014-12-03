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
void BrezilianflaglambdatGu(){
   
  TCanvas *c1 = new TCanvas("c1","upper limit results ",200,10,700,500);
  c1->SetGrid();

const int N=5;
double X[N]={5,10,20,30,40};
double sigma[N];
double lambda[N];
double eff=0.01864555;

double lumi=19145;
double sfactor=1.3746;


double nO[N]={1.92,0.93,0.5,0.326,0.208};
double nE[N]={3.71,1.9,0.988,0.663,0.479};


double nE1sigmaup[N]={5.31,2.53,1.41,0.907,0.666};
double nE1sigmadown[N]={2.77,1.36,0.739,0.494,0.353};

double nE2sigmaup[N]={7.65,3.71,1.95,1.27,0.968};
double nE2sigmadown[N]={2.05,1.11,0.567,0.394,0.265};
double sigma1up[N];
double sigma1down[N];
double sigma2up[N];
double sigma2down[N];
for (int idx=0; idx<N; ++idx){
sigma1up[idx]=nE1sigmaup[idx]-nE[idx];
sigma1down[idx]=nE[idx]-nE1sigmadown[idx];
sigma2up[idx]=nE2sigmaup[idx]-nE[idx];
sigma2down[idx]=nE[idx]-nE2sigmadown[idx];
sigma[idx]=X[idx]/(sfactor*lumi*eff);
lambda[idx]=sqrt(sigma[idx]/67.2);
cout<<lambda[idx]<<endl;
}

   TGraphAsymmErrors *grafexp1sigma=new TGraphAsymmErrors(N);
   grafexp1sigma->SetName("grafexp1sigma");
   grafexp1sigma->SetTitle("Graph");
   grafexp1sigma->SetFillColor(1);
//grafexp1sigma->SetPoint(point number ,x value,y value);
//grafexp1sigma->SetPointError(point number ,0,0,y lower error ,y upper error);
   grafexp1sigma->SetPoint(0,lambda[0],nE[0]);
   grafexp1sigma->SetPointError(0,0,0,sigma1down[0],sigma1up[0]);
   grafexp1sigma->SetPoint(1,lambda[1],nE[1]);
   grafexp1sigma->SetPointError(1,0,0,sigma1down[1],sigma1up[1]);
   grafexp1sigma->SetPoint(2,lambda[2],nE[2]);
   grafexp1sigma->SetPointError(2,0,0,sigma1down[2],sigma1up[2]);
   grafexp1sigma->SetPoint(3,lambda[3],nE[4]);
   grafexp1sigma->SetPointError(3,0,0,sigma1down[3],sigma1up[3]);
   grafexp1sigma->SetPoint(4,lambda[4],nE[4]);
   grafexp1sigma->SetPointError(4,0,0,sigma1down[4],sigma1up[4]);

   grafexp1sigma->SetFillColor(kGreen);
//   grafexp1sigma->
//   grafexp1sigma->
//   grafexp1sigma->
//   grafexp1sigma->

   TGraphAsymmErrors *grafexp2sigma=new TGraphAsymmErrors(N);
   grafexp2sigma->SetPoint(0,lambda[0],nE[0]);
   grafexp2sigma->SetPointError(0,0,0,sigma2down[0],sigma2up[0]);
   grafexp2sigma->SetPoint(1,lambda[1],nE[1]);
   grafexp2sigma->SetPointError(1,0,0,sigma2down[1],sigma2up[1]);
   grafexp2sigma->SetPoint(2,lambda[2],nE[2]);
   grafexp2sigma->SetPointError(2,0,0,sigma2down[2],sigma2up[2]);
   grafexp2sigma->SetPoint(3,lambda[3],nE[4]);
   grafexp2sigma->SetPointError(3,0,0,sigma2down[3],sigma2up[3]);
   grafexp2sigma->SetPoint(4,lambda[4],nE[4]);
   grafexp2sigma->SetPointError(4,0,0,sigma2down[4],sigma2up[4]);

   grafexp2sigma->SetFillColor(kYellow);
   grafexp2sigma->GetXaxis()->SetLabelFont(42);
   grafexp2sigma->GetXaxis()->SetLabelOffset(0.007);
   grafexp2sigma->GetXaxis()->SetLabelSize(0.03);
   grafexp2sigma->GetXaxis()->SetTitleSize(0.045);
   grafexp2sigma->GetXaxis()->SetTitleFont(42);
   grafexp2sigma->GetXaxis()->SetTitleOffset(0.80);
   grafexp2sigma->GetXaxis()->SetRangeUser(0, lambda[4]+0.1*lambda[4]);
   grafexp2sigma->GetYaxis()->SetRangeUser(0, 5);
   grafexp2sigma->GetYaxis()->SetTitle("95% CL Limit on #sigma / #sigma_{pred}");
   grafexp2sigma->GetYaxis()->SetLabelFont(42);
   grafexp2sigma->GetYaxis()->SetLabelOffset(0.007);
   grafexp2sigma->GetYaxis()->SetLabelSize(0.03);
   grafexp2sigma->GetYaxis()->SetTitleSize(0.045);
   grafexp2sigma->GetYaxis()->SetTitleOffset(0.9);
   grafexp2sigma->GetYaxis()->SetTitleFont(42);
   grafexp2sigma->GetXaxis()->SetTitle("#kappa_{tu#gamma}");
   grafexp2sigma->SetTitle("");
   grafexp2sigma->SetMaximum(8);

   TGraph *Expected= new TGraph(N,lambda,nE);

   Expected->SetLineColor(4);
   Expected->SetLineWidth(3);
   Expected->SetMarkerColor(1);
   Expected->SetMarkerStyle(20);
   Expected->SetLineStyle(2);
   Expected->GetYaxis()->SetTitle("Cross Section [pb]");
   Expected->GetXaxis()->SetTitle("#Lambda_{tu#gamma} [GeV]");
   Expected->SetTitle("with ");
   Expected->SetFillColor(10);



   TGraph *Observed= new TGraph(N,lambda,nO);
   Observed->SetLineColor(1);
   Observed->SetLineWidth(3);
   Observed->SetMarkerColor(1);
   Observed->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   Observed->GetYaxis()->SetTitle("Cross Section [pb]");
   Observed->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   Observed->SetTitle("with ");
   Observed->SetFillColor(10);


   grafexp2sigma->Draw("AL3");
   grafexp1sigma->Draw("L3same");
   Expected->Draw("LP");
   Observed->Draw("LP");

    TLegend *leg1 = new TLegend(0.6, 0.6, 0.7, 0.8);
    leg1->SetTextSize(0.03);
    leg1->SetBorderSize(0);
    leg1->SetLineColor(0);
    leg1->SetLineWidth(0);
    leg1->SetFillColor(kWhite);
    leg1->AddEntry(Observed, "95% CL Observed Limit", "L");
    leg1->AddEntry(Expected, "95% CL Expected Limit", "L");
    leg1->AddEntry(grafexp1sigma, "#pm1#sigma Exp.Limit", "F");
    leg1->AddEntry(grafexp2sigma, "#pm2#sigma Exp.Limit", "F");

    leg1->Draw();

    TLine *line1 = new TLine(lambda[0], 1, lambda[4], 1);
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
    line1->Draw("same");
  
    TPaveText *pt = new TPaveText(0.45,0.80,0.7,0.87, "NDC"); // NDC sets coords
    pt->SetLineColor(10);                                              // relative to pad dimensions
    pt->SetFillColor(10); // text is black on white
    pt->SetTextSize(0.03);
    pt->SetTextAlign(12);
    pt->AddText("CMS Preliminary #sqrt{s} = 8 TeV, #int L dt= 19.145 fb^{-1}");
    pt->SetShadowColor(10);
    pt->Draw("same");


}

