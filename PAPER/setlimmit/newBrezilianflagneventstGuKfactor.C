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
void newBrezilianflagneventstGuKfactor(){

 TF1 *fa = new TF1("fa","67.2*x*x",0.001,0.1); 
   
  TCanvas *c1 = new TCanvas("c1","upper limit results ",200,10,700,500);
  c1->SetGrid();

const int N=2;
/*
double X[N]={0.001,0.05};
double nO[N]={0.0298,0.0298};
double nE[N]={0.0408,0.0408};
*/
double X[N]={0.001,0.05};
double nO[N]={0.0217,0.0217};
double nE[N]={0.0408,0.0408};

cout<<"jadid"<<endl;
cout<<"expected tgammau coupling including k-factor"<<1.5*sqrt(nE[0]/(1.375*67.2))<<endl;
cout<<"observed tgammau coupling including k-factor"<<1.5*sqrt(nO[0]/(1.375*67.2))<<endl;
cout<<"ghadimi"<<endl;
cout<<"expected tgammau coupling including k-factor"<<sqrt(nE[0]/(1.375*67.2))<<endl;
cout<<"observed tgammau coupling including k-factor"<<sqrt(nO[0]/(1.375*67.2))<<endl;



double nE1sigmaup[N]={0.0555,0.0555};
double nE1sigmadown[N]={0.0293,0.0293};

double nE2sigmaup[N]={0.0847,0.0847};
double nE2sigmadown[N]={0.0224,0.0224};
double sigma1up[N];
double sigma1down[N];
double sigma2up[N];
double sigma2down[N];
for (int idx=0; idx<N; ++idx){
sigma1up[idx]=nE1sigmaup[idx]-nE[idx];
sigma1down[idx]=nE[idx]-nE1sigmadown[idx];
sigma2up[idx]=nE2sigmaup[idx]-nE[idx];
sigma2down[idx]=nE[idx]-nE2sigmadown[idx];
}

   TGraphAsymmErrors *grafexp1sigma=new TGraphAsymmErrors(N);
   grafexp1sigma->SetName("grafexp1sigma");
   grafexp1sigma->SetTitle("Graph");
   grafexp1sigma->SetFillColor(1);
//grafexp1sigma->SetPoint(point number ,x value,y value);
//grafexp1sigma->SetPointError(point number ,0,0,y lower error ,y upper error);
   grafexp1sigma->SetPoint(0,X[0],nE[0]);
   grafexp1sigma->SetPointError(0,0,0,sigma1down[0],sigma1up[0]);
   grafexp1sigma->SetPoint(1,X[1],nE[1]);
   grafexp1sigma->SetPointError(1,0,0,sigma1down[1],sigma1up[1]);
/*
   grafexp1sigma->SetPoint(2,X[2],nE[2]);
   grafexp1sigma->SetPointError(2,0,0,sigma1down[2],sigma1up[2]);
   grafexp1sigma->SetPoint(3,X[3],nE[3]);
   grafexp1sigma->SetPointError(3,0,0,sigma1down[3],sigma1up[3]);
   grafexp1sigma->SetPoint(4,X[4],nE[4]);
   grafexp1sigma->SetPointError(4,0,0,sigma1down[4],sigma1up[4]);
*/
   grafexp1sigma->SetFillColor(kGreen);
//   grafexp1sigma->
//   grafexp1sigma->
//   grafexp1sigma->
//   grafexp1sigma->

   TGraphAsymmErrors *grafexp2sigma=new TGraphAsymmErrors(N);
   grafexp2sigma->SetPoint(0,X[0],nE[0]);
   grafexp2sigma->SetPointError(0,0,0,sigma2down[0],sigma2up[0]);
   grafexp2sigma->SetPoint(1,X[1],nE[1]);
   grafexp2sigma->SetPointError(1,0,0,sigma2down[1],sigma2up[1]);
/*
   grafexp2sigma->SetPoint(2,X[2],nE[2]);
   grafexp2sigma->SetPointError(2,0,0,sigma2down[2],sigma2up[2]);
   grafexp2sigma->SetPoint(3,X[3],nE[3]);
   grafexp2sigma->SetPointError(3,0,0,sigma2down[3],sigma2up[3]);
   grafexp2sigma->SetPoint(4,X[4],nE[4]);
   grafexp2sigma->SetPointError(4,0,0,sigma2down[4],sigma2up[4]);
*/
   grafexp2sigma->SetFillColor(kYellow);
   grafexp2sigma->GetXaxis()->SetLabelFont(42);
   grafexp2sigma->GetXaxis()->SetLabelOffset(0.007);
   grafexp2sigma->GetXaxis()->SetLabelSize(0.034);
   grafexp2sigma->GetXaxis()->SetTitleSize(0.070);
   grafexp2sigma->GetXaxis()->SetTitleFont(22);
   grafexp2sigma->GetXaxis()->SetTitleOffset(0.60);
   grafexp2sigma->GetXaxis()->SetRangeUser(0, 0.2);
   grafexp2sigma->GetYaxis()->SetRangeUser(0, 0.2);
   grafexp2sigma->GetYaxis()->SetTitle("#sigma_{tu#gamma}*Br(w#rightarrowl#nu) [pb] (95% CL)");
   grafexp2sigma->GetYaxis()->SetLabelFont(42);
   grafexp2sigma->GetYaxis()->SetLabelOffset(0.007);
   grafexp2sigma->GetYaxis()->SetLabelSize(0.034);
   grafexp2sigma->GetYaxis()->SetTitleSize(0.05);
   grafexp2sigma->GetYaxis()->SetTitleOffset(0.8);
   grafexp2sigma->GetYaxis()->SetTitleFont(22);
   grafexp2sigma->GetXaxis()->SetTitle("#kappa_{u}");
   grafexp2sigma->SetTitle("");
 //  grafexp2sigma->SetMaximum(8);

   TGraph *Expected= new TGraph(N,X,nE);

   Expected->SetLineColor(4);
   Expected->SetLineWidth(2);
   Expected->SetMarkerColor(1);
   Expected->SetMarkerStyle(20);
   Expected->SetLineStyle(2);
   Expected->GetYaxis()->SetTitle("Cross Section [pb]");
   Expected->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   Expected->SetTitle("with ");
   Expected->SetFillColor(10);



   TGraph *Observed= new TGraph(N,X,nO);
   Observed->SetLineColor(1);
   Observed->SetLineWidth(2);
   Observed->SetMarkerColor(1);
   Observed->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   Observed->GetYaxis()->SetTitle("Cross Section [pb]");
   Observed->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   Observed->SetTitle("with ");
   Observed->SetFillColor(10);


   grafexp2sigma->Draw("AL3");
   grafexp1sigma->Draw("L3same");
   Expected->Draw("L");
   Observed->Draw("L");

   fa->SetLineColor(2);
   fa->SetLineWidth(3);
   fa->SetMarkerColor(1);
   fa->Draw("same");
  Observed->GetHistogram()->Draw("AXISSAMEY+");
  Observed->GetHistogram()->Draw("AXISSAMEX+");

    TLegend *leg1 = new TLegend(0.2, 0.6, 0.35, 0.8);
    leg1->SetTextSize(0.03);
    leg1->SetBorderSize(0);
    leg1->SetLineColor(0);
    leg1->SetLineWidth(0);
    leg1->SetFillColor(kWhite);
    leg1->AddEntry(fa, "Predicted", "L");
    leg1->AddEntry(Observed, "95% CL Observed Limit", "L");
    leg1->AddEntry(Expected, "95% CL Expected Limit", "L");
    leg1->AddEntry(grafexp1sigma, "#pm1#sigma Exp.Limit", "F");
    leg1->AddEntry(grafexp2sigma, "#pm2#sigma Exp.Limit", "F");

    leg1->Draw();

    TLine *line1 = new TLine(5, 1, 40, 1);
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
//    line1->Draw("same");
  
    TPaveText *pt = new TPaveText(0.1,0.95,0.4,0.95, "NDC"); // NDC sets coords
    pt->SetLineColor(10);                                              // relative to pad dimensions
    pt->SetFillColor(10); // text is black on white
    pt->SetTextSize(0.045);
    pt->SetTextAlign(12);
    pt->AddText("CMS Preliminary, 19.1 fb^{-1}, #sqrt{s} = 8 TeV");
    pt->SetShadowColor(10);
    pt->Draw("same");

}

