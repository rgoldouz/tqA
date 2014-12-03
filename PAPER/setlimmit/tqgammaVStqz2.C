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
#include "TF2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TPaveLabel.h"
#include "TLatex.h"

void tqgammaVStqz2(){

TCanvas *c1 = new TCanvas("c1","upper limit results ",200,10,700,500);
c1->SetLogy(1);
c1->SetLogx(1);
//c1->Range(0,0,1,1);

double a=456.57;
double b=156.25;

double c=1000000/36;
double d=1/0.09;


TF2 *f2 = new TF2("f2","[0]*pow(x,2)+[1]*pow(y,2)",0.0000001,0.05,0.0000001,0.09);
f2->SetParameters(a,b);

Double_t l[1]={1};
f2->SetContour(1,l);
f2->SetFillColor(4);
f2->SetLineColor(4);
f2->SetLineWidth(3);
f2->SetLineStyle(1);
f2->SetMarkerColor(4);
f2->SetNpx(1000);
f2->SetNpy(1000);

TF2 *fhera = new TF2("fhera","[0]*pow(x,2)+[1]*pow(y,2)",0.0000000001,0.02,0.00000000001,0.5);
fhera->SetParameters(c,d);

Double_t lhera[1]={1};
fhera->SetContour(1,lhera);
fhera->SetFillColor(6);
fhera->SetLineColor(6);
fhera->SetLineWidth(3);
fhera->SetLineStyle(1);
fhera->SetMarkerColor(6);
fhera->SetNpx(1000);
fhera->SetNpy(1000);

const int N=3;
double cdfx[4]={0,0.032,0.032,1};
double cdfy[4]={0.037,0.037,0,0};

double dzerox[4]={0.000001,0.02,0.06,1};
double dzeroy[4]={0.032,0.032,0.032,0.032};

//double cdfx[7]={0,0,1,1,0.45,0.45,0};
//double cdfy[7]={0.032,0.1,0.1,0,0,0.032,0.032};
double lepx2[4]={0.0468,0.0468,0.0468,0.0468};
double lepy2[4]={0,0.001,0.0012,0.0013};
double lepx[4]={0,0.001,0.0011,0.0012};
double lepy[4]={0.08,0.08,0.08,0.08};

double herax2[4]={0.0060,0.0060,0.0060,0.0060};
double heray2[4]={0,0.01,0.01,0.01};

double herax[4]={0,0.00005,0.00051,0.00052};
double heray[4]={0.3,0.3,0.3,0.3};

double h1x[4]={0.0064,0.0064,0.0064,0.0064};
double h1y[4]={0.000001,0.02,0.06,1};
double cmsx[4]={0.000001,0.02,0.06,1};
double cmsy[4]={0.0005,0.0005,0.0005,0.0005};

double atlasx[4]={0.000001,0.02,0.06,1};
double atlasy[4]={0.0073,0.0073,0.0073,0.0073};

double mycmsx[4]={0.000161,0.000161,0.000161,0.000161};
double mycmsy[4]={0.000001,0.02,0.06,1};

double mycmscx[4]={0.00182,0.00182,0.00182,0.00182};
double mycmscy[4]={0.000001,0.02,0.06,1};



double cmsphux[4]={0.00001,0.000161,0.000161001,1};
double cmsphuy[4]={0,0,0,0};

double cmsphcx[4]={0.00182,0.00182,0.00182,0.00182};
double cmsphcy[4]={0,0.02,0.06,1};

   TGraph *cdf= new TGraph(4,cdfx,cdfy);
   cdf->SetLineColor(12);
   cdf->SetLineWidth(3);
   cdf->SetMarkerColor(1);
   cdf->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   cdf->GetYaxis()->SetTitle("Cross Section [pb]");
   cdf->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   cdf->SetTitle("with ");
   cdf->SetFillColor(45);
   cdf->GetXaxis()->SetLabelOffset(0.0001);

   TGraph *dzero= new TGraph(4,dzerox,dzeroy);
   dzero->SetLineColor(3);
   dzero->SetLineWidth(3);
   dzero->SetMarkerColor(1);
   dzero->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   dzero->GetYaxis()->SetTitle("Cross Section [pb]");
   dzero->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   dzero->SetTitle("with ");
   dzero->SetFillColor(45);
   dzero->GetXaxis()->SetLabelOffset(0.0001);

   TGraph *cms= new TGraph(4,mycmsx,mycmsy);
   cms->SetLineColor(2);
   cms->SetLineWidth(3);
   cms->SetMarkerColor(1);
   cms->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   cms->GetYaxis()->SetTitle("Cross Section [pb]");
   cms->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   cms->SetTitle("with ");
   cms->SetFillColor(45);

  TGraph *cmsc= new TGraph(4,mycmscx,mycmscy);
   cmsc->SetLineColor(2);
   cmsc->SetLineWidth(3);
   cmsc->SetMarkerColor(1);
   cmsc->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   cmsc->GetYaxis()->SetTitle("Cross Section [pb]");
   cmsc->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   cmsc->SetTitle("with ");
   cmsc->SetFillColor(45);

   TGraph *cmsz= new TGraph(4,cmsx,cmsy);
   cmsz->SetLineColor(1);
   cmsz->SetLineWidth(3);
   cmsz->SetMarkerColor(1);
   cmsz->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   cmsz->GetYaxis()->SetTitle("Cross Section [pb]");
   cmsz->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   cmsz->SetTitle("with ");
   cmsz->SetFillColor(45);

   TGraph *hera= new TGraph(4,herax,heray);
   hera->SetLineColor(6);
   hera->SetLineWidth(3);
   hera->SetMarkerColor(1);
   hera->SetMarkerStyle(20);
//   hera->SetLineStyle(2);
//   hera->SetLineStyle(2);
   hera->GetYaxis()->SetTitle("Cross Section [pb]");
   hera->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   hera->SetTitle("with ");
   hera->SetFillColor(45);

   TGraph *hera2= new TGraph(4,herax2,heray2);
   hera2->SetLineColor(6);
   hera2->SetLineWidth(3);
   hera2->SetMarkerColor(1);
   hera2->SetMarkerStyle(20);
//   hera->SetLineStyle(2);
//   hera->SetLineStyle(2);
   hera2->GetYaxis()->SetTitle("Cross Section [pb]");
   hera2->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   hera2->SetTitle("with ");
   hera2->SetFillColor(45);

   TGraph *h1= new TGraph(4,h1x,h1y);
   h1->SetLineColor(28);
   h1->SetLineWidth(3);
   h1->SetMarkerColor(1);
 //  h1->SetMarkerStyle(20);
//   hera->SetLineStyle(2);
   h1->GetYaxis()->SetTitle("Cross Section [pb]");
   h1->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
//   h1->SetLineStyle(2);
   h1->SetTitle("with ");
   h1->SetFillColor(45);

   TGraph *lep= new TGraph(4,lepx,lepy);
   lep->SetLineColor(4);
   lep->SetLineWidth(3);
   lep->SetMarkerColor(1);
   lep->SetMarkerStyle(20);
//   lep->SetLineStyle(2);
   lep->GetYaxis()->SetTitle("Cross Section [pb]");
   lep->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   lep->SetTitle("with ");
   lep->SetFillColor(39);

   TGraph *lep2= new TGraph(4,lepx2,lepy2);
   lep2->SetLineColor(4);
   lep2->SetLineWidth(3);
   lep2->SetMarkerColor(1);
   lep2->SetMarkerStyle(20);
//   lep->SetLineStyle(2);
   lep2->GetYaxis()->SetTitle("Cross Section [pb]");
   lep2->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   lep2->SetTitle("with ");
   lep2->SetFillColor(39);

   TGraph *atlas= new TGraph(4,atlasx,atlasy);
   atlas->SetLineColor(46);
   atlas->SetLineWidth(3);
   atlas->SetMarkerColor(1);
   atlas->SetMarkerStyle(20);
//   lep->SetLineStyle(2);
   atlas->GetYaxis()->SetTitle("Cross Section [pb]");
   atlas->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   atlas->SetTitle("with ");
   atlas->SetFillColor(39);

   TGraphAsymmErrors *grafphotoncms=new TGraphAsymmErrors(4);
   grafphotoncms->SetPoint(0,cmsphux[0],cmsphuy[0]);
   grafphotoncms->SetPointError(0,0,0,0,0);
   grafphotoncms->SetPoint(1,cmsphux[1],cmsphuy[1]);
   grafphotoncms->SetPointError(1,0,0,0,0);
   grafphotoncms->SetPoint(2,cmsphux[2],cmsphuy[2]);
   grafphotoncms->SetPointError(2,0,0,0,1);
   grafphotoncms->SetPoint(3,cmsphux[3],cmsphuy[3]);
   grafphotoncms->SetPointError(3,0,0,0,1);
//grafhera->SetFillColor(kYellow+3);
grafphotoncms->SetFillColor(kYellow-9);
//grafphotoncms->SetFillColor(kTeal-9);
//   grafphotoncms->SetFillStyle(3004);
   grafphotoncms->GetXaxis()->SetRangeUser(0.00005, 1);
   grafphotoncms->GetYaxis()->SetRangeUser(0.00001, 1);
   grafphotoncms->GetYaxis()->SetTitle("BR(t #rightarrow qZ)");
   grafphotoncms->GetXaxis()->SetTitle("BR(t #rightarrow q#gamma)");
   grafphotoncms->GetXaxis()->SetLabelSize(0.044);
   grafphotoncms->GetYaxis()->SetLabelSize(0.044);
   grafphotoncms->GetXaxis()->SetTitleFont(62);
   grafphotoncms->GetYaxis()->SetTitleFont(62);
   grafphotoncms->GetXaxis()->SetTitleSize(0.055);
   grafphotoncms->GetYaxis()->SetTitleSize(0.055);
   grafphotoncms->GetXaxis()->SetTitleOffset(1);
   grafphotoncms->GetYaxis()->SetTitleOffset(0.8);

   TGraphAsymmErrors *grafphotoncms2=new TGraphAsymmErrors(4);
   grafphotoncms2->SetPoint(0,cmsx[0],cmsy[0]);
   grafphotoncms2->SetPointError(0,0,0,0,1-cmsy[0]);
   grafphotoncms2->SetPoint(1,cmsphux[2],cmsy[1]);
   grafphotoncms2->SetPointError(1,0,0,0,1-cmsy[1]);
   grafphotoncms2->SetPoint(2,cmsphux[2],cmsphuy[2]);
   grafphotoncms2->SetPointError(2,0,0,0,1);
   grafphotoncms2->SetPoint(3,cmsphux[3],cmsphuy[3]);
   grafphotoncms2->SetPointError(3,0,0,0,1);
//grafhera->SetFillColor(kYellow+3);
//grafhera->SetFillColor(kRed-8);
grafphotoncms2->SetFillColor(kYellow-9);
//   grafphotoncms2->SetFillStyle(3005);
   grafphotoncms2->GetXaxis()->SetRangeUser(0.00005, 1);
   grafphotoncms2->GetYaxis()->SetRangeUser(0.00001, 1);
 //  grafphotoncms2->GetYaxis()->SetTitle("BR(t #rightarrow qZ)");
   grafphotoncms2->GetXaxis()->SetTitle("BR(t #rightarrow q#gamma)");
   grafphotoncms2->GetXaxis()->SetLabelSize(0.03);
   grafphotoncms2->GetYaxis()->SetLabelSize(0.03);
   grafphotoncms2->GetXaxis()->SetTitleFont(62);
   grafphotoncms2->GetYaxis()->SetTitleFont(62);
   grafphotoncms2->GetXaxis()->SetTitleSize(0.065);
   grafphotoncms2->GetYaxis()->SetTitleSize(0.045);
   grafphotoncms2->GetXaxis()->SetLabelOffset(0.0001);
   grafphotoncms2->GetXaxis()->SetTitleOffset(0.60);


   TGraphAsymmErrors *grafcms=new TGraphAsymmErrors(4);
   grafcms->SetPoint(0,cmsx[0],cmsy[0]);
   grafcms->SetPointError(0,0,0,0,1-cmsy[0]);
   grafcms->SetPoint(1,cmsx[1],cmsy[1]);
   grafcms->SetPointError(1,0,0,0,1-cmsy[1]);
   grafcms->SetPoint(2,cmsx[2],cmsy[2]);
   grafcms->SetPointError(2,0,0,0,1-cmsy[2]);
   grafcms->SetPoint(3,cmsx[3],cmsy[3]);
   grafcms->SetPointError(3,0,0,0,1-cmsy[3]);
grafcms->SetFillColor(kTeal-9);
//   grafcms->SetFillColor(kCyan-10);
//   grafcms->SetFillColor(kCyan-8);
//  grafcms->SetFillColor(kMagenta-10);
//   grafcms->SetFillStyle(3005);
   grafcms->SetLineColor(1);
   grafcms->GetXaxis()->SetRangeUser(0.00005, 1);
   grafcms->GetYaxis()->SetRangeUser(0.00001, 1);
   grafcms->GetYaxis()->SetTitle("BR(t #rightarrow #gamma u)");
   grafcms->GetXaxis()->SetTitle("#kappa_{tu#gamma}");
   grafcms->GetXaxis()->SetLabelSize(0.03);
   grafcms->GetYaxis()->SetLabelSize(0.03);
   grafcms->GetXaxis()->SetTitleFont(62);
   grafcms->GetYaxis()->SetTitleFont(62);
   grafcms->GetXaxis()->SetTitleSize(0.065);
   grafcms->GetYaxis()->SetTitleSize(0.045);
   grafcms->GetXaxis()->SetLabelOffset(0.0001);
   grafcms->GetXaxis()->SetTitleOffset(0.60);


TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
pad1->SetFillStyle(0);
pad1->SetFrameFillStyle(0);
pad1->SetLeftMargin(0.12);
pad1->SetBottomMargin(0.12);
pad1->SetLogy(1);
pad1->SetLogx(1);
pad1->Draw();
pad1->cd();

   grafphotoncms->Draw("AL3");
  grafphotoncms2->Draw("L3SAME");
   grafcms->Draw("L3SAME");
 //  grafphotoncms2->Draw("L3SAME");
  grafphotoncms->Draw("L3SAME");

//   grafphotoncms->Draw("L3SAME");
   cdf->Draw("same");
   dzero->Draw("same");
   cms->Draw("same");
   cmsc->Draw("same");
   cmsz->Draw("same");
fhera->Draw("cont3SAME");
   hera->Draw("same");
   hera2->Draw("same");
   h1->Draw("same");
f2->Draw("cont3SAME");
   lep->Draw("same");
   lep2->Draw("same");
   atlas->Draw("same");  



  TPaveText *pt = new TPaveText(0.1,0.95,0.4,0.95, "NDC"); // NDC sets coords
    pt->SetLineColor(10);                                              // relative to pad dimensions
    pt->SetFillColor(10); // text is black on white
    pt->SetTextSize(0.045);
    pt->SetTextAlign(12);
    pt->AddText("CMS Preliminary, 19.1 fb^{-1}, #sqrt{s} = 8 TeV");
    pt->SetShadowColor(10);
 //   pt->Draw("same");


   cdf->GetHistogram()->Draw("AXISSAMEY+");
   cdf->GetHistogram()->Draw("AXISSAMEX+");
   grafphotoncms->GetHistogram()->Draw("AXISSAME");


 TLatex latex;
   latex.SetTextSize(0.035);
   latex.SetTextAlign(13);  //align at top
   latex.DrawLatex(0.000014,cdfy[0]+cdfy[0]/1.5,"#color[12]{CDF}");
   latex.DrawLatex(0.000014,dzeroy[0]-dzeroy[0]/10,"#color[3]{D0}");
   latex.DrawLatex(0.000014,lepy[0]+lepy[0]/1.5,"#color[4]{DELPHI}");
//    latex.DrawLatex(0.000014,heray[0]+heray[0]/1.5,"#color[6]{ZEUS}");
   latex.DrawLatex(0.000014,atlasy[0]+atlasy[0]/1.5,"#color[46]{ATLAS}");
   latex.DrawLatex(0.000014,cmsy[0]+cmsy[0]/1.5,"#color[1]{CMS}");

 latex.DrawLatex(mycmsx[0]+mycmsx[0]/8,0.000035,"#splitline{#color[2]{}}{#color[2]{(q=u)}}");
 latex.DrawLatex(mycmscx[0]-mycmscx[0]/1.85,0.000035,"#splitline{#color[2]{}}{#color[2]{(q=c)}}");
 latex.DrawLatex(mycmsx[0]+mycmsx[0]/2,0.0001,"#splitline{#color[2]{     CMS}}{#color[2]{Preliminary}}");


 latex.DrawLatex(0.000014,heray[0]+heray[0]/1.5,"#color[6]{ZEUS (q=u)}");
 latex.DrawLatex(h1x[0]+h1x[0]/8,0.00005,"#splitline{#color[28]{H1}}{#color[28]{(q=u)}}");
 //  latex.DrawLatex(0.03,cmsy[0]+cmsy[0]/1.5,"#color[46]{ #scale[10.2]{#rightarrow}CMS}");


   latex.DrawLatex(0.03,0.55,"#splitline{#color[1]{          95% C.L} }{#color[1]{EXCLUDED REGION}}");

TArrow ar1(0.1,0.1,1,0.1);
   ar1.Draw("same");
pad1->Draw();
//  TText *th1 = new TText(0.005,cdfy[0]+0.005,"#color[2]{Left adjusted}" );
//   th1->SetTextAlign(5); 
//   th1->SetTextSize(0.05);
//   th1->Draw();
 
//   grafexp1sigma->Draw("L3same");
//   cdf->Draw("AL");
//   Observed->Draw("L");

 
//   TGaxis *axis1 = new TGaxis(-8,-0.5,8,-0.5,0.0001,0.5,510,"-");
//   axis1->SetName("axis1");
//   axis1->Draw("same");
//   c1->Update();
/*
    TLegend *leg1 = new TLegend(0.2, 0.6, 0.35, 0.8);
    leg1->SetTextSize(0.03);
    leg1->SetBorderSize(0);
    leg1->SetLineColor(0);
    leg1->SetLineWidth(0);
    leg1->SetFillColor(kWhite);
    leg1->AddEntry(cdf, "Predicted", "L");
//    leg1->AddEntry(Observed, "95% CL Observed Limit", "L");
//    leg1->AddEntry(Expected, "95% CL Expected Limit", "L");
//    leg1->AddEntry(grafexp1sigma, "#pm1#sigma Exp.Limit", "F");
//    leg1->AddEntry(grafexp2sigma, "#pm2#sigma Exp.Limit", "F");

    leg1->Draw();

    TLine *line1 = new TLine(5, 1, 40, 1);
    line1->SetLineColor(2);
    line1->SetLineWidth(2);
//    line1->Draw("same");
  
    TPaveText *pt = new TPaveText(0.15,0.80,0.4,0.87, "NDC"); // NDC sets coords
    pt->SetLineColor(10);                                              // relative to pad dimensions
    pt->SetFillColor(10); // text is black on white
    pt->SetTextSize(0.03);
    pt->SetTextAlign(12);
    pt->AddText("CMS Preliminary #sqrt{s} = 8 TeV, #int L dt= 19.145 fb^{-1}");
    pt->SetShadowColor(10);
    pt->Draw("same");

*/
}

