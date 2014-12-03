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
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TPaveLabel.h"
#include "TLatex.h"
#include "TMath.h"

void limitINFO(){

 TF1 *fa1 = new TF1("fa","0.446*(4/9)*x*x",0.00001,0.6);
 TF1 *fa2 = new TF1("fa","0.479*(4/9)*x*x",0.000001,0.5);  

 const Int_t n = 25;
 Double_t x1[n], x2[n], x3[n], y1[n], y2[n], y3[n];

 for (Int_t i=1;i<n+1;i++) {
     x1[i]  = i*0.02;

     y1[i] = 0.416*(4/9)*x1[i]*x1[i];
     y2[i] = 0.479*(4/9)*x1[i]*x1[i];
     y3[i] = (y1[i]+y2[i])/2;

     x2[i] = abs(y1[i]-y3[i]);
     x3[i] = 0.0010;
   }

 TGraphErrors* error = new TGraphErrors(25, x1, y3, x3, x2);
   error->SetFillColor(4);
   error->SetFillStyle(3010);




TCanvas *c1 = new TCanvas("c1","upper limit results ",200,10,700,500);
c1->SetLogy(1);
//c1->SetLogx(1);
const int N=3;
double cdfx[4]={0,0.4005,0.40051,1};
double cdfy[4]={0.032,0.032,0,0};

//double cdfx[7]={0,0,1,1,0.45,0.45,0};
//double cdfy[7]={0.032,0.1,0.1,0,0,0.032,0.032};
double lepx[4]={0,0.486,0.486001,1};
double lepy[4]={0.0468,0.0468,0,0};
double herax[4]={0,0.174,0.174001,1};
double heray[4]={0.00597,0.00597,0,0};
double h1x[4]={0,0.18,0.18001,1};
double h1y[4]={0.0064,0.0064,0,0};
double cmsx[4]={0,0.0279,0.027901,1};
double cmsy[4]={0.000161,0.000161,0,0};



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

   TGraph *cms= new TGraph(4,cmsx,cmsy);
   cms->SetLineColor(2);
   cms->SetLineWidth(3);
   cms->SetMarkerColor(1);
   cms->SetMarkerStyle(20);
//   Observed->SetLineStyle(2);
   cms->GetYaxis()->SetTitle("Cross Section [pb]");
   cms->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   cms->SetTitle("with ");
   cms->SetFillColor(45);

   TGraph *hera= new TGraph(4,herax,heray);
   hera->SetLineColor(1);
   hera->SetLineWidth(3);
   hera->SetMarkerColor(1);
   hera->SetMarkerStyle(20);
//   hera->SetLineStyle(2);
//   hera->SetLineStyle(2);
   hera->GetYaxis()->SetTitle("Cross Section [pb]");
   hera->GetXaxis()->SetTitle("#Lambda_{T} [GeV]");
   hera->SetTitle("with ");
   hera->SetFillColor(45);

   TGraph *h1= new TGraph(4,h1x,h1y);
   h1->SetLineColor(46);
   h1->SetLineWidth(3);
   h1->SetMarkerColor(1);
   h1->SetMarkerStyle(20);
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

   TGraphAsymmErrors *grafhera=new TGraphAsymmErrors(4);
   grafhera->SetPoint(0,herax[0],heray[0]);
   grafhera->SetPointError(0,0,0,0,0.1-heray[0]);
   grafhera->SetPoint(1,herax[1],heray[1]);
   grafhera->SetPointError(1,0,0,0,0.1-heray[1]);
   grafhera->SetPoint(2,herax[2],heray[2]);
   grafhera->SetPointError(2,0,0,0,0.1-heray[2]);
   grafhera->SetPoint(3,herax[3],heray[3]);
   grafhera->SetPointError(3,0,0,0,0.1-heray[3]);
//grafhera->SetFillColor(kYellow+3);
//grafhera->SetFillColor(kRed-8);
grafhera->SetFillColor(kYellow-8);
 //  grafcdf->SetFillStyle(3010);
   grafhera->GetXaxis()->SetRangeUser(0, 1);
   grafhera->GetYaxis()->SetRangeUser(0, 0.1);


   TGraphAsymmErrors *grafcms=new TGraphAsymmErrors(4);
   grafcms->SetPoint(0,cmsx[0],cmsy[0]);
   grafcms->SetPointError(0,0,0,0,0.1-cmsy[0]);
   grafcms->SetPoint(1,cmsx[1],cmsy[1]);
   grafcms->SetPointError(1,0,0,0,0.1-cmsy[1]);
   grafcms->SetPoint(2,cmsx[2],cmsy[2]);
   grafcms->SetPointError(2,0,0,0,0.1-cmsy[2]);
   grafcms->SetPoint(3,cmsx[3],cmsy[3]);
   grafcms->SetPointError(3,0,0,0,0.1-cmsy[3]);
//grafcdf->SetFillColor(kYellow);
//   grafcms->SetFillColor(kViolet-1);
   grafcms->SetFillColor(kYellow-4);
//   grafcms->SetFillStyle(3004);
   grafcms->SetLineColor(4);
   grafcms->GetXaxis()->SetRangeUser(0, 0.5);
   grafcms->GetYaxis()->SetRangeUser(0, 0.1);
   grafcms->GetYaxis()->SetTitle("BR(t #rightarrow #gamma u)");
   grafcms->GetXaxis()->SetTitle("#kappa_{tu#gamma}");
   grafcms->GetXaxis()->SetLabelSize(0.044);
   grafcms->GetYaxis()->SetLabelSize(0.044);
   grafcms->GetXaxis()->SetTitleFont(62);
   grafcms->GetYaxis()->SetTitleFont(62);
   grafcms->GetXaxis()->SetTitleSize(0.075);
   grafcms->GetYaxis()->SetTitleSize(0.055);
   grafcms->GetXaxis()->SetTitleOffset(0.65);
   grafcms->GetYaxis()->SetTitleOffset(0.80);

TPad *pad1 = new TPad("pad1","pad1",0,0,1,1);
pad1->SetFillStyle(0);
pad1->SetFrameFillStyle(0);
pad1->SetLeftMargin(0.12);
pad1->SetBottomMargin(0.12);
pad1->SetLogy(1);
//pad1->SetLogx(1);
pad1->Draw();
pad1->cd();

   grafcms->Draw("AL3");
   grafhera->Draw("L3SAME");
   fa1->SetLineColor(6);
   fa1->SetLineWidth(3);
   fa1->SetMarkerColor(6);
   fa1->SetLineStyle(9);
  // fa1->SetLineWidth(3504);
//fa1->SetFillStyle(3005);
   fa2->SetLineWidth(3);
   fa2->SetLineColor(0);
   fa2->SetMarkerColor(0);
   fa2->SetLineStyle(9);
 //  fa2->Draw("ACsame");

   fa1->Draw("same");
   cdf->Draw("same");
   cms->Draw("same");
   hera->Draw("same");
   h1->Draw("same");
   lep->Draw("same");

//   error->Draw("A3SAME");


//fa2->SetFillStyle(3005);




  TPaveText *pt = new TPaveText(0.1,0.95,0.4,0.95, "NDC"); // NDC sets coords
    pt->SetLineColor(10);                                              // relative to pad dimensions
    pt->SetFillColor(10); // text is black on white
    pt->SetTextSize(0.045);
    pt->SetTextAlign(12);
    pt->AddText("CMS Preliminary, #sqrt{s} = 8 TeV");
    pt->SetShadowColor(10);
 //   pt->Draw("same");


   cdf->GetHistogram()->Draw("AXISSAMEY+");
   cdf->GetHistogram()->Draw("AXISSAMEX+");
   cdf->GetHistogram()->Draw("AXISSAME");


 TLatex latex;
   latex.SetTextSize(0.035);
   latex.SetTextAlign(13);  //align at top
   latex.DrawLatex(0.02,cdfy[0]+cdfy[0]/3.5,"#color[12]{CDF}");
   latex.DrawLatex(0.02,lepy[0]+lepy[0]/3.5,"#color[4]{DELPHI}");
    latex.DrawLatex(0.02,heray[0]+heray[0]/2,"#color[1]{ZEUS,} #color[46]{H1}");
   latex.DrawLatex(0.033,cmsy[0]-cmsy[0]/50,"#color[2]{CMS Preliminary}");
   latex.DrawLatex(0.25,0.075,"#color[1]{95% C.L EXCLUDED REGION}");
   latex.DrawLatex(0.24,heray[0]+heray[0],"#color[6]{#splitline{         Theory}{(#Lambda = m_{top} = 175 GeV)}}");
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

