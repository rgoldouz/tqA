#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
void plotsystematicsanalyzer() {
 	  	  	
std::vector<TFile*> files;
TString fn="systematichistos.root";
files.push_back(new TFile(fn));

//open files

/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

// load histograms
std::vector<TH1F*> JERhists;
std::vector<TH1F*> JEShists;
std::vector<TH1F*> PUhists;
std::vector<TH1F*> TRhists;
std::vector<TH1F*> MUONhists;
std::vector<TH1F*> PHhists;
std::vector<TH1F*> BTAGhists;
std::vector<TH1F*> MISTAGhists;
std::vector<TH1F*> nominalhists;
JERhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__JER__plus"));
JERhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__JER__minus"));
JEShists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__JES__plus"));
JEShists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__JES__minus"));
PUhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__PU__plus"));
PUhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__PU__minus"));
TRhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__TRIG__plus"));
TRhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__TRIG__minus"));
MUONhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__MUON__plus"));
MUONhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__MUON__minus"));
PHhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__PHOTON__plus"));
PHhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__PHOTON__minus"));
BTAGhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__BTAG__plus"));
BTAGhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__BTAG__minus"));
MISTAGhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__MISSTAG__plus"));
MISTAGhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist__MISSTAG__minus"));
nominalhists.push_back((TH1F*)files[0]->Get("BDT__zgammahist"));
nominalhists.push_back((TH1F*)files[0]->Get("BDT__wjet"));




JEShists[0]->Add(nominalhists[1]);
JEShists[1]->Add(nominalhists[1]);
JERhists[0]->Add(nominalhists[1]);
JERhists[1]->Add(nominalhists[1]);
PUhists[0]->Add(nominalhists[1]);
PUhists[1]->Add(nominalhists[1]);
TRhists[0]->Add(nominalhists[1]);
TRhists[1]->Add(nominalhists[1]);
MUONhists[0]->Add(nominalhists[1]);
MUONhists[1]->Add(nominalhists[1]);
PHhists[0]->Add(nominalhists[1]);
PHhists[1]->Add(nominalhists[1]);
BTAGhists[0]->Add(nominalhists[1]);
BTAGhists[1]->Add(nominalhists[1]);
MISTAGhists[0]->Add(nominalhists[1]);
MISTAGhists[1]->Add(nominalhists[1]);
nominalhists[0]->Add(nominalhists[1]);


// setup the canvas and draw the histograms
TCanvas* convjes = new TCanvas("jes", "jes" , 600, 600);
convjes->cd(0);
//conv->SetLogy(1);
//JEShists[0]->SetMaximum(1.3*JEShists[0]->GetMaximum());
JEShists[0]->SetLineColor(kRed);
JEShists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
JEShists[1]->SetLineColor(kBlue);
JEShists[1]->Draw("same");
TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);

  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( JEShists[0] , "JES UP"                           , "F"); 
  leg->AddEntry( nominalhists[0], "Nominal"                           , "F");
  leg->AddEntry( JEShists[1], "JES Down"              , "F");
  leg->Draw("same");
////////////////////////////////////////////////////////////////////////////
TCanvas* convjer = new TCanvas("jer", "jer" , 600, 600);
convjer->cd(0);
JERhists[0]->SetLineColor(kRed);
JERhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
JERhists[1]->SetLineColor(kBlue);
JERhists[1]->Draw("same");
TLegend* legjer = new TLegend(0.35,0.6,0.85,0.90);

  legjer->SetFillStyle ( 0);
  legjer->SetFillColor ( 0);
  legjer->SetBorderSize( 0);
  legjer->AddEntry( JEShists[0] , "JER UP"                           , "F");
  legjer->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legjer->AddEntry( JEShists[1], "JER Down"              , "F");
  legjer->Draw("same");
/////////////////////////////////////////////////////////////////////////////
TCanvas* convpui = new TCanvas("puil up", "pile up" , 600, 600);
convpui->cd(0);
PUhists[0]->SetLineColor(kRed);
PUhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
PUhists[1]->SetLineColor(kBlue);
PUhists[1]->Draw("same");
TLegend* legpui = new TLegend(0.35,0.6,0.85,0.90);

  legpui->SetFillStyle ( 0);
  legpui->SetFillColor ( 0);
  legpui->SetBorderSize( 0);
  legpui->AddEntry( JEShists[0] , "PU UP"                          , "F");
  legpui->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legpui->AddEntry( JEShists[1], "PU Down"              , "F");
  legpui->Draw("same");
//////////////////////////////////////////////////////////////////

TCanvas* convph = new TCanvas("photon sf", "photon sf" , 600, 600);
convpui->cd(0);
PHhists[0]->SetLineColor(kRed);
PHhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
PHhists[1]->SetLineColor(kBlue);
PHhists[1]->Draw("same");
TLegend* legpho= new TLegend(0.35,0.6,0.85,0.90);

  legpho->SetFillStyle ( 0);
  legpho->SetFillColor ( 0);
  legpho->SetBorderSize( 0);
  legpho->AddEntry( JEShists[0] , "PHOTON_SF UP"                  , "F");
  legpho->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legpho->AddEntry( JEShists[1], "PHOTON_SF  Down"              , "F");
  legpho->Draw("same");
//////////////////////////////////////////////////////////////////
TCanvas* convmu = new TCanvas("muon sf", "muon sf" , 600, 600);
convmu->cd(0);
MUONhists[0]->SetLineColor(kRed);
MUONhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
MUONhists[1]->SetLineColor(kBlue);
MUONhists[1]->Draw("same");
TLegend* legmuo =new TLegend(0.35,0.6,0.85,0.90);

  legmuo->SetFillStyle ( 0);
  legmuo->SetFillColor ( 0);
  legmuo->SetBorderSize( 0);
  legmuo->AddEntry( JEShists[0] , "MUON_SF UP"                  , "F");
  legmuo->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legmuo->AddEntry( JEShists[1], "MUON_SF  Down"              , "F");
  legmuo->Draw("same");
//////////////////////////////////////////////////////////////////
TCanvas* convb = new TCanvas("btag sf", "btag sf" , 600, 600);
convb->cd(0);
BTAGhists[0]->SetLineColor(kRed);
BTAGhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
BTAGhists[1]->SetLineColor(kBlue);
BTAGhists[1]->Draw("same");
TLegend* legbta =new TLegend(0.35,0.6,0.85,0.90);

  legbta->SetFillStyle ( 0);
  legbta->SetFillColor ( 0);
  legbta->SetBorderSize( 0);
  legbta->AddEntry( BTAGhists[0] , "BTAG_SF UP"                  , "F");
  legbta->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legbta->AddEntry( BTAGhists[1], "BTAG_SF  Down"              , "F");
  legbta->Draw("same");
//////////////////////////////////////////////////////////////////
TCanvas* convmistag = new TCanvas("mistag sf", "mistag sf" , 600, 600);
convmistag->cd(0);
MISTAGhists[0]->SetLineColor(kRed);
MISTAGhists[0]->Draw();
nominalhists[0]->SetLineColor(kGreen);
nominalhists[0]->Draw("same");
MISTAGhists[1]->SetLineColor(kBlue);
MISTAGhists[1]->Draw("same");
TLegend* legmis =new TLegend(0.35,0.6,0.85,0.90);

  legmis->SetFillStyle (0);
  legmis->SetFillColor (0);
  legmis->SetBorderSize(0);
  legmis->AddEntry( MISTAGhists[0] , "MISTAG_SF UP"                  , "F");
  legmis->AddEntry( nominalhists[0], "Nominal"                           , "F");
  legmis->AddEntry( MISTAGhists[1], "MISTAG_SF  Down"              , "F");
  legmis->Draw("same");
}
