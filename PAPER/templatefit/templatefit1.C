#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"


using namespace RooFit;

void templatefit1() {
std::vector<string> samples_;
float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02};
samples_.push_back("WPHJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("PHJET200400.root");
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
//samples_.push_back("SIGNAL.root");
//samples_.push_back("REALDATA.root");

std::vector<string> samplesreverse_;
samplesreverse_.push_back("etarev/WPHJET.root");
samplesreverse_.push_back("etarev/ZJET.root");
samplesreverse_.push_back("etarev/PHJET200400.root");
samplesreverse_.push_back("etarev/T-W-CH.root");
samplesreverse_.push_back("etarev/TBAR-W-CH.root");
samplesreverse_.push_back("etarev/T-S-CH.root");
samplesreverse_.push_back("etarev/TBAR-S-CH.root");
samplesreverse_.push_back("etarev/T-T-CH.root");
samplesreverse_.push_back("etarev/TBAR-T-CH.root");
samplesreverse_.push_back("etarev/TTBAR1.root");
samplesreverse_.push_back("etarev/TTBAR2.root");
samplesreverse_.push_back("etarev/TTBAR3.root");
samplesreverse_.push_back("etarev/TTG.root");
samplesreverse_.push_back("etarev/WWG.root");
samplesreverse_.push_back("etarev/WW.root");
samplesreverse_.push_back("etarev/WZ.root");
samplesreverse_.push_back("etarev/ZZ.root");
samplesreverse_.push_back("etarev/ZGAMMA.root");


//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}


std::vector<TFile*> reversefiles;
for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
reversefiles.push_back(new TFile(samplesreverse_[idx].c_str()));
}

TFile* data = new TFile("REALDATA.root");
TFile* datarev=new TFile("etarev/REALDATA.root");
/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<string> variables2_;

  variables2_.push_back("photon_Pt");
//    variables2_.push_back("photon_Eta");
//    variables2_.push_back("muon_Pt");
//    variables2_.push_back("muon_Eta");
//  variables2_.push_back("Jet_Pt");
//  variables2_.push_back("Jet_Eta");
//  variables2_.push_back("Jet_Multiplicity");
//  variables2_.push_back("topmass");
//  variables2_.push_back("WTmass");
//  variables2_.push_back("Cosmuonjet");
//  variables2_.push_back("Cosmuonphoton");
//  variables2_.push_back("Cosphotonjet");
 // variables2_.push_back("Deltaphiphotonjet");
//  variables2_.push_back("Deltaphiphotonmuon");
//  variables2_.push_back("Deltaphimuonjet");
//  variables2_.push_back("DeltaRphotonmuon");
//  variables2_.push_back("DeltaRphotonjet");
//  variables2_.push_back("DeltaRmuonjet");
//  variables2_.push_back("HT");
//  variables2_.push_back("Photonmuonmass");
//  variables2_.push_back("Costopphoton");
//  variables2_.push_back("sigmaetaetaendcap");
//  variables2_.push_back("sigmaetaetabarrel");


// load histograms
std::vector<TH1F*> BGhists;
for(unsigned int idx=0; idx<files.size(); ++idx){
BGhists.push_back((TH1F*)files[idx]->Get((std::string("analyzestep2/").append(variables2_[0]).c_str())));
}

std::vector<TH1F*> revBGhists;
for(unsigned int idx=0; idx<reversefiles.size(); ++idx){
revBGhists.push_back((TH1F*)reversefiles[idx]->Get((std::string("analyzestep2/").append(variables2_[0]).c_str())));
}

TH1F* DATAhists;
DATAhists = (TH1F*)data->Get((std::string("analyzestep2/").append(variables2_[0]).c_str()));
TH1F* SIGhists;
TH1F* revDATAhists;
revDATAhists = (TH1F*)datarev->Get((std::string("analyzestep2/").append(variables2_[0]).c_str()));

float lumi = 5.319;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
BGhists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
revBGhists[idx]->Scale(lumi*scales[idx]);
}

for(unsigned int idx=1; idx<samples_.size(); ++idx){
   BGhists[idx]->Add(BGhists[idx-1]);
 }
for(unsigned int idx=1; idx<samplesreverse_.size(); ++idx){
revBGhists[idx]->Add(revBGhists[idx-1]);
 }
SIGhists=revDATAhists;
SIGhists->Add(revBGhists[samples_.size()-1],-1);

//TCanvas* conv = new TCanvas("BGBG", "BG" , 600, 600);
//conv->cd(0);
//DATAhists->Draw();
//revBGhists[samples_.size()-1]->Draw("same");

//TCanvas* nv = new TCanvas("data", "data" , 600, 600);
//nv->cd(0);
//SIGhists->Draw();


int norm=1;
double scalesig =norm/SIGhists->Integral();
SIGhists->Scale(scalesig);

//double scaledata =norm/DATAhists->Integral();
//DATAhists->Scale(scaledata);

double scalebg =norm/BGhists[samples_.size()-1]->Integral();
BGhists[samples_.size()-1]->Scale(scalebg);

//TCanvas* conv = new TCanvas("signal", "signal" , 600, 600);
//conv->cd(0);
//SIGhists->Draw();

//TCanvas* nv = new TCanvas("data", "data" , 600, 600);
//nv->cd(0);
//DATAhists->Draw();

//TCanvas* con = new TCanvas("BG", "BG" , 600, 600);
//con->cd(0);
//BGhists[samples_.size()-1]->Draw();


//conv->cd(0);
////Ghists->Draw();

    RooRealVar x("x", "x", 50.,500.);
//    RooRealVar x("x", "x", -3,3);
//    RooRealVar x("x", "x", 26.,500.);
//      RooRealVar x("x", "x", 25.,150.);
//      RooRealVar x("x", "x", 0.,400.);

    RooDataHist* template1 = new RooDataHist("template1", "template1", RooArgList(x), SIGhists);
    RooHistPdf modelTemplate1("modelTemplate1", "modelTemplate1", RooArgSet(x), *template1);

    RooDataHist* template2 = new RooDataHist("template2", "template2", RooArgList(x), BGhists[samples_.size()-1]);
    RooHistPdf modelTemplate2("modelTemplate2", "modelTemplate2", RooArgSet(x), *template2);

    RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), DATAhists);

    std::cout << " executing fit..." << std::endl;

    RooRealVar coeffTemplate1("coeffTemplate1", "coeffTemplate1", 0, 2000);
    RooRealVar coeffTemplate2("coeffTemplate2", "coeffTemplate2", 0, 2000);

    RooAddPdf sumTemplate("sumTemplate", "coeff*modelTemplate1 + (1 - coeff)*modelTemplate2",RooArgList(modelTemplate1, modelTemplate2),RooArgList(coeffTemplate1,coeffTemplate2));
    sumTemplate.fitTo(*sumData);

std::cout << "--> coeffTemplate1 = " << coeffTemplate1.getVal() << std::endl;
std::cout << "--> coeffTemplate2 = " << coeffTemplate2.getVal() << std::endl;

RooPlot* frame = x.frame() ;
sumData->plotOn(frame,LineStyle(kDashed)) ;
sumTemplate.plotOn(frame) ;

//modelTemplate2.plotOn(frame) ;

frame->Draw() ;

}
