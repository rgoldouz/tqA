#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

//using namespace RooFit;

void WjetShapeReweight(){
std::vector<string> samples_;
float scales[] = {0.0961,0.0978,1.491,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,0.020,0.0017,0.0055,0.0032,0.00084,0.02};
/*
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
*/
samples_.push_back("WJET.root");

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

std::vector<string> datasamples_;
std::vector<string> datasamplesreverse_;
datasamples_.push_back("REALDATA1.root");
datasamples_.push_back("REALDATA2.root");
datasamples_.push_back("REALDATA3.root");
datasamplesreverse_.push_back("etarev/REALDATA1.root");
datasamplesreverse_.push_back("etarev/REALDATA2.root");
datasamplesreverse_.push_back("etarev/REALDATA3.root");


//open files
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx].c_str()));
}


std::vector<TFile*> reversefiles;
for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
reversefiles.push_back(new TFile(samplesreverse_[idx].c_str()));
}

 std::vector<TFile*> data;
for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
data.push_back(new TFile(datasamples_[idx].c_str()));
}
 
 std::vector<TFile*> datarev;
for(unsigned int idx=0; idx<datasamplesreverse_.size(); ++idx){
datarev.push_back(new TFile(datasamplesreverse_[idx].c_str()));
}
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

std::vector<TH1F*>  DATAhists;
for(unsigned int idx=0; idx<data.size(); ++idx){
DATAhists.push_back((TH1F*)data[idx]->Get((std::string("analyzestep2/").append(variables2_[0]).c_str())));
}
TH1F* SIGhists;
std::vector<TH1F*> revDATAhists;
for(unsigned int idx=0; idx<datarev.size(); ++idx){
 revDATAhists.push_back((TH1F*)datarev[idx]->Get((std::string("analyzestep2/").append(variables2_[0]).c_str())));
 }
float lumi = 1;
/*
for(unsigned int idx=0; idx<samples_.size(); ++idx){
BGhists[idx]->Scale(lumi*scales[idx]);
}
*/
for(unsigned int idx=0; idx<samplesreverse_.size(); ++idx){
revBGhists[idx]->Scale(lumi*scales[idx]);
}

// TH1F* wphjethists = BGhists[0];
/*
for(unsigned int idx=2; idx<samples_.size(); ++idx){
   BGhists[idx]->Add(BGhists[idx-1]);
 }
*/
for(unsigned int idx=1; idx<samplesreverse_.size(); ++idx){
revBGhists[idx]->Add(revBGhists[idx-1]);
 }
for(unsigned int idx=1; idx<datasamplesreverse_.size(); ++idx){
revDATAhists[idx]->Add(revDATAhists[idx-1]);
 }
for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
DATAhists[idx]->Add(DATAhists[idx-1]);
 }
SIGhists=revDATAhists[2];
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

double scalebg =norm/BGhists[0]->Integral();
BGhists[0]->Scale(scalebg);
SIGhists->SetFillColor(kMagenta+2);
SIGhists->Draw();
BGhists[0]->Draw("same");
int bin=1;
double totalsig=0;
double totalbg=0;

TH1F *newwjet = new TH1F("newwjet","newwjet", 25,50., 500.);
for ( bin=0; bin!=22; ++bin){
cout<<"nbin"<<bin<<endl;
//cout<<"data driven pdf"<< SIGhists->GetBinContent(bin)<<endl;
//cout<<"MC pdf         "<<BGhists[0]->GetBinContent(bin)<<endl;
cout<<"********************************************************************"<<endl;
//newwjet->SetBinContent(bin,(SIGhists->GetBinContent(bin)/BGhists[0]->GetBinContent(bin))*BGhists[0]->GetBinContent(bin));
if (bin<3) {cout<<SIGhists->GetBinContent(bin)/BGhists[0]->GetBinContent(bin);}
else {totalsig=totalsig+SIGhists->GetBinContent(bin);
totalbg=totalbg+BGhists[0]->GetBinContent(bin);}
}
cout<<totalsig/totalbg<<endl;
/*
cout<<SIGhists->Integral()<<endl;
cout<<BGhists[0]->Integral()<<endl;
cout<<newwjet->Integral()<<endl;
newwjet->Scale(1/newwjet->Integral());
cout<<newwjet->Integral()<<endl;

for ( bin=1; bin!=7; ++bin){
cout<<newwjet->GetBinContent(bin)/BGhists[0]->GetBinContent(bin)<<endl;
}
*/
//double scalewphjet =norm/wphjethists->Integral();
//wphjethists->Scale(scalewphjet);
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
/*
    RooRealVar x("x", "x", 50.,500.);
//    RooRealVar x("x", "x", -3,3);
//    RooRealVar x("x", "x", 26.,500.);
//      RooRealVar x("x", "x", 25.,150.);
//      RooRealVar x("x", "x", 0.,400.);

    RooDataHist* template1 = new RooDataHist("template1", "template1", RooArgList(x), SIGhists);
    RooHistPdf modelTemplate1("modelTemplate1", "modelTemplate1", RooArgSet(x), *template1);

    RooDataHist* template2 = new RooDataHist("template2", "template2", RooArgList(x), BGhists[samples_.size()-1]);
    RooHistPdf modelTemplate2("modelTemplate2", "modelTemplate2", RooArgSet(x), *template2);


    RooDataHist* template3 = new RooDataHist("template3", "template3", RooArgList(x), wphjethists);
    RooHistPdf modelTemplate3("modelTemplate3", "modelTemplate3", RooArgSet(x), *template3);

    RooDataHist* sumData = new RooDataHist("sumData", "sumData", RooArgList(x), DATAhists[2]);

    std::cout << " executing fit..." << std::endl;

    RooRealVar coeffTemplate1("coeffTemplate1", "coeffTemplate1", 0, 5000);
    RooRealVar coeffTemplate2("coeffTemplate2", "coeffTemplate2",1174.33, 0, 5000);
    RooRealVar coeffTemplate3("coeffTemplate3", "coeffTemplate3", 0, 5000);
    coeffTemplate2.setConstant(kTRUE);

    RooAddPdf sumTemplate("sumTemplate", "coeff*modelTemplate1 + (1 - coeff)*modelTemplate2",RooArgList(modelTemplate1, modelTemplate2,modelTemplate3),RooArgList(coeffTemplate1,coeffTemplate2,coeffTemplate3));
    sumTemplate.fitTo(*sumData);

std::cout << "--> coeffTemplate1 = " << coeffTemplate1.getVal() << std::endl;
std::cout << "--> coeffTemplate2 = " << coeffTemplate2.getVal() << std::endl;
std::cout << "--> coeffTemplate3 = " << coeffTemplate3.getVal() << std::endl;


RooPlot* frame = x.frame() ;
sumData->plotOn(frame,LineStyle(kDashed),LineWidth(3)) ;
 sumTemplate.plotOn(frame,LineWidth(3)) ;

//modelTemplate2.plotOn(frame) ;

frame->Draw() ;
            TH1D *  h = new TH1D("dummy_fit","dummy_fit",100,0,100);
	    //	    h->SetLineColor(kblue);
	    h->SetLineWidth(3);
	    //sumData->SetLineWidth(3);
TLegend* leg = new TLegend(0.35,0.6,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry(sumData , "CMS Data 2012(19.145 /fb)"               , "F");
  leg->AddEntry(h , "result of the fit" , "F");

  leg->Draw("same");
*/

}
