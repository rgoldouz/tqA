#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLine.h"
#include "THStack.h"

void pileup() {
 	  	  	
// list of valid histogram names


// list of input samples
std::vector<string> samples_;
std::vector<string> datasamples_;
std::vector<string> datasamplesreverse_;
float scales[] = {0.628,0.0978,34.01/19.145,6.133/19.145,1.04/19.145,0.32/19.145,0.02/19.145,0.002/19.145,1.17*0.0961,0.0253,0.0224,0.0145,0.0125,0.0160,0.0158,0.0341,0.0341,0.0341,
0.020,0.0017,0.0055,0.0032,0.00084,0.02,0.01139,0.01139,(0.049094905/19.145)};
double lumi=19.76836;

samples_.push_back("WJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("G_Pt_50to80.root");
samples_.push_back("G_Pt_80to120.root");
samples_.push_back("G_Pt_120to170.root");
samples_.push_back("G_Pt_170to300.root");
samples_.push_back("G_Pt_300to470.root");
samples_.push_back("G_Pt_470to800.root");
samples_.push_back("WPHJET.root");
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
samples_.push_back("SINGLE-TOP-PH.root");
samples_.push_back("SINGLE-ANTITOP-PH.root");

samples_.push_back("SIGNALtGu.root");

std::vector<string> samplesreverse_;
samplesreverse_.push_back("etarev/WJET.root");
samplesreverse_.push_back("etarev/ZJET.root");
samplesreverse_.push_back("etarev/G_Pt_50to80.root");
samplesreverse_.push_back("etarev/G_Pt_80to120.root");
samplesreverse_.push_back("etarev/G_Pt_120to170.root");
samplesreverse_.push_back("etarev/G_Pt_170to300.root");
samplesreverse_.push_back("etarev/G_Pt_300to470.root");
samplesreverse_.push_back("etarev/G_Pt_470to800.root");
samplesreverse_.push_back("etarev/WPHJET.root");
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
samplesreverse_.push_back("etarev/SINGLE-TOP-PH.root");
samplesreverse_.push_back("etarev/SINGLE-ANTITOP-PH.root");
samples_.push_back("SIGNALtGu.root");

datasamples_.push_back("REALDATAab.root");
datasamples_.push_back("REALDATAc.root");
datasamples_.push_back("REALDATAd.root");
datasamplesreverse_.push_back("etarev/REALDATA1.root");
datasamplesreverse_.push_back("etarev/REALDATA2.root");
datasamplesreverse_.push_back("etarev/REALDATA3.root");


/////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////

std::vector<TH1F*> hists_Nvertex;
std::vector<TH1F*> hists_NvertexnoPU;

std::vector<TH1F*> hists_photonpt;
std::vector<TH1F*> hists_muonpt;

std::vector<TH1F*> hists_photonptnoPU;
std::vector<TH1F*> hists_muonptnoPU;



//be carefull____________ to make code shorter, histograms of MC samples are from zero to sample size and after that 3 histo for real data and after that 3 files for data in controlregion for wjet sample.
 
for(unsigned int idx=0; idx<samples_.size(); ++idx){


hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(samples_[idx]).c_str(), std::string("Nvertex").append(samples_[idx]).c_str() ,40, -0.5, 39.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(samples_[idx]).c_str(), std::string("NvertexnoPU").append(samples_[idx]).c_str() ,40, -0.5, 39.5));
hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(samples_[idx]).c_str(), std::string("photon_pt").append(samples_[idx]).c_str() , 50,50., 120.));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(samples_[idx]).c_str(), std::string("muon_pt").append(samples_[idx]).c_str() , 30, 26, 140));
hists_photonptnoPU.push_back(new TH1F(  std::string("photon_ptnoPU").append(samples_[idx]).c_str(), std::string("photon_ptnoPU").append(samples_[idx]).c_str() , 50,50., 120.));
hists_muonptnoPU.push_back(new TH1F(  std::string("muon_ptnoPU").append(samples_[idx]).c_str(), std::string("muon_ptnoPU").append(samples_[idx]).c_str() , 30, 26, 140));

}

for(unsigned int idx=0; idx<datasamples_.size(); ++idx){

hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(datasamples_[idx]).c_str(), std::string("Nvertex").append(datasamples_[idx]).c_str() ,40, -0.5, 39.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(datasamples_[idx]).c_str(), std::string("NvertexnoPU").append(datasamples_[idx]).c_str() ,40, -0.5, 39.5));
hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(datasamples_[idx]).c_str(), std::string("photon_pt").append(datasamples_[idx]).c_str() , 50,50., 120.));
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamples_[idx]).c_str(), std::string("muon_pt").append(datasamples_[idx]).c_str() , 30, 26, 140));

hists_photonptnoPU.push_back(new TH1F(  std::string("photon_ptnoPU").append(datasamples_[idx]).c_str(), std::string("photon_ptnoPU").append(datasamples_[idx]).c_str() , 50,50., 120.));
hists_muonptnoPU.push_back(new TH1F(  std::string("muon_ptnoPU").append(datasamples_[idx]).c_str(), std::string("muon_ptnoPU").append(datasamples_[idx]).c_str() , 30, 26, 140));
}



for(unsigned int idx=0; idx<samples_.size(); ++idx){
cout<<samples_[idx]<<endl;
	TFile *input(0);
	TString fname =samples_[idx];
	input = TFile::Open( fname ); 

        std::vector<double> *myNvertex=0;
	std::vector<double> *mypileupSF=0;
	std::vector<double> *myweight=0;
	std::vector<double> *myNpileup=0;
std::vector<double> *myphotonSF=0;
std::vector<double> *mytriggerSF=0;
std::vector<double> *mymuonSF=0;
	std::vector<double> *myptphoton=0;
	std::vector<double> *myptmuon=0;

	TTree* theTree = (TTree*)input->Get("STEP1/FCNC");
//	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");

        theTree->SetBranchAddress( "Nvertex", &myNvertex);
	theTree->SetBranchAddress( "pileupSF", &mypileupSF);
	theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "NPileup", &myNpileup );
   theTree->SetBranchAddress( "photonSF", &myphotonSF);
   theTree->SetBranchAddress( "muonSF", &mymuonSF);
   theTree->SetBranchAddress( "triggerSF", &mytriggerSF);
//	theTree->SetBranchAddress("ptphoton", &myptphoton  );
//   	theTree->SetBranchAddress( "ptmuon", &myptmuon );


	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);
	double finalweight;
	finalweight=(*myweight)[0];


//        if (samples_[idx]=="WPHJET.root") finalweight=finalweight*(1.46604+0.00923692*(*myptphoton)[0]-1.48871e-06*(*myptphoton)[0]*(*myptphoton)[0]);
//        if (samples_[idx]=="ZGAMMA.root") finalweight=finalweight*(1.3292+0.000952237*(*myptphoton)[0]+2.00623e-05*(*myptphoton)[0]*(*myptphoton)[0]-1.41325e-07 *(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]+2.48614e-10*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]*(*myptphoton)[0]);


//	if (finalweight<100) {


	hists_Nvertex[idx] ->Fill( (*myNvertex)[0],finalweight*scales[idx]*lumi);
	hists_NvertexnoPU[idx] ->Fill( (*myNvertex)[0],(finalweight/(*mypileupSF)[0])*scales[idx]*lumi);

//	hists_photonpt[idx] ->Fill( (*myptphoton)[0],finalweight *scales[idx] );
//	hists_muonpt[idx] ->Fill((*myptmuon )[0] ,finalweight*scales[idx] );

//	hists_photonptnoPU[idx] ->Fill( (*myptphoton)[0],(finalweight/(*mypileupSF)[0]) *scales[idx] );
//	hists_muonptnoPU[idx] ->Fill((*myptmuon )[0] ,(finalweight/(*mypileupSF)[0])*scales[idx] );

//	hists_Nvertex[idx] ->Fill( (*myNvertex)[0], (*myphotonSF)[0]* (*mytriggerSF)[0]* (*mymuonSF)[0]* (*mypileupSF)[0]*scales[idx]);
//	hists_NvertexnoPU[idx] ->Fill( (*myNvertex)[0],(*myphotonSF)[0]* (*mytriggerSF)[0]* (*mymuonSF)[0]*scales[idx]);

//	}
	}

        delete myNvertex;
        delete mypileupSF;
	delete myNpileup;
}

for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
	TFile *input(0);
	TString fname =datasamples_[idx];
	input = TFile::Open( fname ); 

        std::vector<double> *myNvertex=0;
	std::vector<double> *myweight=0;
	std::vector<double> *myNpileup=0;
	std::vector<double> *myptphoton=0;
	std::vector<double> *myptmuon=0;

	TTree* theTree = (TTree*)input->Get("STEP1/FCNC");
//	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	cout<<theTree->GetEntries()<<endl;

        theTree->SetBranchAddress( "Nvertex", &myNvertex);
	theTree->SetBranchAddress( "weight", &myweight);
	theTree->SetBranchAddress( "NPileup", &myNpileup );
//	theTree->SetBranchAddress("ptphoton", &myptphoton  );
//   	theTree->SetBranchAddress( "ptmuon", &myptmuon );


	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);


	hists_Nvertex[idx+samples_.size()] ->Fill( (*myNvertex)[0]);
	hists_NvertexnoPU[idx+samples_.size()] ->Fill( (*myNvertex)[0]);
//	hists_photonpt[idx+samples_.size()] ->Fill( (*myptphoton)[0]);
//	hists_muonpt[idx+samples_.size()] ->Fill((*myptmuon )[0] );
//	hists_photonptnoPU[idx+samples_.size()] ->Fill( (*myptphoton)[0]);
//	hists_muonptnoPU[idx+samples_.size()] ->Fill((*myptmuon )[0] );

	}

        delete myNvertex;
	delete input;
	delete myNpileup;
}






//hists_Nvertex[0]->Scale(219.373/hists_Nvertex[0]->Integral());
//hists_NvertexnoPU[0]->Scale(219.373/hists_NvertexnoPU[0]->Integral());


//hists_Nvertex[8]->Scale(1112.2/hists_Nvertex[8]->Integral());
//hists_NvertexnoPU[8]->Scale(1112.2/hists_NvertexnoPU[8]->Integral());

for(unsigned int idx=1; idx<datasamples_.size(); ++idx){
        hists_Nvertex[idx+samples_.size()]->Add(hists_Nvertex[idx+samples_.size()-1]);
        hists_NvertexnoPU[idx+samples_.size()]->Add(hists_NvertexnoPU[idx+samples_.size()-1]);
	hists_photonptnoPU[idx+samples_.size()]->Add(hists_photonpt[idx+samples_.size()-1]); 
	hists_muonptnoPU[idx+samples_.size()]->Add(hists_muonpt[idx+samples_.size()-1]);
	hists_photonpt[idx+samples_.size()]->Add(hists_photonpt[idx+samples_.size()-1]); 
	hists_muonpt[idx+samples_.size()]->Add(hists_muonpt[idx+samples_.size()-1]);
}

/*
TCanvas *c1 = new TCanvas("c1","photonpt",900,700);
c1->cd(0);

THStack *hs1 = new THStack("hs11","photonpt");
hists_photonpt[0]->SetFillColor(kBlue-4);
hs1->Add(hists_photonpt[0]);
hists_photonpt[1]->SetFillColor(kOrange+7);
hs1->Add(hists_photonpt[1]);
hists_photonpt[3]->Add(hists_photonpt[2]);
hists_photonpt[4]->Add(hists_photonpt[3]);
hists_photonpt[5]->Add(hists_photonpt[4]);
hists_photonpt[6]->Add(hists_photonpt[5]);
hists_photonpt[7]->Add(hists_photonpt[6]);
hists_photonpt[7]->SetFillColor(19);
hs1->Add(hists_photonpt[7]);
hists_photonpt[8]->SetFillColor(kGreen+1);
hs1->Add(hists_photonpt[8]);
hists_photonpt[10]->Add(hists_photonpt[9]);
hists_photonpt[6+5]->Add(hists_photonpt[5+5]);
hists_photonpt[7+5]->Add(hists_photonpt[6+5]);
hists_photonpt[8+5]->Add(hists_photonpt[7+5]);
hists_photonpt[9+5]->Add(hists_photonpt[8+5]);
hists_photonpt[9+5]->SetFillColor(kAzure+10);
hs1->Add(hists_photonpt[9+5]);
hists_photonpt[11+5]->Add(hists_photonpt[10+5]);
hists_photonpt[12+5]->Add(hists_photonpt[11+5]);
hists_photonpt[12+5]->SetFillColor(6);
hs1->Add(hists_photonpt[12+5]);
hists_photonpt[13+5]->SetFillColor(32);
hs1->Add(hists_photonpt[13+5]);
hists_photonpt[14+5]->SetFillColor(kSpring-9);
hs1->Add(hists_photonpt[14+5]);
hists_photonpt[16+5]->Add(hists_photonpt[15+5]);
hists_photonpt[17+5]->Add(hists_photonpt[16+5]);
hists_photonpt[17+5]->SetFillColor(kViolet-7);
hs1->Add(hists_photonpt[17+5]);
hists_photonpt[18+5]->SetFillColor(kOrange-2);
hs1->Add(hists_photonpt[18+5]);
hists_photonpt[20+5]->Add(hists_photonpt[19+5]);
hists_photonpt[20+5]->SetFillColor(kYellow+3);
hs1->Add(hists_photonpt[20+5]);

hs1->Draw("hist");
hs1->SetMaximum(1.2*hists_photonpt[samples_.size()+2]->GetMaximum());
hs1->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
hs1->GetYaxis()->SetTitle("NEVENTS");

hists_photonpt[21+5]->SetLineColor(kRed+3);
hists_photonpt[21+5]->SetLineWidth(3);
hists_photonpt[21+5]->Draw("histsame");

hists_photonpt[samples_.size()+2]->SetLineWidth(3.);
hists_photonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_photonpt[samples_.size()+2]->Draw("esame");


TCanvas *c2 = new TCanvas("c2","photonptnoPU",900,700);
c2->cd(0);

THStack *hs2 = new THStack("hs11","photonptnoPU");
hists_photonptnoPU[0]->SetFillColor(kBlue-4);
hs2->Add(hists_photonptnoPU[0]);
hists_photonptnoPU[1]->SetFillColor(kOrange+7);
hs2->Add(hists_photonptnoPU[1]);
hists_photonptnoPU[3]->Add(hists_photonptnoPU[2]);
hists_photonptnoPU[4]->Add(hists_photonptnoPU[3]);
hists_photonptnoPU[5]->Add(hists_photonptnoPU[4]);
hists_photonptnoPU[6]->Add(hists_photonptnoPU[5]);
hists_photonptnoPU[7]->Add(hists_photonptnoPU[6]);
hists_photonptnoPU[7]->SetFillColor(19);
hs2->Add(hists_photonptnoPU[7]);
hists_photonptnoPU[8]->SetFillColor(kGreen+1);
hs2->Add(hists_photonptnoPU[8]);
hists_photonptnoPU[10]->Add(hists_photonptnoPU[9]);
hists_photonptnoPU[6+5]->Add(hists_photonptnoPU[5+5]);
hists_photonptnoPU[7+5]->Add(hists_photonptnoPU[6+5]);
hists_photonptnoPU[8+5]->Add(hists_photonptnoPU[7+5]);
hists_photonptnoPU[9+5]->Add(hists_photonptnoPU[8+5]);
hists_photonptnoPU[9+5]->SetFillColor(kAzure+10);
hs2->Add(hists_photonptnoPU[9+5]);
hists_photonptnoPU[11+5]->Add(hists_photonptnoPU[10+5]);
hists_photonptnoPU[12+5]->Add(hists_photonptnoPU[11+5]);
hists_photonptnoPU[12+5]->SetFillColor(6);
hs2->Add(hists_photonptnoPU[12+5]);
hists_photonptnoPU[13+5]->SetFillColor(32);
hs2->Add(hists_photonptnoPU[13+5]);
hists_photonptnoPU[14+5]->SetFillColor(kSpring-9);
hs2->Add(hists_photonptnoPU[14+5]);
hists_photonptnoPU[16+5]->Add(hists_photonptnoPU[15+5]);
hists_photonptnoPU[17+5]->Add(hists_photonptnoPU[16+5]);
hists_photonptnoPU[17+5]->SetFillColor(kViolet-7);
hs2->Add(hists_photonptnoPU[17+5]);
hists_photonptnoPU[18+5]->SetFillColor(kOrange-2);
hs2->Add(hists_photonptnoPU[18+5]);
hists_photonptnoPU[20+5]->Add(hists_photonptnoPU[19+5]);
hists_photonptnoPU[20+5]->SetFillColor(kYellow+3);
hs2->Add(hists_photonptnoPU[20+5]);

hs2->Draw("hist");
hs2->SetMaximum(1.2*hists_photonptnoPU[samples_.size()+2]->GetMaximum());
hs2->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
hs2->GetYaxis()->SetTitle("NEVENTS");

hists_photonptnoPU[21+5]->SetLineColor(kRed+3);
hists_photonptnoPU[21+5]->SetLineWidth(3);
hists_photonptnoPU[21+5]->Draw("histsame");

hists_photonptnoPU[samples_.size()+2]->SetLineWidth(3.);
hists_photonptnoPU[samples_.size()+2]->SetLineColor(kBlack);
hists_photonptnoPU[samples_.size()+2]->SetMarkerColor(kBlack);
hists_photonptnoPU[samples_.size()+2]->SetMarkerStyle(20.);
hists_photonptnoPU[samples_.size()+2]->Draw("esame");

TCanvas *c3 = new TCanvas("c3","muonpt",900,700);
c3->cd(0);

THStack *hs3 = new THStack("hs11","muonpt");
hists_muonpt[0]->SetFillColor(kBlue-4);
hs3->Add(hists_muonpt[0]);
hists_muonpt[1]->SetFillColor(kOrange+7);
hs3->Add(hists_muonpt[1]);
hists_muonpt[3]->Add(hists_muonpt[2]);
hists_muonpt[4]->Add(hists_muonpt[3]);
hists_muonpt[5]->Add(hists_muonpt[4]);
hists_muonpt[6]->Add(hists_muonpt[5]);
hists_muonpt[7]->Add(hists_muonpt[6]);
hists_muonpt[7]->SetFillColor(19);
hs3->Add(hists_muonpt[7]);
hists_muonpt[8]->SetFillColor(kGreen+1);
hs3->Add(hists_muonpt[8]);
hists_muonpt[10]->Add(hists_muonpt[9]);
hists_muonpt[6+5]->Add(hists_muonpt[5+5]);
hists_muonpt[7+5]->Add(hists_muonpt[6+5]);
hists_muonpt[8+5]->Add(hists_muonpt[7+5]);
hists_muonpt[9+5]->Add(hists_muonpt[8+5]);
hists_muonpt[9+5]->SetFillColor(kAzure+10);
hs3->Add(hists_muonpt[9+5]);
hists_muonpt[11+5]->Add(hists_muonpt[10+5]);
hists_muonpt[12+5]->Add(hists_muonpt[11+5]);
hists_muonpt[12+5]->SetFillColor(6);
hs3->Add(hists_muonpt[12+5]);
hists_muonpt[13+5]->SetFillColor(32);
hs3->Add(hists_muonpt[13+5]);
hists_muonpt[14+5]->SetFillColor(kSpring-9);
hs3->Add(hists_muonpt[14+5]);
hists_muonpt[16+5]->Add(hists_muonpt[15+5]);
hists_muonpt[17+5]->Add(hists_muonpt[16+5]);
hists_muonpt[17+5]->SetFillColor(kViolet-7);
hs3->Add(hists_muonpt[17+5]);
hists_muonpt[18+5]->SetFillColor(kOrange-2);
hs3->Add(hists_muonpt[18+5]);
hists_muonpt[20+5]->Add(hists_muonpt[19+5]);
hists_muonpt[20+5]->SetFillColor(kYellow+3);
hs3->Add(hists_muonpt[20+5]);

hs3->Draw("hist");
hs3->SetMaximum(1.2*hists_muonpt[samples_.size()+2]->GetMaximum());
hs3->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
hs3->GetYaxis()->SetTitle("NEVENTS");

hists_muonpt[21+5]->SetLineColor(kRed+3);
hists_muonpt[21+5]->SetLineWidth(3);
hists_muonpt[21+5]->Draw("histsame");

hists_muonpt[samples_.size()+2]->SetLineWidth(3.);
hists_muonpt[samples_.size()+2]->SetLineColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muonpt[samples_.size()+2]->SetMarkerStyle(20.);
hists_muonpt[samples_.size()+2]->Draw("esame");

TCanvas *c4 = new TCanvas("c4","muonptnoPU",900,700);
c4->cd(0);

THStack *hs4 = new THStack("hs4","muonptnoPU");
hists_muonptnoPU[0]->SetFillColor(kBlue-4);
hs4->Add(hists_muonptnoPU[0]);
hists_muonptnoPU[1]->SetFillColor(kOrange+7);
hs4->Add(hists_muonptnoPU[1]);
hists_muonptnoPU[3]->Add(hists_muonptnoPU[2]);
hists_muonptnoPU[4]->Add(hists_muonptnoPU[3]);
hists_muonptnoPU[5]->Add(hists_muonptnoPU[4]);
hists_muonptnoPU[6]->Add(hists_muonptnoPU[5]);
hists_muonptnoPU[7]->Add(hists_muonptnoPU[6]);
hists_muonptnoPU[7]->SetFillColor(19);
hs4->Add(hists_muonptnoPU[7]);
hists_muonptnoPU[8]->SetFillColor(kGreen+1);
hs4->Add(hists_muonpt[8]);
hists_muonptnoPU[10]->Add(hists_muonptnoPU[9]);
hists_muonptnoPU[6+5]->Add(hists_muonptnoPU[5+5]);
hists_muonptnoPU[7+5]->Add(hists_muonptnoPU[6+5]);
hists_muonptnoPU[8+5]->Add(hists_muonptnoPU[7+5]);
hists_muonptnoPU[9+5]->Add(hists_muonptnoPU[8+5]);
hists_muonptnoPU[9+5]->SetFillColor(kAzure+10);
hs4->Add(hists_muonptnoPU[9+5]);
hists_muonptnoPU[11+5]->Add(hists_muonptnoPU[10+5]);
hists_muonptnoPU[12+5]->Add(hists_muonptnoPU[11+5]);
hists_muonptnoPU[12+5]->SetFillColor(6);
hs4->Add(hists_muonptnoPU[12+5]);
hists_muonptnoPU[13+5]->SetFillColor(32);
hs4->Add(hists_muonptnoPU[13+5]);
hists_muonptnoPU[14+5]->SetFillColor(kSpring-9);
hs4->Add(hists_muonptnoPU[14+5]);
hists_muonptnoPU[16+5]->Add(hists_muonptnoPU[15+5]);
hists_muonptnoPU[17+5]->Add(hists_muonptnoPU[16+5]);
hists_muonptnoPU[17+5]->SetFillColor(kViolet-7);
hs4->Add(hists_muonptnoPU[17+5]);
hists_muonptnoPU[18+5]->SetFillColor(kOrange-2);
hs4->Add(hists_muonptnoPU[18+5]);
hists_muonptnoPU[20+5]->Add(hists_muonptnoPU[19+5]);
hists_muonptnoPU[20+5]->SetFillColor(kYellow+3);
hs4->Add(hists_muonptnoPU[20+5]);

hs4->Draw("hist");
hs4->SetMaximum(1.2*hists_muonptnoPU[samples_.size()+2]->GetMaximum());
hs4->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
hs4->GetYaxis()->SetTitle("NEVENTS");

hists_muonptnoPU[21+5]->SetLineColor(kRed+3);
hists_muonptnoPU[21+5]->SetLineWidth(3);
hists_muonptnoPU[21+5]->Draw("histsame");

hists_muonptnoPU[samples_.size()+2]->SetLineWidth(3.);
hists_muonptnoPU[samples_.size()+2]->SetLineColor(kBlack);
hists_muonptnoPU[samples_.size()+2]->SetMarkerColor(kBlack);
hists_muonptnoPU[samples_.size()+2]->SetMarkerStyle(20.);
hists_muonptnoPU[samples_.size()+2]->Draw("esame");

*/

TLegend* leg = new TLegend(0.65,0.4,0.85,0.90);
  leg->SetFillStyle ( 0);
  leg->SetFillColor ( 0);
  leg->SetBorderSize( 0);
  leg->AddEntry( hists_Nvertex[25], "Single top+#gamma"               , "F");
  leg->AddEntry( hists_Nvertex[23], "Z#gamma"               , "F");
  leg->AddEntry( hists_Nvertex[22], "Diboson"               , "F");
  leg->AddEntry( hists_Nvertex[19], "WW#gamma"               , "F");
  leg->AddEntry( hists_Nvertex[18], "t#bar{t}#gamma"               , "F");
  leg->AddEntry( hists_Nvertex[17], "t#bar{t}"               , "F");
  leg->AddEntry( hists_Nvertex[14], "Single top"               , "F");
  leg->AddEntry( hists_Nvertex[8], "W#gamma"              , "F");
  leg->AddEntry( hists_Nvertex[7], "#gamma+jets"                           , "F");
  leg->AddEntry( hists_Nvertex[1], "Z+jets"                           , "F");
  leg->AddEntry( hists_Nvertex[0], "W+jets"                           , "F");
  leg->AddEntry( hists_Nvertex[26], "Signal(t#gammau)"               , "F");
  leg->AddEntry( hists_Nvertex[samples_.size()+2], "CMS Data 2012(19.145/fb)"               , "PL");


TCanvas *c11 = new TCanvas("c11","Nvertex",900,700);
c11->cd(0);

THStack *hs11 = new THStack("hs11","Nvertex");
hists_Nvertex[0]->SetFillColor(kBlue-4);
hs11->Add(hists_Nvertex[0]);
hists_Nvertex[1]->SetFillColor(kOrange+7);
hs11->Add(hists_Nvertex[1]);
hists_Nvertex[3]->Add(hists_Nvertex[2]);
hists_Nvertex[4]->Add(hists_Nvertex[3]);
hists_Nvertex[5]->Add(hists_Nvertex[4]);
hists_Nvertex[6]->Add(hists_Nvertex[5]);
hists_Nvertex[7]->Add(hists_Nvertex[6]);
hists_Nvertex[7]->SetFillColor(19);
hs11->Add(hists_Nvertex[7]);
hists_Nvertex[8]->SetFillColor(kGreen+1);
hs11->Add(hists_Nvertex[8]);
hists_Nvertex[10]->Add(hists_Nvertex[9]);
hists_Nvertex[6+5]->Add(hists_Nvertex[5+5]);
hists_Nvertex[7+5]->Add(hists_Nvertex[6+5]);
hists_Nvertex[8+5]->Add(hists_Nvertex[7+5]);
hists_Nvertex[9+5]->Add(hists_Nvertex[8+5]);
hists_Nvertex[9+5]->SetFillColor(kAzure+10);
hs11->Add(hists_Nvertex[9+5]);
hists_Nvertex[11+5]->Add(hists_Nvertex[10+5]);
hists_Nvertex[12+5]->Add(hists_Nvertex[11+5]);
hists_Nvertex[12+5]->SetFillColor(6);
hs11->Add(hists_Nvertex[12+5]);
hists_Nvertex[13+5]->SetFillColor(32);
hs11->Add(hists_Nvertex[13+5]);
hists_Nvertex[14+5]->SetFillColor(kSpring-9);
hs11->Add(hists_Nvertex[14+5]);
hists_Nvertex[16+5]->Add(hists_Nvertex[15+5]);
hists_Nvertex[17+5]->Add(hists_Nvertex[16+5]);
hists_Nvertex[17+5]->SetFillColor(kViolet-7);
hs11->Add(hists_Nvertex[17+5]);
hists_Nvertex[18+5]->SetFillColor(kOrange-2);
hs11->Add(hists_Nvertex[18+5]);
hists_Nvertex[20+5]->Add(hists_Nvertex[19+5]);
hists_Nvertex[20+5]->SetFillColor(kYellow+3);
hs11->Add(hists_Nvertex[20+5]);

hs11->Draw("hist");
hs11->SetMaximum(1.2*hists_Nvertex[samples_.size()+2]->GetMaximum());
hs11->GetXaxis()->SetTitle("Nvertex");
hs11->GetYaxis()->SetTitle("NEVENTS");

hists_Nvertex[21+5]->SetLineColor(kRed+3);
hists_Nvertex[21+5]->SetLineWidth(3);
hists_Nvertex[21+5]->Draw("histsame");

hists_Nvertex[samples_.size()+2]->SetLineWidth(3.);
hists_Nvertex[samples_.size()+2]->SetLineColor(kBlack);
hists_Nvertex[samples_.size()+2]->SetMarkerColor(kBlack);
hists_Nvertex[samples_.size()+2]->SetMarkerStyle(20.);
hists_Nvertex[samples_.size()+2]->Draw("esame");
leg->Draw("same");

TCanvas *c111 = new TCanvas("c111","Nvertex",900,700);
c111->cd(0);

THStack *hs111 = new THStack("hs111","NvertexnoPU");
hists_NvertexnoPU[0]->SetFillColor(kBlue-4);
hs111->Add(hists_NvertexnoPU[0]);
hists_NvertexnoPU[1]->SetFillColor(kOrange+7);
hs111->Add(hists_NvertexnoPU[1]);
hists_NvertexnoPU[3]->Add(hists_NvertexnoPU[2]);
hists_NvertexnoPU[4]->Add(hists_NvertexnoPU[3]);
hists_NvertexnoPU[5]->Add(hists_NvertexnoPU[4]);
hists_NvertexnoPU[6]->Add(hists_NvertexnoPU[5]);
hists_NvertexnoPU[7]->Add(hists_NvertexnoPU[6]);
hists_NvertexnoPU[7]->SetFillColor(19);
hs111->Add(hists_NvertexnoPU[7]);
hists_NvertexnoPU[8]->SetFillColor(kGreen+1);
hs111->Add(hists_NvertexnoPU[8]);
hists_NvertexnoPU[10]->Add(hists_NvertexnoPU[9]);
hists_NvertexnoPU[6+5]->Add(hists_NvertexnoPU[5+5]);
hists_NvertexnoPU[7+5]->Add(hists_NvertexnoPU[6+5]);
hists_NvertexnoPU[8+5]->Add(hists_NvertexnoPU[7+5]);
hists_NvertexnoPU[9+5]->Add(hists_NvertexnoPU[8+5]);
hists_NvertexnoPU[9+5]->SetFillColor(kAzure+10);
hs111->Add(hists_NvertexnoPU[9+5]);
hists_NvertexnoPU[11+5]->Add(hists_NvertexnoPU[10+5]);
hists_NvertexnoPU[12+5]->Add(hists_NvertexnoPU[11+5]);
hists_NvertexnoPU[12+5]->SetFillColor(6);
hs111->Add(hists_NvertexnoPU[12+5]);
hists_NvertexnoPU[13+5]->SetFillColor(32);
hs111->Add(hists_NvertexnoPU[13+5]);
hists_NvertexnoPU[14+5]->SetFillColor(kSpring-9);
hs111->Add(hists_NvertexnoPU[14+5]);
hists_NvertexnoPU[16+5]->Add(hists_NvertexnoPU[15+5]);
hists_NvertexnoPU[17+5]->Add(hists_NvertexnoPU[16+5]);
hists_NvertexnoPU[17+5]->SetFillColor(kViolet-7);
hs111->Add(hists_NvertexnoPU[17+5]);
hists_NvertexnoPU[18+5]->SetFillColor(kOrange-2);
hs111->Add(hists_NvertexnoPU[18+5]);
hists_NvertexnoPU[20+5]->Add(hists_NvertexnoPU[19+5]);
hists_NvertexnoPU[20+5]->SetFillColor(kYellow+3);
hs111->Add(hists_NvertexnoPU[20+5]);

hs111->Draw("hist");
hs111->SetMaximum(1.2*hists_NvertexnoPU[samples_.size()+2]->GetMaximum());
hs111->GetXaxis()->SetTitle("Nvertex");
hs111->GetYaxis()->SetTitle("NEVENTS");

hists_Nvertex[21+5]->SetLineColor(kRed+3);
hists_Nvertex[21+5]->SetLineWidth(3);
hists_Nvertex[21+5]->Draw("histsame");

hists_Nvertex[samples_.size()+2]->SetLineWidth(3.);
hists_Nvertex[samples_.size()+2]->SetLineColor(kBlack);
hists_Nvertex[samples_.size()+2]->SetMarkerColor(kBlack);
hists_Nvertex[samples_.size()+2]->SetMarkerStyle(20.);
hists_Nvertex[samples_.size()+2]->Draw("esame");
leg->Draw("same");
/*
for(unsigned int idx=1; idx<samples_.size(); ++idx){
        hists_Nvertex[idx]->Add(hists_Nvertex[idx-1]);
        hists_NvertexnoPU[idx]->Add(hists_NvertexnoPU[idx-1]);
}
hists_Nvertex[samples_.size()-1]->Scale(1/hists_Nvertex[samples_.size()-1]->Integral());
hists_NvertexnoPU[samples_.size()-1]->Scale(1/hists_NvertexnoPU[samples_.size()-1]->Integral());

TH1F *sum_h= new TH1F ( *hists_Nvertex[samples_.size()+2]) ;
sum_h->Scale(1/sum_h->Integral());

TCanvas *c13 = new TCanvas("c13","Nvertexnormalize",900,700);
c13->cd(0);
sum_h->SetMaximum(1.2*sum_h->GetMaximum());
sum_h->SetLineColor(kRed);
sum_h->SetLineWidth(2);
sum_h->GetYaxis()->SetTitle("probability");
sum_h->GetXaxis()->SetTitle("#Delta #phi (MET,#gamma)");


//hists[0]->SetFillColor(5);
sum_h->Draw();
//hists[0]->Draw();

//hists[1]->SetFillColor(1);
//hists[1]->SetFillStyle(3004);
hists_Nvertex[samples_.size()-1]->SetLineWidth(2);
hists_Nvertex[samples_.size()-1]->Draw("sames");
//hists[1]->Draw("sames");

//TLegend* leg = new TLegend(0.55,0.6,0.85,0.80);
//  leg->SetFillStyle ( 0);
//  leg->SetFillColor ( 0);
//  leg->SetBorderSize( 0);
//  leg->AddEntry(sum_h, "DATA"                           , "L");
//  leg->AddEntry( hists_Nvertex[1], "W#gammajets outside top mass window"                           , "L");




TCanvas *c14 = new TCanvas("c14","Nvertexnormalize",900,700);
c14->cd(0);
sum_h->SetMaximum(1.2*sum_h->GetMaximum());
sum_h->SetLineColor(kRed);
sum_h->SetLineWidth(2);
sum_h->GetYaxis()->SetTitle("probability");
sum_h->GetXaxis()->SetTitle("#Delta #phi (MET,#gamma)");


//hists[0]->SetFillColor(5);
sum_h->Draw();
//hists[0]->Draw();

//hists[1]->SetFillColor(1);
//hists[1]->SetFillStyle(3004);
hists_NvertexnoPU[samples_.size()-1]->SetLineWidth(2);
hists_NvertexnoPU[samples_.size()-1]->Draw("sames");
//hists[1]->Draw("sames");


  //leg->Draw("same");
*/
}






