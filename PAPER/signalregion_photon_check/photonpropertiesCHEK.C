#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLine.h"
#include "THStack.h"

#include <cmath> 
void photonpropertiesCHEK() {
 	  	  	
// list of valid histogram names


// list of input samples

std::vector<string> datasamples_;


TH2F *H_ch = new TH2F("H_ch","#sigma_{i#eta i#eta} vs CH isol for DATA",500,0.,0.08,500,0.,200.);
TH2F *H_chbarrel = new TH2F("H_chbarrel","#sigma_{i#eta i#eta} vs CH isol for DATA (Barrel)",500,0.,0.08,500,0.,200.);
TH2F *H_chendcap = new TH2F("H_chendcap","#sigma_{i#eta i#eta} vs CH isol for DATA (Endcaps)",500,0.,0.08,500,0.,200.);


TH2F *H_ph = new TH2F("H_ph","#sigma_{i#eta i#eta} vs PH isol for DATA",500,0.,0.08,500,0.,200.);
TH2F *H_phbarrel = new TH2F("H_phbarrel","#sigma_{i#eta i#eta} vs PH isol for DATA (Barrel)",500,0.,0.08,500,0.,200.);
TH2F *H_phendcap = new TH2F("H_phendcap","#sigma_{i#eta i#eta} vs PH isol for DATA (Endcaps)",500,0.,0.08,500,0.,200.);


TH2F *H_sigmacoswgammadata = new TH2F("H_sigmacoswgamma","#sigma_{i#eta i#eta} vs cosWgamma for DATA ",100,-3,3,100,-1,1);


TH1F *coswgamma= new TH1F("cos","cosw gamma for DATA",20,-1,1);
TH2F *H_etapeakendcapdata = new TH2F("H_etapeak","#sigma_{i#eta i#eta} vs eta for data (ENDCAP)",100,-3,3,100,0.,0.08);

TH1F *eta = new TH1F("eta","eta for DATA IN CR",100,-3,3);
TH1F *phi = new TH1F("phi","phi for DATA IN CR",100,-3,3);

TH1F *etadata = new TH1F("etadata","eta for DATA",100,-3,3);
TH1F *ptdata = new TH1F("ptdata","pt for DATA",100,0,500);
TH1F *HOEdata = new TH1F("HOEdata","HOE for DATA",100,0,1);
TH1F *pfCHdata = new TH1F("pfCHdata","pfCH for DATA",500,0,200);
TH1F *pfNUdata = new TH1F("pfNUdata","pfNU for DATA",100,0,100);
TH1F *pfPHdata = new TH1F("pfPHdata","pfPH for DATA",500,0,200);
TH1F *sigmaietaietadata = new TH1F("sigmaietaietadata","sigmaietaieta for DATA",100,0,0.1);

TH2F *etapt= new TH2F("etapt","#eta vs Pt for DATA",100,-3,3,100,0,500);
TH2F *etaHOE= new TH2F("etaHOE","#eta vs H/E for DATA",100,-3,3,400,0.,0.1);
TH2F *etapfCH= new TH2F("etapfCH","#eta vs CH isol for DATA",100,-3,3,500,0.,200.);
TH2F *etapfNU= new TH2F("etapfNU","#eta vs NU isol for DATA",100,-3,3,500,0.,100.);
TH2F *etapfPH= new TH2F("etapfPH","#eta vs PH isol for DATA",100,-3,3,500,0.,200.);
TH2F *etasigmaietaieta= new TH2F("etasigmaietaieta","#eta vs sigmaietaieta for DATA",100,-3,3,100,0.,0.1);


TH2F *ptHOE= new TH2F("ptHOE","pt vs H/E for DATA",100,0,500,400,0.,0.1);
TH2F *ptpfCH= new TH2F("ptpfCH","pt vs CH isol for DATA",100,0,500,500,0.,200.);
TH2F *ptpfNU= new TH2F("ptpfNU","pt vs NU isol for DATA",100,0,500,500,0.,100.);
TH2F *ptpfPH= new TH2F("ptpfPH","pt vs PH isol for DATA",100,0,500,500,0.,200.);
TH2F *ptsigmaietaieta= new TH2F("ptsigmaietaieta","pt vs sigmaietaieta for DATA",100,0,500,100,0.,0.1);

TH2F *HOEpfCH= new TH2F("HOEpfCH","HOE vs CH isol for DATA",400,0.,0.1,500,0.,200.);
TH2F *HOEpfNU= new TH2F("HOEpfNU","HOE vs NU isol for DATA",400,0.,0.1,500,0.,100.);
TH2F *HOEpfPH= new TH2F("HOEpfPH","HOE vs PH isol for DATA",400,0.,0.1,500,0.,200.);
TH2F *HOEsigmaietaieta= new TH2F("HOEsigmaietaieta","HOE vs sigmaietaieta for DATA",400,0.,0.1,100,0.,0.1);

TH2F *pfCHpfNU= new TH2F("pfCHpfNU","pfCH vs NU isol for DATA",500,0.,200.,500,0.,100.);
TH2F *pfCHpfPH= new TH2F("pfCHpfPH","pfCH vs PH isol for DATA",500,0.,200.,500,0.,200.);
TH2F *pfCHsigmaietaieta= new TH2F("pfCHsigmaietaieta","pfCH vs sigmaietaieta for DATA",500,0.,200.,100,0.,0.1);

TH2F *pfNUpfPH= new TH2F("pfNUpfPH","pfNU vs PH isol for DATA",500,0.,100.,500,0.,200.);
TH2F *pfNUsigmaietaieta= new TH2F("pfNUsigmaietaieta","pfNU vs sigmaietaieta for DATA",500,0.,100.,100,0.,0.1);

TH2F *pfPHsigmaietaieta= new TH2F("pfPHsigmaietaieta","pfPH vs sigmaietaieta for DATA",500,0.,200.,100,0.,0.1);


//datasamples_.push_back("REALDATA1.root");
//datasamples_.push_back("REALDATA2.root");
//datasamples_.push_back("REALDATA3.root");

//datasamples_.push_back("REALDATAab.root");
//datasamples_.push_back("REALDATAc.root");
datasamples_.push_back("REALDATAd.root");


std::vector<TH1F*> hists_photonpt;
std::vector<TH1F*> hists_photoneta;
std::vector<TH1F*> hists_muonpt;
std::vector<TH1F*> hists_muoneta;
std::vector<TH1F*> hists_jetpt;
std::vector<TH1F*> hists_jeteta;
std::vector<TH1F*> hists_topmass;
std::vector<TH1F*> hists_toppt;
std::vector<TH1F*> hists_topeta;
std::vector<TH1F*> hists_topphotonmass;
std::vector<TH1F*> hists_Nvertex;
std::vector<TH1F*> hists_NvertexnoPU;

std::vector<TH1F*> hists_CSV;
std::vector<TH1F*> hists_deltaphimetphoton;
std::vector<TH1F*> hists_deltarphotonmuon;
std::vector<TH1F*> hists_deltarphotonbjet;
std::vector<TH1F*> hists_jetmulti;
std::vector<TH1F*> hists_costopphoton;

std::vector<TH1F*> hists_sigmaietaieta;
std::vector<TH1F*> hists_HOE;




for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
hists_photonpt.push_back(new TH1F(  std::string("photon_pt").append(datasamples_[idx]).c_str(), std::string("photon_pt").append(datasamples_[idx]).c_str() , 10,10., 410.));
hists_photonpt[idx]->Sumw2();
hists_photoneta.push_back(new TH1F(  std::string("photon_eta").append(datasamples_[idx]).c_str(),  std::string("photon_eta").append(datasamples_[idx]).c_str() , 12, -3., 3));
hists_photoneta[idx]->Sumw2();
hists_muonpt.push_back(new TH1F(  std::string("muon_pt").append(datasamples_[idx]).c_str(), std::string("muon_pt").append(datasamples_[idx]).c_str() , 12, 0, 312));
hists_muoneta.push_back(new TH1F(  std::string("muon_eta").append(datasamples_[idx]).c_str(), std::string("muon_eta").append(datasamples_[idx]).c_str() ,  12, -3., 3));
hists_jetpt.push_back(new TH1F(  std::string("jet_pt").append(datasamples_[idx]).c_str(),  std::string("jet_pt").append(datasamples_[idx]).c_str(),  12, 0, 312));
hists_jeteta.push_back(new TH1F(  std::string("bjet_eta").append(datasamples_[idx]).c_str(),  std::string("bjet_eta").append(datasamples_[idx]).c_str() , 12, -3., 3));
hists_topmass.push_back(new TH1F(  std::string("top_mass").append(datasamples_[idx]).c_str(), std::string("top_mass").append(datasamples_[idx]).c_str() , 20, 70, 370.));
hists_toppt.push_back(new TH1F(  std::string("toppt").append(datasamples_[idx]).c_str(), std::string("toppt").append(datasamples_[idx]).c_str() ,35,0., 750.));
hists_topeta.push_back(new TH1F(  std::string("topeta").append(datasamples_[idx]).c_str(), std::string("topeta").append(datasamples_[idx]).c_str() ,30, -3., 3.));
hists_topphotonmass.push_back(new TH1F(  std::string("topphotonmass").append(datasamples_[idx]).c_str(), std::string("topphotonmass").append(datasamples_[idx]).c_str() ,40, 0., 1000.));
hists_Nvertex.push_back(new TH1F(  std::string("Nvertex").append(datasamples_[idx]).c_str(), std::string("Nvertex").append(datasamples_[idx]).c_str() ,60, -0.5, 59.5));
hists_NvertexnoPU.push_back(new TH1F(  std::string("NvertexnoPU").append(datasamples_[idx]).c_str(), std::string("NvertexnoPU").append(datasamples_[idx]).c_str() ,60, -0.5, 59.5));
hists_CSV.push_back(new TH1F(  std::string("CSV").append(datasamples_[idx]).c_str(), std::string("CSV").append(datasamples_[idx]).c_str() ,20, 0, 1));

hists_deltaphimetphoton.push_back(new TH1F(  std::string("deltaphimetphoton").append(datasamples_[idx]).c_str(), std::string("deltaphimetphoton").append(datasamples_[idx]).c_str() ,12, 0, 4));
hists_deltaphimetphoton[idx]->Sumw2();
hists_deltarphotonmuon.push_back(new TH1F(  std::string("deltarphotonmuon").append(datasamples_[idx]).c_str(), std::string("deltarphotonmuon").append(datasamples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonmuon[idx]->Sumw2();
hists_deltarphotonbjet.push_back(new TH1F(  std::string("deltarphotonbjet").append(datasamples_[idx]).c_str(), std::string("deltarphotonbjet").append(datasamples_[idx]).c_str() ,14, 0, 6));
hists_deltarphotonbjet[idx]->Sumw2();
hists_jetmulti.push_back(new TH1F(  std::string("jetmulti").append(datasamples_[idx]).c_str(), std::string("jetmulti").append(datasamples_[idx]).c_str() ,8, -0.5, 7.5));
hists_jetmulti[idx]->Sumw2();
hists_costopphoton.push_back(new TH1F(  std::string("costopphoton").append(datasamples_[idx]).c_str(), std::string("costopphoton").append(datasamples_[idx]).c_str() ,22, -1, 1));
hists_costopphoton[idx]->Sumw2();

hists_sigmaietaieta.push_back(new TH1F(  std::string("sigmaietaieta").append(datasamples_[idx]).c_str(), std::string("sigmaietaieta").append(datasamples_[idx]).c_str() ,32, 0, 0.1));
hists_sigmaietaieta[idx]->Sumw2();
hists_HOE.push_back(new TH1F(  std::string("HOE").append(datasamples_[idx]).c_str(), std::string("HOE").append(datasamples_[idx]).c_str() ,8, 0, 0.05));
hists_HOE[idx]->Sumw2();
}


for(unsigned int idx=0; idx<datasamples_.size(); ++idx){
	TFile *input(0);
	TString fname =datasamples_[idx];
	input = TFile::Open( fname ); 
	std::vector<double> *myptphoton=0;
	std::vector<double> *myetaphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myetamuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *myetajet=0;
	std::vector<double> *mymasstop=0;
        std::vector<double> *mypttop=0;
        std::vector<double> *myetatop=0;
        std::vector<double> *mytopphotonmass=0;
        std::vector<double> *myNvertex=0;
	std::vector<double> *myweight=0;
        std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *myjetmultiplicity=0;
	std::vector<double> *mysigmaietaieta=0;
	std::vector<double> *myHOE=0;
	std::vector<double> *mypfChargedPU=0;
	std::vector<double> *mypfNeutralPU=0;
	std::vector<double> *mypfGammaPU=0;
        std::vector<double> *mycoswphoton=0;
        std::vector<double> *myphiphoton=0;

	TTree* theTree = (TTree*)input->Get("analyzestep2/atq");
	cout<<theTree->GetEntries()<<endl;
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "etaphoton", &myetaphoton );
   	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "etamuon", &myetamuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "etajet", &myetajet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
        theTree->SetBranchAddress( "pttop", &mypttop);
        theTree->SetBranchAddress( "etatop", &myetatop);
        theTree->SetBranchAddress( "topphotonmass", &mytopphotonmass);
        theTree->SetBranchAddress( "Nvertex", &myNvertex);
	theTree->SetBranchAddress( "weight", &myweight);
        theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
	theTree->SetBranchAddress( "sigmaieta", &mysigmaietaieta );
	theTree->SetBranchAddress( "HOE", &myHOE );
	theTree->SetBranchAddress( "pfChargedPU", &mypfChargedPU );
	theTree->SetBranchAddress( "pfNeutralPU", &mypfNeutralPU );
	theTree->SetBranchAddress( "pfGammaPU", &mypfGammaPU );
        theTree->SetBranchAddress( "coswphoton", &mycoswphoton );
        theTree->SetBranchAddress( "phiphoton", &myphiphoton );



	for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
	theTree->GetEntry(ievt);

//if (abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031 && (*mypfChargedPU )[0]<0.7 ) eta->Fill((*myetaphoton )[0]);
//if (abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011 && (*mypfChargedPU )[0]<0.5 ) eta->Fill((*myetaphoton )[0]);

if (abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031 && (*mypfChargedPU )[0]<0.7 ) 
{eta->Fill((*myetaphoton )[0]);
phi->Fill((*myphiphoton )[0]);}
if (abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011 && (*mypfChargedPU )[0]<0.5 ) 
{eta->Fill((*myetaphoton )[0]);
phi->Fill((*myphiphoton )[0]);}
//if (abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031 && (*mypfChargedPU )[0]<0.7 ) coswgamma->Fill((*mycoswphoton )[0]);
//if (abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011 && (*mypfChargedPU )[0]<0.5 ) coswgamma->Fill((*mycoswphoton )[0]);
coswgamma->Fill((*mycoswphoton )[0]);

if ( abs((*myetaphoton )[0]) >1.444) H_etapeakendcapdata->Fill((*myetaphoton )[0],(*mysigmaietaieta)[0]);

 H_ch->Fill((*mysigmaietaieta)[0],(*mypfChargedPU )[0]);
if ( abs((*myetaphoton )[0]) <1.444) H_chbarrel->Fill((*mysigmaietaieta)[0],(*mypfChargedPU )[0]);
if ( abs((*myetaphoton )[0]) >1.444) H_chendcap->Fill((*mysigmaietaieta)[0],(*mypfChargedPU )[0]);

 H_ph->Fill((*mysigmaietaieta)[0],(*mypfGammaPU )[0]);
if ( abs((*myetaphoton )[0]) <1.444) H_phbarrel->Fill((*mysigmaietaieta)[0],(*mypfGammaPU )[0]);
if ( abs((*myetaphoton )[0]) >1.444) H_phendcap->Fill((*mysigmaietaieta)[0],(*mypfGammaPU )[0]);

if (abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] >0.011 && (*mypfChargedPU )[0]<0.5 ) H_sigmacoswgammadata->Fill((*myetaphoton)[0],(*mycoswphoton)[0] );
if (abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] >0.031 && (*mypfChargedPU )[0]<0.7 ) H_sigmacoswgammadata->Fill((*myetaphoton)[0],(*mycoswphoton)[0] );


//if ((abs((*myetaphoton )[0]) >1.444 && (*mysigmaietaieta)[0] < 0.031) || (abs((*myetaphoton )[0]) <1.444 && (*mysigmaietaieta)[0] < 0.011)  ) {
//if ((abs((*myetaphoton )[0]) >1.444 &&  (*mypfChargedPU )[0]<0.7) || (abs((*myetaphoton )[0]) <1.444 && (*mypfChargedPU )[0]<0.5)  ) {
//if (((*myetaphoton )[0]<1.444 && (*mypfNeutralPU )[0]<(0.4+0.04*(*myptphoton)[0])) || ((*myetaphoton )[0]>1.444 && (*mypfNeutralPU )[0]<(1.5+0.04*(*myptphoton)[0]))) {
//if (((*myetaphoton )[0]<1.444 && (*mypfGammaPU )[0]<(0.5+0.005*(*myptphoton)[0])) || ((*myetaphoton )[0]>1.56 && (*mypfGammaPU )[0]<(1+0.005*(*myptphoton)[0]))) {
//if ((*myHOE )[0]<0.05){
etadata->Fill( (*myetaphoton )[0] );
ptdata->Fill( (*myptphoton )[0] );
HOEdata->Fill( (*myHOE )[0] );
pfCHdata->Fill( (*mypfChargedPU )[0] );
pfNUdata->Fill( (*mypfNeutralPU )[0] );
pfPHdata->Fill( (*mypfGammaPU )[0] );
sigmaietaietadata->Fill( (*mysigmaietaieta )[0] );

//}

etapt->Fill((*myetaphoton)[0],(*myptphoton)[0] );
etaHOE->Fill((*myetaphoton)[0],(*myHOE)[0] );
etapfCH->Fill((*myetaphoton)[0],(*mypfChargedPU)[0] );
etapfNU->Fill((*myetaphoton)[0],(*mypfNeutralPU)[0] );
etapfPH->Fill((*myetaphoton)[0],(*mypfGammaPU)[0] );
etasigmaietaieta->Fill((*myetaphoton)[0],(*mysigmaietaieta)[0] );

ptHOE->Fill((*myptphoton)[0],(*myHOE)[0] );
ptpfCH->Fill((*myptphoton)[0],(*mypfChargedPU)[0] );
ptpfNU->Fill((*myptphoton)[0],(*mypfNeutralPU)[0] );
ptpfPH->Fill((*myptphoton)[0],(*mypfGammaPU)[0] );
ptsigmaietaieta->Fill((*myptphoton)[0],(*mysigmaietaieta)[0] );

HOEpfCH->Fill((*myHOE)[0],(*mypfChargedPU)[0] );
HOEpfNU->Fill((*myHOE)[0],(*mypfNeutralPU)[0] );
HOEpfPH->Fill((*myHOE)[0],(*mypfGammaPU)[0] );
HOEsigmaietaieta->Fill((*myHOE)[0],(*mysigmaietaieta)[0] );

pfCHpfNU->Fill((*mypfChargedPU)[0],(*mypfNeutralPU)[0] );
pfCHpfPH->Fill((*mypfChargedPU)[0],(*mypfGammaPU)[0] );
pfCHsigmaietaieta->Fill((*mypfChargedPU)[0],(*mysigmaietaieta)[0] );

pfNUpfPH->Fill((*mypfNeutralPU)[0],(*mypfGammaPU)[0] );
pfNUsigmaietaieta->Fill((*mypfNeutralPU)[0],(*mysigmaietaieta)[0] );

pfPHsigmaietaieta->Fill((*mypfGammaPU)[0],(*mysigmaietaieta)[0] );

	}
	delete myptphoton;
	delete myetaphoton;
	delete myptmuon;
	delete myetamuon;
	delete myptjet;
	delete myetajet;
	delete mymasstop;
        delete mypttop;
        delete myetatop;
        delete mytopphotonmass;
        delete myNvertex;
	delete input;
	delete mycvsdiscriminant;
	delete mydeltaRphotonjet;
	delete mydeltaRphotonmuon;
	delete mycostopphoton;
	delete mydeltaphiphotonmet;
	delete myjetmultiplicity;
	delete mysigmaietaieta;
	delete myHOE;
	delete mypfChargedPU;
	delete mypfNeutralPU;
	delete mypfGammaPU;
}






  TCanvas *c22 = new TCanvas("c22","22222",900,700);
  c22->cd(0);
H_ch->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_ch->GetYaxis()->SetTitle("CH isol");
H_ch->Draw("colz");

  TCanvas *c23 = new TCanvas("c23","23",900,700);
  c23->cd(0);
H_ph->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_ph->GetYaxis()->SetTitle("PH isol");
H_ph->Draw("colz");


  TCanvas *c24 = new TCanvas("c24","24",900,700);
  c24->cd(0);
H_chbarrel->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_chbarrel->GetYaxis()->SetTitle("CH isol");
H_chbarrel->Draw("colz");

  TCanvas *c25 = new TCanvas("c25","25",900,700);
  c25->cd(0);
H_phbarrel->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_phbarrel->GetYaxis()->SetTitle("PH isol");
H_phbarrel->Draw("colz");

  TCanvas *c26 = new TCanvas("c26","26",900,700);
  c26->cd(0);
H_chendcap->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_chendcap->GetYaxis()->SetTitle("CH isol");
H_chendcap->Draw("colz");

  TCanvas *c27 = new TCanvas("c27","27",900,700);
  c27->cd(0);
H_phendcap->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
H_phendcap->GetYaxis()->SetTitle("PH isol");
H_phendcap->Draw("colz");





  TCanvas *c31 = new TCanvas("c31","31",900,700);
  c31->cd(0);
eta->GetXaxis()->SetTitle("#eta ");
eta->GetYaxis()->SetTitle("NEvents");
eta->SetLineWidth(3);
eta->SetLineColor(1);
eta->Draw("");

  TCanvas *c431 = new TCanvas("c431","431",900,700);
  c431->cd(0);
phi->GetXaxis()->SetTitle("#phi ");
phi->GetYaxis()->SetTitle("NEvents");
phi->SetLineWidth(3);
phi->SetLineColor(1);
phi->Draw("");

  TCanvas *c32 = new TCanvas("c32","32",900,700);
  c32->cd(0);
H_etapeakendcapdata->GetXaxis()->SetTitle("#eta ");
H_etapeakendcapdata->GetYaxis()->SetTitle("#sigma_{i#eta i#eta}");
H_etapeakendcapdata->Draw("colz");


  TCanvas *c34 = new TCanvas("c34","34",900,700);
  c34->cd(0);
coswgamma->GetXaxis()->SetTitle("#eta ");
coswgamma->GetYaxis()->SetTitle("NEvents");
coswgamma->Draw("");

  TCanvas *c35 = new TCanvas("c35","35",900,700);
  c35->cd(0);
H_sigmacoswgammadata->GetXaxis()->SetTitle("#eta ");
H_sigmacoswgammadata->GetYaxis()->SetTitle("cos(W,gamma)");
H_sigmacoswgammadata->Draw("colz");

  TCanvas *c36 = new TCanvas("c36","36",900,700);
  c36->cd(0);
etadata->GetXaxis()->SetTitle("#eta ");
etadata->GetYaxis()->SetTitle("NEvents");
etadata->SetLineWidth(3);
etadata->SetLineColor(1);
etadata->Draw("");



  TCanvas *c37 = new TCanvas("c37","37",900,700);
  c37->cd(0);
HOEdata->GetXaxis()->SetTitle("HOE ");
HOEdata->GetYaxis()->SetTitle("NEvents");
HOEdata->SetLineWidth(3);
HOEdata->SetLineColor(1);
HOEdata->Draw("");

  TCanvas *c38 = new TCanvas("c38","38",900,700);
  c38->cd(0);
pfCHdata->GetXaxis()->SetTitle("pfch ");
pfCHdata->GetYaxis()->SetTitle("NEvents");
pfCHdata->SetLineWidth(3);
pfCHdata->SetLineColor(1);
pfCHdata->Draw("");

  TCanvas *c39 = new TCanvas("c39","39",900,700);
  c39->cd(0);
pfNUdata->GetXaxis()->SetTitle("pfnu ");
pfNUdata->GetYaxis()->SetTitle("NEvents");
pfNUdata->SetLineWidth(3);
pfNUdata->SetLineColor(1);
pfNUdata->Draw("");

  TCanvas *c40 = new TCanvas("c40","40",900,700);
  c40->cd(0);
pfPHdata->GetXaxis()->SetTitle("pfphoton ");
pfPHdata->GetYaxis()->SetTitle("NEvents");
pfPHdata->SetLineWidth(3);
pfPHdata->SetLineColor(1);
pfPHdata->Draw("");

  TCanvas *c41 = new TCanvas("c41","41",900,700);
  c41->cd(0);
ptdata->GetXaxis()->SetTitle("Pt ");
ptdata->GetYaxis()->SetTitle("NEvents");
ptdata->SetLineWidth(3);
ptdata->SetLineColor(1);
ptdata->Draw("");

  TCanvas *c42 = new TCanvas("c42","42",900,700);
  c42->cd(0);
sigmaietaietadata->GetXaxis()->SetTitle("#sigma_{i#eta i#eta} ");
sigmaietaietadata->GetYaxis()->SetTitle("NEvents");
sigmaietaietadata->SetLineWidth(3);
sigmaietaietadata->SetLineColor(1);
sigmaietaietadata->Draw("");



  TCanvas *c43 = new TCanvas("c43","43",900,700);
  c43->cd(0);
etapt->GetXaxis()->SetTitle("#eta ");
etapt->GetYaxis()->SetTitle("pt");
etapt->Draw("colz");

  TCanvas *c44 = new TCanvas("c44","44",900,700);
  c44->cd(0);
etaHOE->GetXaxis()->SetTitle("#eta ");
etaHOE->GetYaxis()->SetTitle("HOE");
etaHOE->Draw("colz");

  TCanvas *c440 = new TCanvas("c440","440",900,700);
  c440->cd(0);
etapfCH->GetXaxis()->SetTitle("#eta ");
etapfCH->GetYaxis()->SetTitle("CHiso");
etapfCH->Draw("colz");

  TCanvas *c45 = new TCanvas("c45","45",900,700);
  c45->cd(0);
etapfNU->GetXaxis()->SetTitle("#eta ");
etapfNU->GetYaxis()->SetTitle("NUiso");
etapfNU->Draw("colz");

  TCanvas *c46 = new TCanvas("c46","46",900,700);
  c46->cd(0);
etapfPH->GetXaxis()->SetTitle("#eta ");
etapfPH->GetYaxis()->SetTitle("PHiso");
etapfPH->Draw("colz");

  TCanvas *c47 = new TCanvas("c47","47",900,700);
  c47->cd(0);
etasigmaietaieta->GetXaxis()->SetTitle("#eta ");
etasigmaietaieta->GetYaxis()->SetTitle("sigmaietaieta");
etasigmaietaieta->Draw("colz");


 TCanvas *c48 = new TCanvas("c48","48",900,700);
  c48->cd(0);
ptHOE->GetXaxis()->SetTitle("Pt ");
ptHOE->GetYaxis()->SetTitle("HOE");
ptHOE->Draw("colz");

 TCanvas *c49 = new TCanvas("c49","49",900,700);
  c49->cd(0);
ptpfCH->GetXaxis()->SetTitle("Pt ");
ptpfCH->GetYaxis()->SetTitle("CHiso");
ptpfCH->Draw("colz");

 TCanvas *c50 = new TCanvas("c50","50",900,700);
  c50->cd(0);
ptpfNU->GetXaxis()->SetTitle("Pt ");
ptpfNU->GetYaxis()->SetTitle("Nuiso");
ptpfNU->Draw("colz");

 TCanvas *c51 = new TCanvas("c51","51",900,700);
  c51->cd(0);
ptpfPH->GetXaxis()->SetTitle("Pt ");
ptpfPH->GetYaxis()->SetTitle("PHiso");
ptpfPH->Draw("colz");

 TCanvas *c52 = new TCanvas("c52","52",900,700);
  c52->cd(0);
ptsigmaietaieta->GetXaxis()->SetTitle("Pt ");
ptsigmaietaieta->GetYaxis()->SetTitle("sigma ieta ieta");
ptsigmaietaieta->Draw("colz");


 TCanvas *c54 = new TCanvas("c54","54",900,700);
  c54->cd(0);
HOEpfCH->GetXaxis()->SetTitle("HOE ");
HOEpfCH->GetYaxis()->SetTitle("CHiso");
HOEpfCH->Draw("colz");

 TCanvas *c55 = new TCanvas("c55","55",900,700);
  c55->cd(0);
HOEpfNU->GetXaxis()->SetTitle("HOE ");
HOEpfNU->GetYaxis()->SetTitle("NUiso");
HOEpfNU->Draw("colz");

 TCanvas *c56 = new TCanvas("c56","56",900,700);
  c56->cd(0);
HOEpfPH->GetXaxis()->SetTitle("HOE ");
HOEpfPH->GetYaxis()->SetTitle("PHiso");
HOEpfPH->Draw("colz");

 TCanvas *c57 = new TCanvas("c57","57",900,700);
  c57->cd(0);
HOEsigmaietaieta->GetXaxis()->SetTitle("HOE ");
HOEsigmaietaieta->GetYaxis()->SetTitle("sigma ieta ieta");
HOEsigmaietaieta->Draw("colz");




 TCanvas *c58 = new TCanvas("c58","58",900,700);
  c58->cd(0);
pfCHpfNU->GetXaxis()->SetTitle("CHiso ");
pfCHpfNU->GetYaxis()->SetTitle("NUiso");
pfCHpfNU->Draw("colz");

 TCanvas *c59 = new TCanvas("c59","59",900,700);
  c59->cd(0);
pfCHpfPH->GetXaxis()->SetTitle("CHiso ");
pfCHpfPH->GetYaxis()->SetTitle("PHiso");
pfCHpfPH->Draw("colz");

 TCanvas *c60 = new TCanvas("c60","60",900,700);
  c60->cd(0);
pfCHsigmaietaieta->GetXaxis()->SetTitle("CHiso ");
pfCHsigmaietaieta->GetYaxis()->SetTitle("sigma ieta ieta");
pfCHsigmaietaieta->Draw("colz");

 TCanvas *c61 = new TCanvas("c61","61",900,700);
  c61->cd(0);
pfNUpfPH->GetXaxis()->SetTitle("NUiso ");
pfNUpfPH->GetYaxis()->SetTitle("PHiso");
pfNUpfPH->Draw("colz");

 TCanvas *c62 = new TCanvas("c62","62",900,700);
  c62->cd(0);
pfNUsigmaietaieta->GetXaxis()->SetTitle("NUiso ");
pfNUsigmaietaieta->GetYaxis()->SetTitle("sigmaietaieta");
pfNUsigmaietaieta->Draw("colz");

 TCanvas *c63 = new TCanvas("c63","63",900,700);
  c63->cd(0);
pfPHsigmaietaieta->GetXaxis()->SetTitle("PH iso ");
pfPHsigmaietaieta->GetYaxis()->SetTitle("sigmaietaieta");
pfPHsigmaietaieta->Draw("colz");


}






