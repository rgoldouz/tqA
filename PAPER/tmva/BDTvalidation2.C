/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

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
#include "TBranch.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVAGui.C"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void BDTvalidation2(){
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

std::vector<TString> samples_;
std::vector<TString> datasamples_;
samples_.push_back("WPHJET.root");
samples_.push_back("ZJET.root");
samples_.push_back("PHJET200400.root");
samples_.push_back("WJET.root");
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
samples_.push_back("SIGNAL.root");
samples_.push_back("REALDATA1.root");
samples_.push_back("REALDATA2.root");
samples_.push_back("REALDATA3.root");
//open files


/////////////////////////////////////////////////////////////////////////////////////

std::vector<TString> variables_;
  variables_.push_back("photon_Pt");
  variables_.push_back("muon_Pt");
  variables_.push_back("Jet_Pt");
  variables_.push_back("topmass");
  variables_.push_back("DeltaRphotonjet");
  variables_.push_back("DeltaRphotonmuon");
  variables_.push_back("Costopphoton");
  variables_.push_back("Jet_Multiplicity");
  variables_.push_back("Deltaphiphotonmet");
  variables_.push_back("b_tag_info");


//for(unsigned int i=0; i<variables_.size(); ++i){
std::vector<TFile*> files;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
files.push_back(new TFile(samples_[idx]));
}

std::vector<TFile*> newfiles;
for(unsigned int idx=0; idx<samples_.size(); ++idx){
newfiles.push_back(new TFile(TString("AddBDTbranch/")+samples_[idx], "RECREATE"));
}

for(unsigned int idx=0; idx<samples_.size(); ++idx){

	std::vector<double> *myptphoton=0;
	std::vector<double> *myptmuon=0;
	std::vector<double> *myptjet=0;
	std::vector<double> *mymasstop=0;
	std::vector<double> *mydeltaRphotonjet=0;
	std::vector<double> *mydeltaRphotonmuon=0;
	std::vector<double> *mycostopphoton=0;
	std::vector<double> *mydeltaphiphotonmet=0;
	std::vector<double> *mycvsdiscriminant=0;
	std::vector<double> *myjetmultiplicity=0;
	
	TTree* theTree = (TTree*) files[idx]->Get("analyzestep2/atq");
	theTree->SetBranchAddress("ptphoton", &myptphoton  );
	theTree->SetBranchAddress( "ptmuon", &myptmuon );
	theTree->SetBranchAddress( "ptjet", &myptjet );
	theTree->SetBranchAddress( "masstop", &mymasstop );
	theTree->SetBranchAddress( "deltaRphotonjet", &mydeltaRphotonjet );
	theTree->SetBranchAddress( "deltaRphotonmuon", &mydeltaRphotonmuon );
	theTree->SetBranchAddress( "costopphoton", &mycostopphoton );
	theTree->SetBranchAddress( "jetmultiplicity", &myjetmultiplicity );
	theTree->SetBranchAddress( "deltaphiphotonmet", &mydeltaphiphotonmet );
	theTree->SetBranchAddress( "cvsdiscriminant", &mycvsdiscriminant );

	newfiles[idx]->cd();
	TTree* BDTTree = theTree->CloneTree(0);

for(unsigned int i=0; i<variables_.size(); ++i){

TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
TString dir    = "weights/";
TString prefix = "TMVA";
TString methodName = TString("BDT method");
TString weightfile = dir + prefix + TString("_") + TString("BDT")+ variables_[i] + TString(".weights.xml");
Float_t ptphoton,ptmuon,ptjet,masstop,deltaRphotonjet,deltaRphotonmuon,costopphoton,deltaphiphotonmet,cvsdiscriminant,jetmultiplicity;

cout<<i<<endl;
	if (i==0){
//	reader->AddVariable ("ptphoton", &ptphoton);
	reader->AddVariable ("ptmuon", &ptmuon);
	reader->AddVariable ("ptjet", &ptjet);
	reader->AddVariable ("masstop", &masstop);
	reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
	reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
	reader->AddVariable ("costopphoton", &costopphoton);
	reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
	reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
	reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
	}       

        if (i==1){
        reader->AddVariable ("ptphoton", &ptphoton);
//        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==2){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
//        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==3){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
//        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==4){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
//        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==5){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
//        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==6){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
//        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==7){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
//        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==8){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
//        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }

        if (i==9){
        reader->AddVariable ("ptphoton", &ptphoton);
        reader->AddVariable ("ptmuon", &ptmuon);
        reader->AddVariable ("ptjet", &ptjet);
        reader->AddVariable ("masstop", &masstop);
        reader->AddVariable ("deltaRphotonjet", &deltaRphotonjet);
        reader->AddVariable ("deltaRphotonmuon", &deltaRphotonmuon);
        reader->AddVariable ("costopphoton", &costopphoton);
        reader->AddVariable ("jetmultiplicity", &jetmultiplicity);
        reader->AddVariable ("deltaphiphotonmet", &deltaphiphotonmet);
//        reader->AddVariable ("cvsdiscriminant", &cvsdiscriminant);
        }
	reader->BookMVA( methodName, weightfile );
        Float_t new_v;
        TString newBDT = TString("BDToutput_removed_")+ variables_[i];
      BDTTree->Branch(newBDT, &new_v, newBDT+TString("/F"));
//        BDTTree->Branch("new_v", &new_v, "new_v/F");

		for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
		theTree->GetEntry(ievt);
		ptphoton=(float)(*myptphoton)[0];
		ptmuon=(float)(*myptmuon )[0];
		ptjet=(float)(*myptjet )[0];
		masstop=(float)(*mymasstop )[0];
		deltaRphotonjet=(float)(*mydeltaRphotonjet )[0];
		deltaRphotonmuon=(float)(*mydeltaRphotonmuon )[0];
		costopphoton=(float)(*mycostopphoton )[0];
		jetmultiplicity=(float)(*myjetmultiplicity )[0];
		deltaphiphotonmet=(float)(*mydeltaphiphotonmet )[0];
		cvsdiscriminant=(float)(*mycvsdiscriminant)[0];
		new_v=reader->EvaluateMVA("BDT method");
		BDTTree->Fill();
		}
}
//        theTree->Print();
	  BDTTree->Write();
       	files[idx]->Close();
delete myptphoton;
delete myptmuon;
delete myptjet;
delete mymasstop;
delete mydeltaRphotonjet;
delete mydeltaRphotonmuon;
delete mycostopphoton;
delete mydeltaphiphotonmet;
delete mycvsdiscriminant;
delete myjetmultiplicity;

	
}
}
