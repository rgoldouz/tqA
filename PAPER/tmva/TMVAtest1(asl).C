#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAtest1()
{

  // ---------------------------------------------------------------
  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;
  // --------------------------------------------------------------------------------------------------
  
  // --- Here the preparation phase begins
  
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString TrainName = "TMVA";
  
  TString outfileName( TrainName+".root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is 
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( TrainName, outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  
  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  
  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  
  factory->AddVariable ("ptphoton", 'F');
//  factory->AddVariable ("etaphoton", 'F');
  factory->AddVariable ("ptmuon", 'F');
// factory->AddVariable ("etamuon", 'F');
  factory->AddVariable ("ptjet", 'F');
//  factory->AddVariable ("etajet", 'F');
//  factory->AddVariable ("masstop", 'F');
//  factory->AddVariable ("mtw", 'F');
  factory->AddVariable ("deltaRphotonjet", 'F');
  factory->AddVariable ("deltaRphotonmuon", 'F');
//  factory->AddVariable ("ht", 'F');
//  factory->AddVariable ("photonmuonmass", 'F');
  factory->AddVariable ("costopphoton", 'F');
//  factory->AddVariable ("topphotonmass", 'F');
//factory->AddVariable ("pttop", 'F');
//factory->AddVariable ("etatop", 'F');
  factory->AddVariable ("jetmultiplicity", 'F');
//  factory->AddVariable ("bjetmultiplicity", 'F');
  factory->AddVariable ("deltaphiphotonmet", 'F');
  factory->AddVariable ("cvsdiscriminant", 'F');
//  factory->AddVariable ("leptoncharge", 'F');


  //Load the signal and background event samples from ROOT trees
  
//  TFile *inputSTrain(0);
//  TFile *inputBTrain(0);
  TString sigFileTrain  = "SIGNALtGu.root";
  TString bkgFileTrain = "WPHJET.root";
  TString ttbarFileTrain1 ="TTBAR1.root";
  TString ttbarFileTrain2 ="TTBAR2.root";
  TString ttbarFileTrain3 ="TTBAR3.root";
  TString TTGFileTrain ="TTG.root";
  TString WWFileTrain ="WW.root";
  TString WZFileTrain ="WZ.root";
  TString ZZFileTrain ="ZZ.root";
  TString ZGAMMAFileTrain ="ZGAMMA.root";


//  TString bkgFileTrain = "TTG.root";
  
  TFile *inputSTrain  = TFile::Open( sigFileTrain );
 if (!inputSTrain) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   } 
  TFile *inputBTrain = TFile::Open( bkgFileTrain );
  TFile *inputBttbar1 = TFile::Open( ttbarFileTrain1);
  TFile *inputBttbar2 = TFile::Open( ttbarFileTrain2);
  TFile *inputBttbar3 = TFile::Open( ttbarFileTrain3);
  TFile *inputTTG= TFile::Open(TTGFileTrain );
  TFile *inputWW= TFile::Open(WWFileTrain );
  TFile *inputWZ= TFile::Open(WZFileTrain );
  TFile *inputZZ= TFile::Open(ZZFileTrain );
  TFile *inputZGAMMA= TFile::Open( ZGAMMAFileTrain);

if (!inputBTrain) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }

 
  std::cout << "--- TMVAnalysis    : Accessing Signal Train: " << sigFileTrain << std::endl;
  std::cout << "--- TMVAnalysis    : Accessing Background Train: " << bkgFileTrain << std::endl;
  
  TTree *signalTrain      = (TTree*)inputSTrain->FindObjectAny("atq");
  TTree *backgroundTrain = (TTree*)inputBTrain->FindObjectAny("atq");
  TTree *ttbarbackground1 = (TTree*)inputBttbar1->FindObjectAny("atq");
  TTree *ttbarbackground2 = (TTree*)inputBttbar2->FindObjectAny("atq");
  TTree *ttbarbackground3 = (TTree*)inputBttbar3->FindObjectAny("atq");
  TTree *TTGbackground = (TTree*)inputTTG->FindObjectAny("atq");
  TTree *WWbackground = (TTree*)inputWW->FindObjectAny("atq");
  TTree *WZbackground = (TTree*)inputWZ->FindObjectAny("atq");
  TTree *ZZbackground = (TTree*)inputZZ->FindObjectAny("atq");
  TTree *ZGAMMAbackground = (TTree*)inputZGAMMA->FindObjectAny("atq");

//  TTree *signalTrain      = (TTree*)inputSTrain->FindObjectAny("insidetopmass");
//  TTree *backgroundTrain = (TTree*)inputBTrain->FindObjectAny("insidetopmass");
  
 // factory->AddSignalTree( signalTrain, 1);
//cout<<"reza1";
  
 // factory->AddBackgroundTree( backgroundTrain, 1, );
//cout<<"reza2";

 factory->SetInputTrees(signalTrain,backgroundTrain,1,1);
 factory->AddBackgroundTree( ttbarbackground1, 1.1);
 factory->AddBackgroundTree( ttbarbackground2, 1.1);
// factory->AddBackgroundTree( ttbarbackground3, 0.3);
 factory->AddBackgroundTree( TTGbackground, 1);
 factory->AddBackgroundTree( WWbackground, 1);
 factory->AddBackgroundTree(WZbackground , 1);
 factory->AddBackgroundTree(ZZbackground , 1);
 factory->AddBackgroundTree(ZGAMMAbackground , 1);

  // Set xs-weight
//  factory->SetSignalWeightExpression    ("weight");
  //factory->SetBackgroundWeightExpression("weight");
  
  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycut = "";
 // TCut mycutb = ""; 
  
  Int_t nSignalTrain = signalTrain->GetEntries();
  
  Int_t nBackTrain  = backgroundTrain->GetEntries();
  factory->PrepareTrainingAndTestTree( "","", "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:SplitSeed=88!V" );

  
//factory->PrepareTrainingAndTestTree( "", "",
//				       "nTrain_Signal=100:nTrain_Background=100:nTest_Signal=100:nTest_Background=100:!V" );

//  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
				       //":NSigTrain=:NBkgTrain=:NSigTest=:NBkgTest=:SplitMode=Alternate:!V" );  
				       //":nTrain_Signal=10000:nTest_Signal=2260:nTrain_Background=100000:nTest_Background=100000:SplitMode=Alternate:!V" );  
//				       ":nTrain_Signal=nSignalTrain:nTrain_Background=nBackTrain:!V");
  //":nTrain_Signal=nSignalTrain:nTest_Signal=nSignalTest:nTrain_Background=nBackTrain:nTest_Background=nBackTest:SplitMode=Block:!V");
  //":SplitMode=Alternate:!V");
  //":nTrain_Signal=:nTest_Signal=:nTrain_Background=:nTest_Background=:SplitMode=Block:!V");
  //":SplitMode=Alternate:!V" );  
  
  // ---- Book MVA methods
//  factory->BookMethod( TMVA::Types::kBDT, "BDT",
//		       "!H:!V:NTrees=300:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=-1:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=4" );
//factory->BookMethod( TMVA::Types::kBDT, "BDT",
 //                      "!H:!V:NTrees=400:BoostType=AdaBoost:AdaBoostBeta=0.5:MaxDepth=3:SeparationType=GiniIndex:PruneMethod=NoPruning" );
factory->BookMethod( TMVA::Types::kBDT, "BDT",
                       "!H:!V:NTrees=400:BoostType=AdaBoost:AdaBoostBeta=0.8:MaxDepth=3:SeparationType=GiniIndex:PruneMethod=NoPruning" );


//  factory->BookMethod( TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5");
 // factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher");   
  // ---- Now you can tell the factory to train, test, and evaluate the MVAs
// factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" ); 
// factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001" );
// factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit","H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" ); 

  
// Train MVAs using the set of training events
  factory->TrainAllMethods();
  
   // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  // --------------------------------------------------------------
  
  // Save the output
  outputFile->Close();
  
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  
  delete factory;
  
  // Launch the GUI for the root macros
//  if (!gROOT->IsBatch()) TMVAGui( outfileName );
}

//SetSignalWeightExpression

