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
#include <cmath> 

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

void errorestimation() {

double wjetnominal=302.91;
double wgammajetnominal=985.18;

double errorwjet=0;
double errorwgammajet=0;

errorwjet=sqrt(pow(90.86,2)+pow(wjetnominal-285.79,2)+pow(wjetnominal-299.52,2)+pow(wjetnominal*0.03,2)+pow(wjetnominal*0.02,2)+pow(wjetnominal*0.01,2)+pow(wjetnominal*0.04,2)+pow(wjetnominal*0.01,2));

errorwgammajet=sqrt(pow(159.16,2)+pow(wgammajetnominal-977.60,2)+pow(wgammajetnominal-953.51,2)+pow(wgammajetnominal*0.02,2)+pow(wgammajetnominal*0.01,2)+pow(wgammajetnominal*0.01,2)+pow(wgammajetnominal*0.01,2)+pow(wgammajetnominal*0.04,2));

cout<<"errorwjet     "<<errorwjet<<"               "<<(errorwjet/wjetnominal)*100<<"  %"<<endl;
cout<<"errorwgammajet     "<<errorwgammajet<<"               "<<(errorwgammajet/wgammajetnominal)*100<<"  %"<<endl;

}










