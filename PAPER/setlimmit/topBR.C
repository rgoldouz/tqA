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

void topBR(){

double mt;
double mw=80.419;
double alpha=1/127.9;
double sw=0.234;
double gammat;
double gammaFCNC;
mt=178;
gammat = ((alpha * pow(mt,3))/(16* sw*pow(mw,2)) )*(1-3*pow(mw/mt,4)+2*pow(mw/mt,6));
gammaFCNC= (alpha/2)*mt;
cout<<"mt = 178 gammat=    "<< gammat<<endl;
cout<<"mt = 178  gammaFCNC=   "<< gammaFCNC<<endl;
cout<<"mt = 178  BR=   "<< gammaFCNC/gammat<<endl;

mt=170;
gammat = ((alpha * pow(mt,3))/(16* sw*pow(mw,2)) )*(1-3*pow(mw/mt,4)+2*pow(mw/mt,6));
gammaFCNC= (alpha/2)*mt;
cout<<"mt = 170  gammat=   "<< gammat<<endl;
cout<<"mt = 170  gammaFCNC=   "<< gammaFCNC<<endl;
cout<<"mt = 170  BR=   "<< gammaFCNC/gammat<<endl;

mt=180;
gammat = ((alpha * pow(mt,3))/(16* sw*pow(mw,2)) )*(1-3*pow(mw/mt,4)+2*pow(mw/mt,6));
gammaFCNC= (alpha/2)*mt;
cout<<"mt = 180 gammat=    "<< gammat<<endl;
cout<<"mt = 180 gammaFCNC=    "<< gammaFCNC<<endl;
cout<<"mt = 180  BR=   "<< gammaFCNC/gammat<<endl;

mt=172.5;
gammat = ((alpha * pow(mt,3))/(16* sw*pow(mw,2)) )*(1-3*pow(mw/mt,4)+2*pow(mw/mt,6));
gammaFCNC= (alpha/2)*(0.4444)*mt;
cout<<"mt = 175 gammat=    "<< gammat<<endl;
cout<<"mt = 175 gammaFCNC=    "<< gammaFCNC<<endl;
cout<<"mt = 175  BR=   "<< gammaFCNC/gammat<<endl;


}




