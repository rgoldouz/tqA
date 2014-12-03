#ifndef scalefactorcalc_H 
#define scalefactorcalc_H 

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <math.h>
#include <utility>

//ROOT includes
#include <TH1.h>
#include <TH2.h>
#include "TF1.h"

using namespace std;

class ScaleFactor{
 public:

double photonScaleFactor(double ph_eta);
double photonScaleFactorup(double ph_eta);
double photonScaleFactordown(double ph_eta);

//double elevetoScaleFactor(double ph_eta);
//double elevetoScaleFactorup(double ph_eta);
//double elevetoScaleFactordown(double ph_eta);
 
double muonScaleFactor(double muon_eta);
double muonScaleFactorup(double muon_eta);
double muonScaleFactordown(double muon_eta);

double btagScaleFactor(double b_pt);

double triggerScaleFactor(double mu_eta);
double triggerScaleFactorup(double mu_eta);
double triggerScaleFactordown(double mu_eta);

double wphjetkfactor(double ph_pt);
double zphjetkfactor(double ph_pt);

private:
TF1  *CSVM_SFb_0to2p4;
double phsf;
// double PtBins[]={20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
double jetPt;
double musf;
double trigmusf;
};



#endif
