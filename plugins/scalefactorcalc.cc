#include "myanalysis/Atq/interface/scalefactorcalc.h"
#include <cmath>
#include <iostream>
double ScaleFactor::photonScaleFactor(double ph_eta)
{
if (0.000<=abs(ph_eta) && abs(ph_eta)<0.8) phsf=0.9749;
else if (0.8<=abs(ph_eta) && abs(ph_eta)<1.4442) phsf=0.9771;
else if (1.4442<=abs(ph_eta) && abs(ph_eta)<2.0) phsf=1.01092;
else if (2.0<=abs(ph_eta) && abs(ph_eta)<2.5) phsf=1.020396;
return phsf;
}

double ScaleFactor::photonScaleFactorup(double ph_eta)
{
if (0.000<=abs(ph_eta) && abs(ph_eta)<0.8) phsf=0.9749+0.0119;
else if (0.8<=abs(ph_eta) && abs(ph_eta)<1.4442) phsf=0.9771+0.0119;
else if (1.4442<=abs(ph_eta) && abs(ph_eta)<2.0) phsf=1.0109+0.0286 ;
else if (2.0<=abs(ph_eta) && abs(ph_eta)<2.5) phsf=1.0203+0.0288;
return phsf;
}

double ScaleFactor::photonScaleFactordown(double ph_eta)
{
if (0.000<=abs(ph_eta) && abs(ph_eta)<0.8) phsf=0.9749-0.0119;
else if (0.8<=abs(ph_eta) && abs(ph_eta)<1.4442) phsf=0.9771-0.0119;
else if (1.4442<=abs(ph_eta) && abs(ph_eta)<2.0) phsf=1.0109-0.0286 ;
else if (2.0<=abs(ph_eta) && abs(ph_eta)<2.5) phsf=1.0203-0.0288;
return phsf;
}
//////////////////////////////////////////////////////////////////////
/*
double ScaleFactor::elevetoScaleFactor(double ph_eta)
{
if (abs(ph_eta)<1.442) phsf=0.9933;
else if (1.4442<=abs(ph_eta)) phsf=1.0075;
return phsf;
}

double ScaleFactor::elevetoScaleFactorup(double ph_eta)
{
if (abs(ph_eta)<1.442) phsf=0.9933+0.0066;
else if (1.4442<=abs(ph_eta)) phsf=1.0075+0.02669;
return phsf;
}

double ScaleFactor::elevetoScaleFactordown(double ph_eta)
{
if (abs(ph_eta)<1.442) phsf=0.9933-0.0066;
else if (1.4442<=abs(ph_eta)) phsf=1.0075-0.02669;
return phsf;
}
*/


double ScaleFactor::btagScaleFactor(double b_pt)
{
//PtBins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
CSVM_SFb_0to2p4 = new TF1("CSVM_SFb_0to2p4","0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)))", 20.,800.);
jetPt = b_pt;
if ( jetPt>800 ) jetPt = 800;
else if ( jetPt<20 ) jetPt = 20;
return CSVM_SFb_0to2p4->Eval(jetPt);
}

double ScaleFactor::muonScaleFactor(double muon_eta)
{
if(0.000<=abs(muon_eta) && abs(muon_eta)<0.9) musf=0.9884;
else if (0.9<=abs(muon_eta) && abs(muon_eta)<1.2) musf=0.9806;
else if (1.2<=abs(muon_eta) && abs(muon_eta)<2.1) musf=0.9986;
return musf;
}

double ScaleFactor::muonScaleFactorup(double muon_eta)
{
if(0.000<=abs(muon_eta) && abs(muon_eta)<0.9) musf=0.9884+0.0532;
else if (0.9<=abs(muon_eta) && abs(muon_eta)<1.2) musf=0.9806+0.0527;
else if (1.2<=abs(muon_eta) && abs(muon_eta)<2.1) musf=0.9986+0.0537;
return musf;
}

double ScaleFactor::muonScaleFactordown(double muon_eta)
{
if(0.000<=abs(muon_eta) && abs(muon_eta)<0.9) musf=0.9884-0.0532;
else if (0.9<=abs(muon_eta) && abs(muon_eta)<1.2) musf=0.9806-0.0527;
else if (1.2<=abs(muon_eta) && abs(muon_eta)<2.1) musf=0.9986-0.0537;
return musf;
}


double ScaleFactor::triggerScaleFactor(double mu_eta)
{
if(0.000<=abs(mu_eta) && abs(mu_eta)<0.9) trigmusf=0.9815;
else if (0.9<=abs(mu_eta) && abs(mu_eta)<1.2) trigmusf=0.9651;
else if (1.2<=abs(mu_eta) && abs(mu_eta)<2.1) trigmusf=0.9968;
    return trigmusf;
}

double ScaleFactor::triggerScaleFactorup(double mu_eta)
{
if(0.000<=abs(mu_eta) && abs(mu_eta)<0.9) trigmusf=0.9815+0.00024;
else if (0.9<=abs(mu_eta) && abs(mu_eta)<1.2) trigmusf=0.9651+0.00068;
else if (1.2<=abs(mu_eta) && abs(mu_eta)<2.1) trigmusf=0.9968+0.0053;
    return trigmusf;
}

double ScaleFactor::triggerScaleFactordown(double mu_eta)
{
if(0.000<=abs(mu_eta) && abs(mu_eta)<0.9) trigmusf=0.9815-0.00024;
else if (0.9<=abs(mu_eta) && abs(mu_eta)<1.2) trigmusf=0.9651-0.00068;
else if (1.2<=abs(mu_eta) && abs(mu_eta)<2.1) trigmusf=0.9968-0.0053;
    return trigmusf;
}

double ScaleFactor::wphjetkfactor(double ph_pt){
return 1.46604+0.00923692*ph_pt-1.48871e-06*ph_pt*ph_pt;}
double ScaleFactor::zphjetkfactor(double ph_pt){
return 1.3292+0.000952237*ph_pt+2.00623e-05*ph_pt*ph_pt-1.41325e-07 *ph_pt*ph_pt*ph_pt+2.48614e-10*ph_pt*ph_pt*ph_pt*ph_pt;}



