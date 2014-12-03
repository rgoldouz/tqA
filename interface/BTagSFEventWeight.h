#ifndef BTagSFEventWeight_h
#define BTagSFEventWeight_h

#include <memory>
#include <string>
#include <iostream>
 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TString.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
 
#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "FWCore/Framework/interface/ESHandle.h"

 
class BTagSFEventWeight : public edm::EDProducer {
 
  public:
   explicit BTagSFEventWeight(const edm::ParameterSet&);
   ~BTagSFEventWeight();
   
  private:
   virtual void produce(edm::Event&, const edm::EventSetup&);
 
  private:
   edm::InputTag jets_;
   std::string sysVar_;
   edm::FileInPath filename_;
   double maxPtDB_;
   double maxPt11004_;
   double maxPt2012_;
   double maxPtMisTag_;
   double maxEta_;
   std::map<std::string, TH2F*> effHists_;
   TFile * file_;   
   double effBTagSF2012(double);
   double effBTagSFerr2012(double);
   double effMisTagSF2012(double, double, TString);
   double effBTag    (double, double);
   double effBTagSF  (double, double, bool);
   double effBTagCjet(double, double);
   double effMisTag  (double, double);
   double effMisTagSF(double, double);
   double effBTagEvent(std::vector<double> &, std::vector<double> &);
 };
 
#endif
