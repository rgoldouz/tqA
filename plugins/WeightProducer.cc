// system include files
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TH1.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//
// class decleration
//
class WeightProducer: public edm::EDProducer {
public:
  explicit WeightProducer(const edm::ParameterSet&);
  ~WeightProducer();
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  edm::LumiReWeighting generateWeights(const TH1* data_npu_estimated ) const;
 edm::LumiReWeighting LumiWeights_;
};
WeightProducer::WeightProducer(const edm::ParameterSet& iConfig) 
{
std::string fileNamePU = iConfig.getParameter<std::string> ("FileNamePUDataDistribution");
      edm::FileInPath filePUDataDistr(fileNamePU);
      TFile file(filePUDataDistr.fullPath().c_str(), "READ");
      TH1 *h = 0;
      file.GetObject("pileup", h);
      if (h) {
         h->SetDirectory(0);
      } else {
         std::cerr << "ERROR in WeightProducer: Histogram 'pileup' does not exist in file '"
               << filePUDataDistr.fullPath() << "'\n.";
         std::cerr
               << "See https://twiki.cern.ch/twiki/bin/view/CMS/HamburgWikiAnalysisCalibration#Pile_Up_Reweighting for available input distributions."
               << std::endl;
         exit(1);
      }
      file.Close();
     LumiWeights_= generateWeights(h);
     delete h;
   //register your products
   produces<double> ("");
}
WeightProducer::~WeightProducer() {
}

// ------------ method called to produce the data  ------------
void WeightProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
   double resultWeight = 1.;
   // Optionally, multiply PU weight
   edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
   iEvent.getByLabel("addPileupInfo", puInfo);
   float npu = 0;
      if (puInfo.isValid()) {
         std::vector<PileupSummaryInfo>::const_iterator puIt;
         int n = 0;
         for (puIt = puInfo->begin(); puIt != puInfo->end(); ++puIt, ++n) {
            if (puIt->getBunchCrossing() == 0) {
               npu = puIt->getTrueNumInteractions();
               break;
            }
         }
         resultWeight = LumiWeights_.weight( npu );

      }

   // put weight into the Event
   std::auto_ptr<double> pOut(new double(resultWeight));
   iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
//void
//WeightProducer::beginJob()
//{
//}

// ------------ method called once each job just after ending the event loop  ------------
void WeightProducer::endJob() {
}

// Get weight factor dependent on number of
// added PU interactions
// --------------------------------------------------

edm::LumiReWeighting WeightProducer::generateWeights(const TH1* data_npu_estimated ) const {
  unsigned int nMaxPU = 0;

   nMaxPU = 60;
std::vector< float > TrueDisSummer12;
float npuSummer12_S10[60] = {
      2.560E-06,
      5.239E-06,
      1.420E-05,
      5.005E-05,
      1.001E-04,
      2.705E-04,
      1.999E-03,
      6.097E-03,
      1.046E-02,
      1.383E-02,
      1.685E-02,
      2.055E-02,
      2.572E-02,
      3.262E-02,
      4.121E-02,
      4.977E-02,
      5.539E-02,
      5.725E-02,
      5.607E-02,
      5.312E-02,
      5.008E-02,
      4.763E-02,
      4.558E-02,
      4.363E-02,
      4.159E-02,
      3.933E-02,
      3.681E-02,
      3.406E-02,
      3.116E-02,
      2.818E-02,
      2.519E-02,
      2.226E-02,
      1.946E-02,
      1.682E-02,
      1.437E-02,
      1.215E-02,
      1.016E-02,
      8.400E-03,
      6.873E-03,
      5.564E-03,
      4.457E-03,
      3.533E-03,
      2.772E-03,
      2.154E-03,
      1.656E-03,
      1.261E-03,
      9.513E-04,
      7.107E-04,
      5.259E-04,
      3.856E-04,
      2.801E-04,
      2.017E-04,
      1.439E-04,
      1.017E-04,
      7.126E-05,
      4.948E-05,
      3.405E-05,
      2.322E-05,
      1.570E-05,
      5.005E-06};
    
  std::vector<float> result;
  double s = 0.0;
  for(unsigned int npu = 0; npu < nMaxPU; ++npu) {
    float npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
//    result[npu] = npu_estimated / npuProbs[npu];
//    s += npu_estimated;
result.push_back(npu_estimated);
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for (unsigned int npu = 0; npu < nMaxPU; ++npu) {
//    result[npu] /= s;
TrueDisSummer12.push_back(npuSummer12_S10[npu]);
  }
 edm::LumiReWeighting LW= edm::LumiReWeighting(TrueDisSummer12,result);

  return LW ;
}
#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE( WeightProducer);
