#ifndef bjetcollection_h
#define bjetcollection_h
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"

class bjetcollection : public edm::EDProducer {
    
  public:
    
    explicit bjetcollection (const edm::ParameterSet & iConfig);
    ~bjetcollection();
    
    virtual void beginJob();
    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
    double resolutionFactor(const pat::Jet&, std::vector<double> );
    void scaleJetEnergy(pat::Jet&, double);
  private:

    edm::InputTag theJetLabel_;
    edm::InputTag inputMETs_;
    std::vector<double> resolutionFactorabs_;
    std::vector<double> resolutionFactorup_;
    std::vector<double> resolutionFactordown_;
    std::vector<double> resolutionRanges_;
    std::vector<double> resolutionEnergyRanges_;
    std::vector<double> sigmaMC_;
    edm::FileInPath JECUncSrcFile_;
    std::vector<std::string> allowedTypes_;
    double minJetPt_;
    double maxJetEta_;
    double bdis_;
    std::string JESscaleType_;
    std::string JERscaleType_;
    std::string outputMETs_;
    std::string outputJets_;
    edm::InputTag inputJets_;
    bool ifrealdata_;

//    edm::InputTag inputMETs_;
//  std::vector<std::vector<double>> myparam;

};


#endif
