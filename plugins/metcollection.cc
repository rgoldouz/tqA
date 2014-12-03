#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <vector>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"


class metcollection : public edm::EDProducer {

  public:

    explicit metcollection (const edm::ParameterSet & iConfig);
    ~metcollection();
    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:

    edm::InputTag theJetLabel_;
    double minJetPt_;
};
metcollection::metcollection(const edm::ParameterSet & iConfig) {
  theJetLabel_ = iConfig.getParameter<edm::InputTag>("JetCollection");
  minJetPt_    = iConfig.getParameter<double>("MinJetPt");

produces<std::vector<pat::MET> >();
}

metcollection::~metcollection() {
}


void metcollection::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup) {
  // read in the objects
  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(theJetLabel_, mets);
  std::vector<pat::MET> * mypatmets = new std::vector<pat::MET>(); 

  for (std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it) {
    if (it->pt() > minJetPt_ ) {
      mypatmets->push_back(*it);
    }
  }
  std::auto_ptr<std::vector<pat::MET> > mymets(mypatmets);
  iEvent.put(mymets);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(metcollection);

