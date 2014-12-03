#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
//#include "DataFormats/Common/interface/View.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"




class ElectronSelector : public edm::EDFilter {

  public:

    explicit ElectronSelector(const edm::ParameterSet & iConfig);
    ~ElectronSelector();
  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    edm::InputTag electronSrc_;
    double minElePt_, maxEleEta_;
};


ElectronSelector::ElectronSelector(const edm::ParameterSet & iConfig) {
  electronSrc_   = iConfig.getParameter<edm::InputTag>("ElectronSource");
  minElePt_      = iConfig.getParameter<double>("MinElePt");
  maxEleEta_     = iConfig.getParameter<double>("MaxEleEta");
  produces<std::vector<reco::GsfElectron>>("");
}
ElectronSelector::~ElectronSelector() {
}
bool ElectronSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // electrons
  edm::Handle< std::vector<reco::GsfElectron> > electrons;
  iEvent.getByLabel(electronSrc_, electrons);

  // check which ones to keep
  std::auto_ptr<std::vector<reco::GsfElectron> > prod(new std::vector<reco::GsfElectron>());

  for(unsigned int i = 0; i < electrons->size(); ++i) {
    reco::GsfElectronRef ele(electrons, i);

    if (ele->pt() > minElePt_ & ele->eta()<maxEleEta_) {
prod->push_back(reco::GsfElectron(*ele));
}
}

  iEvent.put(prod);
//bool result = ((prod->size() == 0) ? true :false);
bool result =true;
  return result;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ElectronSelector);


