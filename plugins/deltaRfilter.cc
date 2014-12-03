#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include <TMath.h>
#include <TLorentzVector.h>
#include <memory>

class DeltaRfilter : public edm::EDFilter {
  public:

    explicit DeltaRfilter(const edm::ParameterSet & iConfig);
    ~DeltaRfilter() {}

  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    edm::InputTag topSrc_;
    edm::InputTag muonSrc_;
    edm::InputTag photonSrc_;
    double mindr_;
};

DeltaRfilter::DeltaRfilter(const edm::ParameterSet & iConfig) {
  topSrc_ = iConfig.getParameter<edm::InputTag>("TopSource");
  muonSrc_ = iConfig.getParameter<edm::InputTag>("muonSource");
  photonSrc_= iConfig.getParameter<edm::InputTag>("photonSource");
  mindr_  = iConfig.getParameter<double>("Mindr");
}

bool DeltaRfilter::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  edm::Handle<std::vector< reco::NamedCompositeCandidate >> rectop;
  iEvent.getByLabel(topSrc_, rectop);

  edm::Handle<edm::View<pat::Photon> > Ph    ;
  iEvent.getByLabel(photonSrc_, Ph  );
  
  edm::Handle<std::vector<pat::Muon> > muonHandle;
  iEvent.getByLabel(muonSrc_, muonHandle);

TLorentzVector bjet4vector;
TLorentzVector photon4vector;
TLorentzVector muon4vector;
double drphmu;
double drphjet;

       for(edm::View<pat::Photon>::const_iterator photon=Ph->begin(); photon!=Ph->end(); ++photon )
        {
	photon4vector.SetPxPyPzE(photon->px() ,photon->py() , photon->pz () , photon-> energy());
	}

        for(std::vector<pat::Muon>::const_iterator mu1=muonHandle->begin(); mu1!=muonHandle->end(); ++mu1)
        {
        muon4vector.SetPxPyPzE(mu1->px() ,mu1->py() ,mu1->pz () , mu1-> energy());
	}

	for(std::vector< reco::NamedCompositeCandidate >::const_iterator top1=rectop->begin(); top1!=rectop->end(); ++top1)
        {
	bjet4vector.SetPxPyPzE(top1->daughter("BJet")->px(),top1->daughter("BJet")->py(),top1->daughter("BJet")->pz(), top1->daughter("BJet")->energy());
	}
drphmu=abs(deltaR(photon4vector.Eta(),photon4vector.Phi(),muon4vector.Eta(),muon4vector.Phi()));
drphjet=abs(deltaR(photon4vector.Eta(),photon4vector.Phi(),bjet4vector.Eta(),bjet4vector.Phi()));

return (drphmu>mindr_ && drphjet>mindr_);
}
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DeltaRfilter);
