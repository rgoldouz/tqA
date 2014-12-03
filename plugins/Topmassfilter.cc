#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include <TMath.h>
#include <memory>


class TopFilter : public edm::EDFilter {

  public:

    explicit TopFilter(const edm::ParameterSet & iConfig);
    ~TopFilter() {}

  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);

    edm::InputTag topSrc_;
    double minmass_;
    double maxmass_;
    bool istopwindow_; 
};


TopFilter::TopFilter(const edm::ParameterSet & iConfig) {
  topSrc_ = iConfig.getParameter<edm::InputTag>("TopSource");
  minmass_  = iConfig.getParameter<double>("Minmass");
  maxmass_  = iConfig.getParameter<double>("Maxmass");
  istopwindow_ = iConfig.getParameter<bool>("istopwindow");
}


bool TopFilter::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  edm::Handle<std::vector< reco::NamedCompositeCandidate >> rectop;
  iEvent.getByLabel(topSrc_, rectop);
 double sss;
for(std::vector< reco::NamedCompositeCandidate >::const_iterator top1=rectop->begin(); top1!=rectop->end(); ++top1)
        {
	sss=sqrt( TMath::Power(top1->daughter("W")->energy()+top1->daughter("BJet")->energy(),2)-(TMath::Power(top1->daughter("W")->px()+top1->daughter("BJet")->px(),2)+TMath::Power(top1->daughter("W")->py()+top1->daughter("BJet")->py(),2)+TMath::Power(top1->daughter("W")->pz()+top1->daughter("BJet")->pz(),2))); 
}
if (istopwindow_){
  return (sss > minmass_ & sss < maxmass_);}
else{return (sss < minmass_ || sss > maxmass_);}

}
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TopFilter);

