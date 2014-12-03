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

#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include <vector>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "CommonTools/Utils/interface/EtComparator.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/UserData.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "CommonTools/Utils/interface/EtComparator.h"

#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

class photoncollection : public edm::EDProducer {

  public:

    explicit photoncollection (const edm::ParameterSet & iConfig);
    ~photoncollection();

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
PFIsolationEstimator isolator;

    edm::InputTag theJetLabel_;
  std::vector<edm::InputTag> floatLabels_ ;
  std::vector<std::string>   floatNames_ ;
  bool useUserData_;
  bool useAlternateIsolations_;
  pat::PATUserDataHelper<pat::Photon>      userDataHelper_;
    bool ifrealdata_;
    std::vector<std::string> allowedTypes_;
    std::string photonESscaleType_;


};
photoncollection::photoncollection(const edm::ParameterSet & iConfig){ 
  theJetLabel_ = iConfig.getParameter< edm::InputTag >  ( "photonCollection" );
  useAlternateIsolations_ = iConfig.getParameter< bool >  ( "useAlternateIsolations" );
  useUserData_   = iConfig.exists ("userData");
  floatLabels_   = iConfig.getParameter< std::vector< edm::InputTag >> ( "floatLabels" );
  floatNames_   = iConfig.getParameter< std::vector<std::string> >  ( "floatNames" );
  photonESscaleType_   =iConfig.getParameter<std::string>  ("photonESscaleType"  );
  ifrealdata_ = iConfig.getParameter<bool>("realdata");

  allowedTypes_.push_back(std::string("phes:up"));
  allowedTypes_.push_back(std::string("phes:down"));
  allowedTypes_.push_back(std::string("abs"));


  if ( useUserData_ ) {
    userDataHelper_ = iConfig.getParameter<edm::ParameterSet>("userData");
  }

produces<std::vector<pat::Photon> >();
 isolator.initializePhotonIsolation(kTRUE);
 isolator. setConeSize(0.3);
}
photoncollection::~photoncollection() {
}


void photoncollection::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup ) {
using namespace pat;  
// read in the objects
  edm::Handle<std::vector<pat::Photon> > photons;
  iEvent.getByLabel(theJetLabel_, photons);
  edm::Handle<reco::BeamSpot> bs;
  iEvent.getByLabel("offlineBeamSpot", bs);
  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByLabel("conversions",conversions);
  std::vector<pat::Photon> * myph = new std::vector<pat::Photon>();
  edm::Handle<reco::GsfElectronCollection> gsfElecs;
  iEvent.getByLabel("gsfElectrons", gsfElecs);
  edm::Handle<reco::PFCandidateCollection> pfCandidates;
  iEvent.getByLabel("particleFlow",pfCandidates);
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel("offlinePrimaryVertices", vertices);
  if ( floatLabels_.size()!=floatNames_.size()) {
    std::cout<<"mismatch between supplied floats and names"<<std::endl;
    std::cout<<"floatLabels_.size()="<<floatLabels_.size()<<std::endl;
    std::cout<<"floatNames_.size()="<<floatNames_.size()  <<std::endl;
    return;
  }
  std::vector<double> myFloats;
  myFloats.reserve(floatLabels_.size());
  std::vector<edm::InputTag>::const_iterator tag = floatLabels_.begin();
  for (; tag != floatLabels_.end(); ++tag) {
    edm::Handle<double> floatVal;
    iEvent.getByLabel(*tag,floatVal);
    myFloats.push_back(*floatVal);
  }

  for (std::vector<pat::Photon>::const_iterator it = photons->begin(); it != photons->end(); ++it) {
    pat::Photon aPhoton(*it);
  if (useAlternateIsolations_) {
      reco::VertexRef vtx(vertices,0);
      isolator.fGetIsolation(&(*it), pfCandidates.product(), vtx, vertices);
      const reco::BeamSpot &beamspot = *bs.product();
      bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(it->superCluster(), gsfElecs, conversions, beamspot.position());
      aPhoton.addUserInt("passEleConvVeto",(int)passelectronveto);
      aPhoton.addUserFloat("pfChaIsoAlt",isolator.getIsolationCharged());
      aPhoton.addUserFloat("pfNeuIsoAlt",isolator.getIsolationNeutral());
      aPhoton.addUserFloat("pfGamIsoAlt",  isolator.getIsolationPhoton());
    
    std::vector<std::string>::const_iterator name = floatNames_.begin();
    for (std::vector<double>::const_iterator val = myFloats.begin(); val != myFloats.end(); ++val) {
      aPhoton.addUserFloat(*name,*val);
      ++name;
    }
}
    if ( useUserData_ ) {
      userDataHelper_.add( aPhoton, iEvent,  iEventSetup );
    }
double phScaleFactor=1;
if (!ifrealdata_){
//std::cout<<"bjet pt = "<<(*it).pt()<<std::endl;
  // JER scaled for all possible methods
  if(photonESscaleType_.substr(photonESscaleType_.find(':')+1)=="up"){
  if(aPhoton.eta()<1.4442)phScaleFactor = 1.01;
  else phScaleFactor = 1.03;
}
  else if(photonESscaleType_.substr(photonESscaleType_.find(':')+1)=="down"){
  if(aPhoton.eta()<1.4442) phScaleFactor = 0.99;
  else phScaleFactor = 0.97;
}
}
aPhoton.setP4(aPhoton.p4()* phScaleFactor );


//	if (abs(it->eta()) < 1.4442){
//		if (it->pt() > 50 && it->hadronicOverEm() < 0.05  && passelectronveto && isolator.getIsolationCharged()< 0.7 && isolator.getIsolationNeutral()< (0.4 + 0.04*it->pt()) && isolator.getIsolationPhoton() < (0.5 + 0.005*it->pt()) ) {
//		myph->push_back(*it);
//		}
//	}
//        if (abs(it->eta()) > 1.566 && abs(it->eta())<2.5){
//	{
//		if (it->pt() > 50 && it->hadronicOverEm() < 0.05 && passelectronveto && isolator.getIsolationCharged()< 0.5 && isolator.getIsolationNeutral()< (1.5 + 0.04*it->pt()) && isolator.getIsolationPhoton() < (1 + 0.005*it->pt()) ) {
            myph->push_back(aPhoton);
    //        }
//	}
//  }
}
  std::auto_ptr<std::vector<pat::Photon> > mypho(myph);
  iEvent.put(mypho);
}
void photoncollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(photoncollection);


