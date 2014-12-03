#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <vector>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "myanalysis/Atq/interface/bjetcollection.h"

bjetcollection::bjetcollection(const edm::ParameterSet & iConfig){ 
  theJetLabel_ = iConfig.getParameter<edm::InputTag>("JetCollection");
  inputMETs_  =iConfig.getParameter<edm::InputTag>("inputMETs");
  JESscaleType_   =iConfig.getParameter<std::string>  ("JESscaleType"  );
  JERscaleType_   =iConfig.getParameter<std::string>  ("JERscaleType"  );
  resolutionFactorabs_ = iConfig.getParameter<std::vector<double> > ("resolutionFactors"   );
  resolutionFactorup_ = iConfig.getParameter<std::vector<double> > ("resolutionFactorsup"   );
  resolutionFactordown_ = iConfig.getParameter<std::vector<double> > ("resolutionFactorsdown"   );
  resolutionRanges_ = iConfig.getParameter<std::vector<double> > ("resolutionEtaRanges" );
  JECUncSrcFile_  =iConfig.getParameter<edm::FileInPath>("JECUncSrcFile"); 
  resolutionEnergyRanges_ = iConfig.getParameter<std::vector<double> > ("resolutionEnergyRanges" );
  sigmaMC_ = iConfig.getParameter<std::vector<double> > ("sigmaMC" );
  minJetPt_    = iConfig.getParameter<double>("MinJetPt");
  maxJetEta_   = iConfig.getParameter<double>("MaxJetEta");
  bdis_   = iConfig.getParameter<double>("bdis");
  ifrealdata_ = iConfig.getParameter<bool>("realdata");

/*
  for(int i = 0; i < 5; i++){ //loop over energy for det d
		  //resize vector for eta range of det d
		  myparam.resize(5);
		  
		  for(int j = 0; j < 5; j++){ //loop over eta for det d
		    //fill in parameters vector from python
			myparam[i][j] = sigmaMC_[i*5 + j];
		  }
		}
*/

  allowedTypes_.push_back(std::string("jes:up"));
  allowedTypes_.push_back(std::string("jes:down"));
  allowedTypes_.push_back(std::string("abs"));
  allowedTypes_.push_back(std::string("jer:up"));
  allowedTypes_.push_back(std::string("jer:down"));

  outputJets_ = theJetLabel_.label();
  outputMETs_ = inputMETs_.label(); 
  produces<std::vector<pat::MET> >(""); 
  produces<std::vector<pat::Jet> >("");
}

 void bjetcollection::beginJob()
 { 
   // check if scaleType is ok
   if(std::find(allowedTypes_.begin(), allowedTypes_.end(), JESscaleType_)==allowedTypes_.end()){
     edm::LogError msg("JetEnergyScale"); 
     msg << "Unknown scaleType: " << JESscaleType_ << " allowed types are: \n";
     for(std::vector<std::string>::const_iterator type=allowedTypes_.begin(); type!=allowedTypes_.end(); ++type){
     msg << *type << "\n";
     }
     msg << "Please modify your configuration accordingly \n";
     throw cms::Exception("Configuration Error");
   }             
 }

bjetcollection::~bjetcollection() {
}


void bjetcollection::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup) {
  // read in the objects
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(theJetLabel_, jets);
  // access MET
  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(inputMETs_, mets);

   // keep differences for met rescaling
   double dPx = 0., dPy = 0.;
//std::auto_ptr<std::vector<pat::MET> > pMETs(new std::vector<pat::MET>);
//std::auto_ptr<std::vector<pat::Jet> > mypatbJets(new std::vector<pat::Jet>);

  std::vector<pat::Jet> * mypatbJets = new std::vector<pat::Jet>(); 
  std::vector<pat::MET> * pMETs=new std::vector<pat::MET>;
  for (std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it) {
  pat::Jet scaledJet = *it;
if (it->pt()>10){
double jerScaleFactor=1;
if (!ifrealdata_){
//std::cout<<"bjet pt = "<<(*it).pt()<<std::endl;
  // JER scaled for all possible methods
  if(JERscaleType_.substr(JERscaleType_.find(':')+1)=="up"){
  jerScaleFactor = resolutionFactor(scaledJet,resolutionFactorup_);}
  else if(JERscaleType_.substr(JERscaleType_.find(':')+1)=="down"){
   jerScaleFactor = resolutionFactor(scaledJet,resolutionFactordown_);}
  else jerScaleFactor = resolutionFactor(scaledJet,resolutionFactorabs_);
}

  scaleJetEnergy( scaledJet, jerScaleFactor );

       // get the uncertainty parameters from file, see
       // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
       JetCorrectorParameters* param = new JetCorrectorParameters(JECUncSrcFile_.fullPath(), "");
       // instantiate the jec uncertainty object
       JetCorrectionUncertainty* deltaJEC = new JetCorrectionUncertainty(*param);
       deltaJEC->setJetEta(it->eta()); deltaJEC->setJetPt(it->pt()); 



      if(JESscaleType_.substr(JESscaleType_.find(':')+1)=="up"){
      // JetMET JES uncertainty
      float jetMet = deltaJEC->getUncertainty(true);
      scaleJetEnergy( scaledJet, 1+jetMet );
       }
       else if(JESscaleType_.substr(JESscaleType_.find(':')+1)=="down"){
        // JetMET JES uncertainty
        float jetMet = deltaJEC->getUncertainty(true);
         scaleJetEnergy( scaledJet, 1-jetMet*jetMet );
       }
       dPx    += scaledJet.px() - it->px();
       dPy    += scaledJet.py() - it->py();
       delete deltaJEC;
       delete param;

//    if (scaledJet.pt() > minJetPt_ && fabs(scaledJet.eta()) < maxJetEta_ && scaledJet.bDiscriminator("combinedSecondaryVertexBJetTags")> bdis_) {
      mypatbJets->push_back(scaledJet);
//std::cout<<"bjet pt = "<<(*it).pt()<<std::endl;
//std::cout<<"nbjet  = "<<mypatbJets->size()<<std::endl;
}


//      mypatbJets->push_back(scaledJet);

//    }
  }
   // scale MET accordingly
   pat::MET met = *(mets->begin());
   //std::cout<<"met before: "<<met.pt()<<std::endl;
   double scaledMETPx = met.px() - dPx;
   double scaledMETPy = met.py() - dPy;
met.setP4(reco::MET::LorentzVector(scaledMETPx, scaledMETPy, 0, sqrt(scaledMETPx*scaledMETPx+scaledMETPy*scaledMETPy)));

//   std::cout<<"met after: "<<met.pt()<<std::endl;
   pMETs->push_back( met );
//   std::cout<<"size met : "<<pMETs->size()<<std::endl;

//   std::cout<<"size bjet00 : "<<mypatbJets->size()<<std::endl;

  std::auto_ptr<std::vector<pat::Jet> > myJets(mypatbJets);
 iEvent.put(myJets);
//     std::cout<<"size bjet : "<<mypatbJets->size()<<std::endl;
std::auto_ptr<std::vector<pat::MET> > fmet(pMETs);
  iEvent.put(fmet);
//     std::cout<<"size bjet : "<<mypatbJets->size()<<std::endl;
//     std::cout<"size jet : "<<jets->size()<<std::endl;

}
double bjetcollection::resolutionFactor(const pat::Jet& jet , std::vector<double> JERfactors)
 {
if (jet.pt()<50 || jet.pt()>500) return 1.; 
   std::vector<double> resolutionFactor_=JERfactors;
    int mypt=0;
    int myeta=0; 
   for(unsigned int numberOfJERvariation=0; numberOfJERvariation<resolutionFactor_.size(); ++numberOfJERvariation){
     int etaMin = 2*numberOfJERvariation;
     int etaMax = etaMin+1;
     if(std::abs(jet.eta())>=resolutionRanges_[etaMin] && std::abs(jet.eta())<resolutionRanges_[etaMax])myeta=numberOfJERvariation;}

   for(unsigned int JERenergy=0; JERenergy<resolutionFactor_.size(); ++JERenergy){
     int Emin = 2*JERenergy;
     int Emax =Emin+1; 
     if(std::abs(jet.pt())>=resolutionEnergyRanges_[Emin]&&std::abs(jet.pt())<resolutionEnergyRanges_[Emax]){mypt=JERenergy;}}
     if (resolutionFactor_[myeta]<1) return 1.; 
     double modifiedResolution = std::sqrt(std::pow(resolutionFactor_[myeta],2)-1)*sigmaMC_[mypt*5+myeta];
     double factor=gRandom->Gaus(jet.pt(),modifiedResolution)/jet.pt(); 


/*   if(!jet.genJet()) { return 1.; }
   std::vector<double> resolutionFactor_=JERfactors;
   // check if vectors are filled properly
   if((2*resolutionFactor_.size())!=resolutionRanges_.size()){
     // eta range==infinity
   if(resolutionFactor_.size()==resolutionRanges_.size()&&resolutionRanges_.size()==1&&resolutionRanges_[0]==-1.){
    resolutionRanges_[0]=0;
    resolutionRanges_.push_back(-1.);
    }
    // others
    else{
       edm::LogError msg("JetEnergyResolution");
       msg << "\n resolutionEtaRanges or resolutionFactors in module JetEnergyScale not filled properly.\n";
       msg << "\n resolutionEtaRanges needs a min and max value for each entry in resolutionFactors.\n";
       throw cms::Exception("invalidVectorFilling");
     }
   }
   // calculate eta dependend JER factor
   double modifiedResolution = 1.;
   for(unsigned int numberOfJERvariation=0; numberOfJERvariation<resolutionFactor_.size(); ++numberOfJERvariation){
     int etaMin = 2*numberOfJERvariation;
     int etaMax = etaMin+1;
     if(std::abs(jet.eta())>=resolutionRanges_[etaMin]&&(std::abs(jet.eta())<resolutionRanges_[etaMax]||resolutionRanges_[etaMax]==-1.)){
     modifiedResolution*=resolutionFactor_[numberOfJERvariation];
     // take care of negative scale factors 
     if(resolutionFactor_[numberOfJERvariation]<0){
      edm::LogError msg("JetEnergyResolution");
      msg << "\n chosen scale factor " << resolutionFactor_[numberOfJERvariation] << " is not valid, must be positive.\n";
      throw cms::Exception("negJERscaleFactors");
      }
    }
   }
   // calculate pt smearing factor
//   double factor = (jet.genJet()->pt() - modifiedResolution*(jet.pt() - jet.genJet()->pt()))/jet.pt();
   double factor = (5 - modifiedResolution*(jet.pt() - 7))/jet.pt();
*/
   return (factor<0 ? 0. : factor);
 }
void bjetcollection::scaleJetEnergy(pat::Jet& jet, double factor)
 { jet.scaleEnergy( factor );}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(bjetcollection);

