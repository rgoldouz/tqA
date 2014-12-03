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
#include "myanalysis/Atq/interface/Fjetcollection.h"

/*
class Fjetcollection : public edm::EDProducer {

  public:

    explicit Fjetcollection (const edm::ParameterSet & iConfig);
    ~Fjetcollection();

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:

    edm::InputTag theJetLabel_;
    double minJetPt_;
    double maxJetEta_;
    double bdis_;

};
*/

Fjetcollection::Fjetcollection(const edm::ParameterSet & iConfig) {
  theJetLabel_ = iConfig.getParameter<edm::InputTag>("JetCollection");
  minJetPt_    = iConfig.getParameter<double>("MinJetPt");
  maxJetEta_   = iConfig.getParameter<double>("MaxJetEta");
  bscalefactorType_   =iConfig.getParameter<std::string>  ("bscalefactorType"  );
  bdis_   = iConfig.getParameter<double>("bdis");
  ifrealdata_ = iConfig.getParameter<bool>("realdata");



produces<std::vector<pat::Jet> >();
}

Fjetcollection::~Fjetcollection() {
}


void Fjetcollection::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup) {
  // read in the objects
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(theJetLabel_, jets);
  std::vector<pat::Jet> * mypatbJets = new std::vector<pat::Jet>(); 
  double correctedbtag=0;
  for (std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it) {
  if (ifrealdata_) correctedbtag= it->bDiscriminator("combinedSecondaryVertexBJetTags");
  if (!ifrealdata_ && it->genParticleRefs().size()!=1) correctedbtag= it->bDiscriminator("combinedSecondaryVertexBJetTags");
  if (!ifrealdata_ && it->genParticleRefs().size() ==1) correctedbtag= Bmodification(it->bDiscriminator("combinedSecondaryVertexBJetTags"),abs(it->genParticle(0)->pdgId()), bscalefactorType_);
//std::cout<< it->bDiscriminator("combinedSecondaryVertexBJetTags") <<"                "<< correctedbtag<<std::endl;
//std::cout<< bscalefactorType_<<ifrealdata_<<"        real btag   "<< it->bDiscriminator("combinedSecondaryVertexBJetTags") << "corrected    "<< correctedbtag<<std::endl;

 if (it->pt() > minJetPt_ && fabs(it->eta()) < maxJetEta_ && correctedbtag > bdis_) {
          mypatbJets->push_back(*it);}
  }
  std::auto_ptr<std::vector<pat::Jet> > myJets(mypatbJets);
  iEvent.put(myJets);
}

double Fjetcollection::Bmodification(double bvalue, double id,  std::string sftype){
if (id==0) {return bvalue;}
//std::cout<<"0  results    " <<bvalue<<std::endl;}

else{
const int N=5;
double X[N]={0,0.244,0.679,0.898,1};
double newlight[N]={0,0.2305,0.6117,0.8777,1 };
double newlightup[N]={0,0.2185,0.5921,0.8724,1 };
double newlightdown[N]={0,0.2425,0.6641,0.8981,1 };

double newB[N]={0,0.2237,0.7756,0.9125,1};
double newBup[N]={0,0.1739,0.7342,0.9057,1};
double newBdown[N]={0,0.2856,0.8066,0.9193,1};

double newC[N]={0,0.2414,0.7026,0.9022,1};
double newCup[N]={0,0.2267,0.6751,0.8982,1};
double newCdown[N]={0,0.2521,0.7213,0.9062,1};


int ie=0;

for (int i = 0; i < N; i++) {
    if(bvalue < X[i])  {
    if(i == 0) ie = 0; // less than minimal - back extrapolation with the 1st interval
    else  ie = i-1;
    break;
    }
}
  double y1,y2;
  double x1;
  x1 = X[ie];
  double x2;
  x2= X[ie+1];
  double result;

if (abs(id)==5) {
if (bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="up"){
y1 =newBup[ie];
y2 = newBup[ie+1];}
else if(bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="down"){
y1 =newBdown[ie];
y2 = newBdown[ie+1];}
else{
y1 =newB[ie];
y2 = newB[ie+1];}
} 

else if (abs(id)==4) {
if (bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="up"){
y1 =newCup[ie];
y2 = newCup[ie+1];}
else if(bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="down"){
y1 =newCdown[ie];
y2 = newCdown[ie+1];}
else{
y1 =newC[ie];
y2 = newC[ie+1];}
}

else {
if (bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="up"){
y1 =newlightup[ie];
y2 = newlightup[ie+1];}
else if(bscalefactorType_.substr(bscalefactorType_.find(':')+1)=="down"){
y1 =newlightdown[ie];
y2 = newlightdown[ie+1];}
else{
y1 =newlight[ie];
y2 = newlight[ie+1];}
} 

  result= (y1*(x2-bvalue) + y2*(bvalue-x1))/(x2-x1);
return result;

//std::cout<<"corrected  results    " <<result<<std::endl;

}
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Fjetcollection);

