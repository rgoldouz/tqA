#ifndef Fjetcollection_h
#define Fjetcollection_h
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

class Fjetcollection : public edm::EDProducer {
  public:

    explicit Fjetcollection (const edm::ParameterSet & iConfig);
    ~Fjetcollection();

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:

    edm::InputTag theJetLabel_;
    double minJetPt_;
    double maxJetEta_;
    std::string bscalefactorType_;
    double bdis_;
    bool ifrealdata_;
    double Bmodification(double , double , std::string );

};



#endif
       

