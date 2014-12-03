// system include files
#include <memory>
#include <map>
#include <string>
#include <fstream>
#include <math.h>
// user include files
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "PhysicsTools/PatAlgos/interface/PATPrimaryVertexSelector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "PhysicsTools/SelectorUtils/interface/SimpleCutBasedElectronIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TH1.h"
#include "TH1I.h"
#include "TH2.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <TStyle.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "myanalysis/Atq/interface/scalefactorcalc.h"

//
// class declaration
class Firstanalyse : public edm::EDAnalyzer {
   public:
 /// default constructor
      explicit Firstanalyse(const edm::ParameterSet&);
 /// default destructor
      ~Firstanalyse();

 //     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
 /// everything that needs to be done before the event loop
      virtual void beginJob() ;
 /// everything that needs to be done during the event loop
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
/// everything that needs to be done after the event loop
      virtual void endJob() ;

 /// check if histogram was booked
 bool booked(const std::string histName) const { return hists_.find(histName.c_str())!=hists_.end(); };
 /// fill histogram if it had been booked before
 void fill(const std::string histName, double value) const { if(booked(histName.c_str())) hists_.find(histName.c_str())->second->Fill(value); };
/*
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&); virtuao void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup & iSetup);
*/
      // ----------member data ---------------------------

 // simple map to contain all histograms; 
 // histograms are booked in the beginJob() 
 // method
std::map<std::string, TH1F*> hists_;
std::map<std::string, TH2F*> histstwo_;
edm::InputTag photonsrc_;
edm::InputTag muons_;
edm::InputTag jets_;
edm::InputTag met_;
edm::InputTag tops_;
edm::InputTag pileupw_;
bool ifrealdata_;
std::string mysample_;
std::string mystep_;

TTree* ntuple;
std::string mytreename_;
std::vector<double>  *weight;
std::vector<double>  *Nvertex;
std::vector<double>  *NPileup;
std::vector<double>  *triggerSF;
std::vector<double>  *photonSF;
std::vector<double>  *muonSF;
std::vector<double>  *pileupSF;
 std::vector<double> *phiphoton;

 std::vector<double>  *sigmaieta;
 std::vector<double>  *HOE;
std::vector<double> *ptphoton;
 std::vector<double> *etaphoton;
 std::vector<double>  *pfChargedPU;
 std::vector<double>  *pfNeutralPU;
 std::vector<double>  *pfGammaPU;

std::vector<double>  *fTRun;
 std::vector<double>  *fTEvent;
 std::vector<double>  *fTLumiSection;

     std::vector<double>  *fTPhoR9;
     std::vector<double>  *fTPhoCaloPositionX;
     std::vector<double>  *fTPhoCaloPositionY;
     std::vector<double>  *fTPhoCaloPositionZ;
     std::vector<double>  *fTPhoHoverE;
     std::vector<double>  *fTPhoH1overE;
     std::vector<double>  *fTPhoH2overE;
     std::vector<double>  *fTPhoSigmaIetaIeta;
     std::vector<double>  *fTPhoSigmaEtaEta;

     std::vector<double>  *fTPhoE1x5;
     std::vector<double>  *fTPhoE2x5;
     std::vector<double>  *fTPhoE3x3;
     std::vector<double>  *fTPhoE5x5;
     std::vector<double>  *fTPhomaxEnergyXtal;
     std::vector<double>  *fTPhoIso03HcalDepth1;
     std::vector<double>  *fTPhoIso03HcalDepth2;
     std::vector<double>  *fTPhoIso04HcalDepth1;
     std::vector<double>  *fTPhoIso04HcalDepth2;
     std::vector<double>  *fTPhoIso03nTrksSolid;
     std::vector<double>  *fTPhoIso03nTrksHollow;
     std::vector<double>  *fTPhoIso04nTrksSolid;
     std::vector<double>  *fTPhoIso04nTrksHollow;
	std::vector<double> *Jetmulti;
	std::vector<double> *JetID;
	std::vector<double> *GenJetID;


};

//
// constants, enums and typedefs
//



//vector Int_t LepInd;
//
// static data member definitions
//

//
// constructors and destructor
//


//now do what ever initialization is needed

Firstanalyse::Firstanalyse(const edm::ParameterSet& iConfig):
  hists_(),
  photonsrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonsrc")),
  muons_(iConfig.getUntrackedParameter<edm::InputTag>("muons")),
  jets_ (iConfig.getUntrackedParameter<edm::InputTag>("jets" )),
  met_  (iConfig.getUntrackedParameter<edm::InputTag>("met"  )),
  tops_ (iConfig.getUntrackedParameter<edm::InputTag>("tops"  )),
//  pileupw_(iConfig.getUntrackedParameter<edm::InputTag>("Pileupweight")),
  ifrealdata_(iConfig.getParameter<bool>("realdata")),
  mysample_(iConfig.getParameter<std::string>("sample")),
  mystep_(iConfig.getParameter<std::string>("step")),
  mytreename_(iConfig.getParameter<std::string>("treename"))

{
pileupw_ =iConfig.getParameter<edm::InputTag>("pileupw");

}

Firstanalyse::~Firstanalyse()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
void
Firstanalyse::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
using namespace edm;
using namespace std;
using namespace reco;



double PUW;
edm::Handle<double> WG;
  iEvent.getByLabel(pileupw_, WG);

if (ifrealdata_) {PUW=1;}
else {PUW=*WG;}
double SFtrigger=0; 
double SFmu=0;  
double SFph=0;
double finalSF=0;
double wphjetkf=0;
double zphjetkf=0;
ScaleFactor mySF;

weight= new std::vector<double>;
Nvertex= new std::vector<double>;
triggerSF= new std::vector<double>;
photonSF= new std::vector<double>;
muonSF= new std::vector<double>;
pileupSF= new std::vector<double>;
NPileup= new std::vector<double>;

ptphoton= new std::vector<double>;
etaphoton= new std::vector<double>;
phiphoton= new std::vector<double>;

sigmaieta= new std::vector<double>;
pfChargedPU= new std::vector<double>;
pfNeutralPU= new std::vector<double>;
pfGammaPU= new std::vector<double>;
HOE= new std::vector<double>;

fTRun = new std::vector<double>;
fTEvent = new std::vector<double>;
fTLumiSection = new std::vector<double>;

fTPhoR9 = new std::vector<double>;
fTPhoCaloPositionX = new std::vector<double>;
fTPhoCaloPositionY = new std::vector<double>;
fTPhoCaloPositionZ = new std::vector<double>;
fTPhoHoverE = new std::vector<double>;
fTPhoH1overE = new std::vector<double>;
fTPhoH2overE = new std::vector<double>;
fTPhoSigmaIetaIeta = new std::vector<double>;
fTPhoSigmaEtaEta = new std::vector<double>;

fTPhoE1x5 = new std::vector<double>;
fTPhoE2x5 = new std::vector<double>;
fTPhoE3x3 = new std::vector<double>;
fTPhoE5x5 = new std::vector<double>;
fTPhomaxEnergyXtal = new std::vector<double>;
fTPhoIso03HcalDepth1 = new std::vector<double>;
fTPhoIso03HcalDepth2 = new std::vector<double>;
fTPhoIso04HcalDepth1 = new std::vector<double>;
fTPhoIso04HcalDepth2 = new std::vector<double>;
fTPhoIso03nTrksSolid = new std::vector<double>;
fTPhoIso03nTrksHollow = new std::vector<double>;
fTPhoIso04nTrksSolid = new std::vector<double>;
fTPhoIso04nTrksHollow = new std::vector<double>;
Jetmulti = new std::vector<double>;
JetID = new std::vector<double>;
GenJetID = new std::vector<double>;



fTRun->push_back(iEvent.id().run());
fTEvent->push_back(iEvent.id().event());
fTLumiSection->push_back(iEvent.luminosityBlock());


edm::Handle<edm::View<pat::Photon> > Ph    ;
  iEvent.getByLabel(photonsrc_, Ph  );
edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(muons_, muonHandle);
Handle<std::vector< reco::NamedCompositeCandidate >> mytop;
   iEvent.getByLabel(tops_ , mytop);
for(edm::View<pat::Photon>::const_iterator photon=Ph->begin(); photon!=Ph->end(); ++photon )
        {
        wphjetkf=mySF.wphjetkfactor(photon->pt ());
        zphjetkf=mySF.zphjetkfactor(photon->pt ());
	SFph=mySF.photonScaleFactor(photon->eta ());
        ptphoton->push_back(photon->pt ());
        etaphoton->push_back(photon->eta ());
        phiphoton->push_back(photon->phi ());
        sigmaieta->push_back(photon->sigmaIetaIeta());
        HOE->push_back(photon->hadTowOverEm());
        pfChargedPU->push_back(photon->userFloat("pfChargedPU"));
        pfNeutralPU->push_back(photon->userFloat("pfNeutralPU"));
        pfGammaPU->push_back(photon->userFloat("pfGammaPU"));


    fTPhoR9              ->push_back(photon->r9());
    fTPhoCaloPositionX   ->push_back(photon->caloPosition().X());
    fTPhoCaloPositionY   ->push_back(photon->caloPosition().Y());
    fTPhoCaloPositionZ   ->push_back(photon->caloPosition().Z());
    fTPhoHoverE          ->push_back(photon->hadronicOverEm());
    fTPhoH1overE         ->push_back(photon->hadronicDepth1OverEm());
    fTPhoH2overE         ->push_back(photon->hadronicDepth2OverEm());
    fTPhoSigmaIetaIeta   ->push_back(photon->sigmaIetaIeta());
    fTPhoSigmaEtaEta     ->push_back(photon->sigmaEtaEta());

    fTPhoE1x5 ->push_back(photon->e1x5());
    fTPhoE2x5 ->push_back(photon->e2x5());
    fTPhoE3x3 ->push_back(photon->e3x3());
    fTPhoE5x5 ->push_back(photon->e5x5());
    fTPhomaxEnergyXtal   ->push_back(photon->maxEnergyXtal());
    fTPhoIso03HcalDepth1 ->push_back(photon->hcalDepth1TowerSumEtConeDR03());
    fTPhoIso03HcalDepth2 ->push_back(photon->hcalDepth2TowerSumEtConeDR03());
    fTPhoIso04HcalDepth1 ->push_back(photon->hcalDepth1TowerSumEtConeDR04());
    fTPhoIso04HcalDepth2 ->push_back(photon->hcalDepth2TowerSumEtConeDR04());
    fTPhoIso03nTrksSolid ->push_back(photon->nTrkSolidConeDR03());
    fTPhoIso03nTrksHollow->push_back(photon->nTrkHollowConeDR03());
    fTPhoIso04nTrksSolid ->push_back(photon->nTrkSolidConeDR04());
    fTPhoIso04nTrksHollow->push_back(photon->nTrkHollowConeDR04());


	}

        for(std::vector<pat::Muon>::const_iterator mu1=muonHandle->begin(); mu1!=muonHandle->end(); ++mu1)
        {
	SFmu=mySF.muonScaleFactor(mu1->eta());
	SFtrigger=mySF.triggerScaleFactor(mu1->eta());
	}
finalSF=PUW*SFmu*SFtrigger*SFph;
 photonSF->push_back(SFph);
 muonSF->push_back(SFmu);
 pileupSF->push_back(PUW);
 triggerSF->push_back(SFtrigger);

if (mysample_=="wphjet") {finalSF=finalSF*wphjetkf;}
if (mysample_=="zphjet") {finalSF=finalSF*zphjetkf;}
if (mystep_=="afterb"){
for(std::vector< reco::NamedCompositeCandidate >::const_iterator top1=mytop->begin(); top1!=mytop->end(); ++top1)
        {finalSF=finalSF*1;}}
if (ifrealdata_) {finalSF=1;}
      hists_["recophotonPt" ]->Fill( 5,finalSF );
weight->push_back(finalSF);


edm::Handle<std::vector<pat::Jet> >  myjets;
iEvent.getByLabel(jets_, myjets);
        double btaginfo=0;
        double ctaginfo=0;
Jetmulti->push_back(myjets->size ()) ;

        for (std::vector<pat::Jet>::const_iterator jet1= myjets->begin(); jet1!=myjets->end(); ++jet1)
        {
 if (jet1->genParticleRefs().size()==1){JetID->push_back(jet1->genParticle(0)->pdgId());}
        btaginfo=0;
        ctaginfo=0;
//      for(uint i = 0 ; i < jet1->genParticleRefs().size() ; ++i ){
        if (jet1->genParticleRefs().size()==1){
        if( jet1->genParticle(0)->pdgId()!=5 && jet1->genParticle(0)->pdgId()!=4) btaginfo=btaginfo+1;
        if( jet1->genParticle(0)->pdgId()!=5 && btaginfo==0) ctaginfo=ctaginfo+1;
//cout<<jet1->genParticle(0)->pdgId()<<endl;
// }

        for (int j=0 ; j<15 ;++j){
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.244-0.1+j*0.02) & (btaginfo!=0)) hists_["lighttageffL" ]->Fill(2*j+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.244-0.1+j*0.02) & (btaginfo!=0)) hists_["lighttageffL" ]->Fill(2*j+1.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.244-0.1+j*0.02) && (btaginfo==0) && (ctaginfo!=0)) hists_["ctageffL" ]->Fill(2*j+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.244-0.1+j*0.02) & (btaginfo==0)  && (ctaginfo!=0)) hists_["ctageffL" ]->Fill(2*j+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.244-0.1+j*0.02) & (btaginfo==0) && (ctaginfo==0)) hists_["btageffL" ]->Fill(2*j+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.244-0.1+j*0.02) & (btaginfo==0) && (ctaginfo==0)) hists_["btageffL" ]->Fill(2*j+1.5);

        }
        for (int l=0 ; l<13 ;++l){
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.679-0.1+l*0.02) & (btaginfo!=0)) hists_["lighttageffM" ]->Fill(2*l+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.679-0.1+l*0.02) & (btaginfo!=0)) hists_["lighttageffM" ]->Fill(2*l+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.679-0.1+l*0.02) & (btaginfo==0)&& (ctaginfo!=0)) hists_["ctageffM" ]->Fill(2*l+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.679-0.1+l*0.02) & (btaginfo==0)&& (ctaginfo!=0)) hists_["ctageffM" ]->Fill(2*l+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.679-0.1+l*0.02) & (btaginfo==0)&& (ctaginfo==0)) hists_["btageffM" ]->Fill(2*l+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.679-0.1+l*0.02) & (btaginfo==0)&& (ctaginfo==0)) hists_["btageffM" ]->Fill(2*l+1.5);

       }

        for (int v=0 ; v<15 ;++v){
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo!=0)) hists_["lighttageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo!=0)) hists_["lighttageffT" ]->Fill(2*v+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo==0) && (ctaginfo!=0)) hists_["ctageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo==0)&& (ctaginfo!=0)) hists_["ctageffT" ]->Fill(2*v+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo==0) && (ctaginfo==0)) hists_["btageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo==0)&& (ctaginfo==0)) hists_["btageffT" ]->Fill(2*v+1.5);
}
}
}

/*
        for(edm::View<pat::Photon>::const_iterator photon=Ph->begin(); photon!=Ph->end(); ++photon )
        {
        hists_["recophotonPt" ]->Fill( photon->pt (),finalSF );
//        hists_["recophotoneta" ]->Fill( photon->eta () , PUW);
//        hists_["chargeiso"]->Fill( photon->userFloat("pfChaIsoAlt"));
//        hists_["neutraliso"]->Fill( photon->userFloat("pfNeuIsoAlt"));
//        hists_["gammaiso"]->Fill( photon->userFloat("pfGamIsoAlt"));
//        hists_["rhochargeiso"]->Fill( photon->userFloat("pfChargedPU"));
//        hists_["rhoneutraliso"]->Fill( photon->userFloat("pfNeutralPU"));
//        hists_["rhogammaiso"]->Fill( photon->userFloat("pfGammaPU"));
-------------------------------------------------------
cout<<photon->userFloat("phorho")<<"efferea"<<photon->userFloat("pfChargedEANew")<<endl;
cout<<"hadTowOverEmTightCut"<<photon->userFloat("hadTowOverEmTightCut")<<endl;
cout<<"showerShapeTightCut"<<photon->userFloat("showerShapeTightCut")<<endl;
cout<<"pfChargedPU"<<photon->userFloat("pfChargedPU")<<endl;
cout<<"pfChargedTightCut"<<photon->userFloat("pfChargedTightCut")<<endl;
cout<<"pfNeutralPU"<<photon->userFloat("pfNeutralPU")<<endl;
cout<<"pfNeutralTightCut"<<photon->userFloat("pfNeutralTightCut")<<endl;
cout<<"pfGammaPU"<<photon->userFloat("pfGammaPU")<<endl;
cout<<"pfGammaTightCut"<<photon->userFloat("pfGammaTightCut")<<endl;
cout<<endl;
cout<<endl;
-----------------------------------------------------------
// - photon->userFloat("pfNeutralEANew")*photon->userFloat("phorho")<<"JNJNNFVNNVFJNN"<<photon->userFloat("pfNeutralPU")<<"                    "<<photon->userFloat("reza")<<endl;        }
//cout<<PUW<<endl;
}

        for(std::vector<pat::Muon>::const_iterator mu1=muonHandle->begin(); mu1!=muonHandle->end(); ++mu1)
        {
        hists_["recomuonPt" ]->Fill( mu1->pt () ,PUW);
        hists_["recomuoneta"]->Fill( mu1->eta() ,PUW);
        }

edm::Handle<std::vector<pat::Jet> >  myjets;
   iEvent.getByLabel(jets_, myjets);
//        hists_["jetmulti" ]->Fill( myjets->size () );
	double n=0;
	double m=0;
  for (std::vector<pat::Jet>::const_iterator jet1= myjets->begin(); jet1!=myjets->end(); ++jet1)
        {
if (jet1->eta()<2.4){
         if(jet1->partonFlavour() == 5 || jet1->partonFlavour() == -5){histstwo_["btagtotal"]->Fill(jet1 ->pt (),jet1->eta());
         if (jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.623) histstwo_["btagRECO"]->Fill(jet1 ->pt (),jet1->eta()); 
       }
       else if(jet1->partonFlavour() == 4 || jet1->partonFlavour() == -4){
	histstwo_["ctagtotal"]->Fill(jet1 ->pt (),jet1->eta());
         if (jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.623) histstwo_["ctagRECO"]->Fill(jet1 ->pt (),jet1->eta()); 
      }
       else{
histstwo_["lighttotal"]->Fill(jet1 ->pt (),jet1->eta());
         if (jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.623) histstwo_["lightRECO"]->Fill(jet1 ->pt (),jet1->eta()); 
    }}
		if (jet1 ->pt ()>40 & jet1->eta()<2.5){
		m=m+1; 
                if (jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.244)
                {n=n+1;
		}
		}
		}
  hists_["bjetmulti" ]->Fill(n,PUW); 
 hists_["jetmulti" ]->Fill(m,PUW);
*/
edm::Handle<VertexCollection> vertex;
iEvent.getByLabel("offlinePrimaryVertices", vertex);
hists_["Nvertex" ]->Fill(vertex->size(),PUW); 
hists_["NvertexNOPU" ]->Fill(vertex->size());

 Nvertex->push_back(vertex->size());

edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
iEvent.getByLabel("addPileupInfo", puInfo);
      if (puInfo.isValid()) {
         std::vector<PileupSummaryInfo>::const_iterator puIt;
         int n = 0;
         for (puIt = puInfo->begin(); puIt != puInfo->end(); ++puIt, ++n) {
            if (puIt->getBunchCrossing() == 0) {
            NPileup->push_back(puIt->getTrueNumInteractions());
               break;
            }
         }}
else NPileup->push_back(-1);
GenJetID->push_back(-1);
//  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
//    iEvent.getByLabel("genParticles",mcHandle);
//  for (edm::View<reco::GenParticle>::const_iterator gen_iter = mcHandle->begin(); gen_iter!=mcHandle->end(); ++gen_iter) {
//                if (gen_iter->status()==3 && gen_iter->pdgId()<50) GenJetID->push_back(abs(gen_iter->pdgId()));}
  

/*
TLorentzVector genmuon;
TLorentzVector genneutrino;
TLorentzVector genw;
TLorentzVector genbjet;
TLorentzVector gentop; 
TLorentzVector genph;

  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
    iEvent.getByLabel("genParticles",mcHandle);
  for (edm::View<reco::GenParticle>::const_iterator gen_iter = mcHandle->begin(); gen_iter!=mcHandle->end(); ++gen_iter) {
                if (gen_iter->status()==3) 
                { 
//cout<<"pdjid"<< gen_iter->pdgId()<<endl;

                 switch (abs(gen_iter->pdgId())) 
                 {
                 case 11:
                 genmuon.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                 case 12:
                 genneutrino.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                case 13:
                 genmuon.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                 case 14:
                 genneutrino.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                case 15:
                 genmuon.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                 case 16:
                 genneutrino.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;

                 case 5:
                 genbjet.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
//                 cout<<gen_iter->px()<<endl; 
                case 22:
                 genph.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());

                 break;
                 default:
                 break;
                 }
//cout<<"pdjid"<< gen_iter->pdgId()<<endl;
                }
        }      
genw = genmuon + genneutrino;
gentop = genw + genbjet;

hists_["GENphotonPt"]->Fill(genph.Pt());
hists_["GENphotoneta"]->Fill(genph.Eta());
hists_["GENmuonPt"]->Fill(genmuon.Pt());
hists_["GENmuoneta"]->Fill(genmuon.Eta());
hists_["GENjetPt"]->Fill(genbjet.Pt());
hists_["GENjeteta"]->Fill(genbjet.Eta());
*/
ntuple->Fill();

delete weight;
delete pileupSF;
delete muonSF;
delete photonSF;
delete triggerSF;
delete Nvertex;
delete NPileup;
delete ptphoton;
delete etaphoton;
delete phiphoton;
delete sigmaieta;
delete HOE;
delete pfChargedPU;
delete pfNeutralPU;
delete pfGammaPU;

delete fTRun;
delete fTEvent;
delete fTLumiSection;

delete fTPhoR9;
delete fTPhoCaloPositionX;
delete fTPhoCaloPositionY;
delete fTPhoCaloPositionZ;
delete fTPhoHoverE;
delete fTPhoH1overE;
delete fTPhoH2overE;
delete fTPhoSigmaIetaIeta;
delete fTPhoSigmaEtaEta;

delete fTPhoE1x5;
delete fTPhoE2x5;
delete fTPhoE3x3;
delete fTPhoE5x5;
delete fTPhomaxEnergyXtal;
delete fTPhoIso03HcalDepth1;
delete fTPhoIso03HcalDepth2;
delete fTPhoIso04HcalDepth1;
delete fTPhoIso04HcalDepth2;
delete fTPhoIso03nTrksSolid;
delete fTPhoIso03nTrksHollow;
delete fTPhoIso04nTrksSolid;
delete fTPhoIso04nTrksHollow;
delete Jetmulti;
delete GenJetID;
delete JetID;
}
void
Firstanalyse::beginJob()
{
edm::Service<TFileService> fs;
ntuple=fs->make<TTree>(mytreename_.c_str() ,mytreename_.c_str());
ntuple->Branch("weight","std::vector<double>",&weight);
ntuple->Branch("Nvertex","std::vector<double>",&Nvertex);
ntuple->Branch("triggerSF","std::vector<double>",&triggerSF);
ntuple->Branch("photonSF","std::vector<double>",&photonSF);
ntuple->Branch("muonSF","std::vector<double>",&muonSF);
ntuple->Branch("pileupSF","std::vector<double>",&pileupSF);
ntuple->Branch("NPileup","std::vector<double>",&NPileup);
ntuple->Branch("ptphoton","std::vector<double>",&ptphoton);
ntuple->Branch("etaphoton","std::vector<double>",&etaphoton);
ntuple->Branch("phiphoton","std::vector<double>",&phiphoton);
ntuple->Branch("sigmaieta","std::vector<double>",&sigmaieta);
ntuple->Branch("HOE","std::vector<double>",&HOE);
ntuple->Branch("pfChargedPU","std::vector<double>",&pfChargedPU);
ntuple->Branch("pfNeutralPU","std::vector<double>",&pfNeutralPU);
ntuple->Branch("pfGammaPU","std::vector<double>",&pfGammaPU);
ntuple->Branch("fTRun","std::vector<double>",&fTRun);
ntuple->Branch("fTEvent","std::vector<double>",&fTEvent);
ntuple->Branch("fTLumiSection","std::vector<double>",&fTLumiSection);

ntuple->Branch("fTPhoR9","std::vector<double>",&fTPhoR9);
ntuple->Branch("fTPhoCaloPositionX","std::vector<double>",&fTPhoCaloPositionX);
ntuple->Branch("fTPhoCaloPositionY","std::vector<double>",&fTPhoCaloPositionY);
ntuple->Branch("fTPhoCaloPositionZ","std::vector<double>",&fTPhoCaloPositionZ);
ntuple->Branch("fTPhoHoverE","std::vector<double>",&fTPhoHoverE);
ntuple->Branch("fTPhoH1overE","std::vector<double>",&fTPhoH1overE);
ntuple->Branch("fTPhoH2overE","std::vector<double>",&fTPhoH2overE);
ntuple->Branch("fTPhoSigmaIetaIeta","std::vector<double>",&fTPhoSigmaIetaIeta);
ntuple->Branch("fTPhoSigmaEtaEta","std::vector<double>",&fTPhoSigmaEtaEta);
ntuple->Branch("fTPhoE1x5","std::vector<double>",&fTPhoE1x5);
ntuple->Branch("fTPhoE2x5","std::vector<double>",&fTPhoE2x5);
ntuple->Branch("fTPhoE3x3","std::vector<double>",&fTPhoE3x3);
ntuple->Branch("fTPhoE5x5","std::vector<double>",&fTPhoE5x5);
ntuple->Branch("fTPhomaxEnergyXtal","std::vector<double>",&fTPhomaxEnergyXtal);
ntuple->Branch("fTPhoIso03HcalDepth1","std::vector<double>",&fTPhoIso03HcalDepth1);
ntuple->Branch("fTPhoIso03HcalDepth2","std::vector<double>",&fTPhoIso03HcalDepth2);
ntuple->Branch("fTPhoIso04HcalDepth1","std::vector<double>",&fTPhoIso04HcalDepth1);
ntuple->Branch("fTPhoIso04HcalDepth2","std::vector<double>",&fTPhoIso04HcalDepth2);
ntuple->Branch("fTPhoIso03nTrksSolid","std::vector<double>",&fTPhoIso03nTrksSolid);
ntuple->Branch("fTPhoIso03nTrksHollow","std::vector<double>",&fTPhoIso03nTrksHollow);
ntuple->Branch("fTPhoIso04nTrksSolid","std::vector<double>",&fTPhoIso04nTrksSolid);
ntuple->Branch("fTPhoIso04nTrksHollow","std::vector<double>",&fTPhoIso04nTrksHollow);
ntuple->Branch("Jetmulti","std::vector<double>",&Jetmulti);
ntuple->Branch("GenJetID","std::vector<double>",&GenJetID);
ntuple->Branch("JetID","std::vector<double>",&JetID);


double dpt[]={30,50,70,100,160,2000};
double deta[]={0,0.8,1.6,2.4};

hists_["recophotonPt"] = fs->make<TH1F>("Photon_Pt", "Pt", 25, 0., 500.);
hists_["recophotoneta"] = fs->make<TH1F>("Photon Eta", "Eta", 30, -3., 3.);
hists_["recomuonPt"] = fs->make<TH1F>("Muon Pt", "Pt", 25, 0., 500.);
hists_["recomuoneta"] = fs->make<TH1F>("Muon Eta", "eta", 30, -3., 3.);
hists_["jetmulti"] = fs->make<TH1F>("Jet Multiplicity", "Jet Multiplicity", 11, -0.5, 10.5);
hists_["bjetmulti"] = fs->make<TH1F>("b Jet Multiplicity", "b Jet Multiplicity", 6, -0.5, 5.5);
hists_["chargeiso"] = fs->make<TH1F>("chargeiso", "chargeiso", 15, 0, 5);
hists_["neutraliso"] = fs->make<TH1F>("neutraliso", "neutraliso", 50, 0, 50);
hists_["gammaiso"] = fs->make<TH1F>("gammaiso", "gammaiso", 10, 0, 10);
hists_["rhochargeiso"] = fs->make<TH1F>("rhochargeiso", "rhochargeiso", 15, 0, 5);
hists_["rhoneutraliso"] = fs->make<TH1F>("rhoneutraliso", "rhoneutraliso", 50, 0, 50);
hists_["rhogammaiso"] = fs->make<TH1F>("rhogammaiso", "rhogammaiso", 10, 0, 10);
hists_["Nvertex"] = fs->make<TH1F>("Nvertex", "Nvertex",60,-0.5, 59.5);
hists_["NvertexNOPU"] = fs->make<TH1F>("NvertexNOPU", "NvertexNOPU",60,-0.5, 59.5);

hists_["GENphotonPt"] = fs->make<TH1F>("GENphoton_Pt", "Pt", 25,50., 500.);
hists_["GENphotoneta"] = fs->make<TH1F>("GENphoton_Eta", "eta", 30, -3., 3.);
hists_["GENmuonPt"] = fs->make<TH1F>("GENmuon_Pt", "Pt", 25, 26., 400.);
hists_["GENmuoneta"] = fs->make<TH1F>("GENmuon_Eta", "eta", 30, -3., 3.);
hists_["GENjetPt"] = fs->make<TH1F>("GENJet_Pt", "Pt", 25, 30., 500.);
hists_["GENjeteta"] = fs->make<TH1F>("GENJet_Eta", "Eta", 30, -3., 3.);
hists_["GENtopmass"] = fs->make<TH1F>("GENtopmass", "topmass", 50, 0., 400.);
hists_["GENWTmass"] = fs->make<TH1F>("GENWTmass", "WTmass", 25, 25., 150.);

histstwo_["btagtotal"] = fs->make<TH2F>("btagtotal", "btagtotal", 5, dpt,1,0,2.4);
histstwo_["ctagtotal"] = fs->make<TH2F>("ctagtotal", "ctagtotal", 5, dpt,1,0,2.4);
histstwo_["lighttotal"] = fs->make<TH2F>("lighttotal", "lighttotal", 5, dpt, 3, deta );

histstwo_["btagRECO"] = fs->make<TH2F>("btagRECO", "btagRECO", 5, dpt,1,0,2.4);
histstwo_["ctagRECO"] = fs->make<TH2F>("ctagRECO", "ctagRECO", 5, dpt,1,0,2.4);
histstwo_["lightRECO"] = fs->make<TH2F>("lightRECO", "lightRECO", 5, dpt, 3, deta );

hists_["btageffL"] = fs->make<TH1F>("btageffL", "btageffL", 30,0., 30.);
hists_["btageffM"] = fs->make<TH1F>("btageffM", "btageffM", 30,0., 30.);
hists_["btageffT"] = fs->make<TH1F>("btageffT", "btageffT", 30,0., 30.);
hists_["ctageffL"] = fs->make<TH1F>("ctageffL", "ctageffL", 30,0., 30.);
hists_["ctageffM"] = fs->make<TH1F>("ctageffM", "ctageffM", 30,0., 30.);
hists_["ctageffT"] = fs->make<TH1F>("ctageffT", "ctageffT", 30,0., 30.);
hists_["lighttageffL"] = fs->make<TH1F>("lighttageffL", "lighttageffL", 30,0., 30.);
hists_["lighttageffM"] = fs->make<TH1F>("lighttageffM", "lighttageffM", 30,0., 30.);
hists_["lighttageffT"] = fs->make<TH1F>("lighttageffT", "lighttageffT", 30,0., 30.);

}
void
Firstanalyse::endJob()
{
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(Firstanalyse);


