// Author:  Reza Goldouzian,

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
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
 #include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
 #include "DataFormats/DetId/interface/DetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

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
////
// class declaration
//

class PatTopSelectionAnalyzer : public edm::EDAnalyzer {
   public:
 /// default constructor
      explicit PatTopSelectionAnalyzer(const edm::ParameterSet&);
 /// default destructor
      ~PatTopSelectionAnalyzer();

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
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup & iSetup);
*/
      // ----------member data ---------------------------

 // simple map to contain all histograms; 
 // histograms are booked in the beginJob() 
 // method
std::map<std::string, TH1F*> hists_; 
std::map<std::string, TH2F*> histstwo_;
std::string mytreename_;
std::string mysample_;
edm::InputTag photonsrc_;
edm::InputTag muons_;
edm::InputTag jets_;
edm::InputTag bjets_;
edm::InputTag met_;
edm::InputTag tops_;
TTree* ntuple;
 edm::InputTag rHInputProducerB_;
 edm::InputTag rHInputProducerE_;
edm::InputTag pileupw_;
bool ifrealdata_;

std::vector<double> *ptphoton;
 std::vector<double> *ptmuon;
 std::vector<double> *etaphoton;
 std::vector<double> *etamuon;
 std::vector<double> *ptjet;
 std::vector<double> *etajet;
 std::vector<double> *phimuon;
 std::vector<double> *phijet;
 std::vector<double> *phiphoton;
 std::vector<double> *masstop;
 std::vector<double> *mtw;
 std::vector<double> *cosmuonjet;
 std::vector<double> *cosmuonphoton;
 std::vector<double> *cosphotonjet;
 std::vector<double> *deltaphiphotonjet;
 std::vector<double> *deltaphiphotonmuon;
 std::vector<double> *deltaphimuonjet;
 std::vector<double> *deltaRphotonjet;
 std::vector<double> *deltaRphotonmuon;
 std::vector<double> *deltaRmuonjet;
 std::vector<double> *ht;
 std::vector<double> *photonmuonmass;
 std::vector<double> *costopphoton;
 std::vector<double> *coswphoton;

 std::vector<double> *topphotonmass; 
 std::vector<double> *lbphmass;
 std::vector<double> *pttop;
 std::vector<double> *etatop; 
 std::vector<double> *jetmultiplicity; 
 std::vector<double> *Genjetmultiplicity;
 std::vector<double> *bjetmultiplicity; 
 std::vector<double> *deltaphiphotonmet; 
 std::vector<double>  *cvsdiscriminant;
 std::vector<double>  *leptoncharge;
 std::vector<double>  *sigmaieta;
 std::vector<double>  *HOE;
 std::vector<double>  *Nvertex;
 std::vector<double>  *jetmatchinginfo;
//no matching =0 , light jet matching =1 , c jet matching =2 , b jet matching=3
 std::vector<double>  *MET;
 std::vector<double>  *METphi;
 std::vector<double>  *pfChargedPU;
 std::vector<double>  *pfNeutralPU;
 std::vector<double>  *pfGammaPU;
std::vector<double>  *NPileup;


 std::vector<double>  *btagSF;
 std::vector<double>  *btagSFup;
 std::vector<double>  *btagSFdown;

 std::vector<double>  *triggerSF;
 std::vector<double>  *triggerSFup;
 std::vector<double>  *triggerSFdown;

 std::vector<double>  *photonSF;
 std::vector<double>  *photonSFup;
 std::vector<double>  *photonSFdown;

 std::vector<double>  *muonSF;
 std::vector<double>  *muonSFup;
 std::vector<double>  *muonSFdown;

 std::vector<double>  *pileupSF;
 std::vector<double>  *pileupSFup;
 std::vector<double>  *pileupSFdown;

 std::vector<double>  *mistagSFup;
 std::vector<double>  *mistagSFdown;

 std::vector<double>  *weight;
 std::vector<double>  *soluinfo;

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

     std::vector<double>  *swissCross;

     std::vector<double>  *DeltaPhiPhotonTop;
     std::vector<double>  *CosWBjet;
     std::vector<double>  *SumJetPt;
     std::vector<double>  *DeltaRMuonNeutrino;
     std::vector<double>  *DeltaRMuonJet;
     std::vector<double>  *PhotonMuonPtRatio;
     std::vector<double>  *PtTopPhoton;
     std::vector<double>  *PzTopPhoton;
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
 
PatTopSelectionAnalyzer::PatTopSelectionAnalyzer(const edm::ParameterSet& iConfig):
  hists_(),
  photonsrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonsrc")),
  muons_(iConfig.getUntrackedParameter<edm::InputTag>("muons")),
  jets_ (iConfig.getUntrackedParameter<edm::InputTag>("jets" )),
  bjets_ (iConfig.getUntrackedParameter<edm::InputTag>("mybjets" )),
  met_  (iConfig.getUntrackedParameter<edm::InputTag>("met"  )),
  tops_ (iConfig.getUntrackedParameter<edm::InputTag>("tops"  )),
  ifrealdata_(iConfig.getParameter<bool>("realdata")),
  mytreename_(iConfig.getParameter<std::string>("treename")),
  mysample_(iConfig.getParameter<std::string>("sample"))

{
pileupw_ =iConfig.getParameter<edm::InputTag>("pileupw");
}



PatTopSelectionAnalyzer::~PatTopSelectionAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PatTopSelectionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
using namespace edm;
using namespace std;
using namespace reco;


//cout<<iEvent.id().event()<<endl;
double PUW;
edm::Handle<double> WG;
  iEvent.getByLabel(pileupw_, WG);
edm::Handle<double> PUWup;
  iEvent.getByLabel("pileupweightup", PUWup);
edm::Handle<double> PUWdown;
  iEvent.getByLabel("pileupweightdown", PUWdown);

edm::Handle<double> bsf;
  iEvent.getByLabel("BSF", bsf);
edm::Handle<double> bsfup;
  iEvent.getByLabel("BSFup", bsfup);
edm::Handle<double> bsfdown;
  iEvent.getByLabel("BSFdown", bsfdown);
edm::Handle<double> mistagup;
  iEvent.getByLabel("mistagup", mistagup);
edm::Handle<double> mistagdown;
  iEvent.getByLabel("mistagdown", mistagdown);


if (ifrealdata_) {PUW=1;}
else {PUW=*WG;}


ptphoton= new std::vector<double>; 
etaphoton= new std::vector<double>;
phiphoton= new std::vector<double>; 
ptmuon= new std::vector<double>; 
etamuon= new std::vector<double>; 
phimuon= new std::vector<double>; 
ptjet= new std::vector<double>; 
etajet= new std::vector<double>;
phijet= new std::vector<double>;
masstop=new std::vector<double>;
mtw=new std::vector<double>;
cosmuonjet=new std::vector<double>;
cosmuonphoton=new std::vector<double>;
cosphotonjet=new std::vector<double>;
deltaphiphotonjet= new std::vector<double>;
deltaphiphotonmuon= new std::vector<double>;
deltaphimuonjet= new std::vector<double>;
deltaRphotonjet= new std::vector<double>;
deltaRphotonmuon= new std::vector<double>;
deltaRmuonjet= new std::vector<double>;
ht= new std::vector<double>;
photonmuonmass= new std::vector<double>;
lbphmass= new std::vector<double>;
costopphoton= new std::vector<double>;
coswphoton= new std::vector<double>;

MET = new std::vector<double>;
METphi=  new std::vector<double>;


topphotonmass= new std::vector<double>;
lbphmass= new std::vector<double>;
pttop= new std::vector<double>;
etatop= new std::vector<double>;
jetmultiplicity= new std::vector<double>;
Genjetmultiplicity= new std::vector<double>;
bjetmultiplicity= new std::vector<double>;
deltaphiphotonmet= new std::vector<double>;
cvsdiscriminant= new std::vector<double>;
leptoncharge= new std::vector<double>;
sigmaieta= new std::vector<double>;
HOE= new std::vector<double>;
Nvertex= new std::vector<double>;
jetmatchinginfo= new std::vector<double>;
pfChargedPU= new std::vector<double>;
pfNeutralPU= new std::vector<double>;
pfGammaPU= new std::vector<double>;
NPileup= new std::vector<double>;

btagSF= new std::vector<double>;
btagSFup =new std::vector<double>;
btagSFdown =new std::vector<double>;
mistagSFup=new std::vector<double>;
mistagSFdown =new std::vector<double>;

triggerSF= new std::vector<double>;
triggerSFup= new std::vector<double>;
triggerSFdown= new std::vector<double>;

photonSF= new std::vector<double>;
photonSFup= new std::vector<double>;
photonSFdown= new std::vector<double>;

muonSF= new std::vector<double>;
muonSFup= new std::vector<double>;
muonSFdown= new std::vector<double>;

pileupSF= new std::vector<double>;
pileupSFup= new std::vector<double>;
pileupSFdown= new std::vector<double>;

weight= new std::vector<double>;
soluinfo= new std::vector<double>;

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
swissCross = new std::vector<double>;

DeltaPhiPhotonTop = new std::vector<double>;
CosWBjet = new std::vector<double>;
SumJetPt = new std::vector<double>;
DeltaRMuonNeutrino = new std::vector<double>;
DeltaRMuonJet = new std::vector<double>;
PhotonMuonPtRatio = new std::vector<double>;
PtTopPhoton = new std::vector<double>;
PzTopPhoton = new std::vector<double>;



TLorentzVector top4v;
TLorentzVector w4v;
TLorentzVector jet4v;
TLorentzVector bjet4v;

TLorentzVector muon4v;
TLorentzVector nu4v;
TLorentzVector photon4v;
TLorentzVector photon4vboost;
TLorentzVector jet4vboost;
TLorentzVector muon4vboost;
TLorentzVector w4vboost;

TVector3 topboost;
double wphjetkfac;
double zphjetkfac;


fTRun->push_back(iEvent.id().run());
fTEvent->push_back(iEvent.id().event());
fTLumiSection->push_back(iEvent.luminosityBlock());

/*
TLorentzVector genmuon;
TLorentzVector genneutrino;
TLorentzVector genw;
TLorentzVector genbjet;
TLorentzVector gentop; 
TVector3 topboost;

bool runonmc=false; 
math::XYZPoint vertexPosition;
vertexPosition.SetCoordinates(0,0,0);
edm::Handle<vector<reco::Vertex> > hpv;
iEvent.getByLabel( "offlinePrimaryVertices", hpv );
vector<reco::Vertex> goodVertices;
for (unsigned int i = 0; i < hpv->size(); i++) {
	af ( (*hpv)[i].ndof() > 4 && ( fabs((*hpv)[i].z()) <= 24. ) &&( fabs((*hpv)[i].position().rho()) <= 2.0 ) )
goodVertices("push_back((*hpv)[i]);
std::vector<double>  *btagSF;
 std::vector<double>  *triggerSF;
 std::vector<double>  *photonSF;
 std::vector<douweight;ble>  *muonSF;
 std::vector<double>  *pileupSF;
 std::vector<double>  *weight;
if (goodVertices.size()>0) {
vertexPosition = goodVertices[0].position();
}
}
*/
//sssssssssssssssssssssssss
//getting HLT menu
/*
 edm::Handle<edm::TriggerResults> hHLT;
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),hHLT);
  const edm::TriggerNames & triggerNames = iEvent.triggerNames(*hHLT);
  for (std::vector<std::string>::const_iterator it=triggerNames.triggerNames().begin(); it != triggerNames.triggerNames().end(); it++) {
    std::cout << *it << std::endl;
  }
*/
/*
  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
    iEvent.getByLabel("genParticles",mcHandle);
  for (edm::View<reco::GenParticle>::const_iterator gen_iter = mcHandle->begin(); gen_iter!=mcHandle->end(); ++gen_iter) {
		if (gen_iter->status()==3) 
		{ 
		 switch (abs(gen_iter->pdgId())) 
		 {
                 case 13:
		 genmuon.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
		 break;
                 case 14:
		 genneutrino.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
                 break;
                 case 5:
                 genbjet.SetPxPyPzE(gen_iter->px() ,gen_iter->py() ,gen_iter->pz () , gen_iter->energy());
//                 cout<<gen_iter->px()<<endl; 
                 break;
                 default:
                 break;
                 }
//cout<<"pdjid"<< gen_iter->pdgId()<<endl;
                }
	}    	


genw = genmuon + genneutrino;
gentop = genw + genbjet;
hists_["gentopmass1" ]->Fill(sqrt( TMath::Power(genw.Energy()+genbjet.Energy(),2)-(TMath::Power(genw.Px()+genbjet.Px(),2)+TMath::Power(genw.Py()+genbjet.Py(),2)+TMath::Power(genw.Pz()+genbjet.Pz(),2))));
hists_["gentopmass2" ]->Fill(sqrt( TMath::Power(gentop.Energy(),2)-(TMath::Power(gentop.Px(),2)+TMath::Power(gentop.Py(),2)+TMath::Power(gentop.Pz(),2))));
*/
edm::Handle<edm::View<pat::Photon> > Ph    ;
  iEvent.getByLabel(photonsrc_, Ph  );

edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(muons_, muonHandle);

Handle<std::vector< reco::NamedCompositeCandidate >> mytop;
   iEvent.getByLabel(tops_ , mytop);

edm::Handle<double> solutionINFO;
  iEvent.getByLabel(tops_, solutionINFO);
soluinfo->push_back(*solutionINFO);
ScaleFactor mysf;

////getting different scale factors
for(edm::View<pat::Photon>::const_iterator photon=Ph->begin(); photon!=Ph->end(); ++photon )
        {
	photonSF->push_back(mysf.photonScaleFactor(photon->eta ()));
        photonSFup->push_back(mysf.photonScaleFactorup(photon->eta ()));
        photonSFdown->push_back(mysf.photonScaleFactordown(photon->eta ()));


        wphjetkfac=mysf.wphjetkfactor(photon->pt ());
        zphjetkfac=mysf.zphjetkfactor(photon->pt ());
	}
for(std::vector<pat::Muon>::const_iterator mu1=muonHandle->begin(); mu1!=muonHandle->end(); ++mu1)
        {
        muonSF->push_back(mysf.muonScaleFactor(mu1->eta()));
        muonSFup->push_back(mysf.muonScaleFactorup(mu1->eta()));
        muonSFdown->push_back(mysf.muonScaleFactordown(mu1->eta()));

        triggerSF->push_back(mysf.triggerScaleFactor(mu1->eta()));
        triggerSFup->push_back(mysf.triggerScaleFactorup(mu1->eta()));
        triggerSFdown->push_back(mysf.triggerScaleFactordown(mu1->eta()));

	}
//for(std::vector< reco::NamedCompositeCandidate >::const_iterator top1=mytop->begin(); top1!=mytop->end(); ++top1)
//        {
//        btagSF->push_back(mysf.btagScaleFactor(top1->daughter("BJet")->pt()));
//	}
        btagSF->push_back(*bsf);
        btagSFup->push_back(*bsfup);
        btagSFdown->push_back(*bsfdown);
        mistagSFup ->push_back(*mistagup);
        mistagSFdown->push_back(*mistagdown);


        pileupSF->push_back(PUW);
        pileupSFup->push_back(*PUWup);
        pileupSFdown->push_back(*PUWdown);

        weight->push_back(PUW*photonSF->at(0)*muonSF->at(0)*triggerSF->at(0));

if (ifrealdata_) {weight->at(0)=1;}
if (mysample_=="wphjet") {weight->at(0)=weight->at(0)*wphjetkfac;}
if (mysample_=="zphjet") {weight->at(0)=weight->at(0)*zphjetkfac;}


	for(edm::View<pat::Photon>::const_iterator photon=Ph->begin(); photon!=Ph->end(); ++photon )
	{
        ptphoton->push_back(photon->pt ());
        etaphoton->push_back(photon->eta ());
        phiphoton->push_back(photon->phi ());

	hists_["recophotonPt" ]->Fill( photon->pt () ,weight->at(0) );
	hists_["recophotoneta" ]->Fill( photon->eta () ,weight->at(0) );
        hists_["recophotonphi" ]->Fill( photon->phi() ,weight->at(0) );

	photon4v.SetPxPyPzE(photon->px() ,photon->py() , photon->pz () , photon-> energy());
        if (abs(photon->eta()) < 1.4442){
        hists_["sigmaetaetabarrel" ]->Fill( photon->sigmaIetaIeta(),weight->at(0));}
        
        if (abs(photon->eta()) > 1.566 && abs(photon->eta())<2.5){
        hists_["sigmaetaetaendcap" ]->Fill( photon->sigmaIetaIeta(),weight->at(0));}
	
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
/*
  Handle<EcalRecHitCollection> pRecHitsB;
 iEvent.getByLabel( "reducedEcalRecHitsEB" , pRecHitsB);

 Handle<EcalRecHitCollection> pRecHitsE;
 iEvent.getByLabel("reducedEcalRecHitsEE", pRecHitsE);


	edm::InputTag i1("reducedEcalRecHitsEB");
	edm::InputTag i2("reducedEcalRecHitsEE");	
  EcalClusterLazyTools lazyTool(iEvent, iSetup, i1 , i2);
 */
//swissCross->push_back(EcalTools::swissCross(seedId, *myRecHits,0.));


	for(std::vector<pat::Muon>::const_iterator mu1=muonHandle->begin(); mu1!=muonHandle->end(); ++mu1)
	{
        ptmuon->push_back( mu1->pt () );
        etamuon->push_back(mu1->eta() );
        phimuon->push_back(mu1->phi() );
	muon4v.SetPxPyPzE(mu1->px() ,mu1->py() ,mu1->pz () , mu1-> energy());
	hists_["recomuonPt" ]->Fill( mu1->pt () ,weight->at(0) );
	hists_["recomuoneta"]->Fill( mu1->eta() ,weight->at(0) );
        hists_["recomuonphi"]->Fill( mu1->phi() ,weight->at(0) );

        hists_["muoncharge"]->Fill(mu1->pdgId()/abs(mu1->pdgId()) ,weight->at(0) );
//        leptoncharge->push_back(mu1->pdgId()/abs(mu1->pdgId())) ;
        leptoncharge->push_back(mu1->charge()) ;
        
	}

edm::Handle<std::vector<pat::Jet> >  myjets;
   iEvent.getByLabel(jets_, myjets);
//	hists_["jetmulti" ]->Fill( myjets->size () );
        double hj=0;
	double ii=0;
        double largestbdis=0;
        double btaginfo=0;
        double ctaginfo=0;
	double taginfo=0;
	for (std::vector<pat::Jet>::const_iterator jet1= myjets->begin(); jet1!=myjets->end(); ++jet1)
	{
//cout<<jet1->pt ()<<endl;

        btaginfo=0;
        ctaginfo=0;
//	for(uint i = 0 ; i < jet1->genParticleRefs().size() ; ++i ){
	if (jet1->genParticleRefs().size()==1){
	if( jet1->genParticle(0)->pdgId()!=5 && jet1->genParticle(0)->pdgId()!=4){
	btaginfo=btaginfo+1; 
	}
        if( jet1->genParticle(0)->pdgId()!=5 && btaginfo==0) {
	ctaginfo=ctaginfo+1;	
        }
        if( jet1->genParticle(0)->pdgId()==5 ) {
        }


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

        for (int v=0 ; v<10 ;++v){
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo!=0)) hists_["lighttageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo!=0)) hists_["lighttageffT" ]->Fill(2*v+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo==0) && (ctaginfo!=0)) hists_["ctageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo==0)&& (ctaginfo!=0)) hists_["ctageffT" ]->Fill(2*v+1.5);

        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")<0.898-0.1+v*0.02) & (btaginfo==0) && (ctaginfo==0)) hists_["btageffT" ]->Fill(2*v+0.5);
        if ((jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>0.898-0.1+v*0.02) & (btaginfo==0)&& (ctaginfo==0)) hists_["btageffT" ]->Fill(2*v+1.5);
}
}
	hj=hj+jet1 ->energy();
        ii=ii+1;
        if (jet1->bDiscriminator("combinedSecondaryVertexBJetTags")>largestbdis) {
	largestbdis=jet1->bDiscriminator("combinedSecondaryVertexBJetTags");
        if (jet1->genParticleRefs().size()==1) taginfo=jet1->genParticle(0)->pdgId();
	}
//		ptjet->push_back(jet1 ->pt () );
//		etajet->push_back(jet1->eta() );
//		hists_["recojetPt" ]->Fill( jet1->pt () );
//		hists_["recojeteta"]->Fill( jet1->eta() );
//		jet4v.SetPxPyPzE(jet1->px() ,jet1->py() ,jet1->pz () , jet1->energy());
	}
//cout<<" ----------------------------------------------------"<<endl;

//edm::Handle<reco::GenJetCollection> genJets;
//iEvent.getByLabel("ak5GenJets", genJets);
//double jetnum=0;
//for (reco::GenJetCollection::const_iterator genjet1= genJets->begin(); genjet1!=genJets->end(); ++genjet1)
//    {
//if (genjet1->pt()>30) {
//jetnum=jetnum+1;}}

Genjetmultiplicity->push_back(0);

hists_["jetmulti" ]->Fill(ii,weight->at(0) );
jetmultiplicity->push_back(ii);

cvsdiscriminant->push_back(largestbdis);
hists_["bdisc"]->Fill(largestbdis,weight->at(0) );

jetmatchinginfo->push_back(taginfo);

edm::Handle<std::vector<pat::Jet> >  btagjet;
iEvent.getByLabel(bjets_, btagjet);

hists_["bjetmulti" ]->Fill(btagjet->size () ,weight->at(0));
bjetmultiplicity->push_back(btagjet->size ()) ;


edm::Handle<std::vector<pat::MET> >  met;
   iEvent.getByLabel("patMETsPF",met);
	for (std::vector<pat::MET>::const_iterator mymet=met->begin(); mymet!=met->end(); ++mymet)
	{
 hists_["met" ]->Fill(mymet->energy());
// if(mymet->energy()<30)cout<<"ei vaie man"<<endl;
	nu4v.SetPxPyPzE(mymet->px() ,mymet->py() ,mymet->pz () ,mymet->energy());
	}

   double npu = 0;
edm::Handle<VertexCollection> vertex;
   iEvent.getByLabel("offlinePrimaryVertices", vertex);
npu=vertex->size();

Nvertex->push_back(npu);

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

//double mtop=1000;
//double mtopn;
	for(std::vector< reco::NamedCompositeCandidate >::const_iterator top1=mytop->begin(); top1!=mytop->end(); ++top1)
	{
//std::cout<<"top mass is: "<< top1->mass()<<"w mass is "<<top1->daughter("W")->mass()<< "nu pz is "<<top1->daughter("RecoNu")->pz()<<std::endl;


//	mtopn=sqrt( TMath::Power(top1->daughter("W")->energy()+top1->daughter("BJet")->energy(),2)-(TMath::Power(top1->daughter("W")->px()+top1->daughter("BJet")->px(),2)+TMath::Power(top1->daughter("W")->px()+top1->daughter("BJet")->px(),2)+TMath::Power(top1->daughter("W")->pz()+top1->daughter("BJet")->pz(),2)));
//		if (abs(mtop-172.5)>abs(mtopn-172.5)){
		//hists_["topmass" ]->Fill(sqrt( TMath::Power(top1->energy(),2)-(TMath::Power(top1->px(),2)+TMath::Power(top1->py(),2)+TMath::Power(top1->pz(),2) )));
		//masstop->push_back(sqrt( TMath::Power(top1->energy(),2)-(TMath::Power(top1->px(),2)+TMath::Power(top1->py(),2)+TMath::Power(top1->pz(),2) )) );
		hists_["topmass" ]->Fill(sqrt( TMath::Power(top1->daughter("W")->energy()+top1->daughter("BJet")->energy(),2)-(TMath::Power(top1->daughter("W")->px()+top1->daughter("BJet")->px(),2)+TMath::Power(top1->daughter("W")->py()+top1->daughter("BJet")->py(),2)+TMath::Power(top1->daughter("W")->pz()+top1->daughter("BJet")->pz(),2))),weight->at(0) );
masstop->push_back(sqrt( TMath::Power(top1->daughter("W")->energy()+top1->daughter("BJet")->energy(),2)-(TMath::Power(top1->daughter("W")->px()+top1->daughter("BJet")->px(),2)+TMath::Power(top1->daughter("W")->py()+top1->daughter("BJet")->py(),2)+TMath::Power(top1->daughter("W")->pz()+top1->daughter("BJet")->pz(),2))));
		top4v.SetPxPyPzE(top1->px(), top1->py() ,top1->pz () ,top1->energy());
		topboost=top4v.BoostVector();
		w4v.SetPxPyPzE(top1->daughter("W")->px(), top1->daughter("W")->py() ,top1->daughter("W")->pz () ,top1->daughter("W")->energy());
		bjet4v.SetPxPyPzE(top1->daughter("BJet")->px(),top1->daughter("BJet")->py(),top1->daughter("BJet")->pz(),top1->daughter("BJet")->energy());
//		mtop=mtopn;
//		}
pttop->push_back(top1->pt());
etatop->push_back(top1->eta());
hists_["topPt"]->Fill(top1->pt(),weight->at(0) );
hists_["topEta"]->Fill(top1->eta(),weight->at(0) );
	}

              ptjet->push_back(bjet4v.Pt());
              etajet->push_back(bjet4v.Eta()) ;
              phijet->push_back(bjet4v.Phi()) ;

              hists_["recojetPt" ]->Fill( bjet4v.Pt() ,weight->at(0) );
              hists_["recojeteta"]->Fill( bjet4v.Eta() ,weight->at(0) );
              hists_["recojetphi"]->Fill( bjet4v.Phi() ,weight->at(0) );		



hists_["WTmass" ]->Fill(sqrt(TMath::Power(muon4v.Et()+nu4v.Et(),2) - (TMath::Power(muon4v.Px()+nu4v.Px(),2)+TMath::Power(muon4v.Py()+nu4v.Py(),2))),weight->at(0) );
//cout<<muon4v.Et()<<"muon ET, Pt"<<muon4v.Pt()<<endl;
//cout<<nu4v.Et()<<"met ET, Pt"<<nu4v.Pt()<<endl;

MET->push_back(nu4v.Et());
METphi->push_back(nu4v.Phi());


mtw->push_back(sqrt(TMath::Power(muon4v.Et()+nu4v.Et(),2) - (TMath::Power(muon4v.Px()+nu4v.Px(),2)+TMath::Power(muon4v.Py()+nu4v.Py(),2))));

hists_["Cosmuonphoton"]->Fill(TMath::Cos((photon4v.Vect()).Angle(muon4v.Vect())),weight->at(0) );
cosmuonphoton->push_back(TMath::Cos((photon4v.Vect()).Angle(muon4v.Vect())));

hists_["Deltaphiphotonjet"]->Fill(abs(deltaPhi(photon4v.Phi(),bjet4v.Phi())),weight->at(0) );
deltaphiphotonjet->push_back(abs(deltaPhi(photon4v.Phi(),bjet4v.Phi())));

hists_["Deltaphiphotonmuon"]->Fill(abs(deltaPhi(photon4v.Phi(),muon4v.Phi())),weight->at(0) );
deltaphiphotonmuon->push_back(abs(deltaPhi(photon4v.Phi(),muon4v.Phi())));

hists_["Deltaphimuonjet"]->Fill(abs(deltaPhi(bjet4v.Phi(),muon4v.Phi())),weight->at(0) );
deltaphimuonjet->push_back(abs(deltaPhi(bjet4v.Phi(),muon4v.Phi())));

hists_["DeltaRphotonjet"]->Fill(abs(deltaR(photon4v.Eta(),photon4v.Phi(),bjet4v.Eta(),bjet4v.Phi())),weight->at(0) );
deltaRphotonjet->push_back(abs(deltaR(photon4v.Eta(),photon4v.Phi(),bjet4v.Eta(),bjet4v.Phi())));

hists_["DeltaRphotonmuon"]->Fill(abs(deltaR(photon4v.Eta(),photon4v.Phi(),muon4v.Eta(),muon4v.Phi())),weight->at(0) );
deltaRphotonmuon->push_back(abs(deltaR(photon4v.Eta(),photon4v.Phi(),muon4v.Eta(),muon4v.Phi())));

hists_["DeltaRmuonjet"]->Fill(abs(deltaR(bjet4v.Eta(),bjet4v.Phi(),muon4v.Eta(),muon4v.Phi())),weight->at(0));
deltaRmuonjet->push_back(abs(deltaR(bjet4v.Eta(),bjet4v.Phi(),muon4v.Eta(),muon4v.Phi())));

hists_["HT"]->Fill(hj,weight->at(0) );
ht->push_back(hj);

hists_["Photonmuonmass"]->Fill(sqrt(TMath::Power(muon4v.E()+photon4v.E(),2) - (TMath::Power(muon4v.Px()+photon4v.Px(),2)+TMath::Power(muon4v.Py()+photon4v.Py(),2)+TMath::Power(muon4v.Px()+photon4v.Pz(),2))),weight->at(0) );
photonmuonmass->push_back(sqrt(TMath::Power(muon4v.E()+photon4v.E(),2) - (TMath::Power(muon4v.Px()+photon4v.Px(),2)+TMath::Power(muon4v.Py()+photon4v.Py(),2)+TMath::Power(muon4v.Pz()+photon4v.Pz(),2))));

hists_["Topphotonmass"]->Fill(sqrt(TMath::Power(top4v.E()+photon4v.E(),2) - (TMath::Power(top4v.Px()+photon4v.Px(),2)+TMath::Power(top4v.Py()+photon4v.Py(),2)+TMath::Power(top4v.Pz()+photon4v.Pz(),2))),weight->at(0) );
topphotonmass->push_back(sqrt(TMath::Power(top4v.E()+photon4v.E(),2) - (TMath::Power(top4v.Px()+photon4v.Px(),2)+TMath::Power(top4v.Py()+photon4v.Py(),2)+TMath::Power(top4v.Pz()+photon4v.Pz(),2))));

lbphmass->push_back(sqrt(TMath::Power(bjet4v.E()+muon4v.E()+photon4v.E(),2) - (TMath::Power(bjet4v.Px()+muon4v.Px()+photon4v.Px(),2)+TMath::Power(bjet4v.Py()+muon4v.Py()+photon4v.Py(),2)+TMath::Power(bjet4v.Pz()+muon4v.Pz()+photon4v.Pz(),2))));

hists_["Costopphoton"]->Fill(TMath::Cos((photon4v.Vect()).Angle(top4v.Vect())),weight->at(0) );
costopphoton->push_back(TMath::Cos((photon4v.Vect()).Angle(top4v.Vect())));

hists_["Deltaphiphotonmet"]->Fill(abs(deltaPhi(photon4v.Phi(),nu4v.Phi())),weight->at(0) );
deltaphiphotonmet->push_back(abs(deltaPhi(photon4v.Phi(),nu4v.Phi())));

hists_["Coswphoton"]->Fill(TMath::Cos((photon4v.Vect()).Angle(w4v.Vect())),weight->at(0) );
coswphoton->push_back(TMath::Cos((photon4v.Vect()).Angle(w4v.Vect())));

DeltaPhiPhotonTop->push_back(abs(deltaPhi(photon4v.Phi(),top4v.Phi())));
SumJetPt->push_back(hj-sqrt(TMath::Power(bjet4v.Px(),2)+ TMath::Power(bjet4v.Py(),2)));
CosWBjet->push_back(TMath::Cos((bjet4v.Vect()).Angle(w4v.Vect())));
DeltaRMuonNeutrino->push_back(abs(deltaR(nu4v.Eta(),nu4v.Phi(),muon4v.Eta(),muon4v.Phi())));
DeltaRMuonJet->push_back(abs(deltaR(muon4v.Eta(),muon4v.Phi(),bjet4v.Eta(),bjet4v.Phi())));
PhotonMuonPtRatio->push_back(sqrt(TMath::Power(photon4v.Px(),2)+TMath::Power(photon4v.Py(),2))/sqrt(TMath::Power(muon4v.Px(),2)+TMath::Power(muon4v.Py(),2)));
PtTopPhoton->push_back((top4v+photon4v).Pt());
PzTopPhoton->push_back(abs((top4v+photon4v).Pz()));



muon4vboost=muon4v;
jet4vboost=bjet4v;
photon4vboost=photon4v;
w4vboost=w4v;
muon4vboost.Boost(-topboost);
jet4vboost.Boost(-topboost);
photon4vboost.Boost(-topboost);
w4vboost.Boost(-topboost);
//CosWBjet->push_back(9);


hists_["Cosmuonjet"]->Fill(TMath::Cos((muon4vboost.Vect()).Angle(jet4vboost.Vect())),weight->at(0) );
hists_["Cosphotonjet"] ->Fill(TMath::Cos((photon4vboost.Vect()).Angle(jet4vboost.Vect())),weight->at(0) );
cosmuonjet->push_back(TMath::Cos((muon4vboost.Vect()).Angle(jet4vboost.Vect())));
cosphotonjet->push_back(TMath::Cos((photon4vboost.Vect()).Angle(jet4vboost.Vect())));


ntuple->Fill();

delete ptphoton;
delete etaphoton;
delete phiphoton;
delete ptmuon;
delete etamuon;
delete phimuon;
delete ptjet;
delete etajet;
delete phijet;
delete masstop;
delete mtw;
delete cosmuonjet;
delete cosmuonphoton;
delete cosphotonjet;
delete deltaphiphotonjet;
delete deltaphiphotonmuon;
delete deltaphimuonjet;
delete deltaRphotonjet;
delete deltaRphotonmuon;
delete deltaRmuonjet;
delete ht;
delete photonmuonmass;
delete costopphoton;
delete coswphoton;
delete lbphmass;
delete topphotonmass;
delete pttop;
delete etatop;
delete jetmultiplicity;
delete Genjetmultiplicity;
delete bjetmultiplicity;
delete deltaphiphotonmet;
delete cvsdiscriminant;
delete leptoncharge;
delete sigmaieta;
delete HOE; 
delete Nvertex;
delete jetmatchinginfo;
delete MET;
delete METphi;
delete NPileup;

delete pfChargedPU;
delete pfNeutralPU;
delete pfGammaPU;

delete btagSF;
delete btagSFup;
delete btagSFdown;
delete mistagSFup;
delete mistagSFdown;

delete triggerSF;
delete triggerSFup;
delete triggerSFdown;

delete photonSF;
delete photonSFup;
delete photonSFdown;



delete muonSF;
delete muonSFup;
delete muonSFdown;

delete pileupSF;
delete pileupSFup;
delete pileupSFdown;

delete weight;
delete soluinfo;

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
delete swissCross;

delete DeltaPhiPhotonTop;
delete CosWBjet;
delete SumJetPt;
delete DeltaRMuonNeutrino;
delete DeltaRMuonJet;
delete PhotonMuonPtRatio;
delete PtTopPhoton;
delete PzTopPhoton;

}

// ------------ method called once each job just before starting event loop  ------------
void 
PatTopSelectionAnalyzer::beginJob()
{

  // register to the TFileService
edm::Service<TFileService> fs;
ntuple=fs->make<TTree>(mytreename_.c_str() ,mytreename_.c_str());
ntuple->Branch("ptphoton","std::vector<double>",&ptphoton);
ntuple->Branch("ptmuon","std::vector<double>",&ptmuon);
ntuple->Branch("etaphoton","std::vector<double>",&etaphoton);
ntuple->Branch("etamuon","std::vector<double>",&etamuon);
ntuple->Branch("ptjet","std::vector<double>",&ptjet);
ntuple->Branch("etajet","std::vector<double>",&etajet);
ntuple->Branch("phiphoton","std::vector<double>",&phiphoton);
ntuple->Branch("phijet","std::vector<double>",&phijet);
ntuple->Branch("phimuon","std::vector<double>",&phimuon);
ntuple->Branch("masstop","std::vector<double>",&masstop);
ntuple->Branch("mtw","std::vector<double>",&mtw);
ntuple->Branch("cosmuonjet","std::vector<double>",&cosmuonjet);
ntuple->Branch("cosmuonphoton","std::vector<double>",&cosmuonphoton);
ntuple->Branch("cosphotonjet","std::vector<double>",&cosphotonjet);
ntuple->Branch("deltaphiphotonjet","std::vector<double>",&deltaphiphotonjet);
ntuple->Branch("deltaphiphotonmuon","std::vector<double>",&deltaphiphotonmuon);
ntuple->Branch("deltaphimuonjet","std::vector<double>",&deltaphimuonjet);
ntuple->Branch("deltaRphotonjet","std::vector<double>",&deltaRphotonjet);
ntuple->Branch("deltaRphotonmuon","std::vector<double>",&deltaRphotonmuon);
ntuple->Branch("deltaRmuonjet","std::vector<double>",&deltaRmuonjet);
ntuple->Branch("ht","std::vector<double>",&ht);
ntuple->Branch("photonmuonmass","std::vector<double>",&photonmuonmass);
ntuple->Branch("costopphoton","std::vector<double>",&costopphoton);
ntuple->Branch("coswphoton","std::vector<double>",&coswphoton);

ntuple->Branch("topphotonmass","std::vector<double>",&topphotonmass);
ntuple->Branch("lbphmass","std::vector<double>",&lbphmass);
ntuple->Branch("pttop","std::vector<double>",&pttop);
ntuple->Branch("etatop","std::vector<double>",&etatop);
ntuple->Branch("jetmultiplicity","std::vector<double>",&jetmultiplicity);
ntuple->Branch("Genjetmultiplicity","std::vector<double>",&Genjetmultiplicity);
ntuple->Branch("bjetmultiplicity","std::vector<double>",&bjetmultiplicity);
ntuple->Branch("deltaphiphotonmet","std::vector<double>",&deltaphiphotonmet);
ntuple->Branch("cvsdiscriminant","std::vector<double>",&cvsdiscriminant);
ntuple->Branch("leptoncharge","std::vector<double>",&leptoncharge);
ntuple->Branch("sigmaieta","std::vector<double>",&sigmaieta);
ntuple->Branch("HOE","std::vector<double>",&HOE);
ntuple->Branch("Nvertex","std::vector<double>",&Nvertex);
ntuple->Branch("jetmatchinginfo","std::vector<double>",&jetmatchinginfo);
ntuple->Branch("MET","std::vector<double>",&MET);
ntuple->Branch("METphi","std::vector<double>",&METphi);
ntuple->Branch("pfChargedPU","std::vector<double>",&pfChargedPU);
ntuple->Branch("pfNeutralPU","std::vector<double>",&pfNeutralPU);
ntuple->Branch("pfGammaPU","std::vector<double>",&pfGammaPU);

ntuple->Branch("btagSF","std::vector<double>",&btagSF);
ntuple->Branch("btagSFup","std::vector<double>",&btagSFup);
ntuple->Branch("btagSFdown","std::vector<double>",&btagSFdown);
ntuple->Branch("mistagSFup","std::vector<double>",&mistagSFup);
ntuple->Branch("mistagSFdown","std::vector<double>",&mistagSFdown);

ntuple->Branch("triggerSF","std::vector<double>",&triggerSF);
ntuple->Branch("triggerSFup","std::vector<double>",&triggerSFup);
ntuple->Branch("triggerSFdown","std::vector<double>",&triggerSFdown);

ntuple->Branch("photonSF","std::vector<double>",&photonSF);
ntuple->Branch("photonSFup","std::vector<double>",&photonSFup);
ntuple->Branch("photonSFdown","std::vector<double>",&photonSFdown);


ntuple->Branch("muonSF","std::vector<double>",&muonSF);
ntuple->Branch("muonSFup","std::vector<double>",&muonSFup);
ntuple->Branch("muonSFdown","std::vector<double>",&muonSFdown);

ntuple->Branch("pileupSF","std::vector<double>",&pileupSF);
ntuple->Branch("pileupSFup","std::vector<double>",&pileupSFup);
ntuple->Branch("pileupSFdown","std::vector<double>",&pileupSFdown);

ntuple->Branch("weight","std::vector<double>",&weight);
ntuple->Branch("NPileup","std::vector<double>",&NPileup);
ntuple->Branch("soluinfo","std::vector<double>",&soluinfo);

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
ntuple->Branch("swissCross","std::vector<double>",&swissCross);
ntuple->Branch("DeltaPhiPhotonTop","std::vector<double>",&DeltaPhiPhotonTop);
ntuple->Branch("CosWBjet","std::vector<double>",&CosWBjet);
ntuple->Branch("SumJetPt","std::vector<double>",&SumJetPt);
ntuple->Branch("DeltaRMuonNeutrino","std::vector<double>",&DeltaRMuonNeutrino);
ntuple->Branch("DeltaRMuonJet","std::vector<double>",&DeltaRMuonJet);
ntuple->Branch("PhotonMuonPtRatio","std::vector<double>",&PhotonMuonPtRatio);
ntuple->Branch("PtTopPhoton","std::vector<double>",&PtTopPhoton);
ntuple->Branch("PzTopPhoton","std::vector<double>",&PzTopPhoton);

 // book histograms:
hists_["recophotonPt"] = fs->make<TH1F>("photon_Pt", "Pt", 25,50., 500.);
hists_["recophotoneta"] = fs->make<TH1F>("photon_Eta", "eta", 30, -3., 3.);
hists_["recophotonphi"] = fs->make<TH1F>("photon_phi", "phi", 30, -3.5, 3.5);
hists_["recomuonPt"] = fs->make<TH1F>("muon_Pt", "Pt", 25, 26., 400.);
hists_["recomuoneta"] = fs->make<TH1F>("muon_Eta", "eta", 30, -3., 3.);
hists_["recomuonphi"] = fs->make<TH1F>("muon_phi", "phi", 30, -3.5, 3.5);
hists_["recojetPt"] = fs->make<TH1F>("Jet_Pt", "Pt", 25, 30., 500.);
hists_["recojeteta"] = fs->make<TH1F>("Jet_Eta", "Eta", 30, -3., 3.);
hists_["recojetphi"] = fs->make<TH1F>("Jet_phi", "phi", 30, -3.5, 3.5);
hists_["jetmulti"] = fs->make<TH1F>("Jet_Multiplicity", "Jet_Multiplicity", 11, -0.5, 10.5);
hists_["topmass"] = fs->make<TH1F>("topmass", "topmass", 50, 0., 400.);
hists_["WTmass"] = fs->make<TH1F>("WTmass", "WTmass", 25, 25., 150.);
hists_["Cosmuonjet"] = fs->make<TH1F>("Cosmuonjet" , "Cosmuonjet",40 ,-1 ,1);
hists_["Cosmuonphoton"] = fs->make<TH1F>("Cosmuonphoton", "Cosmuonphoton",40 ,-1 ,1);
hists_["Cosphotonjet"] = fs->make<TH1F>("Cosphotonjet", "Cosphotonjet",40 ,-1 ,1);
hists_["Deltaphiphotonjet"] = fs->make<TH1F>("Deltaphiphotonjet", "Deltaphiphotonjet",20,0,3.5);
hists_["Deltaphiphotonmuon"] = fs->make<TH1F>("Deltaphiphotonmuon","Deltaphiphotonmuon",20,0,3.5);
hists_["Deltaphimuonjet"] = fs->make<TH1F>("Deltaphimuonjet","Deltaphimuonjet",20,0,3.5);
hists_["DeltaRphotonmuon"] = fs->make<TH1F>("DeltaRphotonmuon","DeltaRphotonmuon",40,0,7);
hists_["DeltaRphotonjet"] = fs->make<TH1F>("DeltaRphotonjet","DeltaRphotonjet",40,0,7);
hists_["DeltaRmuonjet"] = fs->make<TH1F>("DeltaRmuonjet","DeltaRmuonjet",40,0,7);
hists_["HT"] = fs->make<TH1F>("HT","HT",100,0,2000);
hists_["Photonmuonmass"] = fs->make<TH1F>("Photonmuonmass","Photonmuonmass", 50, 0., 400.);
hists_["Costopphoton"] = fs->make<TH1F>("Costopphoton","Costopphoton",40 ,-1 ,1);
hists_["Coswphoton"] = fs->make<TH1F>("Coswphoton","Coswphoton",40 ,-1 ,1);

hists_["sigmaetaetabarrel" ]= fs->make<TH1F>("sigmaetaetabarrel","sigmaetaetabarrel",15,0,0.03);
hists_["sigmaetaetaendcap" ]= fs->make<TH1F>("sigmaetaetaendcap","sigmaetaetaendcap",35,0,0.07);
hists_["topPt"] = fs->make<TH1F>("top_Pt", "Pt", 35,0., 750.);
hists_["topEta"]= fs->make<TH1F>("top_Eta", "eta", 30, -3., 3.);
hists_["bdisc"]= fs->make<TH1F>("b_tag_info","bdisc",20,0,1);
hists_["Topphotonmass"]= fs->make<TH1F>("top_photon_mass","top_photon_mass", 50, 0., 700.);
hists_["bjetmulti"] = fs->make<TH1F>("bJet_Multiplicity", "bJet_Multiplicity", 3, -0.5, 2.5);
hists_["muoncharge"] = fs->make<TH1F>("muoncharge", "muoncharge", 5, -2.5, 2.5);
hists_["Deltaphiphotonmet"] = fs->make<TH1F>("Deltaphiphotonmet","Deltaphiphotonmet",20,0,3.5);
hists_["met"] = fs->make<TH1F>("met", "met", 60,0., 750.);
hists_["btageffL"] = fs->make<TH1F>("btageffL", "btageffL", 30,0., 30.);
hists_["btageffM"] = fs->make<TH1F>("btageffM", "btageffM", 30,0., 30.);
hists_["btageffT"] = fs->make<TH1F>("btageffT", "btageffT", 30,0., 30.);
hists_["ctageffL"] = fs->make<TH1F>("ctageffL", "ctageffL", 30,0., 30.);
hists_["ctageffM"] = fs->make<TH1F>("ctageffM", "ctageffM", 30,0., 30.);
hists_["ctageffT"] = fs->make<TH1F>("ctageffT", "ctageffT", 30,0., 30.);
hists_["lighttageffL"] = fs->make<TH1F>("lighttageffL", "lighttageffL", 30,0., 30.);
hists_["lighttageffM"] = fs->make<TH1F>("lighttageffM", "lighttageffM", 30,0., 30.);
hists_["lighttageffT"] = fs->make<TH1F>("lighttageffT", "lighttageffT", 30,0., 30.);




//hists_["gentopmass1"] = fs->make<TH1F>("gentopmass1", "gentopmass1", 50, 0., 400.);
//hists_["gentopmass2"] = fs->make<TH1F>("gentopmass2", "gentopmass2", 50, 0., 400.);
//hists_["deltaenergy"] = fs->make<TH1F>("deltaenergy", "deltaenergy", 80, -150., 150.);

//histstwo_["pxjetcomparison"] = fs->make<TH2F>("pxjetcomparison", "pxjetcomparison", 500, -400., 400.,500, -400., 400.);
//histstwo_["pyjetcomparison"] = fs->make<TH2F>("pyjetcomparison", "pyjetcomparison", 500, -400., 400.,500, -400., 400.);
//histstwo_["pzjetcomparison"] = fs->make<TH2F>("pzjetcomparison", "pzjetcomparison", 500, -400., 400.,500, -400., 400.);

}


// ------------ method called once each job just after ending the event loop  ------------
void 
PatTopSelectionAnalyzer::endJob() 
{
}
/*
// ------------ method called when starting to processes a run  ------------
void 
Analysis::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Analysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{

}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analysis::endLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup & iSetup)
{
 using namespace edm;

  Handle<edm::MergeableCounter> numEventsStruct;
  lumi.getByLabel("eventCountProducer", numEventsStruct);

  if (numEventsStruct.isValid()) {
    totNumEvents_+=numEventsStruct->value;
    std::cout << "In endLuminosityBlock; adding "
	      << numEventsStruct->value << " to total events" << std::endl;
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
*/
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatTopSelectionAnalyzer);
