#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/PatAlgos/plugins/PATJetProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackProbabilityTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackCountingTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Selector.h"
#include "myanalysis/Atq/interface/TopProducer.h"
#include <vector>
#include <memory>
#include "DataFormats/Math/interface/LorentzVector.h"


//using namespace pat;


TopProducer::TopProducer(const edm::ParameterSet& iConfig) 
{
  // initialize the configurables
  muonsSrc_                 = iConfig.getParameter<edm::InputTag>	      ( "muonsSource" );
  jetsSrc_                 = iConfig.getParameter<edm::InputTag>	      ( "jetsSource" );
  bjetsSrc_                 = iConfig.getParameter<edm::InputTag>              ( "bjetsSource" );
  METsSrc_                 = iConfig.getParameter<edm::InputTag>	      ( "METsSource" );
  
  useNegativeDeltaSolutions_ = iConfig.getUntrackedParameter<bool> ("useNegativeDeltaSolutions",false); 
  usePositiveDeltaSolutions_ = iConfig.getUntrackedParameter<bool> ("usePositiveDeltaSolutions",true); 

  usePzMinusSolutions_ = iConfig.getUntrackedParameter<bool> ("usePzMinusSolutions",false); 
  usePzPlusSolutions_ = iConfig.getUntrackedParameter<bool> ("usePzPlusSolutions",false); 
  usePzAbsValMinimumSolutions_ = iConfig.getUntrackedParameter<bool> ("usePzAbsValMinimumSolutions",true); 

  useMetForNegativeSolutions_=iConfig.getUntrackedParameter<bool> ("useMetForNegativeSolutions",false);

  usePxMinusSolutions_ = iConfig.getUntrackedParameter<bool> ("usePxMinusSolutions",false); 
  usePxPlusSolutions_ = iConfig.getUntrackedParameter<bool> ("usePxPlusSolutions",false); 

//produces<std::vector<pat::Muon> >();
//produces<std::vector<pat::Jet> >();
//produces<std::vector<pat::MET> >();



//produces<std::vector< pat::TopLeptonic > >();
produces<std::vector< reco::NamedCompositeCandidate > >();
produces<double> (""); 
}

void TopProducer::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup){



edm::Handle<edm::View<pat::Muon> > muons;
iEvent.getByLabel(muonsSrc_,muons);


edm::Handle<edm::View<pat::Jet> > jets;
iEvent.getByLabel(jetsSrc_,jets);

edm::Handle<edm::View<pat::Jet> > bjets;
iEvent.getByLabel(bjetsSrc_,bjets);

edm::Handle<edm::View<pat::MET> > mets;
iEvent.getByLabel(METsSrc_,mets);



//std::cout<<muons->size()<<"muon size"<<jets->size()<<"jetsize"<<mets->size()<<"met size"<<mets->size()<<std::endl;

 std::vector< reco::NamedCompositeCandidate > * TopCandidates = new std::vector<reco::NamedCompositeCandidate>();

// for(size_t i = 0; i < muons->size(); ++i){
//   for(size_t j = 0; j < jets->size(); ++j){
//     for(size_t m = 0; m < mets->size(); ++m){
//std::cout<<muons->at(i).pt()<<"muonpt"<<jets->at(j).pt()<<"jetpt"<<mets->at(m).pt()<<"metpt"<<std::endl;
  if (bjets->size()==0)
{     
       reco::NamedCompositeCandidate Top,W,Nu;
       double BDISCRIMINATOR=0;
	int aaa=0;
       for(size_t j = 0; j < jets->size(); ++j){
        if (jets->at(j).bDiscriminator("combinedSecondaryVertexBJetTags")>BDISCRIMINATOR) 
	{BDISCRIMINATOR=jets->at(j).bDiscriminator("combinedSecondaryVertexBJetTags");
	 aaa=j;
	}
	}
//std::cout<<jets->size()<<"--------------------------------"<<aaa<<std::endl;

//       Top.addDaughter(muons->at(i),"Lepton");
       Top.addDaughter(muons->at(0),"Muon");
//       Top.addDaughter(jets->at(aaa),"BJet");
       Top.addDaughter(mets->at(0),"MET");
 


       std::vector<math::XYZTLorentzVector> NuMomenta = Nu4Momentum(muons->at(0),mets->at(0)).first;
       //       W.setP4(muons->at(i).p4()+mets->at(m).p4());
       if(NuMomenta.size()>0){
	 Nu.setP4((math::XYZTLorentzVector)NuMomenta.at(0));
	 W.setP4(Nu.p4()+(muons->at(0).p4()));
/*
	 for(size_t j = 0; j < jets->size(); ++j){
        if (abs(sqrt(( TMath::Power(W.energy()+jets->at(j).energy(),2))-(TMath::Power(W.px()+jets->at(j).px(),2)+TMath::Power(W.py()+jets->at(j).py(),2)+TMath::Power(W.pz()+jets->at(j).pz(),2)))-172.5)<BDISCRIMINATOR )
      {BDISCRIMINATOR=abs(sqrt(( TMath::Power(W.energy()+jets->at(j).energy(),2))-(TMath::Power(W.px()+jets->at(j).px(),2)+TMath::Power(W.py()+jets->at(j).py(),2)+TMath::Power(W.pz()+jets->at(j).pz(),2)))-172.5);
       aaa=j;
      }}
*/
       Top.addDaughter(jets->at(aaa),"BJet");

	 Top.setP4(W.p4()+jets->at(aaa).p4());
       }
//       else{W.setP4(muons->at(0).p4()+mets->at(0).p4());
 //      Top.setP4(W.p4()+jets->at(aaa).p4());

//}
 //std::cout<< "W mass" <<sqrt(( TMath::Power(W.energy(),2))-(TMath::Power(W.px(),2)+TMath::Power(W.py(),2)+TMath::Power(W.pz(),2)))<<std::endl;
//if (Top.mass()>10 && Top.mass()<130){std::cout<<"top mass is: "<< Top.mass()<<"muonpt"<<muons->at(i).pt()<<"tfgfgf"<<jets->at(j).pt()<<"jetpt"<<mets->at(m).pt()<<"metpt"<<std::endl;}
//            if(NuMomenta.size()>0)std::cout << "top mass is: "<< Top.mass()<< "RecoNuPt "<< Nu.pt()  <<std::endl;       
//     else{std::cout<< "no neutrino solution given!" << Top.mass() <<std::endl;}
       
//       std::cout<<"wmass is"<<W.mass()<<"nuotrinomass"<<Nu.mass()<< std::endl; 
       Top.addDaughter(W,"W");
       Top.addDaughter(Nu,"RecoNu");

       TopCandidates->push_back(Top);
//    std::cout << Top.daughter("RecoNu")->pz() << std::endl;
//    std::cout << "'''''''''''''''''''''''''''''''''''''''''" << std::endl;
//if (Top.mass()<130){std::cout<<"top mass is: "<< Top.mass()<<"w mass is "<<W.mass()<< "nu pz is "<<Nu.pz()<<std::endl; 
//    std::cout << "11111111111111111111111111111111111111111111" << std::endl;}

     }
 else{
       reco::NamedCompositeCandidate Top,W,Nu;

//       Top.addDaughter(muons->at(i),"Lepton");
       Top.addDaughter(muons->at(0),"Muon");
       Top.addDaughter(bjets->at(0),"BJet");
       Top.addDaughter(mets->at(0),"MET");



       std::vector<math::XYZTLorentzVector> NuMomenta = Nu4Momentum(muons->at(0),mets->at(0)).first;

       //       W.setP4(muons->at(i).p4()+mets->at(m).p4());
       if(NuMomenta.size()>0){
         Nu.setP4((math::XYZTLorentzVector)NuMomenta.at(0));
         W.setP4(Nu.p4()+(muons->at(0).p4()));
         Top.setP4(W.p4()+bjets->at(0).p4());
       }
//       else{W.setP4(muons->at(0).p4()+mets->at(0).p4());
//       Top.setP4(W.p4()+bjets->at(0).p4());

//}

//if (Top.mass()>10 && Top.mass()<130){std::cout<<"top mass is: "<< Top.mass()<<"muonpt"<<muons->at(i).pt()<<"tfgfgf"<<jets->at(j).pt()<<"jetpt"<<mets->at(m).pt()<<"metpt"<<std::endl;}
//            if(NuMomenta.size()>0)std::cout << "top mass is: "<< Top.mass()<< "RecoNuPt "<< Nu.pt()  <<std::endl;       
//     else{std::cout<< "no neutrino solution given!" << Top.mass() <<std::endl;}

//       std::cout<<"wmass is"<<W.mass()<<"nuotrinomass"<<Nu.mass()<< std::endl; 
       Top.addDaughter(W,"W");
       Top.addDaughter(Nu,"RecoNu");
       TopCandidates->push_back(Top);
//    std::cout << Top.daughter("RecoNu")->pz() << std::endl;
//    std::cout << "'''''''''''''''''''''''''''''''''''''''''" << std::endl;
//if (Top.mass()<130){std::cout<<"top mass is: "<< Top.mass()<<"w mass is "<<W.mass()<< std::endl; }
//if (Top.mass()<130){std::cout<<"top mass is: "<< Top.mass()<<"w mass is "<<W.mass()<< "nu pz is "<<Nu.pz()<<std::endl;
// std::cout << "22222222222222222222222222222222222222222" << std::endl;}

     }

 
 
 std::auto_ptr< std::vector< reco::NamedCompositeCandidate > > newTopCandidate(TopCandidates);
 std::auto_ptr<double> soluinfo(new double(Nu4Momentum(muons->at(0),mets->at(0)).second));

 
////////

//iEvent.put(newTopLeptonic);

iEvent.put(newTopCandidate);
iEvent.put(soluinfo);

//    std::cout << TopCandidates[0].daughter("W")->energy() << std::endl;
//    std::cout << "'''''''''''''''''''''''''''''''''''''''''" << std::endl;

}

TopProducer::~TopProducer(){;}

std::pair<std::vector<math::XYZTLorentzVector> ,double> TopProducer::Nu4Momentum(const reco::Candidate & Lepton,const reco::Candidate & MET){

  double  mW = 80.38;

  std::vector<math::XYZTLorentzVector> result;
  
  //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );
  double Nuinfo;  
  double MisET2 = (MET.px()*MET.px() + MET.py()*MET.py());
  double mu = (mW*mW)/2 + MET.px()*Lepton.px() + MET.py()*Lepton.py();
  double a  = (mu*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
  double a2 = TMath::Power(a,2);
  double b  = (TMath::Power(Lepton.energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.energy(),2) - TMath::Power(Lepton.pz(),2));
  double pz1(0),pz2(0),pznu(0);
//  int nNuSol(0);

  math::XYZTLorentzVector p4nu_rec;
  math::XYZTLorentzVector p4W_rec;
 math::XYZTLorentzVector p4b_rec;
  math::XYZTLorentzVector p4Top_rec;
  math::XYZTLorentzVector p4lep_rec;    

  p4lep_rec.SetPxPyPzE(Lepton.px(),Lepton.py(),Lepton.pz(),Lepton.energy());
  
  math::XYZTLorentzVector p40_rec(0,0,0,0);

  if(a2-b > 0 ){
    if(!usePositiveDeltaSolutions_)
      {
std::cout << " test1 " << std::endl;
	result.push_back(p40_rec);
	return std::make_pair(result,0);
      }
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
//    nNuSol = 2;     
  
  
//    std::cout << " test2 " << std::endl;

    if(usePzPlusSolutions_)pznu = pz1;    
    if(usePzMinusSolutions_)pznu = pz2;
    if(usePzAbsValMinimumSolutions_){
      pznu = pz1;
      if(fabs(pz1)>fabs(pz2)) pznu = pz2;
    }
    

  double Enu = sqrt(MisET2 + pznu*pznu);
  
  p4nu_rec.SetPxPyPzE(MET.px(), MET.py(), pznu, Enu);
//    p4nu_rec.SetPxPyPzE(0, 0, 0, 0);
  
  result.push_back(p4nu_rec);
  Nuinfo=0;
  }
  else{
pznu=a;
double Enuu= sqrt(MisET2 + pznu*pznu);
  p4nu_rec.SetPxPyPzE(MET.px(), MET.py(), pznu, Enuu);
//  p4nu_rec.SetPxPyPzE(0, 0, 0, 0);
  Nuinfo=1;

  result.push_back(p4nu_rec);
//    std::cout << " test3 " << std::endl;

}
/*
    if(!useNegativeDeltaSolutions_){
std::cout << " test3 " << std::endl;
      result.push_back(p40_rec);
      return result;
    }
    //    double xprime = sqrt(mW;


    double ptlep = Lepton.pt(),pxlep=Lepton.px(),pylep=Lepton.py(),metpx=MET.px(),metpy=MET.py();

    double EquationA = 1;
    double EquationB = -3*pylep*mW/(ptlep);
    double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
    double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

    std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

    std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);

    
    double deltaMin = 14000*14000;
    double zeroValue = -mW*mW/(4*pxlep); 
    double minPx=0;
    double minPy=0;

    //    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl; 
      
    if(usePxMinusSolutions_){
std::cout << " test4 " << std::endl;
      for( int i =0; i< (int)solutions.size();++i){
      if(solutions[i]<0 ) continue;
      double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep); 
      double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
      double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 

      //      std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl; 

      if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
      minPx=p_x;
      minPy=p_y;}
      //     std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    } 
    
    if(usePxPlusSolutions_){
std::cout << " test5 " << std::endl;
      for( int i =0; i< (int)solutions2.size();++i){
	if(solutions2[i]<0 ) continue;
	double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep); 
	double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
	double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 
	//  std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
	if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
	  minPx=p_x;
	  minPy=p_y;
	}
	//	std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    }
    
    double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
    double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);
   std::cout << " test888888888888888888888 " << std::endl;
 
    if(deltaMin==14000*14000)return result;    
    //    else std::cout << " test " << std::endl;
std::cout << " test9999999999999999999999 " << std::endl;

    if(delta2ZeroValue < deltaMin){
      deltaMin = delta2ZeroValue;
      minPx=zeroValue;
      minPy=pyZeroValue;}

    //    std::cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
    ///    ////Y part   

    double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
    double a_Minimum  = (mu_Minimum*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
    pznu = a_Minimum;
    
    if(!useMetForNegativeSolutions_){
std::cout << " test6 " << std::endl;
      double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
      p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
    }
    else{
std::cout << " test7" << std::endl;

      pznu = a;
      double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
      p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
    }
    result.push_back(p4nu_rec);}
*/
  return std::make_pair(result,Nuinfo);    
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopProducer);

