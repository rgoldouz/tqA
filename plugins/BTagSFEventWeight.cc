#include "FWCore/Utilities/interface/EDMException.h"
#include "myanalysis/Atq/interface/BTagSFEventWeight.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
 

BTagSFEventWeight::BTagSFEventWeight(const edm::ParameterSet& cfg):
jets_                   ( cfg.getParameter<edm::InputTag>    ( "jets"   ) ),
sysVar_                 ( cfg.getParameter<std::string>      ("sysVar"  ) ),
filename_               ( cfg.getParameter<edm::FileInPath>  ("filename") )
{
produces<double>();
   
   // set the edges of the last histo bin
   maxPtDB_     = 240.;
   maxPtMisTag_ = 520.;
   maxEta_      = 2.4;
   maxPt2012_  = 800.; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript
 
   // load TFile Service
   edm::Service<TFileService> fs;
 file_ = new TFile((TString)filename_.fullPath());
 effHists_["EffBJetsTaggedPtEta"] = (TH2F*) file_->Get("btagRECO")->Clone();
 effHists_["EffCJetsTaggedPtEta"] = (TH2F*) file_->Get("ctagRECO")->Clone();
 effHists_["EffLJetsTaggedPtEta"] = (TH2F*) file_->Get("lightRECO")->Clone();

}
BTagSFEventWeight::~BTagSFEventWeight()
 {
 }
 
 void BTagSFEventWeight::produce(edm::Event& evt, const edm::EventSetup& setup)
 {
edm::Handle<edm::View< pat::Jet > > jets;
   evt.getByLabel(jets_, jets);
  double pt, eta;
  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0);
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0);
     for(edm::View<pat::Jet>::const_iterator jet = jets->begin();jet != jets->end(); ++jet) {
       pt  = jet->pt();
       eta = std::abs(jet->eta());
       if(jet->partonFlavour() == 5 || jet->partonFlavour() == -5){
         oneMinusBEffies               .push_back(1.- effBTag(pt, eta));
         oneMinusBEffies_scaled        .push_back(1.- (effBTag(pt, eta) * effBTagSF(pt, eta, false)));
//   std::cout<<"  eff B-tag   "<<1.- effBTag(pt, eta)<<std::endl;   
       }
   
       else if(jet->partonFlavour() == 4 || jet->partonFlavour() == -4){
         oneMinusBMistags               .push_back(1.- effBTagCjet(pt, eta));
         oneMinusBMistags_scaled        .push_back(1.-(effBTagCjet(pt, eta) * effBTagSF(pt, eta, true))); // ATTENTION: btag SF used for c-jets to with 2x the uncertainty!
//   std::cout<<"  eff c-tag   "<<1.- effBTagCjet(pt, eta)<<std::endl;

       }
  
       else{
         oneMinusBMistags               .push_back(1.- effMisTag(pt, eta));
         oneMinusBMistags_scaled        .push_back(1.-(effMisTag(pt, eta) * effMisTagSF(pt, eta)));
//   std::cout<<"  eff light-tag   "<<1.- effMisTag(pt, eta)<<std::endl;

       }
    }
      
    double effBTagEvent_unscaled = effBTagEvent( oneMinusBEffies, oneMinusBMistags );
    double effBTagEvent_scaled   = effBTagEvent( oneMinusBEffies_scaled, oneMinusBMistags_scaled );
    double effBTagEventSF = effBTagEvent_scaled / effBTagEvent_unscaled;
   
   std::auto_ptr<double> bTagSFEventWeight(new double);
   *bTagSFEventWeight = effBTagEventSF;    
   evt.put(bTagSFEventWeight);  
 }

// b tag eff. from MC as a function of jet pt, eta
 double BTagSFEventWeight::effBTag(double jetPt, double jetEta)
 {
   double result = -1111.;
  // if histo file exists, take value from there; else return a default value
   if(filename_.location()) {
   TH2F* his = effHists_.find("EffBJetsTaggedPtEta")->second;
     // ensure that pt is in accepted range of BTV DB
     if(jetPt  >= maxPt2012_) jetPt  = maxPt2012_-1.;
     if(jetEta >= maxEta_    ) jetEta = maxEta_-0.1;
     result = his->GetBinContent( his->FindBin(jetPt, jetEta) );
   }
   else {result = 0.7; std::cout<< "WARNING!!! b tag eff. is ALWAYS 0.7!!! CHECK!!!"<<std::endl; }
   return result;
 }

// b tag eff. SF as a function of jet pt, eta
 double BTagSFEventWeight::effBTagSF(double jetPt, double jetEta, bool isCjet)
 {
  double result = -1111., error = -1111.;
  result = effBTagSF2012(jetPt);
  error = effBTagSFerr2012(jetPt);
  if(isCjet) error*=2;
 if(sysVar_ == "bTagSFUp")   result += error;
 else if(sysVar_ == "bTagSFDown") result -= error;
 return result;
}

double BTagSFEventWeight::effBTagSF2012(double x)
{
//https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
return (0.939158+(0.000158694*x))+(-2.53962e-07*(x*x));
}

 double BTagSFEventWeight::effBTagSFerr2012(double x)
 {
   // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript (dataset=ABCD); x = jetPt
   // pt binning
   double pt[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
  // corresponding SFb uncertainties
  double SFb_errorCSVM[] = {
     0.0415694,
     0.023429,
     0.0261074,
     0.0239251,
     0.0232416,
     0.0197251,
     0.0217319,
     0.0198108,
     0.0193,
     0.0276144,
     0.0205839,
     0.026915,
     0.0312739,
     0.0415054,
     0.0740561,
     0.0598311 };
 
   int iBin = -1;
   for(int i=0; i<15; i++) {
    if (x>pt[i] && x<pt[i+1]) {
      iBin =i;
       break;
     }
   }
   double factor = 1.;
   if(iBin<0){
     // outside the quoted pt range: use twice the error
    factor=2;
     // if pt>800: use SFb(800)
     if(x>800) iBin=14;
     // if pt<20: use SFb(20)
     if(x<20 ) iBin=0;
   } 
   return factor * SFb_errorCSVM[iBin];
 }

 // b tag eff. from MC for c jets as a function of jet pt, eta;
 // as first step: take average of b and mis eff.
 double BTagSFEventWeight::effBTagCjet(double jetPt, double jetEta)
 {
   double result = -1111.;
  // if histo file exists, take value from there; else return a default value
   if(filename_.location()) {
     TH2F* his = effHists_.find("EffCJetsTaggedPtEta")->second;
     // ensure that pt is in accepted range
    if(jetPt >= maxPt2012_) jetPt = maxPt2012_-1.;
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
     result = his->GetBinContent( his->FindBin(jetPt,jetEta) );
   }
   else {result = 0.35; std::cout<< "WARNING!!! c tag eff. is ALWAYS 0.35!!! CHECK!!!"<<std::endl; }
   return result;
 }

 
 // mistag eff. from MC as a function of jet pt, eta
 double BTagSFEventWeight::effMisTag(double jetPt, double jetEta)
 {
   double result = -1111.;
  // if histo file exists, take value from there; else return a default value
   if(filename_.location()) {
    TH2F* his = effHists_.find("EffLJetsTaggedPtEta")->second;
    // ensure that pt is in accepted range
    if(jetPt >= maxPtMisTag_) jetPt = maxPtMisTag_-1.;
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
    result = his->GetBinContent( his->FindBin(jetPt, jetEta) );
   }
   else {result = 0.01; std::cout<< "WARNING!!! b tag eff. is ALWAYS 0.01!!! CHECK!!!"<<std::endl; }
   return result;
 }

double BTagSFEventWeight::effMisTagSF(double jetPt, double jetEta)
 {
   double result = -1111., error = -1111.;
    if(jetEta >= maxEta_) jetEta = maxEta_-0.1;
     TString                       meanminmax = "mean";
     if(sysVar_ == "misTagSFUp"  ) meanminmax = "max";
     if(sysVar_ == "misTagSFDown") meanminmax = "min";
     result = effMisTagSF2012(jetPt, jetEta, meanminmax);
     return result;
}

double BTagSFEventWeight::effMisTagSF2012(double x, double jetEta, TString meanminmax)
 {
   // function from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript; dataset=ABCD 
   // x = jetPt; meanminmax = "mean" -> central value; = "min" -> down variation; = "max" -> up variation
     double val =0, valMean=0, valMin=0, valMax=0;
     bool outOfRange = false;
     if(x<20) {x=20; outOfRange = true;}
     if(jetEta>=0. && jetEta <=0.8){
       if(x>1000.) {x=1000.; outOfRange = true;}
       valMean= ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
       valMin= ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
       valMax= ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
     }
     else if(jetEta>0.8 && jetEta <=1.6){
       if(x>1000.) {x=1000.; outOfRange = true;}
       valMean= ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
       valMin= ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
       valMax= ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
    }
     else if(jetEta>1.6 && jetEta <=2.4){
       if(x>850.) {x=850.; outOfRange = true;}
       valMean= ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
       valMin= ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
       valMax= ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
     }
       if( meanminmax == "mean") val= valMean;
       if( meanminmax == "min" ) val= outOfRange? valMin-(valMean-valMin) : valMin; // if outOfRange -> 2x the uncertainty
       if( meanminmax == "max" ) val= outOfRange? valMax+(valMax-valMean) : valMax; // if outOfRange -> 2x the uncertainty
     if(val!=0) return val;
     else { 
     return 1.; 
   }
}

//btag SF for 0 tag+ 1 tag catagory, for doing this I follow the recommended method for >=2 tag and consider 1- >=2tag
double BTagSFEventWeight::effBTagEvent(std::vector<double> &oneMinusBEffies,std::vector<double> &oneMinusBMistags)
 {
   double bTaggingEfficiency = 1.;
   double tmp = 1.;

   for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
     tmp *= (*eff);}
   for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
     tmp *= (*mis);}
   bTaggingEfficiency -= tmp;
   for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
     tmp = 1.-(*eff);
     for(std::vector<double>::const_iterator eff2 =oneMinusBEffies.begin(); eff2 != oneMinusBEffies.end(); ++eff2){
       if(eff != eff2) tmp *= (*eff2);
     }
     for(std::vector<double>::const_iterator mis =oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
       tmp *= (*mis);
     }
     bTaggingEfficiency -= tmp;
   }
   for(std::vector<double>::const_iterator mis =oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
    tmp = 1.-(*mis);
    for(std::vector<double>::const_iterator eff =oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
       tmp *= (*eff);
     }
     for(std::vector<double>::const_iterator mis2 =oneMinusBMistags.begin(); mis2 != oneMinusBMistags.end(); ++mis2){
       if(mis != mis2) tmp *= (*mis2);
     }
    bTaggingEfficiency -= tmp;
   }
   return 1-bTaggingEfficiency;
 }
#include "FWCore/Framework/interface/MakerMacros.h"
 DEFINE_FWK_MODULE( BTagSFEventWeight ); 
