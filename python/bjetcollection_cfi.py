import FWCore.ParameterSet.Config as cms
selectedbJets = cms.EDProducer("bjetcollection",
  JetCollection = cms.InputTag("patJetsAK5PF"),
  inputMETs            = cms.InputTag("patMETsPF"),
  MinJetPt      = cms.double(30),
  MaxJetEta     = cms.double(2.5),
#  bdis      = cms.double(0.244)
   bdis      = cms.double(0.679),
   JESscaleType            = cms.string("abs"), #abs or rel(*eta) or jes:up / jes:down (pt-dependend)
   JERscaleType            = cms.string("abs"), #abs or rel(*eta) or jer:up / jer:down (pt-dependend)
   resolutionFactors    = cms.vdouble(1.052,1.057,1.096,1.134,1.288), # list the different JER factors here: (JER1, JER2)
   resolutionFactorsup    = cms.vdouble(1.126,1.125,1.176,1.256,1.57), # list the different JER factors here: (JER1, JER2)
   resolutionFactorsdown    = cms.vdouble(0.979,0.99,1.017,1.014,1.008), # list the different JER factors here: (JER1, JER2)
   resolutionEnergyRanges= cms.vdouble(50,70,70,100,100,200,200,300,300,500),
   sigmaMC = cms.vdouble(0.1156,0.1288,0.1455,0.1266,0.1346,0.1106,0.1113,0.1166,0.1035,0.1133,0.0921,0.0939,0.0970,0.0814,0.0878,0.0785,0.0737,0.0754,0.0624,0.0668,0.0648,0.0704,0.0685,0.0545,0.0652),
   resolutionEtaRanges  = cms.vdouble(0,0.5,0.5,1.1,1.1,1.7,1.7,2.3,2.3,5),  # list the |eta| ranges for the different JER factors here (etaMin1, etaMax1, etaMin2, etaMax2), etaMax=-1: means |eta|<infinity---https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    realdata   = cms.bool(False),

   JECUncSrcFile        = cms.FileInPath("myanalysis/Atq/plugins/Summer13_V4_MC_Uncertainty_AK5PF.txt")

#   bdis      = cms.double(0.898)

)

