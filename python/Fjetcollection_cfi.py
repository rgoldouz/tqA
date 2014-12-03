import FWCore.ParameterSet.Config as cms
FselectedbJets = cms.EDProducer("Fjetcollection",
  JetCollection = cms.InputTag("patJetsAK5PFPt30"),
  MinJetPt      = cms.double(30),
  MaxJetEta     = cms.double(2.5),
  bscalefactorType            = cms.string("abs"), #abs or or bsf:up / bsf:down
  realdata   = cms.bool(False),
 
#  bdis      = cms.double(0.244)
   bdis      = cms.double(0.679)
#   bdis      = cms.double(0.898)

)

