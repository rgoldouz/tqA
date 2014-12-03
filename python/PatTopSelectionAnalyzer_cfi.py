import FWCore.ParameterSet.Config as cms

analyzePatTopSelection = cms.EDAnalyzer("PatTopSelectionAnalyzer",
    photonsrc = cms.untracked.InputTag("patPhotons"),
    muons = cms.untracked.InputTag("selectedPatMuonsPF"),                                             
    jets  = cms.untracked.InputTag("selectedPatJetsPF"),
    mybjets  = cms.untracked.InputTag("selectedPatJetsPF"),
    met   = cms.untracked.InputTag("patMETs"),
    tops = cms.untracked.InputTag("tops"),
    pileupw = cms.InputTag("pileupweight") ,
    realdata   = cms.bool(False),
    treename =cms.string("atq"),
    sample=cms.string("atq")
 

)
