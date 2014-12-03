import FWCore.ParameterSet.Config as cms

firstselection = cms.EDAnalyzer("Firstanalyse",
    photonsrc = cms.untracked.InputTag("patPhotons"),
    muons = cms.untracked.InputTag("selectedPatMuonsPF"),
    jets  = cms.untracked.InputTag("goodJets"),
    met   = cms.untracked.InputTag("patMETsPF"),
    tops = cms.untracked.InputTag("tops"),
    pileupw = cms.InputTag("pileupweight") ,
    realdata   = cms.bool(False),
    sample=cms.string("atq"),
    treename =cms.string("FCNC"),
    step=cms.string("step")	


)

