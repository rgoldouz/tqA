import FWCore.ParameterSet.Config as cms
selectedmets = cms.EDProducer("metcollection",
  JetCollection = cms.InputTag("patMETsPF"),
  MinJetPt      = cms.double(30),
)

