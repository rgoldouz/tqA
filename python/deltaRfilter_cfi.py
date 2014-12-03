import FWCore.ParameterSet.Config as cms

deltaRfilter=cms.EDFilter(
"DeltaRfilter",
  TopSource = cms.InputTag("top"),
  muonSource=cms.InputTag("tightMuons"),
  photonSource=cms.InputTag("Finalphoton"),
  Mindr = cms.double(0.7)
)
 

