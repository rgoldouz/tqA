import FWCore.ParameterSet.Config as cms
topFilter = cms.EDFilter(
  "TopFilter",
  TopSource = cms.InputTag("top"),
  Minmass = cms.double(130.0),
  Maxmass = cms.double(220.0),
  istopwindow= cms.bool(True)

)
