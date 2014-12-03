import FWCore.ParameterSet.Config as cms

electronSelector = cms.EDFilter(
  "ElectronSelector",
  ElectronSource    = cms.InputTag('patElectronsIDIso'),
  MinElePt       = cms.double(20),
  MaxEleEta      = cms.double(2.5)

)

