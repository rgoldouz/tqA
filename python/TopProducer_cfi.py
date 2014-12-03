import FWCore.ParameterSet.Config as cms
recoTops = cms.EDProducer("TopProducer",
                                  electronsSource = cms.InputTag("topElectrons"),
                                  muonsSource = cms.InputTag("topMuons"),
                                  jetsSource = cms.InputTag("goodJets"),
                                  bjetsSource = cms.InputTag("bJets"),
                                  METsSource = cms.InputTag("preselectedMETs"),
                                )
