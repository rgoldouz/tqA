import FWCore.ParameterSet.Config as cms

weightProducer = cms.EDProducer('WeightProducer',
FileNamePUDataDistribution = cms.string("myanalysis/Atq/plugins/MyJAN2013DataPileupHistogram.root")
)
