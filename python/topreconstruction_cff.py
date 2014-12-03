import FWCore.ParameterSet.Config as cms
from myanalysis.Atq.TopProducer_cfi import *
top = recoTops.clone(muonsSource ='tightMuons' , bjetsSource='bJets' ,jetsSource='goodJets', METsSource='mymet')

topreconstruction = cms.Sequence(
top 
)
