#Author:  Reza Goldouzian,

import FWCore.ParameterSet.Config as cms
from myanalysis.Atq.eostools import *

#Define the process
process = cms.Process("Aqt")
process.load("FWCore.MessageService.MessageLogger_cfi")

#Define the input sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()

#    fileNames = cms.untracked.vstring("root://eoscms//eos/cms/store/user/rgoldouz/REALDATA/susypat_1_1_USb.root")
                            )
#Files = listFiles("/eos/cms/store/user/mojtaba/REALDATARunB/"
#)
#for file in Files:
#    process.source.fileNames.append( "root://eoscms/" + file )

f = open('realdata.txt', 'r')
inputs = (x for x in f.readlines())
files = [file.strip() for file in inputs]
for file in files:
	process.source.fileNames.append( "root://eoscms//eos/cms/store/user/mojtaba/REALDATARunB/" + file )
f.closed


#event number cout
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
## restrict the number of events for testing
process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(10)
)

## Define the TFileService
process.TFileService =cms.Service("TFileService",
      fileName =cms.string("TEST.root")
)

## ----------------------------------------------------------------
## Apply object selection according to TopPAG reference selection
## for ICHEP 2010. This will result in 5 additional collections:
##
## * goodJets
## * tightphoton
## * vetoMuons
## * looseMuons
## * tightMuons
## * met
## Have a look ont the cff file to learn more about the exact
## selection citeria.
## ----------------------------------------------------------------
process.load("myanalysis.Atq.topObjectSelection_cff")
process.topObjectProduction = cms.Path(
    process.topObjectSelection
)

## Trigger bit (HLT_mu9)
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerstep  = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_Mu17_v3"])
## Vertex requirement
process.vertexstep  = cms.EDFilter("VertexSelector", src = cms.InputTag("offlinePrimaryVertices"), cut = cms.string("!isFake && ndof > 4 && abs(z) < 15 && position.Rho < 2"), filter = cms.bool(True))
## Exact one tight muon
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
process.tightmuonstep = countPatMuons.clone(src = 'tightMuons', minNumber = 1, maxNumber = 1)
## Exact one loose muon
process.loosemuonstep = countPatMuons.clone(src = 'looseMuons', minNumber = 1, maxNumber = 1)
## Veto on additional muons
process.vetoloosmuonstep  = countPatMuons.clone(src = 'vetoMuons' , maxNumber = 1)
## Exact one tight phton
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *
process.tightphotonstep  = countPatPhotons.clone(src = 'tightphoton', minNumber = 1 , maxNumber = 1)
## Different jet multiplicity selections
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
process.goodjetstep = countPatJets.clone(src = 'goodJets'   , minNumber = 1)
process.goodbjetstep = countPatJets.clone(src = 'bJets'   , minNumber = 1)

##met selection
from PhysicsTools.PatAlgos.selectionLayer1.metCountFilter_cfi import *
process.metstep= countPatMET.clone(src= 'mymet', minNumber = 1)

##------------------------------------------------------------------

from myanalysis.Atq.PatTopSelectionAnalyzer_cfi import *
#process.analyzestep  = analyzePatTopSelection.clone(photonsrc='tightphoton' ,muons='tightMuons', jets='goodJets')
process.analyzestep  = analyzePatTopSelection.clone(photonsrc='tightphoton' ,muons='tightMuons', jets='bjets')


## Switch output report on
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Define loose event selection path
process.looseEventSelection = cms.Path(
    
    process.vertexstep   *
    process.tightmuonstep  *
    process.vetoloosmuonstep   *
    process.metstep   *
    process.tightphotonstep  *
#    process.goodjetstep  *
    process.goodbjetstep  *
    process.analyzestep

    )


















