import FWCore.ParameterSet.Config as cms
from myanalysis.Atq.eostools import *


#Define the process
process = cms.Process("Aqt")
process.load("FWCore.MessageService.MessageLogger_cfi")

#Define the input sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()

#    fileNames = cms.untracked.vstring("root://eoscms//eos/cms/store/user/rgoldouz/REALDATA/susypat_1_1_USb.root")
#    fileNames = cms.untracked.vstring("root://eoscms//eos/cms/store/user/rgoldouz/pat_LL_005/pattuples_LL_005_2.root")

                            )

f = open('qcd.txt', 'r')
inputs = (x for x in f.readlines())
files = [file.strip() for file in inputs]
for file in files:
	process.source.fileNames.append( "root://eoscms//eos/cms/store/user/mojtaba/QCD/" + file )
f.closed

Filest = listFiles("/eos/cms/store/user/rgoldouz/QCD/"
)
for file in Filest:
    process.source.fileNames.append( "root://eoscms/" + file )

#event number cout
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
## restrict the number of events for testing
process.maxEvents = cms.untracked.PSet( 
input = cms.untracked.int32(-1)
)

## Define the TFileService
process.TFileService =cms.Service("TFileService",
      fileName =cms.string("QCD.root")
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
process.triggerstep  = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_IsoMu24_eta2p1_v*"])
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
process.tightphotonstep  = countPatPhotons.clone(src = 'Finalphoton', minNumber = 1 , maxNumber = 1)

## Different jet multiplicity selections
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
process.goodjetstep = countPatJets.clone(src = 'goodJets'   , minNumber = 1 )
process.goodbjetstep = countPatJets.clone(src = 'bJets'   , minNumber =0  ,maxNumber = 1)

##met selection
from PhysicsTools.PatAlgos.selectionLayer1.metCountFilter_cfi import *
process.metstep= countPatMET.clone(src= 'mymet', minNumber = 1)

##electron veto
from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *
process.myelecveto=countPatElectrons.clone(src= 'myelectron', maxNumber =0) 

##------------------------------------------------------------------
process.load("myanalysis.Atq.topreconstruction_cff")

from myanalysis.Atq.Topmassfilter_cfi import *
process.topmassdisinwindow=topFilter.clone(istopwindow=True)
process.topmassdisoutwindow=topFilter.clone(istopwindow=False)

#topmassdis= cms.Sequence(
#topmassdiscriminator
#)

from myanalysis.Atq.PatTopSelectionAnalyzer_cfi import *
from myanalysis.Atq.Firstanalyse_cfi import *

process.analyzestep0  = firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' ,realdata=True )
process.analyzestep1  = firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' )
process.analyzestep2  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='patJetsAK5PFPt30',met='patMETsPF' ,tops='top')
process.Signalregion  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='patJetsAK5PFPt30',met='patMETsPF' ,tops='top')
process.Sidebandregion  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='patJetsAK5PFPt30',met='patMETsPF' ,tops='top')



## Switch output report on
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


#Define loose event selection path
process.signalreg = cms.Path(
    process.triggerstep *
    process.vertexstep   *
    process.tightmuonstep  *
    process.vetoloosmuonstep   *
    process.myelecveto   *
    process.metstep   *
    process.tightphotonstep  *
#   process.analyzestep0   *
#   process.analyzestep1   *
    process.goodjetstep   *
    process.goodbjetstep  *
    process.topreconstruction * 
    process.analyzestep2 *
    process.topmassdisinwindow  *
    process.Signalregion
    )


process.sidebandreg = cms.Path(
    process.triggerstep *
    process.vertexstep   *
    process.tightmuonstep  *
    process.vetoloosmuonstep   *
    process.myelecveto   *
    process.metstep   *
    process.tightphotonstep  *
#   process.analyzestep0   *
#   process.analyzestep1   *
    process.goodjetstep   *
    process.goodbjetstep  *
    process.topreconstruction *
#    process.analyzestep2 *
    process.topmassdisoutwindow *
    process.Sidebandregion
    )

#process.schedule = cms.Schedule(process.topObjectProduction,process.initialEventSelection,process.Signalre,process.Sidebandre) 


