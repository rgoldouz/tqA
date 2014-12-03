#Author:  Reza Goldouzian,

import FWCore.ParameterSet.Config as cms
#from myanalysis.Atq.eostools import *

#Define the process
process = cms.Process("Aqt")
process.load("FWCore.MessageService.MessageLogger_cfi")

#Define the input sample
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
#    fileNames = cms.untracked.vstring("root://eoscms//eos/cms/store/user/rgoldouz/REALDATA/susypat_1_1_USb.root")
#File2 = listFiles("/eos/cms/store/user/rgoldouz/MuEG_Run2012A_recover_06Aug2012_v1_AOD/")
#for file in File2:
#       process.source.fileNames.append( "root://eoscms/" + file )


#File1 = listFiles("/eos/cms/store/user/rgoldouz/MuEG_Run2012A_13Jul2012_v1_AOD/")
#for file in File1:
#        process.source.fileNames.append( "root://eoscms/" + file )

d = open('runD.txt', 'r')
dinputs = (x for x in d.readlines())
dfiles = [file.strip() for file in dinputs]
for file in dfiles:
        process.source.fileNames.append( " rfio:///dpm/particles.ipm.ac.ir/home/cms/store/user/goldouzian/SingleMu_Run2012D_PromptReco_v1_AOD/" + file )
d.closed


#File2 = listFiles("/eos/cms/store/user/rgoldouz/MuEG_Run2012A_recover_06Aug2012_v1_AOD/")
#for file in File2:
#	process.source.fileNames.append( "root://eoscms/" + file )




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
      fileName =cms.string("REALDATA3.root")
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
process.YJET.realdata=True
process.goodJets.realdata=True
process.bJets.realdata=True


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

process.STEP1=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons',realdata=True  )
process.STEP2=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' ,realdata=True)
process.STEP3=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' ,realdata=True)
process.STEP4=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons',realdata=True )
process.STEP5=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons',realdata=True )
process.STEP6=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' , tops='top', step='afterb',realdata=True)
process.STEP7=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' , tops='top', step='afterb',realdata=True)
process.STEP8=firstselection.clone(photonsrc='Finalphoton' ,muons='tightMuons' , tops='top', step='afterb',realdata=True)


process.firstanalyzestep = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='goodJets',mybjets='bJets',met='patMETsPF' ,tops='top',realdata=True)

process.analyzestep2  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='goodJets',mybjets='bJets',met='patMETsPF' ,tops='top',realdata=True)
process.Signalregion  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='goodJets',mybjets='bJets',met='patMETsPF' ,tops='top',treename='insidetopmass',realdata=True)
process.Sidebandregion  = analyzePatTopSelection.clone(photonsrc='Finalphoton' ,muons='tightMuons', jets='goodJets',mybjets='bJets',met='patMETsPF' ,tops='top',treename='outsidetopmass',realdata=True)


from myanalysis.Atq.deltaRfilter_cfi import *
process.deltarconstraint =deltaRfilter.clone()


## Switch output report on
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_P_V42_AN2::All' 
#process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")
#process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
#process.p = cms.Path(process.ecalLaserCorrFilter) 

#Define loose event selection path
process.signalreg = cms.Path(
#    process.ecalLaserCorrFilter *
    process.triggerstep *
    process.vertexstep   *
    process.tightmuonstep  *
    process.tightphotonstep  *
    process.STEP1 *
    process.vetoloosmuonstep   *
    process.STEP2 *
    process.myelecveto   *
    process.STEP3 *
    process.metstep   *
    process.STEP4 *
#   process.analyzestep0   *
#   process.analyzestep1   *
    process.goodjetstep   *
    process.STEP5 *
    process.goodbjetstep  *
    process.topreconstruction * 
    process.STEP6 *
    process.firstanalyzestep  *
    process.deltarconstraint  *
    process.STEP7 *
    process.analyzestep2 
#    process.topmassdisinwindow  *
#    process.STEP8 *
#    process.Signalregion
    )


#process.sidebandreg = cms.Path(
#    process.triggerstep *
#    process.vertexstep   *
#    process.tightmuonstep  *
#    process.tightphotonstep  *
#    process.STEP1 *
#    process.vetoloosmuonstep   *
#    process.STEP2 *
#    process.myelecveto   *
#    process.STEP3 *
#    process.metstep   *
#    process.STEP4 *
#   process.analyzestep0   *
#   process.analyzestep1   *
#    process.goodjetstep   *
#    process.STEP4 *
#    process.goodbjetstep  *
#    process.topreconstruction *
#    process.STEP5 *
#    process.firstanalyzestep  *
#    process.deltarconstraint  *
#    process.STEP6 *
#    process.analyzestep2 *
#    process.topmassdisoutwindow *
#    process.STEP7 *
#    process.Sidebandregion
#    )

#process.schedule = cms.Schedule(process.topObjectProduction,process.initialEventSelection,process.Signalre,process.Sidebandre) 










