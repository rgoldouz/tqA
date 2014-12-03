import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.muonCleaner_cfi import *
looseMuons = cleanPatMuons.clone(
src = 'patMuonsPFID',
    preselection =
    'pt > 20. &'
    'abs(eta) < 2.1 &'
    'pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5*pfIsolationR04.sumPUPt)/ pt < 0.2'
)
#---------------------------------
tightMuons = cleanPatMuons.clone(
    src = 'looseMuons',
    preselection =
    'pt > 26. &'
    'abs(eta) < 2.1 &'
   'pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5*pfIsolationR04.sumPUPt)/ pt <0.12'
#    'pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5*pfIsolationR04.sumPUPt)/ pt > 0.3'
)
#--------------------------------
vetoMuons =cleanPatMuons.clone(
src = 'patMuonsPFID',
    preselection =
    'pt > 10. &'
    'abs(eta) < 2.5 &'
    'pfIsolationR04.sumChargedHadronPt + max(0., pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt - 0.5*pfIsolationR04.sumPUPt)/ pt < 0.2'
)
#----------------------------------
from PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi import *
from myanalysis.Atq.photoncollection_cfi import *
tightphoton = myphoton.clone()
mytightphoton = myselectedphoton.clone()
#photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && userInt("passEleConvVeto") > 0 && hadTowOverEm < userFloat("hadTowOverEmTightCut") && sigmaIetaIeta > userFloat("showerShapeTightCut")+0.05*userFloat("showerShapeTightCut")  && userFloat("pfChargedPU") < userFloat("pfChargedTightCut") && userFloat("pfNeutralPU") <userFloat("pfNeutralTightCut") && userFloat("pfGammaPU") <userFloat("pfGammaTightCut")'    
#                              )

#photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && userInt("passEleConvVeto") > 0 && hadTowOverEm < userFloat("hadTowOverEmTightCut") && sigmaIetaIeta > userFloat("showerShapeTightCut")+0.1*userFloat("showerShapeTightCut")  && userFloat("pfChargedPU") < userFloat("pfChargedTightCut") && userFloat("pfNeutralPU") <userFloat("pfNeutralTightCut") && userFloat("pfGammaPU") <userFloat("pfGammaTightCut")'    
#                              )

#photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && userInt("passEleConvVeto") > 0 && hadTowOverEm < userFloat("hadTowOverEmTightCut") && sigmaIetaIeta > userFloat("showerShapeTightCut")  && userFloat("pfChargedPU") < userFloat("pfChargedTightCut") && userFloat("pfNeutralPU") <userFloat("pfNeutralTightCut") && userFloat("pfGammaPU") <userFloat("pfGammaTightCut")'

#photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) ' 
photonIDCutTight  = cms.string('et > 50.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && userInt("passEleConvVeto") > 0 && hadTowOverEm < userFloat("hadTowOverEmTightCut") && sigmaIetaIeta < userFloat("NEWshowerShape")  && userFloat("pfChargedPU") < userFloat("pfChargedTightCut") && userFloat("pfNeutralPU") <userFloat("pfNeutralTightCut") && userFloat("pfGammaPU") <userFloat("pfGammaTightCut")'
                              )

Finalphoton=selectedPatPhotons.clone(
 src ='mytightphoton',
cut =photonIDCutTight
)

#----------------------------------
#from PhysicsTools.PatAlgos.selectionLayer1.metSelector_cfi import *
#mymet = selectedPatMET.clone(
# src = 'goodJets',
#  src = 'patMETsPF',

# cut =
#    'et > 30. '
#)
from myanalysis.Atq.metcollection_cfi import *
mymet =selectedmets.clone()
mymet.JetCollection='YJET' 

#----------------------------------------
from myanalysis.Atq.bjetcollection_cfi import *
from myanalysis.Atq.Fjetcollection_cfi import *
YJET=selectedbJets.clone()
bJets=FselectedbJets.clone()
bJets.JetCollection='YJET'

#---------------------------------------------
goodJets=FselectedbJets.clone()
goodJets.JetCollection='YJET'
goodJets.bdis=0

#---------------------------------------
from myanalysis.Atq.puWeightProducer_cfi import *
pileupweight=weightProducer.clone()

pileupweightup=weightProducer.clone()
#pileupweightup.FileNamePUDataDistribution='myanalysis/Atq/plugins/MyDataPileupHistogramup.root'
pileupweightup.FileNamePUDataDistribution='myanalysis/Atq/plugins/MyJAN2013DataPileupHistogramup.root'

pileupweightdown=weightProducer.clone()
#pileupweightdown.FileNamePUDataDistribution='myanalysis/Atq/plugins/MyDataPileupHistogramdown.root'
pileupweightdown.FileNamePUDataDistribution='myanalysis/Atq/plugins/MyJAN2013DataPileupHistogramdown.root'

#--------------------------------------

from myanalysis.Atq.BTagSFEventWeight_cfi import *
BSF=bTagSFEventWeight.clone()
BSFup=bTagSFEventWeight.clone()
BSFup.sysVar='bTagSFUp'
BSFdown=bTagSFEventWeight.clone()
BSFdown.sysVar='bTagSFDown'
mistagup=bTagSFEventWeight.clone()
mistagup.sysVar='misTagSFUp'
mistagdown=bTagSFEventWeight.clone()
mistagdown.sysVar='misTagSFDown'


#--------------------------------------
from myanalysis.Atq.elecSelector_cfi import *
myelectron=electronSelector.clone()

#-----------------------------------------------------
#finding rho factor for each event
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
mykt6PFJets = kt4PFJets.clone( rParam = 0.6 )
mykt6PFJets.doRhoFastjet = True
mykt6PFJets.doAreaFastjet = True
mykt6PFJets.voronoiRfact = 0.9

from myanalysis.Atq.jesChange_cfi import *
topObjectSelection = cms.Sequence(
#    newra2PFchsJets  * 
    mykt6PFJets  *
    YJET *
    goodJets   *
    bJets*
    tightphoton *
    mytightphoton *
    Finalphoton  * 
    looseMuons *
    tightMuons *
    vetoMuons  *
    mymet  *
    myelectron  *
    pileupweight *
    pileupweightup *
    pileupweightdown *
    BSF *
    BSFdown *
    BSFup *
    mistagup *
    mistagdown )

