import FWCore.ParameterSet.Config as cms
 
bTagSFEventWeight = cms.EDProducer("BTagSFEventWeight",
	jets  = cms.InputTag("patJetsAK5PFPt30"), ## jet collection (after jet selection, before b-tagging)
        sysVar   = cms.string(""),                  ## bTagSFUp, bTagSFDown, misTagSFUp, misTagSFDown possible;
	filename = cms.FileInPath("myanalysis/Atq/plugins/MC_B_C_lighte_jet_eff.root")               ## if filename != "", the efficiencies are read from histos
 )
