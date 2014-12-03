import FWCore.ParameterSet.Config as cms
myphoton = cms.EDProducer("photoncollection",
photonCollection = cms.InputTag("patPhotonsRA2"),
    floatLabels    = cms.VInputTag(cms.InputTag("mykt6PFJets","rho")),
    floatNames     = cms.vstring("phorho"),
    realdata   = cms.bool(False),
    photonESscaleType            = cms.string("abs"), #abs or rel(*eta) or jes:up / jes:down (pt-dependend)

    useAlternateIsolations = cms.bool(True),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring('reza', 'NEWshowerShape'),
        userFunctions = cms.vstring('1364',
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.031""",#showerShape

)


))
myselectedphoton = cms.EDProducer("photoncollection",
photonCollection = cms.InputTag("tightphoton"),
    realdata   = cms.bool(False),
    photonESscaleType            = cms.string("abs"), 
    useAlternateIsolations = cms.bool(False),
    floatLabels    = cms.VInputTag(),
    floatNames     = cms.vstring(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(
                                         'pfChargedPURel','pfNeutralPURel','pfGammaPURel',
                                         'pfChargedPU','pfNeutralPU','pfGammaPU'
                                         ),
        userFunctions = cms.vstring(
            'max((userFloat("pfChaIsoAlt") - userFloat("pfChargedEANew")*userFloat("phorho"))/et,0.)',
            'max((userFloat("pfNeuIsoAlt") - userFloat("pfNeutralEANew")*userFloat("phrho"))/et,0.)',
            'max((userFloat("pfGamIsoAlt")   - userFloat("pfGammaEANew")  *userFloat("phorho"))/et,0.)',
            'max((userFloat("pfChaIsoAlt") - userFloat("pfChargedEANew")*userFloat("phorho")),0.)',
            """max((userFloat("pfNeuIsoAlt") - userFloat("pfNeutralEANew")*userFloat("phorho")),0.)""",
            'max((userFloat("pfGamIsoAlt")   - userFloat("pfGammaEANew")  *userFloat("phorho")),0.)'
        )
    )
)
                              
