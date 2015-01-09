import FWCore.ParameterSet.Config as cms

tauGenJets = cms.EDProducer(
    "TauGenJetProducer",
    GenParticles =  cms.InputTag('genParticles'),
    PrunedGenParticles =  cms.InputTag('prunedGenParticles'), #mf
    isMiniAOD = cms.bool( True ), #mf
    includeNeutrinos = cms.bool( False ),
    verbose = cms.untracked.bool( False )
    )
