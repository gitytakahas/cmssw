import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Source
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_4_4_2_patch10/RelValProdTTbar/AODSIM/START44_V7_special_120119-v1/0088/0ECC7FA5-7A43-E111-9B36-002618943861.root',
       '/store/relval/CMSSW_4_4_2_patch10/RelValProdTTbar/AODSIM/START44_V7_special_120119-v1/0084/1AFC5FFF-3143-E111-9E01-001BFCDBD1BC.root'
    )
)
## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")

## Test JEC from test instances of the global DB
#process.load("PhysicsTools.PatAlgos.patTestJEC_cfi")

## Test JEC from local sqlite file
#process.load("PhysicsTools.PatAlgos.patTestJEC_local_cfi")

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
                               )

process.outpath = cms.EndPath(process.out)
