import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
      '/store/relval/CMSSW_7_0_9_patch2/RelValTTbar_13/MINIAODSIM/PLS170_V6AN2-v1/00000/C0AE81D4-8C59-E411-B48B-0025905964B2.root',
       '/store/relval/CMSSW_7_0_9_patch2/RelValTTbar_13/MINIAODSIM/PLS170_V6AN2-v1/00000/CACC80D5-8C59-E411-B6CB-0025905A60CA.root' 
] )


secFiles.extend( [
   ] )
