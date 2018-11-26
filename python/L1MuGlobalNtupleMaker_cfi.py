import FWCore.ParameterSet.Config as cms

L1MuGlobalNtupleMaker = cms.EDAnalyzer('L1MuGlobalNtupleMaker',
                                  RunningOnData      = cms.bool(False),
                                  PileUpScenario     = cms.int(200),
                                  PileupSrc          = cms.InputTag("slimmedAddPileupInfo"),
                                  triggerbits        = cms.InputTag("TriggerResults","","HLT"),
                                 )
