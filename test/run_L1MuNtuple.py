
import FWCore.ParameterSet.Config as cms

# General config options
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '103X_upgrade2023_realistic_v2', #default value 
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('reEmulation',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run re-emulation")

options.register('doPhase2Emul',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run the phase 2 re-emulation")

options.register('redoPrimitives',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Redo the primitives")

options.register('runOnMC',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Set to True when running on MC")

options.parseArguments()


from Configuration.ProcessModifiers.convertHGCalDigisSim_cff import convertHGCalDigisSim
from Configuration.StandardSequences.Eras import eras


if options.doPhase2Emul :
    print "Using track trigger"
    process = cms.Process('L1',eras.Phase2_trigger,convertHGCalDigisSim)
else :
    process = cms.Process('L1',eras.Phase2_timing,convertHGCalDigisSim)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')


process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)
process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")



# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

#Import the emulation configuration
from L1Trigger.L1MuGlobalNtuples.customL1Emu_cff import *
customL1Emu(process, options)

print "Using GlobalTag", options.globalTag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# input source
process.source = cms.Source("PoolSource",
                             fileNames = cms.untracked.vstring ("root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRSpring19DR/Mu_FlatPt2to100-pythia8-gun/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v1/60000/E0D5C6A5-B855-D14F-9124-0B2C9B28D0EA.root"),
                             inputCommands = cms.untracked.vstring("keep *", 
                                                                   "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
                                                                   "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
                                                                   "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
                                                                   "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
                                                                   "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
                                                                   "drop l1tEMTFTrack2016s_simEmtfDigis__HLT"
                                                                   )
                             )
 


#Import the ntuplizer
process.load("L1Trigger.L1MuGlobalNtuples.L1MuGlobalNtupleMaker_cfi")

# we don't have emtfDigis yet, use unpacked input payloads of GMT
#process.L1MuGlobalNtupleMaker.emtfMuon = cms.InputTag("gmtStage2Digis","EMTF") 

process.ntuplizer = cms.Path(process.L1MuGlobalNtupleMaker)

# Output file
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("L1MuPhase2Ntuple_output.root")
                                   )


# print log
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



if options.reEmulation :
    if options.doPhase2Emul : 
        print "Using Phase-2 emulation"
        process.L1simulation_step = cms.Path(process.phase2_SimL1Emulator)
    else : 
        print "Using Phase-1 emulation"
        process.L1simulation_step = cms.Path(process.SimL1Emulator)


process.endjob_step = cms.EndPath(process.endOfProcess)


# # --- save emulated objects in a root file 
# process.outprova = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         #dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
#         #filterName = cms.untracked.string('')
#     ),
#     outputCommands = cms.untracked.vstring("keep *",
#                                            "drop *_hgcalTriggerPrimitiveDigiProducer_*_*",  
#                                            "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
#                                            "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
#                                            "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
#                                            "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
#                                            "drop l1tEMTFTrack2016s_simEmtfDigis__HLT"),
#     fileName = cms.untracked.string('prova.root'),
#     splitLevel = cms.untracked.int32(0)
# )
# process.outprova_step = cms.EndPath(process.outprova)



# Schedule definition
process.schedule = cms.Schedule()

if options.redoPrimitives :
    print "Regenerating trigger primitives"
    process.schedule.append(process.redoPrimitives_step)

if options.doPhase2Emul :
    process.schedule.append(process.L1TrackTrigger_step)

if options.reEmulation :
    process.schedule.append(process.L1simulation_step)

process.schedule.append(process.ntuplizer)
process.schedule.append(process.endjob_step)

# --- save emulated objects in a root file 
# process.schedule.append(process.outprova_step)
