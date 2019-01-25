
import FWCore.ParameterSet.Config as cms

def customL1Emu(process, options):

    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
    process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
    process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')

    process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')


    if options.redoPrimitives :

        process.redoPrimitives_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis, process.hgcalTriggerPrimitives)

    if options.doTTrigger :

        process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')
        process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)
        process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

