from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/samples_MC_phase2/'

config.section_('JobType')
config.JobType.psetName = '../run_L1MuNtuple.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['L1MuPhase2Ntuple_output.root']
config.JobType.pyCfgParams = ['doPhase2Emul=True']

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_CH_CERN'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    
    config.General.requestName = 'L1MuPhase2Ntuples_SingleNu_PU200'
    config.Data.unitsPerJob = 3
    config.Data.inputDataset = '/SingleNeutrino/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'L1MuPhase2Ntuples_SingleMuPt2to100_PU140'
    config.Data.unitsPerJob = 3
    config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU140_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'L1MuPhase2Ntuples_SingleMuPt2to100_PU200'
    config.Data.unitsPerJob = 3
    config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
