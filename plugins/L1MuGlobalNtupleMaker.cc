//ROOT includes
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
  
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

typedef math::XYZTLorentzVector LorentzVector;
 
#include "L1MuGlobalNtupleMaker.h"
 
// constructors and destructor
L1MuGlobalNtupleMaker::L1MuGlobalNtupleMaker(const edm::ParameterSet& iConfig) :
  _runningOnData(iConfig.getParameter<bool>("RunningOnData")),
  _PU_scenario(iConfig.getParameter<int>("PileUpScenario")),
  _PileupSrc(iConfig.getParameter<edm::InputTag>("PileupSrc")),
  _triggerBits(consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerbits")))
{
  _pileupSummaryToken = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));

  _h_Global_Info = fs->make<TH1F>("h_Global_Info", "General information about the sample", 8, 0., 8.);
  _Nevents_processed = 0;

  create_trees();
}

L1MuGlobalNtupleMaker::~L1MuGlobalNtupleMaker()
{
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//
void L1MuGlobalNtupleMaker::create_trees()
{
  _mytree = fs->make<TTree>("mytree", "Tree containing L1 info");

  //Save run number info when running on data
  if(_runningOnData){
    _mytree->Branch("run_number",&_run_number);
  }
  else{
    _mytree->Branch("NTruePU",&_NTruePU);
  }


}

void L1MuGlobalNtupleMaker::endJob() 
{
  _h_Global_Info->Fill(0.5,_Nevents_processed);
  _h_Global_Info->Fill(1.5,_PU_scenario);

  _h_Global_Info->GetXaxis()->SetBinLabel(1,"Events processed");
  _h_Global_Info->GetXaxis()->SetBinLabel(2,"PU scenario in the MC sample");
}

// ------------ method called for each event  ------------
void L1MuGlobalNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  _Nevents_processed++;

  //Retrieve the run number
  if(_RunningOnData_){
    _run_number = iEvent.id().run();
  }

  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//
  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByLabel(_PileupSrc, PupInfo);
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	_NTruePU = PVI->getTrueNumInteractions();
      }
    }
  }


}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MuGlobalNtupleMaker);
