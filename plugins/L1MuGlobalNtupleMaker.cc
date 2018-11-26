//ROOT includes
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
  
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
  _RunningOnData(iConfig.getParameter<bool>("RunningOnData")),
  _PU_scenario(iConfig.getParameter<int>("PileUpScenario")),
  _PileupSrc(iConfig.getParameter<edm::InputTag>("PileupSrc")),
  _maxL1Muons(iConfig.getParameter<int>("maxL1Muons")),
  _maxBMTFMuons(iConfig.getParameter<int>("maxBMTFMuons")),
  _maxOMTFMuons(iConfig.getParameter<int>("maxOMTFMuons")),
  _maxEMTFMuons(iConfig.getParameter<int>("maxEMTFMuons")),
  _maxDTPrimitives(iConfig.getParameter<int>("maxDTPrimitives")),
  _L1MuonToken(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("L1muon"))),
  _bmtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("bmtfMuon"))),
  _omtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("omtfMuon"))),
  _emtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("emtfMuon"))),
  _KbmtfMuonToken(consumes<L1MuKBMTrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("KbmtfMuon"))),
  _bmtfPhInputToken = consumes<L1MuDTChambPhContainer>(iConfig.getUntrackedParameter<edm::InputTag>("bmtfInputPhMuon")),
  _bmtfThInputToken = consumes<L1MuDTChambThContainer>(iConfig.getUntrackedParameter<edm::InputTag>("bmtfInputThMuon"))
{
  _pileupSummaryToken = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(_PileupSrc));

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
  if(_RunningOnData){
    _mytree->Branch("run_number",&_run_number);
  }
  else{
    _mytree->Branch("NTruePU",&_NTruePU);
  }

  //L1 muons
  _mytree->Branch("l1mu_et",&_l1mu_et);
  _mytree->Branch("l1mu_eta",&_l1mu_eta);
  _mytree->Branch("l1mu_phi",&_l1mu_phi);
  _mytree->Branch("l1mu_etaAtVtx",&_l1mu_etaAtVtx);
  _mytree->Branch("l1mu_phiAtVtx",&_l1mu_phiAtVtx);
  _mytree->Branch("l1mu_hwet",&_l1mu_hwet);
  _mytree->Branch("l1mu_hweta",&_l1mu_hweta);
  _mytree->Branch("l1mu_hwphi",&_l1mu_hwphi);

  _mytree->Branch("l1mu_charge",&_l1mu_charge);
  _mytree->Branch("l1mu_hwIso",&_l1mu_hwIso);
  _mytree->Branch("l1mu_hwQual",&_l1mu_hwQual);
  _mytree->Branch("l1mu_bx",&_l1mu_bx);

  _mytree->Branch("l1mu_Nmuons",&_l1mu_Nmuons);

  //BMTF muons
  _mytree->Branch("bmtfmu_hwpt",&_bmtfmu_hwpt);
  _mytree->Branch("bmtfmu_hweta",&_bmtfmu_hweta);
  _mytree->Branch("bmtfmu_hwphi",&_bmtfmu_hwphi);
  _mytree->Branch("bmtfmu_glbphi",&_bmtfmu_glbphi);
  _mytree->Branch("bmtfmu_hwsign",&_bmtfmu_hwsign);
  _mytree->Branch("bmtfmu_hwqual",&_bmtfmu_hwqual);
  _mytree->Branch("bmtfmu_link",&_bmtfmu_link);
  _mytree->Branch("bmtfmu_processor",&_bmtfmu_processor);
  _mytree->Branch("bmtfmu_bx",&_bmtfmu_bx);
  _mytree->Branch("bmtfmu_wheel",&_bmtfmu_wheel);

  _mytree->Branch("bmtfmu_Nmuons",&_bmtfmu_Nmuons);

  //OMTF muons
  _mytree->Branch("omtfmu_hwpt",&_omtfmu_hwpt);
  _mytree->Branch("omtfmu_hweta",&_omtfmu_hweta);
  _mytree->Branch("omtfmu_hwphi",&_omtfmu_hwphi);
  _mytree->Branch("omtfmu_glbphi",&_omtfmu_glbphi);
  _mytree->Branch("omtfmu_hwsign",&_omtfmu_hwsign);
  _mytree->Branch("omtfmu_hwqual",&_omtfmu_hwqual);
  _mytree->Branch("omtfmu_link",&_omtfmu_link);
  _mytree->Branch("omtfmu_processor",&_omtfmu_processor);
  _mytree->Branch("omtfmu_bx",&_omtfmu_bx);
  _mytree->Branch("omtfmu_wheel",&_omtfmu_wheel);

  _mytree->Branch("omtfmu_Nmuons",&_omtfmu_Nmuons);

  //EMTF muons
  _mytree->Branch("emtfmu_hwpt",&_emtfmu_hwpt);
  _mytree->Branch("emtfmu_hweta",&_emtfmu_hweta);
  _mytree->Branch("emtfmu_hwphi",&_emtfmu_hwphi);
  _mytree->Branch("emtfmu_glbphi",&_emtfmu_glbphi);
  _mytree->Branch("emtfmu_hwsign",&_emtfmu_hwsign);
  _mytree->Branch("emtfmu_hwqual",&_emtfmu_hwqual);
  _mytree->Branch("emtfmu_link",&_emtfmu_link);
  _mytree->Branch("emtfmu_processor",&_emtfmu_processor);
  _mytree->Branch("emtfmu_bx",&_emtfmu_bx);
  _mytree->Branch("emtfmu_wheel",&_emtfmu_wheel);

  _mytree->Branch("emtfmu_Nmuons",&_emtfmu_Nmuons);

  //Kalman BMTF muons
  _mytree->Branch("Kbmtfmu_pt",&_Kbmtfmu_pt);
  _mytree->Branch("Kbmtfmu_eta",&_Kbmtfmu_eta);
  _mytree->Branch("Kbmtfmu_phi",&_Kbmtfmu_phi);
  _mytree->Branch("Kbmtfmu_dxy",&_Kbmtfmu_dxy);
  _mytree->Branch("Kbmtfmu_approxChi2",&_Kbmtfmu_approxChi2);
  _mytree->Branch("Kbmtfmu_sector",&_Kbmtfmu_sector);
  _mytree->Branch("Kbmtfmu_wheel",&_Kbmtfmu_wheel);
  _mytree->Branch("Kbmtfmu_qual",&_Kbmtfmu_qual);
  _mytree->Branch("Kbmtfmu_bx",&_Kbmtfmu_bx);

  _mytree->Branch("Kbmtfmu_Nmuons",&_Kbmtfmu_Nmuons);

  //BMTF inputs
  _mytree->Branch("Inputbmtf_phiBX",&_Inputbmtf_phiBX);
  _mytree->Branch("Inputbmtf_phiWheel",&_Inputbmtf_phiWheel);
  _mytree->Branch("Inputbmtf_phiSector",&_Inputbmtf_phiSector);
  _mytree->Branch("Inputbmtf_phiStation",&_Inputbmtf_phiStation);
  _mytree->Branch("Inputbmtf_phiAngle",&_Inputbmtf_phiAngle);
  _mytree->Branch("Inputbmtf_phiBandAngle",&_Inputbmtf_phiBandAngle);

  _mytree->Branch("Inputbmtf_phiNprimitives",&_Inputbmtf_phiNprimitives);

  _mytree->Branch("Inputbmtf_thetaBX",&_Inputbmtf_thetaBX);
  _mytree->Branch("Inputbmtf_thetaWheel",&_Inputbmtf_thetaWheel);
  _mytree->Branch("Inputbmtf_thetaSector",&_Inputbmtf_thetaSector);
  _mytree->Branch("Inputbmtf_thetaStation",&_Inputbmtf_thetaStation);
  _mytree->Branch("Inputbmtf_thetaAngle",&_Inputbmtf_thetaAngle);

  _mytree->Branch("Inputbmtf_thetaNprimitives",&_Inputbmtf_thetaNprimitives);
}

void L1MuGlobalNtupleMaker::beginJob()
{
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
  _Nevents_processed++;

  edm::Handle<l1t::MuonBxCollection> L1muon; 
  edm::Handle<l1t::RegionalMuonCandBxCollection> bmtfMuon;
  edm::Handle<l1t::RegionalMuonCandBxCollection> omtfMuon;
  edm::Handle<l1t::RegionalMuonCandBxCollection> emtfMuon;
  edm::Handle<L1MuKBMTrackCollection> KbmtfMuon;
  edm::Handle<L1MuDTChambPhContainer > L1MuDTChambPhContainer;
  edm::Handle<L1MuDTChambThContainer > L1MuDTChambThContainer;

  iEvent.getByToken(_L1MuonToken,    L1muon);
  iEvent.getByToken(_bmtfMuonToken,  bmtfMuon);
  iEvent.getByToken(_omtfMuonToken,  omtfMuon);
  iEvent.getByToken(_emtfMuonToken,  emtfMuon);
  iEvent.getByToken(_KbmtfMuonToken, KbmtfMuon);
  iEvent.getByToken(_bmtfPhInputToken, L1MuDTChambPhContainer);
  iEvent.getByToken(_bmtfThInputToken, L1MuDTChambThContainer);

  if(L1muon.isValid()){
    SetL1Muons(L1muon, _maxL1Muons);
  } else {
    edm::LogWarning("MissingProduct") << "L1Upgrade Muons not found. Branch will not be filled" << std::endl;
  }
  if(bmtfMuon.isValid()){
    SetBMTFMuons(bmtfMuon, _maxBMTFMuons);
  } else {
    edm::LogWarning("MissingProduct") << "BMTF Muons not found. Branch will not be filled" << std::endl;
  }
  if(omtfMuon.isValid()){
    SetOMTFMuons(omtfMuon, _maxOMTFMuons);
  } else {
    edm::LogWarning("MissingProduct") << "OMTF Muons not found. Branch will not be filled" << std::endl;
  }
  if(emtfMuon.isValid()){
    SetEMTFMuons(emtfMuon, _maxEMTFMuons);
  } else {
    edm::LogWarning("MissingProduct") << "EMTF Muons not found. Branch will not be filled" << std::endl;
  }
  if(KbmtfMuon.isValid()){
    SetKBMTFMuons(KbmtfMuon, _maxKBMTFMuons);
  } else {
    edm::LogWarning("MissingProduct") << "Kalman BMTF Muons not found. Branch will not be filled" << std::endl;
  }
  if(L1MuDTChambPhContainer.isValid()){
    SetBMTFPhiInputs(L1MuDTChambPhContainer, _maxDTPrimitives);
  } else {
    edm::LogWarning("MissingProduct") << "BMTF Phi inputs not found. Branch will not be filled" << std::endl;
  }
  if(L1MuDTChambThContainer.isValid()){
    SetBMTFThetaInputs(L1MuDTChambThContainer, _maxDTPrimitives);
  } else {
    edm::LogWarning("MissingProduct") << "BMTF Theta inputs not found. Branch will not be filled" << std::endl;
  }

  //Retrieve the run number
  if(_RunningOnData){
    _run_number = iEvent.id().run();
  }

  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//
  if(!_RunningOnData){
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

  _mytree->Fill();

}

void L1MuGlobalNtupleMaker::SetL1Muons(const edm::Handle<l1t::MuonBxCollection> muon, int maxL1Upgrade)
{
  _l1mu_et.clear();
  _l1mu_eta.clear();
  _l1mu_phi.clear();
  _l1mu_etaAtVtx.clear();
  _l1mu_phiAtVtx.clear();
  _l1mu_hwet.clear();
  _l1mu_hweta.clear();
  _l1mu_hwphi.clear();

  _l1mu_charge.clear();
  _l1mu_hwIso.clear();
  _l1mu_hwQual.clear();
  _l1mu_bx.clear();

  _l1mu_Nmuons = 0;

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::MuonBxCollection::const_iterator it=muon->begin(ibx); it!=muon->end(ibx) && _l1mu_Nmuons < maxL1Upgrade; it++){
      if (it->pt() > 0){
	_l1mu_et.push_back(it->et());
	_l1mu_eta.push_back(it->eta());
	_l1mu_phi.push_back(it->phi());
	_l1mu_etaAtVtx.push_back(it->etaAtVtx());
	_l1mu_phiAtVtx.push_back(it->phiAtVtx());
	_l1mu_hwet.push_back(it->hwPt());
	_l1mu_hweta.push_back(it->hwEta());
	_l1mu_hwphi.push_back(it->hwPhi());

	_l1mu_charge.push_back(it->charge());
	_l1mu_hwIso.push_back(it->hwIso());
	_l1mu_hwQual.push_back(it->hwQual());

	_l1mu_bx.push_back(ibx);

        _l1mu_Nmuons++;

      }
    }
  }

}

void L1MuGlobalNtupleMaker::SetBMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxBMTFUpgrade)
{
  _bmtfmu_hwpt.clear();
  _bmtfmu_hweta.clear();
  _bmtfmu_hwphi.clear();
  _bmtfmu_glbphi.clear();
  _bmtfmu_hwsign.clear();
  _bmtfmu_hwqual.clear();
  _bmtfmu_link.clear();
  _bmtfmu_processor.clear();
  _bmtfmu_bx.clear();
  _bmtfmu_wheel.clear();

  _bmtfmu_Nmuons = 0;

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::RegionalMuonCandBxCollection::const_iterator it = muon->begin(ibx); it != muon->end(ibx) && _bmtfmu_Nmuons < maxBMTFUpgrade; ++it){
      if (it->hwPt() > 0) {

	_bmtfmu_hwpt.push_back(it->hwPt());
	_bmtfmu_hweta.push_back(it->hwEta());
	_bmtfmu_hwphi.push_back(it->hwPhi());
	_bmtfmu_glbphi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(),it->trackFinderType(),it->processor()));
	_bmtfmu_hwsign.push_back(it->hwSign());
	_bmtfmu_hwqual.push_back(it->hwQual());
	_bmtfmu_link.push_back(it->link());
	_bmtfmu_processor.push_back(it->processor());
	_bmtfmu_bx.push_back(ibx);

	std::map<int, int>  trAdd;
        trAdd = it->trackAddress();
        int wheel = pow(-1,trAdd[0]) * trAdd[1];
	_bmtfmu_wheel.push_back(wheel);

	_bmtfmu_Nmuons++;

      }
    }
  }

}

void L1MuGlobalNtupleMaker::SetOMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxOMTFUpgrade)
{
  _omtfmu_hwpt.clear();
  _omtfmu_hweta.clear();
  _omtfmu_hwphi.clear();
  _omtfmu_glbphi.clear();
  _omtfmu_hwsign.clear();
  _omtfmu_hwqual.clear();
  _omtfmu_link.clear();
  _omtfmu_processor.clear();
  _omtfmu_bx.clear();
  _omtfmu_wheel.clear();

  _omtfmu_Nmuons = 0;

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::RegionalMuonCandBxCollection::const_iterator it = muon->begin(ibx); it != muon->end(ibx) && _omtfmu_Nmuons < maxOMTFUpgrade; ++it){
      if (it->hwPt() > 0) {

	_omtfmu_hwpt.push_back(it->hwPt());
	_omtfmu_hweta.push_back(it->hwEta());
	_omtfmu_hwphi.push_back(it->hwPhi());
	_omtfmu_glbphi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(),it->trackFinderType(),it->processor()));
	_omtfmu_hwsign.push_back(it->hwSign());
	_omtfmu_hwqual.push_back(it->hwQual());
	_omtfmu_link.push_back(it->link());
	_omtfmu_processor.push_back(it->processor());
	_omtfmu_bx.push_back(ibx);

	std::map<int, int>  trAdd;
        trAdd = it->trackAddress();
        int wheel = pow(-1,trAdd[0]) * trAdd[1];
	_omtfmu_wheel.push_back(wheel);

	_omtfmu_Nmuons++;

      }
    }
  }

}

void L1MuGlobalNtupleMaker::SetEMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxEMTFUpgrade)
{
  _emtfmu_hwpt.clear();
  _emtfmu_hweta.clear();
  _emtfmu_hwphi.clear();
  _emtfmu_glbphi.clear();
  _emtfmu_hwsign.clear();
  _emtfmu_hwqual.clear();
  _emtfmu_link.clear();
  _emtfmu_processor.clear();
  _emtfmu_bx.clear();
  _emtfmu_wheel.clear();

  _emtfmu_Nmuons = 0;

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::RegionalMuonCandBxCollection::const_iterator it = muon->begin(ibx); it != muon->end(ibx) && _emtfmu_Nmuons < maxEMTFUpgrade; ++it){
      if (it->hwPt() > 0) {

	_emtfmu_hwpt.push_back(it->hwPt());
	_emtfmu_hweta.push_back(it->hwEta());
	_emtfmu_hwphi.push_back(it->hwPhi());
	_emtfmu_glbphi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(),it->trackFinderType(),it->processor()));
	_emtfmu_hwsign.push_back(it->hwSign());
	_emtfmu_hwqual.push_back(it->hwQual());
	_emtfmu_link.push_back(it->link());
	_emtfmu_processor.push_back(it->processor());
	_emtfmu_bx.push_back(ibx);

	std::map<int, int>  trAdd;
        trAdd = it->trackAddress();
        int wheel = pow(-1,trAdd[0]) * trAdd[1];
	_emtfmu_wheel.push_back(wheel);

	_emtfmu_Nmuons++;

      }
    }
  }
}

void L1MuGlobalNtupleMaker::SetKBMTFMuons(const edm::Handle<L1MuKBMTrackCollection> muon, int maxKBMTFUpgrade)
{
  _Kbmtfmu_pt.clear();
  _Kbmtfmu_eta.clear();
  _Kbmtfmu_phi.clear();
  _Kbmtfmu_dxy.clear();
  _Kbmtfmu_approxChi2.clear();
  _Kbmtfmu_sector.clear();
  _Kbmtfmu_wheel.clear();
  _Kbmtfmu_qual.clear();
  _Kbmtfmu_bx.clear();

  _Kbmtfmu_Nmuons = 0;

    for (L1MuKBMTrackCollection::const_iterator it = muon->begin(); it != muon->end() && _Kbmtfmu_Nmuons < maxKBMTFUpgrade; ++it){
      if (it->pt() > 0) {

	_Kbmtfmu_pt.push_back(it->pt());
	_Kbmtfmu_eta.push_back(it->eta());
	_Kbmtfmu_phi.push_back(it->phi());
	_Kbmtfmu_dxy.push_back(it->dxy());
	_Kbmtfmu_approxChi2.push_back(it->approxChi2());
	_Kbmtfmu_sector.push_back(it->sector());
	_Kbmtfmu_wheel.push_back(it->wheel());
	_Kbmtfmu_qual.push_back(it->quality());

	_Kbmtfmu_bx.push_back(it->bx());

	_Kbmtfmu_Nmuons++;

    }
  }
}

void L1MuGlobalNtupleMaker::SetBMTFPhiInputs(const edm::Handle<L1MuDTChambPhContainer > L1MuDTChambPhContainer, int maxDTPrimitives)
{
  L1MuDTChambPhContainer::Phi_Container const *PhContainer = L1MuDTChambPhContainer->getContainer();

  _Inputbmtf_phiBX.clear();
  _Inputbmtf_phiWheel.clear();
  _Inputbmtf_phiSector.clear();
  _Inputbmtf_phiStation.clear();
  _Inputbmtf_phiAngle.clear();
  _Inputbmtf_phiBandAngle.clear();

  _Inputbmtf_phiNprimitives = PhContainer->size();

  for( L1MuDTChambPhContainer::Phi_Container::const_iterator DTPhDigiItr =  PhContainer->begin(); DTPhDigiItr != PhContainer->end(); ++DTPhDigiItr) {
    if((unsigned int) iphtr>maxDTPrimitives-1) continue;
    _Inputbmtf_phiBX.push_back(DTPhDigiItr->bxNum());
    _Inputbmtf_phiWheel.push_back(DTPhDigiItr->whNum());
    _Inputbmtf_phiSector.push_back(DTPhDigiItr->scNum());
    _Inputbmtf_phiStation.push_back(DTPhDigiItr->stNum());
    _Inputbmtf_phiAng.push_back(DTPhDigiItr->phi());
    _Inputbmtf_phiBandAng.push_back(DTPhDigiItr->phiB());
  }
}

void L1MuGlobalNtupleMaker::SetBMTFThetaInputs(const edm::Handle<L1MuDTChambThContainer > L1MuDTChambThContainer, int maxDTPrimitives)
{
  L1MuDTChambThContainer::The_Container const *ThContainer = L1MuDTChambThContainer->getContainer();

  _Inputbmtf_thetaBX.clear();
  _Inputbmtf_thetaWheel.clear();
  _Inputbmtf_thetaSector.clear();
  _Inputbmtf_thetaStation.clear();
  _Inputbmtf_thetaAngle.clear();

  _Inputbmtf_thetaNprimitives = ThContainer->size();

  for( L1MuDTChambThContainer::The_Container::const_iterator DTThDigiItr =  ThContainer->begin(); DTThDigiItr != ThContainer->end();++DTThDigiItr){
    if((unsigned int) ithtr>maxDTPrimitives-1) continue;
    _Inputbmtf_thetaBX.push_back(DTThDigiItr->bxNum());
    _Inputbmtf_thetaWheel.push_back(DTThDigiItr->whNum());
    _Inputbmtf_thetaSector.push_back(DTThDigiItr->scNum());
    _Inputbmtf_thetaStation.push_back(DTThDigiItr->stNum());

    ostringstream  ss1, ss2; 
    ss1.clear(); ss2.clear();
    ss1<<"9"; ss2<<"9";

    for(int j=0; j<7; j++){
      ss1<<DTThDigiItr->position(j);
      ss2<<DTThDigiItr->code(j) ;
    }
    _Inputbmtf_thetaAng.push_back(stoi(ss1.str()));

  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MuGlobalNtupleMaker);
