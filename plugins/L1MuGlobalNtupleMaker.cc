//ROOT includes
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkGlbMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkGlbMuonParticleFwd.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
  
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef vector<TTTrack<edm::Ref<edm::DetSetVector<Phase2TrackerDigi>,Phase2TrackerDigi,edm::refhelper::FindForDetSetVector<Phase2TrackerDigi> > > > TTTracksCollection;

#include "L1MuGlobalNtupleMaker.h"
 
// constructors and destructor
L1MuGlobalNtupleMaker::L1MuGlobalNtupleMaker(const edm::ParameterSet& iConfig) :
  _RunningOnData(iConfig.getParameter<bool>("RunningOnData")),
  _PU_scenario(iConfig.getParameter<int>("PileUpScenario")),
  _PileupSrc(iConfig.getParameter<edm::InputTag>("PileupSrc")),
  _maxGenMuons(iConfig.getParameter<int>("maxGenMuons")),
  _maxL1Muons(iConfig.getParameter<int>("maxL1Muons")),
  _maxBMTFMuons(iConfig.getParameter<int>("maxBMTFMuons")),
  _maxOMTFMuons(iConfig.getParameter<int>("maxOMTFMuons")),
  _maxEMTFMuons(iConfig.getParameter<int>("maxEMTFMuons")),
  _maxKBMTFMuons(iConfig.getParameter<int>("maxKBMTFMuons")),
  _maxDTPrimitives(iConfig.getParameter<int>("maxDTPrimitives")),
  _maxTkMuons(iConfig.getParameter<int>("maxTkMuons")),
  _maxTkGlbMuons(iConfig.getParameter<int>("maxTkGlbMuons")),
  _maxTTTracks(iConfig.getParameter<int>("maxTTTracks")),
  _maxTrkG4Parts(iConfig.getParameter<int>("maxTrkG4Parts")),
  _genParticleToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticle"))),
  _L1MuonToken(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("L1muon"))),
  _bmtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("bmtfMuon"))),
  _omtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("omtfMuon"))),
  _emtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("emtfMuon"))),
  _KbmtfMuonToken(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("KbmtfMuon"))),
  _bmtfPhInputToken(consumes<L1MuDTChambPhContainer>(iConfig.getParameter<edm::InputTag>("bmtfInputPhMuon"))),
  _bmtfThInputToken(consumes<L1MuDTChambThContainer>(iConfig.getParameter<edm::InputTag>("bmtfInputThMuon"))),
  _TkMuonToken(consumes<l1t::L1TkMuonParticleCollection>(iConfig.getParameter<edm::InputTag>("tkMuon"))),
  _TkGlbMuonToken(consumes<l1t::L1TkGlbMuonParticleCollection>(iConfig.getParameter<edm::InputTag>("tkGlbMuon"))),
  _TTTracksToken(consumes<TTTracksCollection>(iConfig.getParameter<edm::InputTag>("tttracks"))),
  _TrkG4PartsToken(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trkG4Parts")))
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

  //Generated muons
  _mytree->Branch("genmu_pt",&_genmu_pt);
  _mytree->Branch("genmu_eta",&_genmu_eta);
  _mytree->Branch("genmu_phi",&_genmu_phi);

  _mytree->Branch("genmu_Nmuons",&_genmu_Nmuons);

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
  _mytree->Branch("Kbmtfmu_hwpt",&_Kbmtfmu_hwpt);
  _mytree->Branch("Kbmtfmu_hweta",&_Kbmtfmu_hweta);
  _mytree->Branch("Kbmtfmu_hwphi",&_Kbmtfmu_hwphi);
  _mytree->Branch("Kbmtfmu_glbphi",&_Kbmtfmu_glbphi);
  _mytree->Branch("Kbmtfmu_hwsign",&_Kbmtfmu_hwsign);
  _mytree->Branch("Kbmtfmu_hwqual",&_Kbmtfmu_hwqual);
  _mytree->Branch("Kbmtfmu_link",&_Kbmtfmu_link);
  _mytree->Branch("Kbmtfmu_processor",&_Kbmtfmu_processor);
  _mytree->Branch("Kbmtfmu_bx",&_Kbmtfmu_bx);
  _mytree->Branch("Kbmtfmu_wheel",&_Kbmtfmu_wheel);

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

  //Track muons
  _mytree->Branch("tkmu_pt",&_tkmu_pt);
  _mytree->Branch("tkmu_eta",&_tkmu_eta);
  _mytree->Branch("tkmu_phi",&_tkmu_phi);
  _mytree->Branch("tkmu_charge",&_tkmu_charge);
  _mytree->Branch("tkmu_tkiso",&_tkmu_tkiso);

  _mytree->Branch("tkmu_Nmuons",&_tkmu_Nmuons);

  //Track global muons
  _mytree->Branch("tkglbmu_pt",&_tkglbmu_pt);
  _mytree->Branch("tkglbmu_eta",&_tkglbmu_eta);
  _mytree->Branch("tkglbmu_phi",&_tkglbmu_phi);
  _mytree->Branch("tkglbmu_charge",&_tkglbmu_charge);
  _mytree->Branch("tkglbmu_tkiso",&_tkglbmu_tkiso);

  _mytree->Branch("tkglbmu_Nmuons",&_tkglbmu_Nmuons);

  //TTTracks
  _mytree->Branch("tttracks_pt",&_tttracks_pt);
  _mytree->Branch("tttracks_eta",&_tttracks_eta);
  _mytree->Branch("tttracks_phi",&_tttracks_phi);
  _mytree->Branch("tttracks_chi2",&_tttracks_chi2);

  _mytree->Branch("tttracks_Nmuons",&_tttracks_Nmuons);

  //TrkG4Parts
  _mytree->Branch("trkG4Parts_pt",&_trkG4Parts_pt);
  _mytree->Branch("trkG4Parts_eta",&_trkG4Parts_eta);
  _mytree->Branch("trkG4Parts_phi",&_trkG4Parts_phi);
  _mytree->Branch("trkG4Parts_pdgId",&_trkG4Parts_pdgId);

  _mytree->Branch("trkG4Parts_Nmuons",&_trkG4Parts_Nmuons);

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

  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<l1t::MuonBxCollection> L1muon; 
  edm::Handle<l1t::RegionalMuonCandBxCollection> bmtfMuon;
  edm::Handle<l1t::RegionalMuonCandBxCollection> omtfMuon;
  edm::Handle<l1t::RegionalMuonCandBxCollection> emtfMuon;
  edm::Handle<l1t::RegionalMuonCandBxCollection> KbmtfMuon;
  edm::Handle<L1MuDTChambPhContainer > L1MuDTChambPhContainer;
  edm::Handle<L1MuDTChambThContainer > L1MuDTChambThContainer;
  edm::Handle<l1t::L1TkMuonParticleCollection> TkMuon;
  edm::Handle<l1t::L1TkGlbMuonParticleCollection> TkGlbMuon;
  edm::Handle<TTTracksCollection> Tttrack;
  edm::Handle<TrackingParticleCollection> TrkG4Part;

  iEvent.getByToken(_genParticleToken, genParticles);
  iEvent.getByToken(_L1MuonToken,      L1muon);
  iEvent.getByToken(_bmtfMuonToken,    bmtfMuon);
  iEvent.getByToken(_omtfMuonToken,    omtfMuon);
  iEvent.getByToken(_emtfMuonToken,    emtfMuon);
  iEvent.getByToken(_KbmtfMuonToken,   KbmtfMuon);
  iEvent.getByToken(_bmtfPhInputToken, L1MuDTChambPhContainer);
  iEvent.getByToken(_bmtfThInputToken, L1MuDTChambThContainer);
  iEvent.getByToken(_TkMuonToken,      TkMuon);
  iEvent.getByToken(_TkGlbMuonToken,   TkGlbMuon);
  iEvent.getByToken(_TTTracksToken,    Tttrack);
  iEvent.getByToken(_TrkG4PartsToken,  TrkG4Part);

  if(genParticles.isValid()){
    SetGenMuons(genParticles, _maxGenMuons);
  } else {
    edm::LogWarning("MissingProduct") << "Generated Muons not found. Branch will not be filled" << std::endl;
  }
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
  if(TkMuon.isValid()){
    SetTkMuons(TkMuon, _maxTkMuons);
  } else {
    edm::LogWarning("MissingProduct") << "L1 Phase2 TkMuons not found. Branch will not be filled" << std::endl;
  }
  if(TkGlbMuon.isValid()){
    SetTkGlbMuons(TkGlbMuon, _maxTkGlbMuons);
  } else {
    edm::LogWarning("MissingProduct") << "L1 Phase2 TkGlbMuons not found. Branch will not be filled" << std::endl;
  }
  if(Tttrack.isValid()){
    SetTTTracks(Tttrack, _maxTTTracks);
  } else{
    edm::LogWarning("MissingProduct") << "L1 Phase2 Tttrack not found. Branch will not be filled" << std::endl;
  }
  if(TrkG4Part.isValid()){
    SetTrkG4Parts(TrkG4Part, _maxTrkG4Parts);
  } else{
    edm::LogWarning("MissingProduct") << "L1 Phase2 TrkG4Part not found. Branch will not be filled" << std::endl;
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
  _NTruePU = -1;
  if(!_RunningOnData){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByLabel(_PileupSrc, PupInfo);
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	_NTruePU = PVI->getTrueNumInteractions();
	break;
      }
    }
  }

  _mytree->Fill();

}

void L1MuGlobalNtupleMaker::SetGenMuons(const edm::Handle<reco::GenParticleCollection> genParticles, int maxGenMu)
{
  _genmu_pt.clear();
  _genmu_eta.clear();
  _genmu_phi.clear();

  _genmu_Nmuons = 0;

  for(size_t i = 0; i < genParticles->size() && _genmu_Nmuons < maxGenMu ; ++ i) {
    const reco::GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    if(fabs(id) != 13) continue;

    _genmu_pt.push_back(p.pt());
    _genmu_eta.push_back(p.eta());
    _genmu_phi.push_back(p.phi());

    _genmu_Nmuons++;
  }
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

void L1MuGlobalNtupleMaker::SetKBMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxKBMTFUpgrade)
{
  _Kbmtfmu_hwpt.clear();
  _Kbmtfmu_hweta.clear();
  _Kbmtfmu_hwphi.clear();
  _Kbmtfmu_glbphi.clear();
  _Kbmtfmu_hwsign.clear();
  _Kbmtfmu_hwqual.clear();
  _Kbmtfmu_link.clear();
  _Kbmtfmu_processor.clear();
  _Kbmtfmu_bx.clear();
  _Kbmtfmu_wheel.clear();

  _Kbmtfmu_Nmuons = 0;

  for (int ibx = muon->getFirstBX(); ibx <= muon->getLastBX(); ++ibx) {
    for (l1t::RegionalMuonCandBxCollection::const_iterator it = muon->begin(ibx); it != muon->end(ibx) && _Kbmtfmu_Nmuons < maxKBMTFUpgrade; ++it){

      if (it->hwPt() > 0) {
	_Kbmtfmu_hwpt.push_back(it->hwPt());
	_Kbmtfmu_hweta.push_back(it->hwEta());
	_Kbmtfmu_hwphi.push_back(it->hwPhi());
	_Kbmtfmu_glbphi.push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(),it->trackFinderType(),it->processor()));
	_Kbmtfmu_hwsign.push_back(it->hwSign());
	_Kbmtfmu_hwqual.push_back(it->hwQual());
	_Kbmtfmu_link.push_back(it->link());
	_Kbmtfmu_processor.push_back(it->processor());
	_Kbmtfmu_bx.push_back(ibx);

        std::map<int, int>  trAdd;
        trAdd = it->trackAddress();
        int wheel = pow(-1,trAdd[0]) * trAdd[1];
	_Kbmtfmu_wheel.push_back(wheel);

	_Kbmtfmu_Nmuons++;

      }
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

  int iphtr = 0;

  for( L1MuDTChambPhContainer::Phi_Container::const_iterator DTPhDigiItr =  PhContainer->begin(); DTPhDigiItr != PhContainer->end(); ++DTPhDigiItr) {
    if( iphtr > maxDTPrimitives-1) continue;
    iphtr += 1;
    _Inputbmtf_phiBX.push_back(DTPhDigiItr->bxNum());
    _Inputbmtf_phiWheel.push_back(DTPhDigiItr->whNum());
    _Inputbmtf_phiSector.push_back(DTPhDigiItr->scNum());
    _Inputbmtf_phiStation.push_back(DTPhDigiItr->stNum());
    _Inputbmtf_phiAngle.push_back(DTPhDigiItr->phi());
    _Inputbmtf_phiBandAngle.push_back(DTPhDigiItr->phiB());
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

  int ithtr = 0;

  for( L1MuDTChambThContainer::The_Container::const_iterator DTThDigiItr =  ThContainer->begin(); DTThDigiItr != ThContainer->end();++DTThDigiItr){
    if(ithtr > maxDTPrimitives-1) continue;
    ithtr += 1;
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
    _Inputbmtf_thetaAngle.push_back(stoi(ss1.str()));

  }
}

void L1MuGlobalNtupleMaker::SetTkMuons(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, int maxTkMuons)
{
  _tkmu_pt.clear();
  _tkmu_eta.clear();
  _tkmu_phi.clear();
  _tkmu_charge.clear();
  _tkmu_tkiso.clear();

  _tkmu_Nmuons = 0;

  for (l1t::L1TkMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && _tkmu_Nmuons < maxTkMuons; it++){

    if (it->pt() > 0){
      _tkmu_pt.push_back(it->pt());
      _tkmu_eta.push_back(it->eta());
      _tkmu_phi.push_back(it->phi());
      _tkmu_charge.push_back(it->charge());
      _tkmu_tkiso.push_back(it->getTrkIsol());

      _tkmu_Nmuons++;

    }
  }

}

void L1MuGlobalNtupleMaker::SetTkGlbMuons(const edm::Handle<l1t::L1TkGlbMuonParticleCollection> muon, int maxTkGlbMuons)
{
  _tkglbmu_pt.clear();
  _tkglbmu_eta.clear();
  _tkglbmu_phi.clear();
  _tkglbmu_charge.clear();
  _tkglbmu_tkiso.clear();

  _tkglbmu_Nmuons = 0;

  for (l1t::L1TkGlbMuonParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && _tkglbmu_Nmuons < maxTkGlbMuons; it++){

    if (it->pt() > 0){
      _tkglbmu_pt.push_back(it->pt());
      _tkglbmu_eta.push_back(it->eta());
      _tkglbmu_phi.push_back(it->phi());
      _tkglbmu_charge.push_back(it->charge());
      _tkglbmu_tkiso.push_back(it->getTrkIsol());

      _tkglbmu_Nmuons++;

    }
  }

}

void L1MuGlobalNtupleMaker::SetTTTracks(const edm::Handle<TTTracksCollection> muon, int maxTTTracks)
{
  _tttracks_pt.clear();
  _tttracks_eta.clear();
  _tttracks_phi.clear();
  _tttracks_chi2.clear();

  _tttracks_Nmuons = 0;

  for (TTTracksCollection::const_iterator it=muon->begin(); it!=muon->end() && _tttracks_Nmuons < maxTTTracks; it++){

    if (it->getMomentum().perp() > 0){
      _tttracks_pt.push_back(it->getMomentum().perp()); //particle pT
      _tttracks_eta.push_back(it->getMomentum().eta());
      _tttracks_phi.push_back(it->getMomentum().phi());
      _tttracks_chi2.push_back(it->getChi2());

      _tttracks_Nmuons++;

    }
  }

}

void L1MuGlobalNtupleMaker::SetTrkG4Parts(const edm::Handle<TrackingParticleCollection> muon, int maxTrkG4Parts)
{
  _trkG4Parts_pt.clear();
  _trkG4Parts_eta.clear();
  _trkG4Parts_phi.clear();
  _trkG4Parts_pdgId.clear();

  _trkG4Parts_Nmuons = 0;

  for (TrackingParticleCollection::const_iterator it=muon->begin(); it!=muon->end() && _trkG4Parts_Nmuons < maxTrkG4Parts; it++){

    if (it->pt() > 0){
      _trkG4Parts_pt.push_back(it->pt());
      _trkG4Parts_eta.push_back(it->eta());
      _trkG4Parts_phi.push_back(it->phi());
      _trkG4Parts_pdgId.push_back(it->pdgId());

      _trkG4Parts_Nmuons++;

    }
  }


}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MuGlobalNtupleMaker);
