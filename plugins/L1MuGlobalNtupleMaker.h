
//---------- class declaration----------

class L1MuGlobalNtupleMaker : public edm::EDAnalyzer {
public:
  explicit L1MuGlobalNtupleMaker(const edm::ParameterSet&);
  ~L1MuGlobalNtupleMaker();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  void SetGenMuons(const edm::Handle<reco::GenParticleCollection>  genParticles, int maxGenMu);
  void SetL1Muons(const edm::Handle<l1t::MuonBxCollection> muon, int maxL1Upgrade);
  void SetBMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxBMTFUpgrade);
  void SetOMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxOMTFUpgrade);
  void SetEMTFMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection> muon, int maxEMTFUpgrade);
  void SetKBMTFMuons(const edm::Handle<L1MuKBMTrackCollection> muon, int maxKBMTFUpgrade);
  void SetBMTFPhiInputs(const edm::Handle<L1MuDTChambPhContainer > L1MuDTChambPhContainer, int maxDTPrimitives);
  void SetBMTFThetaInputs(const edm::Handle<L1MuDTChambThContainer > L1MuDTChambThContainer, int maxDTPrimitives);
  void SetTkMuons(const edm::Handle<l1t::L1TkMuonParticleCollection> muon, int maxTkMuons);
  void SetTkGlbMuons(const edm::Handle<l1t::L1TkGlbMuonParticleCollection> muon, int maxTkGlbMuons);

  bool _RunningOnData;
  int _PU_scenario;
  const edm::InputTag _PileupSrc;

  int _maxGenMuons;
  int _maxL1Muons;
  int _maxBMTFMuons;
  int _maxOMTFMuons;
  int _maxEMTFMuons;
  int _maxKBMTFMuons;
  int _maxDTPrimitives;
  int _maxTkMuons;
  int _maxTkGlbMuons;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* _h_Global_Info;

  //Counters
  int _Nevents_processed;

  //#########################################
  //TTree and TTree variables
  //#########################################
  TTree *_mytree;

  int _run_number;

  //MC truth
  float _NTruePU;
  std::vector<float> _genmu_pt;
  std::vector<float> _genmu_eta;
  std::vector<float> _genmu_phi;
  short int _genmu_Nmuons;

  //L1 muon
  std::vector<float> _l1mu_et;
  std::vector<float> _l1mu_eta;
  std::vector<float> _l1mu_phi;
  std::vector<float> _l1mu_etaAtVtx;
  std::vector<float> _l1mu_phiAtVtx;
  std::vector<short int> _l1mu_hwet;
  std::vector<short int> _l1mu_hweta;
  std::vector<short int> _l1mu_hwphi;
  std::vector<short int> _l1mu_charge;
  std::vector<unsigned short int> _l1mu_hwIso;
  std::vector<unsigned short int> _l1mu_hwQual;
  std::vector<short int> _l1mu_bx;

  short int _l1mu_Nmuons;

  //BMTF muon
  std::vector<short int> _bmtfmu_hwpt;
  std::vector<short int> _bmtfmu_hweta;
  std::vector<short int> _bmtfmu_hwphi;
  std::vector<short int> _bmtfmu_glbphi;
  std::vector<short int> _bmtfmu_hwsign;
  std::vector<short int> _bmtfmu_hwqual;
  std::vector<short int> _bmtfmu_link;
  std::vector<short int> _bmtfmu_processor;
  std::vector<short int> _bmtfmu_bx;
  std::vector<short int> _bmtfmu_wheel;

  short int _bmtfmu_Nmuons;

  //BMTF input
  std::vector<short int> _Inputbmtf_phiBX;
  std::vector<short int> _Inputbmtf_phiWheel;
  std::vector<short int> _Inputbmtf_phiSector;
  std::vector<short int> _Inputbmtf_phiStation;
  std::vector<short int> _Inputbmtf_phiAngle;
  std::vector<short int> _Inputbmtf_phiBandAngle;

  short int _Inputbmtf_phiNprimitives;

  std::vector<short int> _Inputbmtf_thetaBX;
  std::vector<short int> _Inputbmtf_thetaWheel;
  std::vector<short int> _Inputbmtf_thetaSector;
  std::vector<short int> _Inputbmtf_thetaStation;
  std::vector<short int> _Inputbmtf_thetaAngle;

  short int _Inputbmtf_thetaNprimitives;

  //OMTF muon
  std::vector<short int> _omtfmu_hwpt;
  std::vector<short int> _omtfmu_hweta;
  std::vector<short int> _omtfmu_hwphi;
  std::vector<short int> _omtfmu_glbphi;
  std::vector<short int> _omtfmu_hwsign;
  std::vector<short int> _omtfmu_hwqual;
  std::vector<short int> _omtfmu_link;
  std::vector<short int> _omtfmu_processor;
  std::vector<short int> _omtfmu_bx;
  std::vector<short int> _omtfmu_wheel;

  short int _omtfmu_Nmuons;

  //EMTF muon
  std::vector<short int> _emtfmu_hwpt;
  std::vector<short int> _emtfmu_hweta;
  std::vector<short int> _emtfmu_hwphi;
  std::vector<short int> _emtfmu_glbphi;
  std::vector<short int> _emtfmu_hwsign;
  std::vector<short int> _emtfmu_hwqual;
  std::vector<short int> _emtfmu_link;
  std::vector<short int> _emtfmu_processor;
  std::vector<short int> _emtfmu_bx;
  std::vector<short int> _emtfmu_wheel;

  short int _emtfmu_Nmuons;

  //KBMTF muon
  std::vector<float> _Kbmtfmu_pt;
  std::vector<float> _Kbmtfmu_eta;
  std::vector<float> _Kbmtfmu_phi;
  std::vector<int>   _Kbmtfmu_dxy;
  std::vector<int>   _Kbmtfmu_approxChi2;
  std::vector<short int> _Kbmtfmu_sector;
  std::vector<short int> _Kbmtfmu_wheel;
  std::vector<short int> _Kbmtfmu_qual;
  std::vector<short int> _Kbmtfmu_bx;

  short int _Kbmtfmu_Nmuons;

  //Tk muon
  std::vector<float> _tkmu_pt;
  std::vector<float> _tkmu_eta;
  std::vector<float> _tkmu_phi;
  std::vector<float> _tkmu_charge;
  std::vector<float> _tkmu_tkiso;

  short int _tkmu_Nmuons;

  //Tk Glb muon
  std::vector<float> _tkglbmu_pt;
  std::vector<float> _tkglbmu_eta;
  std::vector<float> _tkglbmu_phi;
  std::vector<float> _tkglbmu_charge;
  std::vector<float> _tkglbmu_tkiso;

  short int _tkglbmu_Nmuons;

  //Tokens
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > _pileupSummaryToken;
  edm::EDGetTokenT<reco::GenParticleCollection> _genParticleToken;
  edm::EDGetTokenT<l1t::MuonBxCollection> _L1MuonToken;
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> _bmtfMuonToken;
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> _omtfMuonToken;
  edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> _emtfMuonToken;
  edm::EDGetTokenT<L1MuKBMTrackCollection> _KbmtfMuonToken;
  edm::EDGetTokenT<L1MuDTChambPhContainer> _bmtfPhInputToken;
  edm::EDGetTokenT<L1MuDTChambThContainer> _bmtfThInputToken;
  edm::EDGetTokenT<l1t::L1TkMuonParticleCollection> _TkMuonToken;
  edm::EDGetTokenT<l1t::L1TkGlbMuonParticleCollection> _TkGlbMuonToken;

};
