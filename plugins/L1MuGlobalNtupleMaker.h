
//---------- class declaration----------

class L1MuGlobalNtupleMaker : public edm::EDAnalyzer {
public:
  explicit L1MuGlobalNtupleMaker(const edm::ParameterSet&);
  ~L1MuGlobalNtupleMaker();

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  bool _RunningOnData_;
  const edm::InputTag _PileupSrc;

  edm::Service<TFileService> fs;

  void create_trees();

  // ----------member data ---------------------------
  TH1F* _h_Global_Info;

  //Counters
  int _Nevents_processed;

  //TTree and TTree variables
  TTree *_mytree;

  int _run_number;

  //MC truth
  float _NTruePU;

  //Tokens
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > _pileupSummaryToken;
  edm::EDGetTokenT<edm::TriggerResults> _triggerBits;

};
