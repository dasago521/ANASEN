class PhysicalEventProcessor: public TPad{
public:

  PhysicalEventProcessor();   // Constructor
  ~PhysicalEventProcessor();  // Destructor

  // Members
  float DataFraction;
  
  // Methods. 
  void ExtractData();
  bool OpenROOTFile(string file_name);
  //  void PrintParticleReactionRelations();
  void PrintEntry(int entry);
  void ResetChain();
  void ReconstructReaction(int hit, int ParticleIndex);
  void SetBeamEnergyCuts(int NKbCuts, float* KbCutLowerLimit, float* KbCutUpperLimit);
  void SetBeamEnergyCuts(int NKbCuts, float InitBeamEne, float FinalBeamEne);
  void SetBeamParticle(string BeamName, float mb, float Kb_after_win, float win_pos);
  void SetBranches();
  bool SetEnergyLossFile(string ParticleName, string ELossFile); 
  void SetHistParameters(string ParamName, int bins, float min, float max, string ParticleName="");
  void SetLightParticle(string LightPName, float ml);
  void SetMaxHitsPerEvent(int MaxHits);
  void SetParticleCoinc(string Particle1Name, string Particle2Name);
  void SetParamDirectory(string ParamDirectory);
  bool SetPIDCuts(string ParticleName, string File, string CutFormat);
  void SetPossibleReaction(string ReactionName, float mh);
  void SetPossibleReaction(string ReactionName, float mh, int NEexc, float* ListEexc);
  bool SetPTCuts(string ParticleName, string File, string CutFormat);
  void SetTargetParticle(string TargetName, float mt);
  void SetTime1Cuts(int NT1Cuts, float* T1CutLowerLimit, float* T1CutUpperLimit);
  //  void SetTime1HistParams(int Bins, float T1Min, float T1Max);
  void SetTreeName(TString TreeName);
  void ShowBeamEnergySlider(int r);
  void ShowEvTh();
  void ShowEexcSpectra();
  void ShowHistKbT1();
  void ShowHistT1TOF();

private:

  void ExecuteEvent(Int_t event, Int_t px, Int_t py);
  TCanvas* Disp;
  TPad* PadL;
  TPad* PadR;

  // Particle related parameters (energy in MeV, mass in MeV/c^2, distance in cm)
  EnergyLoss* BeamInGas;
  string BeamName;
  bool **CoincMatrix;
  bool **IsParticleInReaction;
  float Kb_after_win;      // Beam kinetic energy after Kapton window.
  EnergyLoss** LightInGas;
  string* LightParticleName;
  float* LightParticleMass;
  float** ListEexc;        // List of energies of excitation.
  float mb, mt, *mh, *ml;  // Masses. Must be in MeV/c^2
  int NLightParticles;
  int NReactions;
  int* NEexc;               // Number of excited states.
  string* ReactionName;
  string TargetName;
  float win_pos;         // Window z-position with respect to forward silicon detectors.
  float z_228Th;
  

  // Histograms and graphs
  TH1F*** HistEexc;
  TH2F*** HistEvTh;
  TH2F* HistKbZr;
  TH2F** HistKbT1;
  TH2F* HistTimes;
  TH2F** HistT1TOF;

  // Cuts and gates
  int NKbCuts;        // For cuts in the beam kinetic energy.
  float* KbCutLowerLimit;
  float* KbCutUpperLimit;
  bool ApplyEnergyCuts;
  int NT1Cuts;       // For cuts in the RF time (Time1).
  float* T1CutLowerLimit;
  float* T1CutUpperLimit;
  bool ApplyTime1Cuts;
  TCutG*** PIDGraphCut;
  TCutG*** PTGraphCut;

  // Miscelaneous
  TChain* Chain;
  int CurrentCut;
  int CurrentReaction;
  unsigned int Entries;
  int HP_Kh_bins;
  float HP_Kh_max, HP_Kh_min;
  int *HP_Kl_bins;
  float *HP_Kl_max, *HP_Kl_min;
  int HP_T1_bins;
  float HP_T1_max, HP_T1_min;
  int HP_TOF_bins;
  float HP_TOF_max, HP_TOF_min;
  KinematicCurve*** KineCurve;
  int MaxHits;
  int MaxReactions;
  int MaxLightParticles;
  int NumPCWires;
  string ParamDirectory;
  int PreviousCut;
  TSlider *Slider;
  TString TreeName;

  // Arrays to store the TBranches' data. Similar arrays are used 
  // in a TreeSelector. These are used in ExtractData().
  Int_t NHits;
  Float_t FinalE[30]; 
  Float_t FinalX[30]; 
  Float_t FinalY[30]; 
  Float_t FinalZ[30]; 
  Float_t PC_Vcor[30]; 
  Float_t PC_Z[30]; 
  Int_t PC_wire[30]; 
  Int_t Time1[30]; 
  Int_t Time2[30]; 

  // Switches.
  bool BeamReady;
  bool EnergyLossReady;
  bool FillFinalEvDist;
  bool MassesReady;
  bool ParticleCoincRequested;
  bool* ReactionReady;
  bool TargetReady;

  ClassDef(PhysicalEventProcessor,1);
};

