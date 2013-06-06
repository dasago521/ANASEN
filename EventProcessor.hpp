class EventProcessor{
public:

  EventProcessor();   // Constructor
  ~EventProcessor();  // Destructor

  // Members
  bool DataExtracted;
  float DataFraction;

  // Methods. SX3 and PC exclusive methods start with 'SX3_' and 'PC_', respectively.
  void ConvertEVTtoROOT(string files_list="data_files.list", int BufferWords=13328);
  void CreateDetectors(int NumQ3Dets, int NumSX3Dets, int NumPCWires, int NumCsIDets, int MaxMultiplicity, string CombinationsSX3andPC="", string CombinationsSX3andCsI="");
  void CreatePhysicalEventsFile(string file_name);
  void Enable3DHits(int Max3DPoints);
  void ExtractData(string HitType="any", bool EmptyChamber=0);
  void GlobalGeometryTest(float SourcePos, int FwdDet, int R1Det, int R2Det, string FileELoss);
  void LoadDetectorsParameters();
  void LoadParamsSelectedDetectors(string DetType, int InitDet, int FinalDet);
  bool OpenROOTFile(string file_name);
  void PC_CalibratePosition(float pos_source_DS, float pos_source_US, string AlphaELossFile, int InitWire, int FinalWire);
  void PC_GetPulserOffsets(float* Voltages, int InitWire, int FinalWire);
  void PC_SetThresholds(float DisplayLimit, int InitWire, int FinalWire);
  void PC_ShowPosSummary();
  void ResetChain();
  void SetBranches();
  void SetChannelMap(string ASICS_map, string CAEN_map);
  void SetChannelMap(string ASICS_map, string CAEN_map, string SX3_geometry_map);
  void SetParamDirectory(string ParamDirectory);
  void SetPeakFindingParameters(int NPeaks, float PeakSigma, float PeakMinAmplitude, float NormalChi2Value=1);
  void SetTime1Cuts(int NT1Cuts, float* T1CutLowerLimit, float* T1CutUpperLimit);
  void SetTreeName(TString TreeName);
  void Show3DHits();
  void ShowCalESummary();
  void ShowMultiplicity();
  void ShowSX3CsICorrelation();
  void ShowRawESummary();
  void ShowTimeSummary();
  void SX3_CalibrateEnergy(int* back_ref, float* SourcePeaks, float MinV, float MaxV, int InitDet=0, int FinalDet=-1, bool GetRelCoeffs=1);
  void SX3_CalibrateFrontBack(float MaxDiff, int Cuts, float *EnergyCuts, float Uncertainty, int InitDet=0, int FinalDet=-1);
  void SX3_CalibratePosition(float MinE, float MaxE, int InitDet=0, int FinalDet=-1);
  void SX3_ForceBchPosRange();
  void SX3_GetPulserOffsets(float* Voltages, int InitDet=0, int FinalDet=-1);
  void SX3_GetPulserOffsetsOneChannel(float* Voltages, int Det_ID, int Chan);
  void SX3_MakeBorderHitsFile(string OutputFile, float Min_Energy, float Max_Energy);
  void SX3_MakeFilteredFile(string FilteredFilesDir, int SelectedDet);
  void SX3_SetThresholds(float DisplayLimit, int InitDet=0, int FinalDet=-1);
  void SX3_ShowGraphsQbQf(int DetID, int Cuts, float *EnergyCuts, float Uncertainty);
  void SX3_ShowGraphsQdQu(int DetID, int Cuts, float *EnergyCuts, float Uncertainty);
  void SX3_ShowHistQb(float MinQ, float MaxQ, int InitDet, int FinalDet);
  void SX3_ShowPosSummary();
  void SX3_ShowPvE(int Det, char GetEnergyFrom, int HistEBins, float HistMinE, float HistMaxE);
  void SX3_TestFrontBackMatch(char GetEnergyFrom, int BinsQ, float MinQ, float MaxQ, int InitDet, int FinalDet);
  void SX3_TestPosCal(char GetEnergyFrom, float MinE, float MaxE, int InitDet=0, int FinalDet=-1);
  void WritePhysicalEventsFile();

private:  
  // Pointers to NTuple and detector objects
  ASICS_Signal* ASICS;
  CAEN_Signal* CAEN;
  Physical_Event* Particle;
  PC_Detector** PC;
  SX3_Detector** SX3;

  
  // Detector related parameters
  bool** Good_SX3_CsI_comb;
  bool** Good_SX3_PC_comb;
  int NumCsIDets;
  int NumPCWires;
  int NumQ3Dets;
  int NumSX3Dets;

  // Histograms and graphs
  TH2F* Corr_SX3_CsI;
  TH2F* HistCsIRawESummary;
  TH2F* HistEvD;
  TH2F* HistPCPosSummary;
  TH2F* HistSX3PosSummary;
  TH2F* HistTimes;
  TH2F** PID_PC;
  TH2F** PID_CsI_SX3;
  TH2F* RF_SX3ID;
  TH2F* RF_SiT;
  TGraph2D* Hits3D;
  TGraph2D* Universe;
  

  // Cuts and gates
  int NT1Cuts;       // For cuts in the RF time (Time1).
  float* T1CutLowerLimit;
  float* T1CutUpperLimit;
  bool ApplyTime1Cuts;

  // Miscelaneous
  TChain* Chain;
  unsigned int Entries;
  ChannelMap* Map;
  int Max3DPoints;
  float NormalChi2Value;
  int NPeaks; 
  string ParamDirectory;
  TFile* PE_ROOT_File;
  TTree* PE_Tree;
  float PeakSigma;
  float PeakMinAmplitude;
  string ROOTFileName;
  TString TreeName;
  float z_228Th;

  // Arrays to store the TBranches' data. Similar arrays are used 
  // in a TreeSelector. These are used in ExtractData().
  Int_t ASICS_Nhits;
  Int_t ASICS_MB[500]; 
  Int_t ASICS_Chip[500]; 
  Int_t ASICS_Chan[500];
  Int_t ASICS_Energy[500];
  Int_t ASICS_Time[500]; 
  Int_t CAEN_Nhits;
  Int_t CAEN_Module[500];    
  Int_t CAEN_Channel[500]; 
  Int_t CAEN_Value[500]; 

  // Switches.
  bool BeamReady;
  bool EnergyLossReady;
  bool FillFinalEvDist;
  bool MassesReady;
  bool Save3DHits;
  bool WritePEFile;

  ClassDef(EventProcessor,1);
};

