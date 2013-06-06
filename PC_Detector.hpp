class PC_Detector{
public:
  PC_Detector(string Name, int MaxMultiplicity);
  ~PC_Detector();

  //== Members ===========================================================
  // Histograms 
  TH1I* HistMult;
  TH1I** HistRawE;
  TH2F* HistRawESummary;

  // Miscelaneous
  int InternalEvents;
  bool IsOn;
  string Name;            // Name of the detector.
  
  // Detector's quantities obtained in each event.
  float EventVd;
  float EventVu;
  float EventFinalE;
  float EventFinalX;
  float EventFinalY;
  float EventFinalZ;
  float EventTime1;
  float EventTime2;  
  float EventZ;
  float EventE;
  int EventType;

  // Channel threshold raw energy.
  unsigned short* Threshold;
  // Offsets from the pulser data (used in q -> q-q0).
  float* q0;
  // Voltage calibration coefficients (replacement of only the pulser offset correction)
  float* oV;
  float* mV;
  // Position calibration coefficients.
  int PosPolyOrder;
  float *Cu, *Cd;
  float DetLength;
  // Energy calibration coefficients.
  float Elin, Eshift;


  // Detector's quantities in one event
  unsigned short mult;    // Multiplicity.
  unsigned short* RawE;   // Ponter to array of raw energy.
  short* Channel;         // Pointer to arrays of channel number (0=downstream, 1=upstream).

  
  //== Methods ===========================================================
  void CreateArrays();
  void CreateHistograms();
  void DeleteArrays();
  void DeleteHistograms();
  float GetCorrectedV();
  void GetEventVs();
  void GetEventZ();
  void GetPositionCoeffs(float pos_source_DS, float pos_source_US, string AlphaELossFile);
  void GetPulserOffsets(int NPeaks, float PeakSigma, float MinPeakAmplitude, float NormalChi2Value, float* Volts, int OnlyThisChannel=-1);
  void LoadCoefficients(string file);
  void PrintEvent();
  void Reset();
  void SaveEvent(int event);
  void SetBeamParameter(float Kb_after_win/*MeV*/, float win_pos/*cm*/);
  void SetPIDCut(string cut_file, string cut_name);
  void SetReactionMasses(float mb, float mt, float mh, float ml);
  void SetThresholds(float DisplayLimit);
  void WriteCoefficients(string file);
  
private:
  bool ArraysCreated;
  bool HistogramsCreated;
  int NumChan;            // Number of channels in one wire (2, up- and down-stream).
  int MaxMultiplicity;    // Maximum multiplicity for each wire in one event.

  // Arrays for the internal events
  float* FinalE;
  float* FinalX;
  float* FinalY;
  float* FinalZ;
  float* Time1;
  float* Time2;
  float* Vd;
  float* Vu;


};
