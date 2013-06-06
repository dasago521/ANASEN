class SX3_Detector{
public:
  SX3_Detector(string Name, int MaxMultiplicity);
  ~SX3_Detector();

  string Name;        // Name of the detector.
  bool IsOn;          // Simple variable that tells whether the detector could be use or not.
  bool ForceFullHits; // Used for calibrations.
  bool ForceBchPosRange;
  bool HistDiffQCreated;
  bool HistPvECreated;
  bool HistRFTimeCreated;
  char GetEnergyFrom;
  int BackRef;

  // Histograms
  TH2F*** HistBF;
  TH2F*   HistCalESummary;
  TH2F*** HistDiffQ;
  TH1I*   HistMult;
  TH1F*** HistP; 
  TH1F*** HistPEdge;
  TH2F*** HistPvE;
  TH2F*** HistPvEEdge;
  TH1F**  HistQb;
  TH1F**  HistQf;
  TH1I**  HistRawE;
  TH2F*   HistRawESummary;
  TH2I**  HistRawETime;
  TH1F**  HistRFTime;
  TH2F*   HistTimeSummary; 

  // Miscelaneous
  bool ArraysCreated;
  int InternalEvents;
  int MaxMultiplicity;
  unsigned short mult;

  //Pointers to raw data arrays for one event.
  unsigned short* RawE;
  unsigned short* Time;
  short* Channel;

  // Detector's quantities obtained in each event.
  float EventQd, EventQu, EventQb;
  float EventP, EventZ, EventE, EventT;
  int EventBch, EventStp;
  int EventType;
  int EventBorderHit;

  // Channel threshold raw energy.
  unsigned short* Threshold;
  // Pulser offset (makes all channels start from the same point.
  float *q0;
  // Front-back calibration coefficients.
  float *mU, *mD, *oF;
  // Position calibration coefficients.
  int PosPolyOrder;
  float **Cu, **Cd;
  float DetLength;
  //Energy calibration coefficients.
  float* mB, *oB;
  float Elin, Eshift;
 

  //== Methods ===========================================================
  void CreateArrays();
  void CreateHistDiffQ(int BinsQ, float MinQ, float MaxQ);
  void CreateHistPvE(int BinsE, float MinE, float MaxE);
  void CreateHistRFTime(int BinsT, float MinT, float MaxT);
  void DeleteArrays();
  void DeleteHistDiffQ();
  void DeleteHistPvE();
  void DeleteHistRFTime();
  void FillPvE(int strip=-1);
  void FillHistDiffQ();  
  void FillHistPvE();  
  void GetEnergyCoeffs(int NPeaks, float PeakSigma, float MinPeakAmplitude, float* Th228Peaks, float MinQ, float MaxQ);
  void GetEventEP();
  void GetEventQs(string HitType="normal");
  void GetFrontCoeffs(int strip, float MaxQDifference, int Cuts, float *EnergyCuts, float Uncertainty);
  void GetPositionPoly(int strip, float MinE, float MaxE);
  void GetPulserOffsets(int NPeaks, float PeakSigma, float MinPeakAmplitude, float NormalChi2Value, float* Volts, int OnlyThisChannel=-1);
  float GetQbCal();
  float GetQdCal();
  float GetQuCal();
  void GetRelativeCoeffs(int back_ref, float PeakSigma, float MinPeakAmplitude, float MinQ, float MaxQ);
  void LoadCoefficients(string file);
  void Reset();
  void SaveEventQs(int event);
  void SetPeakFindingParameters(int NPeaks, float PeakSigma, float MinPeakAmplitude, float NormalChi2Value=1);
  void SetThresholds(float DisplayLimit);
  void ShowHistDiffQ(int strip);
  void ShowHistPvE(int strip, float MinE, float MaxE, string CanvasExtension="");
  void ShowGraphsQbQf(int Cuts, float *EnergyCuts, float Uncertainty);
  void ShowGraphsQdQu(int Cuts, float *EnergyCuts, float Uncertainty);
  void WriteCoefficients(string file);


private:

  //Number of channels, usually 12: 8 for the front-resistive strips and 4
  //back channels. 
  int NumChan;

  TRandom3* rand;

  // Peak finding parameters.
  int NPeaks;
  float PeakSigma, MinPeakAmplitude, NormalChi2Value;

  // Arrays to store calibration data.
  float* Qd;
  float* Qu;
  float* Qb;
  int* BackSeg;
  int* FrontStp;
  int* TypeOfEvent;
  int* BorderHit;


};


