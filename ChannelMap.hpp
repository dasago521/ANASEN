class ChannelMap {
public:
  ChannelMap(string ASICS_map, string CAEN_map);

  void GetSX3FwdWorldCoordinates(Int_t DetID, Int_t FrontStp, Float_t DetZ, Float_t& WX, Float_t& WY, Float_t& WZ);
  void GetSX3WorldCoordinates(Int_t Det_ID, Int_t FrontStp, Float_t Det_Z, Float_t& WX, Float_t& WY, Float_t& WZ);
  void IdentifyASICS(Int_t MB, Int_t Chip, Int_t Chan, string& Det_Type, Int_t& Det_ID, Int_t& Det_Ch);
  void IdentifyCAEN(Int_t Module, Int_t Channel, string& Det_Type, Int_t& Det_ID, Int_t& Det_Ch);
  void InitWorldCoordinates(string WorldCoordinatesFilename);

  bool WorldCoordinatesLoaded;

private:
  int Total_ASICS_Channels;
  int Total_CAEN_Channels;
  // Used for world coordinates of SX3 detectors.
  int *SX3_ID;
  float *ZOffset, *XAt0, *XAt4, *YAt0, *YAt4;   
  // Useful random number.
  TRandom3* rand;                               

  Int_t *ASICS_MB,*ASICS_Chip,*ASICS_Chan;
  string* ASICS_Det_Type;
  Int_t *ASICS_Det_ID, *ASICS_Det_Ch;
  Int_t *CAEN_Module, *CAEN_Chan;
  string* CAEN_Det_Type;
  Int_t *CAEN_Det_ID, *CAEN_Det_Ch;
};
