/*******************************************************************
Code: NTupleType.hpp

Description: Contains classes corresponding to different n-tuples
             which are used as the branches of a TTree.

Author: Daniel Santiago-Gonzalez
2013-Jan
*******************************************************************/

/////////////////////////////////////////////////////////////////////////////////
//  This class is for ASICs signals (or hits)
/////////////////////////////////////////////////////////////////////////////////
class ASICS_Signal {

public:
  ASICS_Signal(int Max_Hits){
    this->Max_Hits = Max_Hits;
    MB = new Int_t[Max_Hits];
    Chip = new Int_t[Max_Hits];
    Chan = new Int_t[Max_Hits];
    Energy = new Int_t[Max_Hits];
    Time = new Int_t[Max_Hits];
  };
  
  Int_t Max_Hits;
  Int_t Nhits;
  Int_t* MB;
  Int_t* Chip;
  Int_t* Chan;
  Int_t* Energy;
  Int_t* Time;


  void Reset(){
    Nhits = 0;
    for (Int_t i=0; i<Max_Hits;i++) {
      MB[i] = -1;
      Chip[i] = -1;
      Chan[i] = -1;
      Energy[i] = 0;
      Time[i]  = 0;
    } 	
  };
};





/////////////////////////////////////////////////////////////////////////////////
//  This class is for CAEN ADC and TDC triggers.
/////////////////////////////////////////////////////////////////////////////////
class CAEN_Signal{

public:

  CAEN_Signal(int Max_Hits){
    this->Max_Hits = Max_Hits;
    Module = new Int_t[Max_Hits];
    Channel = new Int_t[Max_Hits];
    Value = new Int_t[Max_Hits];
    ErrorFlag = 0;
  };

  Int_t Max_Hits;
  Int_t Nhits;
  Int_t* Module;
  Int_t* Channel;
  Int_t* Value;
  Int_t ErrorFlag;

  void Reset(){
    ErrorFlag = 0;
    Nhits = 0;
    for (Int_t i=0; i<Max_Hits;i++) {
      Module[i] = -1;
      Channel[i] = -1;
      Value[i] = 0;
    } 	
  };


};



/////////////////////////////////////////////////////////////////////////////////
//  This class is for the physical events.
/////////////////////////////////////////////////////////////////////////////////
class Physical_Event{

public:

  Physical_Event(int Max_Hits){
    this->Max_Hits = Max_Hits;
    FinalE = new Float_t[Max_Hits];
    FinalX = new Float_t[Max_Hits];
    FinalY = new Float_t[Max_Hits];
    FinalZ = new Float_t[Max_Hits];
    PC_V = new Float_t[Max_Hits];
    PC_Z = new Float_t[Max_Hits];
    PC_wire = new Int_t[Max_Hits];
    Time1 = new Int_t[Max_Hits];
    Time2 = new Int_t[Max_Hits];
  };

  Int_t Max_Hits;
  Int_t NHits;
  Float_t* FinalE;
  Float_t* FinalX;
  Float_t* FinalY;
  Float_t* FinalZ;
  Float_t* PC_V;
  Float_t* PC_Z;
  Int_t* PC_wire;
  Int_t* Time1;
  Int_t* Time2;

  void Reset(){
    NHits = 0;
    for (Int_t i=0; i<Max_Hits;i++) {
      FinalE[i] = 0;
      FinalX[i] = 0;
      FinalY[i] = 0;
      FinalZ[i] = -100;
      PC_Z[i] = -100;
      PC_V[i] = 0;
      PC_wire[i] = -1;
      Time1[i] = 0;
      Time2[i] = 0;
    } 	
    return;
  };
};
