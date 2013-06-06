#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>

#include <iostream>
#include <fstream>
#include <string>

#include "ChannelMap.hpp"

using namespace std;

ChannelMap::ChannelMap(string ASICS_map, string CAEN_map)
{

  rand = new TRandom3();   // Useful random number.

  cout << "> Loading channel maps ..." << endl;
  ifstream MapFile;
  string line;
  int line_counter = 0;
  Total_ASICS_Channels = Total_CAEN_Channels = 0;
  WorldCoordinatesLoaded = 0;

  //===============================================================================================
  // ASICs section.
  MapFile.open(ASICS_map.c_str());
  
  if (!MapFile.is_open())
    cout << ">\tERROR: The ASICs channel map file \"" << ASICS_map << "\" couldn't be opened." << endl;
  else {
    
    cout << ">\tThe ASICs channel map file \"" << ASICS_map << "\" was opened successfully." << endl;

    line_counter = 0;
    do{
      getline(MapFile, line);
      // Only count non-empty lines.
      if (!line.empty())          
	line_counter++; 
    } while (!MapFile.eof());
    MapFile.close();

    Total_ASICS_Channels = line_counter-1;

    ASICS_MB = new Int_t[Total_ASICS_Channels];
    ASICS_Chip = new Int_t[Total_ASICS_Channels];
    ASICS_Chan = new Int_t[Total_ASICS_Channels];
    ASICS_Det_Type = new string[Total_ASICS_Channels];
    ASICS_Det_ID = new Int_t[Total_ASICS_Channels];
    ASICS_Det_Ch = new Int_t[Total_ASICS_Channels];

    MapFile.open(ASICS_map.c_str());
    getline(MapFile,line);    // The 1st line is for comments;
    for (int i=0; i<Total_ASICS_Channels; i++) {
      MapFile >> ASICS_MB[i] >> ASICS_Chip[i] >> ASICS_Chan[i] >> ASICS_Det_Type[i] >> ASICS_Det_ID[i] >> ASICS_Det_Ch[i]; 
      getline(MapFile, line);  // The last column is for comments
    }
    MapFile.close();

    cout << ">\t" << Total_ASICS_Channels << " ASCIs channels loaded." << endl;
    
  }


  //===============================================================================================
  // CAEN section.
  MapFile.open(CAEN_map.c_str());
  
  if (!MapFile.is_open())
    cout << ">\tERROR: The CAEN channel map file \"" << CAEN_map << "\" couldn't be opened." << endl;
  else {
    
    cout << ">\tThe CAEN channel map file \"" << CAEN_map << "\" was opened successfully." << endl;
    
    line_counter = 0;
    do{
      getline(MapFile, line);
      // Only count non-empty lines.
      if (!line.empty())          
	line_counter++; 
    } while (!MapFile.eof());
    MapFile.close();

    Total_CAEN_Channels = line_counter-1;

    CAEN_Module = new Int_t[Total_CAEN_Channels];
    CAEN_Chan = new Int_t[Total_CAEN_Channels];
    CAEN_Det_Type = new string[Total_CAEN_Channels];
    CAEN_Det_ID = new Int_t[Total_CAEN_Channels];
    CAEN_Det_Ch = new Int_t[Total_CAEN_Channels];

    MapFile.open(CAEN_map.c_str());
    getline(MapFile, line);    // The 1st line is for column description;
    for (int i=0; i<Total_CAEN_Channels; i++) {
      MapFile >> CAEN_Module[i] >> CAEN_Chan[i] >> CAEN_Det_Type[i] >> CAEN_Det_ID[i] >> CAEN_Det_Ch[i]; 
      getline(MapFile, line);  // The last column is for comments
    }
    MapFile.close();

    cout << ">\t" << Total_CAEN_Channels << " CAEN channels loaded." << endl;
    

  }

}




//---------------------------------------------------------------------------------------------------------
// Description: This method gets the lab (world) coordinates of a hit in a SX3 detector. 
//              It assumes the file with the geometrical parameters has been loaded.
// - Det_Z is the hit z-position obtained from the SX3 detector, where z is between 0
//   and Det_Lentgth, in cm. The convention is that a low (high) z corresponds to a hit
//   on the down (up) stream side of the detector.
// - WZ is the world z-coordinate (in cm) which is 0 at the position of the forward detectors.
// - In case a SX3 is used as forward detector use GetSX3FwdWorldCoordinates() below.
void ChannelMap::GetSX3WorldCoordinates(Int_t Det_ID, Int_t FrontStp, Float_t Det_Z, Float_t& WX, Float_t& WY, Float_t& WZ)
{
  rand->SetSeed();
  float unit_random = rand->Uniform();

  // Case where the lines in the text file with the world coordinates correspond to
  // the same DetID, i.e. line 0 (neglecting the 1st line with the column descriptions)
  // has data for SX3 detector 0.
  if (Det_ID == SX3_ID[Det_ID]) {
    WZ = ZOffset[Det_ID] + Det_Z; 
    WX = (XAt4[Det_ID] - XAt0[Det_ID])*(FrontStp + unit_random)/4.0 + XAt0[Det_ID];
    WY = (YAt4[Det_ID] - YAt0[Det_ID])*(FrontStp + unit_random)/4.0 + YAt0[Det_ID];
  }
  else {
    // Here one could look for the DetID in a for-loop, but since I use the case
    // above I won't spend time typing that.
  }
  return;
};

// Special case to get the SX3 hit global position before the Q3 were installed as 
// the forward detectors.
void ChannelMap::GetSX3FwdWorldCoordinates(Int_t DetID, Int_t FrontStp, Float_t DetZ, Float_t& WX, Float_t& WY, Float_t& WZ)
{
  // Note that in case of Forward location of X3 XAt0 is a radius offset and XAt4 is angle
  // in degrees with respect to the horizontal axis (beam right is positive).
  float DetAngle = XAt4[DetID]*TMath::Pi()/180;
  float r  = XAt0[DetID] - DetZ;
  float DetWidth = 4; // cm
  // This is the width coordinate, w, which is perpendicular to the detectors' z-axis (across its lenght,
  // going through front strips 1 and 2). It tells us how far from this axis the hit was detected.
  rand->SetSeed();   
  float w = (FrontStp - 2 + rand->Uniform())*DetWidth/4;

  rand->SetSeed();   
  WZ = rand->Uniform(-0.05,0.05);  // This small interval is to account for the detector's thickness ~0.1 cm (probably not necessary)
  WX = r*cos(DetAngle) + w*sin(DetAngle);
  WY = r*sin(DetAngle) - w*cos(DetAngle);
  return;
};



//---------------------------------------------------------------------------------------------------------
//Description: From the channel map provided in the constructor this function returns the
//             detector number and detector channel number for a given motherboard, chip and 
//             chip channel.
void ChannelMap::IdentifyASICS(Int_t MB, Int_t Chip, Int_t Chan, string& Det_Type, Int_t& Det_ID, Int_t& Det_Ch)
{
  for (Int_t i=0; i<Total_ASICS_Channels; i++) {
    if (ASICS_MB[i] == MB && ASICS_Chip[i] == Chip && ASICS_Chan[i] == Chan) {
      Det_Type = ASICS_Det_Type[i];
      Det_ID = ASICS_Det_ID[i];
      Det_Ch = ASICS_Det_Ch[i];
      return;
    }	 
  }
}


//---------------------------------------------------------------------------------------------------------
//Description: From the channel map provided in the constructor this function returns the
//             detector number and detector channel number for a given CAEN module and 
//             channel.
void ChannelMap::IdentifyCAEN(Int_t Module, Int_t Channel, string& Det_Type, Int_t& Det_ID, Int_t& Det_Ch)
{
  for (Int_t i=0; i<Total_CAEN_Channels; i++) {
    if (CAEN_Module[i] == Module && CAEN_Chan[i] == Channel) {
      Det_Type = CAEN_Det_Type[i];
      Det_ID = CAEN_Det_ID[i];
      Det_Ch = CAEN_Det_Ch[i];
      return;
    }	 
  }
}



//---------------------------------------------------------------------------------------------------------

void ChannelMap::InitWorldCoordinates(string WorldCoordinatesFilename)
{
  ifstream WorldCoordFile;
  int line_counter, NumDets;
  string line;

  //Some comments
  int x = 0;

  //Yeah babe
  
  WorldCoordFile.open(WorldCoordinatesFilename.c_str());
  
  if (!WorldCoordFile.is_open()) {
    cout << ">\tERROR: The file \"" << WorldCoordinatesFilename << "\" couldn't be opened." << endl;
    WorldCoordinatesLoaded = 0;
  }
  else {
    cout << ">\tThe file \"" << WorldCoordinatesFilename << "\" was opened successfully." << endl;

    // Count the number of lines with data in the text file.
    line_counter = 0;
    do{
      getline(WorldCoordFile, line);
      // Only count non-empty lines.
      if (!line.empty())          
	line_counter++; 
    } while (!WorldCoordFile.eof());
    WorldCoordFile.close();

    NumDets = line_counter-1;
    SX3_ID = new int[NumDets];
    ZOffset = new float[NumDets];
    XAt0 = new float[NumDets];
    XAt4 = new float[NumDets];
    YAt0 = new float[NumDets];
    YAt4 = new float[NumDets];

    WorldCoordFile.open(WorldCoordinatesFilename.c_str());
    getline(WorldCoordFile, line);    // The 1st line is for column description;
    for (int i=0; i<NumDets; i++) {
      WorldCoordFile >> SX3_ID[i] >> ZOffset[i] >> XAt0[i] >> XAt4[i] >> YAt0[i] >> YAt4[i];
      getline(WorldCoordFile, line);  // The last column is for comments
    }
    WorldCoordFile.close();
    WorldCoordinatesLoaded = 1;   // Ready to use the world coordinates.
    cout << ">\tWorld coordinates loaded for " << NumDets << " SX3 detectors." << endl;
  }

  cout << "Test by oz simpel cout" << endl;
  
  return;
}


