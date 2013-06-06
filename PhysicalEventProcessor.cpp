//Methods for the PhysicalEventProcessor class.

//C and C++ libraries.
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>

//ROOT libraries
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSlider.h>
#include <TStopwatch.h>
#include <TTree.h>

// Useful libraries
#include "../include/EnergyLoss.cpp"
#include "../include/KinematicCurve.cpp"

//PhysicalEventProcessor headers.
#include "PhysicalEventProcessor.hpp"

using namespace std;


///////////////////////////////////////////////////////////////////////////////////
// Constructor. 
///////////////////////////////////////////////////////////////////////////////////
PhysicalEventProcessor::PhysicalEventProcessor()
{
  //-----------------------------------------------------------------
  // Initilize parameters
  ApplyEnergyCuts = 0;
  ApplyTime1Cuts = 0;
  BeamReady = 0;
  CurrentCut = 0;
  DataFraction = 1.0;
  EnergyLossReady = 0;
  MaxHits = 10;
  MaxLightParticles = 10;
  MaxReactions = 10;
  NLightParticles = 0;
  NReactions = 0;
  NumPCWires = 19;
  ParticleCoincRequested = 0;
  PreviousCut = 0;
  TargetReady = 0;

  //-----------------------------------------------------------------
  // Initialize histograms and some histogram related parameters
  float Eres = 0.1;
  HistKbZr = new TH2F("HistKbZr","",400,0,40, (int)(50/Eres),40,90);
  HistKbT1 = new TH2F*[MaxLightParticles];
  HistT1TOF = new TH2F*[MaxLightParticles];
  HistEexc = new TH1F**[MaxReactions];
  HistEvTh = new TH2F**[MaxReactions];
  HP_Kl_bins = new int[MaxLightParticles];
  HP_Kl_max = new float[MaxLightParticles];
  HP_Kl_min = new float[MaxLightParticles];
  for (int p=0; p<MaxLightParticles; p++) {
    HP_Kl_bins[p] = 440;
    HP_Kl_max[p] = 22;
    HP_Kl_min[p] = 0;
  }
  HP_Kh_bins = 1000;
  HP_Kh_max = 100;
  HP_Kh_min = 0;
  HP_T1_bins = 1024;
  HP_T1_max = 4096;
  HP_T1_min = 0;
  HP_TOF_bins = 400;
  HP_TOF_max = 40;
  HP_TOF_min = 0;
  //-----------------------------------------------------------------
  // Initialize matrices.
  IsParticleInReaction = new bool*[MaxLightParticles];
  for (int p=0; p<MaxLightParticles; p++) {
    IsParticleInReaction[p] = new bool[MaxReactions];
    for (int r=0; r<MaxReactions; r++) 
      IsParticleInReaction[p][r] = 0;
  }
  CoincMatrix = new bool*[MaxLightParticles];
  for (int p1=0; p1<MaxLightParticles; p1++) {
    CoincMatrix[p1] = new bool[MaxLightParticles];
    for (int p2=0; p2<MaxLightParticles; p2++) {
      CoincMatrix[p1][p2] = 0;
    }
  }
  ReactionReady = new bool[MaxReactions];
  mh = new float[MaxReactions];
  ml = new float[MaxReactions];
  NEexc = new int[MaxReactions];
  ReactionName = new string[MaxReactions];
  for (int r=0; r<MaxReactions; r++) {
    ReactionReady[r] = 0;
    mh[r] = ml[r] = 0;
    NEexc[r] = 0;
    ReactionName[r] = "";
  }
  ListEexc = new float*[MaxReactions]; 
  LightParticleName = new string[MaxLightParticles];
  LightParticleMass = new float[MaxLightParticles];
  LightInGas = new EnergyLoss*[MaxLightParticles];

  //-----------------------------------------------------------------
  // Initialize the graphical cut pointers to 0.
  PIDGraphCut = new TCutG**[MaxLightParticles];
  for (int p=0; p<MaxLightParticles; p++) {
    PIDGraphCut[p] = new TCutG*[NumPCWires];
    for (int w=0; w<NumPCWires; w++) 
      PIDGraphCut[p][w] = 0;
  } 
  PTGraphCut = new TCutG**[MaxLightParticles];
  for (int p=0; p<MaxLightParticles; p++) {
    PTGraphCut[p] = new TCutG*[NumPCWires];
    for (int w=0; w<NumPCWires; w++) 
      PTGraphCut[p][w] = 0;
  }  
}


///////////////////////////////////////////////////////////////////////////////////
// Destructor. 
// Still need some work.
///////////////////////////////////////////////////////////////////////////////////
PhysicalEventProcessor::~PhysicalEventProcessor()
{
  cout << "Destructor starts" << endl;

  cout << "Destructor ends" << endl;
}  


////////////////////////////////////////////////////////////////////////////
//  Once the ROOT file has been opened and the histograms and data arrays
//  have been created by using CreateHistograms() and CreateArrays(), 
//  respectively, one can use this function to read the data in the ROOT
//  file and store it in the arrays for futher computations and use the
//  histograms for data visualization.
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ExtractData()
{
  cout << "> Extracting data ..." << endl;
  Long64_t TotEntries = Chain->GetEntries();
  int w;
  bool GoodTime1 = 0;
  TStopwatch ProcessTime;
  bool* LightPDetected = new bool[NLightParticles];
  int* LightPHitNumber = new int[NLightParticles];

  cout << ">\tTotal entries: " << TotEntries;
  if (DataFraction<1.0)
    cout << " (taking only " << DataFraction*100 << "% of total entries)";
  cout << "\n";

  // Index for the arrays qb and qu. It will increment if a valid hit is found.
  int event=0;
  ProcessTime.Start();
  
  for(Long64_t n=0; n<TotEntries; n++){
    Chain->GetEntry(n);
    //    PrintEntry(n);
    
    if (n == TMath::Nint(0.10*TotEntries))  cout << ">\t10% through the data" << endl;
    if (n == TMath::Nint(0.25*TotEntries))  cout << ">\t25% through the data" << endl;
    if (n == TMath::Nint(0.50*TotEntries))  cout << ">\t50% through the data" << endl;
    if (n == TMath::Nint(0.75*TotEntries))  cout << ">\t75% through the data" << endl;
    if (n == TMath::Nint(0.90*TotEntries))  cout << ">\t90% through the data" << endl;
    if (n>TMath::Nint(DataFraction*TotEntries)) break;

    if (NHits>MaxHits) 
      continue;

    // Start the event by assuming no particle was detected.
    for (int p=0; p<NLightParticles; p++) {
      LightPDetected[p] = 0;
      LightPHitNumber[p] = -1;
    }
      
    // Sort and identify particle type
    for (int h=0; h<NHits; h++) {
      // Check if there is a time cut.
      if (ApplyTime1Cuts) {
	GoodTime1 = 0;
	for (int c=0; c<NT1Cuts; c++) 
	  if (Time1[h]>=T1CutLowerLimit[c] && Time1[h]<=T1CutUpperLimit[c]) {
	    GoodTime1 = 1;
	    break;
	  }
	if (!GoodTime1)
	  continue;
      }      
      w = PC_wire[h];
      // Identify what type of particle we have in this hit.
      for (int p=0; p<NLightParticles; p++) {
	// Check if the PID cut has been loaded.
	if (PIDGraphCut[p][w]!=0) {
	  // Check if the point is inside the PID cut (it is better to use graphical
	  // cuts that do not overlap).
	  if (PIDGraphCut[p][w]->IsInside(FinalE[h], PC_Vcor[h])) {
	    LightPDetected[p] = 1;
	    LightPHitNumber[p] = h;
	    // Check if the particle could have punched throught the Si detecotrs.
	    //if (PTGraphCut[p][w]->IsInside(FinalE[h], PC_Vcor[h])
	    //	PossiblePunghThrough = 1;
	  }
	}
      }
    }// end for(h)      
      
    // Particle coincidences set by the user have priority over single particle reconstruction.
    if (ParticleCoincRequested) {
      for (int p1=0; p1<NLightParticles; p1++)
	for (int p2=p1; p2<NLightParticles; p2++)
	  if (CoincMatrix[p1][p2] && LightPDetected[p1] && LightPDetected[p2]) {
	    ReconstructReaction(LightPHitNumber[p1], p1);
	    ReconstructReaction(LightPHitNumber[p2], p2);
	  }          
    }
    else {
      // Single particle reaction reconstruction
      for (int p=0; p<NLightParticles; p++)
	if (LightPDetected[p])
	  ReconstructReaction(LightPHitNumber[p], p);
    }



  }// end for(n)

  ProcessTime.Stop();
  cout << ">\tTime processing the event loop: " << ProcessTime.RealTime() << " s" << endl;
  
  return;
  
}




////////////////////////////////////////////////////////////////////////////
//  Simple function that opens a ROOT file (read mode).  The method
//  SetTreeName() must be called before using this method.
////////////////////////////////////////////////////////////////////////////
bool PhysicalEventProcessor::OpenROOTFile(string filename)
{
  TFile* ROOT_File = new TFile(filename.c_str());
  if(ROOT_File->IsOpen()){
    cout << "> ROOT file \"" << filename << "\" opened." << endl; 
    Chain->Add(filename.c_str());
    ROOT_File->Close();
    delete ROOT_File;
    return 1;
  }
  else{
    cout << "> ERROR: ROOT file \"" << filename << "\" could not be opened." << endl;
    delete ROOT_File;
    return 0;
  }
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::PrintEntry(int entry)
{
  cout << "Entry=" << entry << "   NHits=" << NHits << endl;
  for (int h=0; h<NHits; h++) {
    cout << FinalE[h] << "  " << FinalX[h] << "  " << FinalY[h] << "  " << FinalZ[h] 
	 << "  " << PC_Vcor[h] << "  " << PC_Z[h] << "  " << Time1[h] << "  " << Time2[h] << endl;
  }

}



////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ReconstructReaction(int hit, int ParticleIndex)
{
  float Q0, Q, Eexc, path_l, path_b, theta, theta_deg, K_l0, K_b0, z0, TOF_beam, TOF_light, TOF_total;
  float xf = FinalX[hit];
  float yf = FinalY[hit];
  float zf = FinalZ[hit];
  float K_lf = FinalE[hit];
  float zPC = PC_Z[hit];
  float rho_SX3 = sqrt(xf*xf + yf*yf);
  float rho_PC = 3;  // Distance from the anode wire to the z-axis [cm].

  // Check in which possible reactions this light particle is found.
  for (int r=0; r<NReactions; r++) {
    if (IsParticleInReaction[ParticleIndex][r]) {    
      // Before doing anything check if this function is ready to be used.
      if (!BeamReady || !TargetReady || !ReactionReady[r])
	return;
      // Compute light particle scattering angle.
      if (zPC>zf) 
	theta = atan((rho_SX3 - rho_PC)/(zPC - zf));      
      else 
	theta = TMath::Pi() - atan((rho_SX3 - rho_PC)/fabs(zPC - zf));
      // Compute the reaction's z-coordinate (origin of reaction).
      z0 = zf + rho_SX3/tan(theta);  // Works for all values of theta in [0,pi] rad.
      // The reaction z must have a reasonable value within ANASEN.
      if (z0>-5 && z0<50) {
	// Energy reconstruction and TOF calculation.
	// For the beam
	path_b = win_pos - z0;
	K_b0 = BeamInGas->GetFinalEnergy(Kb_after_win, path_b, 0.1);
	TOF_beam = BeamInGas->GetTimeOfFlight(K_b0, path_b, 0.1);
	// For the light particle
	path_l = sqrt( pow(rho_SX3,2) + pow(z0-zf,2) );
	K_l0 = LightInGas[ParticleIndex]->GetInitialEnergy(K_lf, path_l, 0.1);
	TOF_light = LightInGas[ParticleIndex]->GetTimeOfFlight(K_l0, path_l, 0.1);
	
	//	HistKbZr->Fill(z0, K_b0);
	HistT1TOF[ParticleIndex]->Fill(TOF_beam + TOF_light, Time1[hit]);
	HistKbT1[ParticleIndex]->Fill(Time1[hit], K_b0);
	// Check whether the user wants the beam energy within a certain range.
	if (ApplyEnergyCuts) {
	  for (int c=0; c<NKbCuts; c++) {
	    if (K_b0>=KbCutLowerLimit[c] && K_b0<=KbCutUpperLimit[c]) {
	      Q0 = mb + mt - mh[r] - ml[r];
	      Q = K_l0*(1+ml[r]/mh[r]) - K_b0*(1-mb/mh[r]) - 2*sqrt(mb*ml[r]*K_b0*K_l0)*cos(theta)/mh[r];
	      Eexc = Q0 - Q;
	      HistEexc[r][c]->Fill(Eexc);    
	      theta_deg = theta*180/TMath::Pi();
	      HistEvTh[r][c]->Fill(theta_deg, K_l0);
	    }
	  }
	}
      }
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ResetChain()
{
  Chain->Reset();
  delete Chain;
  Chain = new TChain(TreeName);  
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetBeamEnergyCuts(int NKbCuts, float* KbCutLowerLimit, float* KbCutUpperLimit)
{
  this->NKbCuts = NKbCuts;
  this->KbCutLowerLimit = new float[NKbCuts];
  this->KbCutUpperLimit = new float[NKbCuts];
  for (int c=0; c<NKbCuts; c++) {
    this->KbCutLowerLimit[c] = KbCutLowerLimit[c];
    this->KbCutUpperLimit[c] = KbCutUpperLimit[c];
  }

  for (int r=0; r<MaxReactions; r++) {
    HistEexc[r] = new TH1F*[NKbCuts];
    HistEvTh[r] = new TH2F*[NKbCuts];
  }
  ApplyEnergyCuts = 1;
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetBeamEnergyCuts(int NKbCuts, float InitBeamEne, float FinalBeamEne)
{
  this->NKbCuts = NKbCuts;
  this->KbCutLowerLimit = new float[NKbCuts];
  this->KbCutUpperLimit = new float[NKbCuts];
  for (int c=0; c<NKbCuts; c++) {
    this->KbCutLowerLimit[c] = InitBeamEne + c*(FinalBeamEne - InitBeamEne)/NKbCuts;
    this->KbCutUpperLimit[c] = InitBeamEne + (c+1)*(FinalBeamEne - InitBeamEne)/NKbCuts;
  }

  for (int r=0; r<MaxReactions; r++) {
    HistEexc[r] = new TH1F*[NKbCuts];
    HistEvTh[r] = new TH2F*[NKbCuts];
  }
  ApplyEnergyCuts = 1;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Set the beam's name, mass (in MeV/c^2), kinetic energy (in MeV) after the
// Kapton window and the window position (distance) from the forward silicon
// detectors (in cm).
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetBeamParticle(string BeamName, float mb, float Kb_after_win /*MeV*/, float win_pos/*cm*/)
{
  this->BeamName = BeamName;
  this->mb = mb;
  this->Kb_after_win = Kb_after_win;
  this->win_pos = win_pos;
  BeamReady = 1;;
  return;
}



////////////////////////////////////////////////////////////////////////////
// This method is meant to be used after all the ROOT files have been opened.  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetBranches()
{
  // This procedure to is similar to the one used in a TreeSelector. 
  Chain->SetBranchAddress("NHits", &NHits);
  Chain->SetBranchAddress("FinalE", FinalE);
  Chain->SetBranchAddress("FinalX", FinalX);
  Chain->SetBranchAddress("FinalY", FinalY);
  Chain->SetBranchAddress("FinalZ", FinalZ);
  Chain->SetBranchAddress("PC_Vcor", PC_Vcor);
  Chain->SetBranchAddress("PC_Z", PC_Z);
  Chain->SetBranchAddress("PC_wire", PC_wire);
  Chain->SetBranchAddress("Time1", Time1);
  Chain->SetBranchAddress("Time2", Time2);
  return;
}


////////////////////////////////////////////////////////////////////////////
// Create new EnergyLoss objects for the corresponding particle type. As 
// with many other file names in this class, only the actual name is expected
// in this function's arguments, then the parameter directory will be
// inserted before the name. This function returns 1 if the file exist or 0
// otherwise.
////////////////////////////////////////////////////////////////////////////
bool PhysicalEventProcessor::SetEnergyLossFile(string ParticleName, string ELossFile)
{
  bool Status = 0;
  int ParticleIndex = -1;
  ELossFile = ParamDirectory + ELossFile;
  // Check if the particle name is "beam" or if it corresponds to a light particle name.
  if (ParticleName==BeamName) {
    BeamInGas = new EnergyLoss(ELossFile, mb);
    if (BeamInGas->GoodELossFile)
      Status = 1;
  }
  else {
    for (int p=0; p<NLightParticles; p++) 
      if (LightParticleName[p]==ParticleName)
	ParticleIndex = p;
    
    if (ParticleIndex==-1) {
      cout << "> WARNING: Particle '" << ParticleName << "' was not found in the particles' list."  
	   << " Energy loss file was not loaded." << endl;
    }
    else {
      LightInGas[ParticleIndex] = new EnergyLoss(ELossFile, LightParticleMass[ParticleIndex]);
      if (LightInGas[ParticleIndex]->GoodELossFile)
	Status = 1;
    }
  }
  if (Status==0)
    cout << "> WARNING: Energy loss file " << ELossFile << " not found." << endl;
  return Status;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetHistParameters(string ParamName, int bins, float min, float max, string Option)
{
  if (ParamName=="T1" || ParamName=="Time1") {
    HP_T1_bins = bins;
    HP_T1_max = max;
    HP_T1_min = min;
  } 
  else if (ParamName=="TOF") {
    HP_TOF_bins = bins;
    HP_TOF_max = max;
    HP_TOF_min = min;
  }
  else if (ParamName=="Kh") {
    HP_Kh_bins = bins;
    HP_Kh_max = max;
    HP_Kh_min = min;
  }
  else if (ParamName=="Kl") {
    for (int p=0; p<NLightParticles; p++) {
      if (LightParticleName[p]==Option) {
	HP_Kl_bins[p] = bins;
	HP_Kl_max[p] = max;
	HP_Kl_min[p] = min;
      }
    }
  } 
  else if (ParamName=="Ex") {
    for (int p=0; p<NLightParticles; p++) {
      if (LightParticleName[p]==Option) {
	HP_Kl_bins[p] = bins;
	HP_Kl_max[p] = max;
	HP_Kl_min[p] = min;
      }
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetLightParticle(string LightPName, float ml) 
{
  float Eres = 0.1;
  LightParticleName[NLightParticles] = LightPName;
  LightParticleMass[NLightParticles] = ml;
  HistKbT1[NLightParticles] = new TH2F(Form("HistKbT1_%d",NLightParticles),
				       Form("Kb-T1 for %s",LightPName.c_str()),
				       HP_T1_bins,HP_T1_min,HP_T1_max, HP_Kh_bins,HP_Kh_min,HP_Kh_max);
  HistT1TOF[NLightParticles] = new TH2F(Form("HistT1TOF_%d",NLightParticles),
				       Form("T1-TOF for %s",LightPName.c_str()),
				       HP_TOF_bins,HP_TOF_min,HP_TOF_max, HP_T1_bins,HP_T1_min,HP_T1_max);
  HistT1TOF[NLightParticles]->GetXaxis()->SetTitle(Form("TOF_{%s} + TOF_{beam} [ns]",LightPName.c_str()));
  NLightParticles++;
  return;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetMaxHitsPerEvent(int MaxHits)
{
  this->MaxHits = MaxHits;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Set the directory where the parameters such as calibration coefficients
//  can be found.
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetParamDirectory(string ParamDirectory)
{
  this->ParamDirectory = ParamDirectory;
  return;
}

 
////////////////////////////////////////////////////////////////////////////
 
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetParticleCoinc(string Particle1Name, string Particle2Name)
{
  for (int p1=0; p1<NLightParticles; p1++) 
    if (LightParticleName[p1]==Particle1Name)
      for (int p2=0; p2<NLightParticles; p2++) 
	if (LightParticleName[p2]==Particle2Name) {
	  CoincMatrix[p1][p2] = 1;
	  CoincMatrix[p2][p1] = 1;
	}
  ParticleCoincRequested = 1;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Loads the PID graphical cuts from a root file in the ParamDirectory for
// a selected particle and with a specific string format.
////////////////////////////////////////////////////////////////////////////
bool PhysicalEventProcessor::SetPIDCuts(string LightParticleName, string File, string CutFormat)
{
  int ParticleIndex = -1;
  bool FileStatus = 0;
  string CutName;

  for (int p=0; p<NLightParticles; p++) 
    if (this->LightParticleName[p]==LightParticleName)
      ParticleIndex = p;
  if (ParticleIndex==-1) {
    cout << "> WARNING: Particle '" << LightParticleName << "' was not found in the particles' list." 
	 << " PID cuts were not set." << endl;
    return 0;
  }
  
  File = ParamDirectory + File;
  TFile FCuts(File.c_str());
  if (FCuts.IsOpen()) {
    FileStatus = 1;
    //cout << "> Getting PID cuts from: " << File <<endl;
    for (int w=0; w<NumPCWires; w++) {
      CutName = CutFormat + Form("%d",w);
      // Reset the pointer to zero.
      PIDGraphCut[ParticleIndex][w] = 0;
      // Then try to get it from the file with the graphical cuts.
      PIDGraphCut[ParticleIndex][w] = (TCutG*)FCuts.Get(CutName.c_str());
      // If the cut was not found send a warning message.
      if (PIDGraphCut[ParticleIndex][w]==0) 
	cout << "> WARNING: Graphical cut " << CutName << " does not exist."<< endl;
    }
    FCuts.Close();
  }
  else 
    cout << "> WARNING: File with PID cuts " << File << " does not exist." <<endl;
  return FileStatus;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetPossibleReaction(string ReactionName, float mh)
{
  bool FoundLightPart = 0;
  int PosTarget, Kl_bins;
  string LightPName = "";
  float Kl_max, Kl_min;
  float Eres = 0.1;
  this->ReactionName[NReactions] = ReactionName;
  this->mh[NReactions] = mh;
  // Determine the which of the light particles was involved in this reaction.
  // 1st find the target name position
  PosTarget = ReactionName.find(TargetName);
  // 2nd look for the light particle name after the target name position.
  for (int p=0; p<NLightParticles; p++)
    if (ReactionName.find(LightParticleName[p], PosTarget+1) < ReactionName.length()) {
      LightPName = LightParticleName[p];
      ml[NReactions] = LightParticleMass[p];
      FoundLightPart = 1;
      IsParticleInReaction[p][NReactions] = 1;
      Kl_bins = HP_Kl_bins[p];
      Kl_max = HP_Kl_max[p];
      Kl_min = HP_Kl_min[p];
    }
  if (FoundLightPart) {
    // Create the histograms for this reaction.
    for (int c=0; c<NKbCuts; c++) {
      HistEexc[NReactions][c] =	new TH1F(Form("HistEexc_%d_%d",NReactions,c), 
					 Form("%s | K_{br}#in[%.1f,%.1f] MeV", ReactionName.c_str(),KbCutLowerLimit[c],KbCutUpperLimit[c]),
					 (int)(25/Eres),-10,15);
      HistEexc[NReactions][c]->GetXaxis()->SetTitle("Excitation energy [MeV]");
      HistEvTh[NReactions][c] =	new TH2F(Form("HistEvTh_%d_%d",NReactions,c),
					 Form("%s | K_{br}#in[%.1f,%.1f] MeV", ReactionName.c_str(),KbCutLowerLimit[c],KbCutUpperLimit[c]),
					 180,0,180, Kl_bins,Kl_min,Kl_max);
      HistEvTh[NReactions][c]->GetXaxis()->SetTitle("Lab. scattering angle [deg]");
      HistEvTh[NReactions][c]->GetYaxis()->SetTitle(Form("Reconstructed E_{%s} [MeV]",LightPName.c_str()));
    }
    
    ReactionReady[NReactions] = 1;
    cout << "> Reaction " << ReactionName << " ready!" << endl;
    NReactions++;
  }
  else 
    cout << "> WARNING: Could not find light particle for reaction " << ReactionName << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetPossibleReaction(string ReactionName, float mh, int NEexc, float* ListEexc)
{
  bool FoundLightPart = 0;
  string LightPName = "";
  int PosTarget, Kl_bins;
  float Kl_max, Kl_min;
  float Eres = 0.1;
  this->ReactionName[NReactions] = ReactionName;
  this->mh[NReactions] = mh;
  // Determine the which of the light particles was involved in this reaction.
  // 1st find the target name position
  PosTarget = ReactionName.find(TargetName);
  // 2nd look for the light particle name after the target name position.
  for (int p=0; p<NLightParticles; p++)
    if (ReactionName.find(LightParticleName[p], PosTarget+1) < ReactionName.length()) {
      LightPName = LightParticleName[p];
      ml[NReactions] = LightParticleMass[p];
      FoundLightPart = 1;
      IsParticleInReaction[p][NReactions] = 1;
      Kl_bins = HP_Kl_bins[p];
      Kl_max = HP_Kl_max[p];
      Kl_min = HP_Kl_min[p];
    }
  if (FoundLightPart) {
    // Create the histograms for this reaction.
    for (int c=0; c<NKbCuts; c++) {
      HistEexc[NReactions][c] =	new TH1F(Form("HistEexc_%d_%d",NReactions,c), 
					 Form("%s | K_{br}#in[%.1f,%.1f] MeV", ReactionName.c_str(),KbCutLowerLimit[c],KbCutUpperLimit[c]),
					 (int)(25/Eres),-10,15);
      HistEexc[NReactions][c]->GetXaxis()->SetTitle("Excitation energy [MeV]");
      HistEvTh[NReactions][c] =	new TH2F(Form("HistEvTh_%d_%d",NReactions,c),
					 Form("%s | K_{br}#in[%.1f,%.1f] MeV", ReactionName.c_str(),KbCutLowerLimit[c],KbCutUpperLimit[c]),
					 180,0,180, Kl_bins,Kl_min,Kl_max);
      HistEvTh[NReactions][c]->GetXaxis()->SetTitle("Lab. scattering angle [deg]");
      HistEvTh[NReactions][c]->GetYaxis()->SetTitle(Form("Reconstructed E_{%s} [MeV]",LightPName.c_str()));      
    }
    // Get the list of excitation energies.
    this->NEexc[NReactions] = NEexc;
    this->ListEexc[NReactions] = new float[NEexc];
    for (int e=0; e<NEexc; e++)
      this->ListEexc[NReactions][e] = ListEexc[e];
    // This reaction is ready to be used.
    ReactionReady[NReactions] = 1;
    cout << "> Reaction " << ReactionName << " with " << NEexc << " excitation energies ready!" << endl;
    NReactions++;
  }
  else 
    cout << "> WARNING: Could not find light particle for reaction " << ReactionName << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Loads the PT (punch through) graphical cuts from a root file in the 
// ParamDirectory for a selected particle and with a specific name format.
////////////////////////////////////////////////////////////////////////////
bool PhysicalEventProcessor::SetPTCuts(string LightParticleName, string File, string CutFormat)
{
  int ParticleIndex = -1;
  bool FileStatus = 0;
  string CutName;

  for (int p=0; p<NLightParticles; p++) 
    if (this->LightParticleName[p]==LightParticleName)
      ParticleIndex = p;
  if (ParticleIndex==-1) {
    cout << "> WARNING: Particle '" << LightParticleName << "' was not found in the particles' list." 
	 << " PT cuts were not set." << endl;
    return 0;
  }
  
  File = ParamDirectory + File;
  TFile FCuts(File.c_str());
  if (FCuts.IsOpen()) {
    FileStatus = 1;
    for (int w=0; w<NumPCWires; w++) {
      CutName = CutFormat + Form("%d",w);
      // Reset the pointer to zero.
      PTGraphCut[ParticleIndex][w] = 0;
      // Then try to get it from the file with the graphical cuts.
      PTGraphCut[ParticleIndex][w] = (TCutG*)FCuts.Get(CutName.c_str());
      // If the cut was not found send a warning message.
      if (PTGraphCut[ParticleIndex][w]==0) 
	cout << "> WARNING: Graphical cut " << CutName << " does not exist."<< endl;
    }
    FCuts.Close();
  }
  else 
    cout << "> WARNING: File with PT cuts " << File << " does not exist." <<endl;
  return FileStatus;
}


////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetTargetParticle(string TargetName, float mt)
{
  this->TargetName = TargetName;
  this->mt = mt;
  TargetReady = 1;
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetTime1Cuts(int NT1Cuts, float* T1CutLowerLimit, float* T1CutUpperLimit)
{
  this->NT1Cuts = NT1Cuts;
  this->T1CutLowerLimit = new float[NT1Cuts];
  this->T1CutUpperLimit = new float[NT1Cuts];
  for (int c=0; c<NT1Cuts; c++) {
    this->T1CutLowerLimit[c] = T1CutLowerLimit[c];
    this->T1CutUpperLimit[c] = T1CutUpperLimit[c];
  }
  ApplyTime1Cuts = 1;
  return;
}



////////////////////////////////////////////////////////////////////////////
//  Set the name of the TTree used in the root files.
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::SetTreeName(TString TreeName)
{
  this->TreeName = TreeName;
  Chain = new TChain(TreeName);  
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ShowBeamEnergySlider(int r)
{
  float Kb_average;
  CurrentReaction = r;
  // Create the pointers for the kinematic curve objects for each reaction,
  // beam energy cut and energy of excitation.
  KineCurve = new KinematicCurve**[NKbCuts];
  for (int c=0; c<NKbCuts; c++) {
    Kb_average = (KbCutLowerLimit[c] + KbCutUpperLimit[c])/2;
    KineCurve[c] = new KinematicCurve*[NEexc[r]];
    if (NEexc[r]==0) {
      KineCurve[c][0] = new KinematicCurve();
      KineCurve[c][0]->SetReaction(Kb_average, mb, mt, mh[r], ml[r], 0);
      KineCurve[c][0]->GetCurve(180);
      KineCurve[c][0]->SetLineStyle(2);
      KineCurve[c][0]->SetLineWidth(2);
      KineCurve[c][0]->SetLineColor(kGray+2);
    }
    else {
      for (int e=0; e<NEexc[r]; e++) {
	KineCurve[c][e] = new KinematicCurve();
	KineCurve[c][e]->SetReaction(Kb_average, mb, mt, mh[r], ml[r], ListEexc[r][e]);
	KineCurve[c][e]->GetCurve(180);
	KineCurve[c][e]->SetLineStyle(2);
	KineCurve[c][e]->SetLineWidth(2);
	KineCurve[c][e]->SetLineColor(kGray+2);
      }    
    }
  }      
  // Create canvas and pads.
  Disp = new TCanvas(Form("Disp.EvThR%d",r),Form("%s EvTh",ReactionName[r].c_str()),1500,750);
  Disp->cd();
  PadL = new TPad("PadL","PadL",0.0,0.0,0.6,0.93);
  PadL->Draw();
  Disp->cd();
  PadR = new TPad("PadR","PadR",0.6,0.0,1.0,0.93);
  PadR->Draw();
  // Left pad
  PadL->cd();
  HistEvTh[r][0]->Draw("colz");
  if (NEexc[r]==0) 
    KineCurve[0][0]->Draw("l");
  else {
    for (int e=0; e<NEexc[r]; e++) 
      KineCurve[0][e]->Draw("l same");
  }
  // Right pad
  PadR->cd();
  HistEexc[r][0]->Draw();
  // Slider
  Disp->cd();
  Slider = new TSlider("Slider","Kb",0.1,0.95,0.9,1.0); 
  Slider->SetRange(0,1.0/NKbCuts);
  Slider->SetObject(this); 

  return;
}

// Overriding the ExecuteEvent function of TPad.
void PhysicalEventProcessor::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  float SliderPos = (Slider->GetMinimum() + Slider->GetMaximum())/2;
  for (int c=0; c<NKbCuts; c++) {
    if (SliderPos>1.0*c/NKbCuts && SliderPos<1.0*(c+1)/NKbCuts) {
      CurrentCut = c;
      break;
    }
  }
  if (CurrentCut!=PreviousCut) {
    PadL->cd();
    HistEvTh[CurrentReaction][CurrentCut]->Draw("colz");
    if (NEexc[CurrentReaction]==0) 
      KineCurve[CurrentCut][0]->Draw("l");
    else {
      for (int e=0; e<NEexc[CurrentReaction]; e++) 
	KineCurve[CurrentCut][e]->Draw("l same");
    }
    PadR->cd();
    HistEexc[0][CurrentCut]->Draw();
  }
  //  cout << SliderPos << "  " << CurrentCut << endl;
  
  PreviousCut = CurrentCut;
  //int SelectedCut
  return;
}

////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ShowEvTh()
{
  TCanvas** Display;
  KinematicCurve**** KC;
  float Kb_average;
  if (ApplyEnergyCuts && NReactions>0) {
    Display = new TCanvas*[NReactions];
    KC = new KinematicCurve***[NReactions];
    for (int r=0; r<NReactions; r++) {
      // Create the pointers for the kinematic curve objects for each reaction,
      // beam energy cut and energy of excitation.
      KC[r] = new KinematicCurve**[NKbCuts];
      for (int c=0; c<NKbCuts; c++) {
	Kb_average = (KbCutLowerLimit[c] + KbCutUpperLimit[c])/2;
	KC[r][c] = new KinematicCurve*[NEexc[r]];
	if (NEexc[r]==0) {
	  KC[r][c][0] = new KinematicCurve();
	  KC[r][c][0]->SetReaction(Kb_average, mb, mt, mh[r], ml[r], 0);
	}
	else {
	  for (int e=0; e<NEexc[r]; e++) {
	    KC[r][c][e] = new KinematicCurve();
	    KC[r][c][e]->SetReaction(Kb_average, mb, mt, mh[r], ml[r], ListEexc[r][e]);
	  }    
	}
      }      
      Display[r] = new TCanvas(Form("Disp.EvThR%d",r),Form("%s EvTh",ReactionName[r].c_str()),1400,700);
      if (NKbCuts<=3)
	Display[r]->Divide(NKbCuts);
      else
	Display[r]->Divide((int)ceil(NKbCuts/3.0),3);
      for (int c=0; c<NKbCuts; c++) {
	Display[r]->cd(c+1);
	HistEvTh[r][c]->Draw("colz");
	if (NEexc[r]==0) {
	  KC[r][c][0]->GetCurve(180);
	  KC[r][c][0]->SetLineStyle(2);
	  KC[r][c][0]->SetLineWidth(2);
	  KC[r][c][0]->SetLineColor(kGray+2);
	  KC[r][c][0]->Draw("l same");
	}
	else {
	  for (int e=0; e<NEexc[r]; e++) {
	    KC[r][c][e]->GetCurve(180);
	    KC[r][c][e]->SetLineStyle(2);
	    KC[r][c][e]->SetLineWidth(2);
	    KC[r][c][e]->SetLineColor(kGray+2);
	    KC[r][c][e]->Draw("l same");
	  }
	}
      }
    }// end for(r)
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ShowEexcSpectra()
{
  TCanvas** Display;
  if (ApplyEnergyCuts && NReactions>0) {
    Display = new TCanvas*[NReactions];
    for (int r=0; r<NReactions; r++) {
      Display[r] = new TCanvas(Form("Disp.EexcR%d",r),Form("%s Eexc",ReactionName[r].c_str()),1400,700);
      if (NKbCuts<=3)
	Display[r]->Divide(NKbCuts);
      else
	Display[r]->Divide((int)ceil(NKbCuts/3.0),3);
      for (int c=0; c<NKbCuts; c++) {
	Display[r]->cd(c+1);
	HistEexc[r][c]->Draw();
      }
    }
  }
  return;
}



////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ShowHistKbT1()
{
  TCanvas** Display;
  Display = new TCanvas*[NLightParticles];
  for (int p=0; p<NLightParticles; p++) {
    Display[p] = new TCanvas(Form("Disp.KbT1_%d",p),Form("Kb-T1 %s",LightParticleName[p].c_str()),1400,700);
    HistKbT1[p]->Draw("colz");
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PhysicalEventProcessor::ShowHistT1TOF()
{
  TCanvas** Display;
  Display = new TCanvas*[NLightParticles];
  for (int p=0; p<NLightParticles; p++) {
    Display[p] = new TCanvas(Form("Disp.T1TOF_%d",p),Form("T1-TOF %s",LightParticleName[p].c_str()),1400,700);
    HistT1TOF[p]->Draw("colz");
  }
  return;
}





