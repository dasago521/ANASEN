//Methods for the ProportionalCounter class.

//C and C++ libraries.
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>

//ROOT libraties
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TObject.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TPolyLine3D.h>
#include <TCutG.h>

using namespace std;

//Detector's library
#include "PC_Detector.hpp"

////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
PC_Detector::PC_Detector(string Name, int MaxMultiplicity)
{
  this->Name = Name;
  this->MaxMultiplicity = MaxMultiplicity;
  NumChan = 2;

  //Arrays to store the event data.
  Channel = new short[MaxMultiplicity];
  RawE = new unsigned short[MaxMultiplicity];
  //Initialize them.
  mult = 0;
  for (int m=0; m<MaxMultiplicity; m++) {
    Channel[m] = -1;
    RawE[m] = 0;
  }

  ArraysCreated = 0;
  HistogramsCreated = 0;
  IsOn = 1;
  InternalEvents = 0;
  
  // Arrays for calibration coefficients.
  Threshold = new unsigned short[NumChan];
  q0 = new float[NumChan];
  oV = new float[NumChan];
  mV = new float[NumChan];
  PosPolyOrder = 1;
  Cd = new float[PosPolyOrder+1];
  Cu = new float[PosPolyOrder+1];
  //Initialize them.
  for (int ch=0; ch<NumChan; ch++) {
    Threshold[ch] = 175;
    q0[ch] = 0;
    oV[ch] = 0;
    mV[ch] = 1;
  }
  for (int order=0; order<PosPolyOrder+1; order++) {
    Cd[order] = 1;
    Cu[order] = 1;
  }    
  Elin = 1;
  Eshift = 0;

}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
PC_Detector::~PC_Detector()
{
  if (ArraysCreated)
    DeleteArrays();
  if (HistogramsCreated)
    DeleteHistograms();
  delete[] Channel;
  delete[] RawE;
  delete[] Threshold;
  delete[] q0;
  delete[] oV;
  delete[] mV;
  delete[] Cd;
  delete[] Cu;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PC_Detector::CreateArrays()
{
  FinalE = new float[InternalEvents];
  FinalX = new float[InternalEvents];
  FinalY = new float[InternalEvents];
  FinalZ = new float[InternalEvents];
  Time1 = new float[InternalEvents];
  Time2 = new float[InternalEvents];
  Vd = new float[InternalEvents];
  Vu = new float[InternalEvents];
  
  for (int e=0; e<InternalEvents; e++) {
    Vd[e] = Vu[e] = 0;
    FinalE[e] = 0;
    FinalX[e] = FinalY[e] = 0;
    FinalZ[e] = -100;
    Time1[e] = Time2[e] = 0;
  }
  ArraysCreated = 1;
  cout << "> Arrays to store data of internal events created. Dimension = "
       << InternalEvents << "." << endl;
  return;
}



////////////////////////////////////////////////////////////////////////////
//  Define histograms used in this class.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::CreateHistograms()
{
  // Histogram for the channel multiplicity.
  HistMult = new TH1I(Form("%s.HistMult",Name.c_str()),Form("Multiplicity in %s",Name.c_str()),
		      MaxMultiplicity,-0.5,MaxMultiplicity-0.5);
  HistMult->GetXaxis()->SetTitle("Channels fired per event");
  HistMult->GetXaxis()->CenterTitle();

  //Histograms for raw signals.
  HistRawE = new TH1I*[NumChan];
  for (int ch=0; ch<NumChan; ch++) {
    HistRawE[ch] = new TH1I(Form("%s.HistRawE%d",Name.c_str(),ch),
			    Form("%s | Raw energy (ch=%d)",Name.c_str(),ch),1024,0,4096);
  }
  HistRawESummary = new TH2F(Form("%s.HistRawESummary",Name.c_str()), Form("%s RawE summary",Name.c_str()),
			     NumChan,-0.5,NumChan-0.5, 1024,0,4096);
  HistRawESummary->GetYaxis()->SetTitle("Raw energy [arb. units]");
  HistRawESummary->GetYaxis()->CenterTitle();
  HistRawESummary->GetXaxis()->SetTitle("Channel");
  HistRawESummary->GetXaxis()->CenterTitle(); 
  HistogramsCreated = 1;
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void PC_Detector::DeleteArrays()
{
  delete[] FinalE;
  delete[] FinalX;
  delete[] FinalY;
  delete[] FinalZ;
  delete[] Time1;
  delete[] Time2;
  delete[] Vd;
  delete[] Vu;
  ArraysCreated = 0;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void PC_Detector::DeleteHistograms()
{
  delete HistMult;
  for (int ch=0; ch<NumChan; ch++)
    delete HistRawE[ch];
  delete[] HistRawE;
  delete HistRawESummary;
  HistogramsCreated = 0;
  return;
}



////////////////////////////////////////////////////////////////////////////
// Calculate the path length corrected voltage (~ energy).
////////////////////////////////////////////////////////////////////////////
float PC_Detector::GetCorrectedV()
{ 
  float Vcor=0, theta, rho_SX3, rho_PC=3;
  rho_SX3 = sqrt( pow(EventFinalX,2) + pow(EventFinalY,2) );
  if (EventZ > EventFinalZ) 
    theta = atan((rho_SX3-rho_PC)/(EventZ-EventFinalZ));
  else
    theta = TMath::Pi() + atan((rho_SX3-rho_PC)/(EventZ-EventFinalZ));  // Here atan()<0.
  Vcor = (EventVd + EventVu)*sin(theta);
  //  cout << EventZ << " " << EventFinalZ << " " << rho_SX3 << " " << (EventVd + EventVu) << "  " << Vcor << "  " << theta*180/TMath::Pi() << endl;
  return Vcor;
}


////////////////////////////////////////////////////////////////////////////
// In this function the raw signals and the calibrated voltages are sorted
// into a more accessible way. To get the voltages the coefficients from 
// the pulser data must have been obtained before calling this function.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::GetEventVs()
{ 
  Int_t B[2] = {0x1,0x2};
  Int_t HitKey = 0x0;
  float Vd = 0, Vu = 0;
  int chan;
  EventType = 0;
  EventVd = EventVu = 0;

  for (int mi=0; mi<mult; mi++) {
    chan = Channel[mi];
    HitKey = HitKey | B[chan];
    if (chan==0) 
      Vd = mV[0]*RawE[mi] + oV[0];
    else if (chan==1) 
      Vu = mV[1]*RawE[mi] + oV[1];
  }

  // We only consider hits where we have up-stream and down-stream signals
  // simultaneously.
  if (HitKey==0x3) {
    EventVd = Vd;
    EventVu = Vu;
    EventType = 1;
  }

  return;
}


////////////////////////////////////////////////////////////////////////////
// Calculate the acutal z-position in cm.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::GetEventZ()
{ 
  float Vtot=0, PosUS=-100, PosDS=-100;

  Vtot = EventVd + EventVu;
  if (Vtot>0) {
    PosUS = 0;
    PosDS = 0;
    for (int order=0; order<=PosPolyOrder; order++) {
      PosDS += Cd[order]*pow(EventVd/Vtot, order);
      PosUS += Cu[order]*pow(EventVu/Vtot, order);
    }
  }
  EventZ = (PosDS + PosUS)/2;  
  return;
}


////////////////////////////////////////////////////////////////////////////
// This method assumes two needle sources of 210Po have been placed close
// to the up-stream and down-stream ends of the proportional counter and 
// that the chamber has been filled with gas. The decay alpha particles
// released by this sources are mono-energetic. The PC up- and down-stream
// signals can then be calibrated with respect to the calulated position
// at which the alpha particle crossed the PC wire.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::GetPositionCoeffs(float pos_source_DS, float pos_source_US, string AlphaELossFile)
{
  EnergyLoss* AlphaInGas;
  TCanvas* Display;
  TF1* PolyDS;
  TF1* PolyUS;
  TGraph* GraphUS;
  TGraph* GraphDS;
  TH2F* BackgroundDS;
  TH2F* BackgroundUS;
  int points = 0;
  float Vtot, rho_SX3, rho_PC=3, path_l, K_lr, K_lf, DecayPoint;
  float MinimumDifference=0.030, Expected_z_PC, Po210AlphaEnergy = 5.304;
  // Alpha energy loss file, for energy reconstruction. Exit this function if the file
  // does not exist.
  AlphaInGas = new EnergyLoss(AlphaELossFile);
  if (!AlphaInGas->GoodELossFile) 
    return;

  BackgroundDS = new TH2F("BackgroundDS",Form("Source at %.1f cm",pos_source_DS),200,0,1,400,0,40);
  BackgroundDS->GetXaxis()->SetTitle("V_{d}/V_{tot}");
  BackgroundDS->GetYaxis()->SetTitle("Expected PC position [cm]");
  BackgroundDS->SetStats(0);
  BackgroundUS = new TH2F("BackgroundUS",Form("Source at %.1f cm",pos_source_US),200,0,1,400,0,40);
  BackgroundUS->GetXaxis()->SetTitle("V_{u}/V_{tot}");
  BackgroundUS->GetYaxis()->SetTitle("Expected PC position [cm]");
  BackgroundUS->SetStats(0);
  PolyDS = new TF1("PolyDS",Form("pol%d",PosPolyOrder),0.05,0.95);
  PolyUS = new TF1("PolyUS",Form("pol%d",PosPolyOrder),0.05,0.95);
  GraphUS = new TGraph(); 
  GraphDS = new TGraph(); 
 
  cout << "> Finding position coefficients ..." << endl;

  // Loop over the internal events to fill the graphs and histograms.
  for (int n=0; n<InternalEvents; n++) {
    EventType = 1;
    EventVd = Vd[n];
    EventVu = Vu[n];
    EventFinalE = FinalE[n];
    EventFinalX = FinalX[n];
    EventFinalY = FinalY[n];
    EventFinalZ = FinalZ[n];
    EventTime1 = Time1[n];
    EventTime2 = Time2[n];
    Vtot = Vd[n] + Vu[n];
    rho_SX3 = sqrt( pow(FinalX[n],2) + pow(FinalY[n],2) );
    K_lf = FinalE[n];
    // Determine the decay point.
    if (Vd[n]>Vu[n]) 
      DecayPoint = pos_source_DS;
    else
      DecayPoint = pos_source_US;
    // Calculate the energy loss for the path traveled.
    path_l = sqrt( pow(rho_SX3,2) + pow(DecayPoint-FinalZ[n],2) );
    K_lr = AlphaInGas->GetInitialEnergy(K_lf, path_l, 0.1);
    // Check if the reconstructed energy matches the alpha decay energy.
    if (fabs(K_lr - Po210AlphaEnergy)<MinimumDifference) {
      Expected_z_PC = DecayPoint - rho_PC*(DecayPoint - FinalZ[n])/rho_SX3;
      if (Vd[n]/Vtot>0.15 && Vd[n]/Vtot<0.9 && Vu[n]/Vtot>0.15 && Vu[n]/Vtot<0.9) {
	GraphDS->SetPoint(points, Vd[n]/Vtot, Expected_z_PC);
	GraphUS->SetPoint(points, Vu[n]/Vtot, Expected_z_PC);
	points++;
      }
    }
  }
  cout << "> Showing graphs with " << points << " points." << endl;

  // Show the graphs and fits.   
  Display = new TCanvas(Form("D.%s.PosCoeffs",Name.c_str()),
			Form("%s pos. coeffs.",Name.c_str()),0,0,1500,700);
  Display->Divide(2);

  Display->cd(1);
  BackgroundDS->Draw();
  GraphDS->Draw("p");
  GraphDS->Fit(PolyDS,"RQROB");
  for (int order=0; order<=PosPolyOrder; order++)
    Cd[order] = PolyDS->GetParameter(order);

  Display->cd(2);
  BackgroundUS->Draw();
  GraphUS->Draw("p");
  GraphUS->Fit(PolyUS,"RQROB");
  for (int order=0; order<=PosPolyOrder; order++)
    Cu[order] = PolyUS->GetParameter(order);
 
  // Print the obtaied coefficients.
  cout << "> Parameters found:\n>\tOrder" << setw(6) <<"Cd" << setw(13) <<"Cu" << endl;
  for (int order=0; order<=PosPolyOrder; order++) 
    cout << ">\t" << order << setw(13) << Cd[order] << setw(13) << Cu[order] << endl;

  Display->Update();
  Display->WaitPrimitive();

  delete AlphaInGas;
  delete BackgroundUS;
  delete BackgroundDS;
  delete PolyDS;
  delete PolyUS;
  delete GraphUS;
  delete GraphDS;
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Function that looks for peaks in the RawE spectra and matches the 
// position of such peaks with the set voltages given by the *Volts array.
// The first 4 arguments are related to the peak-finding function of a
// TSpectrum and the last argument is for when the user wants to look
// at an individual channel, usually when there are some complications 
// between *Volts and the peaks found. The default value of OnlyThisChannel
// is -1, meaning that all the channels will be shown (one by one).
////////////////////////////////////////////////////////////////////////////
void PC_Detector::GetPulserOffsets(int NPeaks, float PeakSigma, float MinPeakAmplitude,
				   float NormalChi2Value, float* Volts, int OnlyThisChannel)
{
  TGraph** FitGraph;
  FitGraph = new TGraph*[NumChan];
  
  TCanvas* Display = new TCanvas("PulserOffsets","Pulser data",0,0,1300,600);
  Display->Divide(2);
  TF1* linear;
  Float_t chi2[NumChan];
  Float_t slope[NumChan];
  Float_t offset[NumChan];

  cout << "> Coefficients found:\n"
       << "> ch\tmV (linear)\t oV (offset)" << endl;
  
  for (Int_t ch=0; ch<NumChan; ch++) {

    if (OnlyThisChannel>-1 && OnlyThisChannel!=ch)
      continue;
    
    // Initialize some parameters.
    chi2[ch]=999;
    slope[ch]=1;
    offset[ch]=0;

    Display->cd(1);
    HistRawE[ch]->Draw();
    TSpectrum *s = new TSpectrum();
    Int_t PeaksFound = s->Search(HistRawE[ch],PeakSigma," ",MinPeakAmplitude);
    HistRawE[ch]->SetTitle(Form("Peaks found = %d",PeaksFound));
    Display->Update();
    
    Float_t *xpeaks = s->GetPositionX();
    Float_t Temp=0;
    //Reorder the peaks' x-coordinate from lower to higher x-value.
    for(Int_t i=0; i<PeaksFound; i++) {
      for(Int_t j=i; j<PeaksFound; j++) {
	if (xpeaks[j] < xpeaks[i]) {
	  Temp = xpeaks[i];
	  xpeaks[i] = xpeaks[j];
	  xpeaks[j] = Temp;  
	}
      }
    }
    
    // This section depends on how your data looks like.
    FitGraph[ch] = new TGraph(PeaksFound,&xpeaks[0],&Volts[0]);
        
    FitGraph[ch]->GetXaxis()->SetTitle("Measured signal [arb. units]");
    FitGraph[ch]->GetXaxis()->SetLimits(0,4096);
    FitGraph[ch]->GetYaxis()->SetTitle("Pulser voltage [V]");
    FitGraph[ch]->GetYaxis()->SetLimits(0,10);
    
    linear = new TF1("linear","pol1",0,4096);
    FitGraph[ch]->Fit(linear,"qrROB");
    chi2[ch] = linear->GetChisquare();
    slope[ch] = linear->GetParameter(1);
    offset[ch] = linear->GetParameter(0);
    
    FitGraph[ch]->SetTitle(Form("Ch=%d  |  Chi^2=%f", ch, chi2[ch]));
    
    Display->cd(2);
    FitGraph[ch]->Draw("AP*");
    Display->Update();
    if(chi2[ch]<NormalChi2Value){
      oV[ch] = offset[ch];
      mV[ch] = slope[ch];
      q0[ch] = -offset[ch]/slope[ch];
      cout << "> " << ch << "\t" << mV[ch] << "  " << oV[ch] << endl;
    }

    else 
      cout << "> Warning: ch=" << ch << " has a high chi^2 value " << chi2[ch] << endl;

    Display->WaitPrimitive();
    delete s;
  }// end for(ch)


  for (Int_t ch=0; ch<NumChan; ch++)
    delete FitGraph[ch];
  delete[] FitGraph;
  delete Display;
  return; 
}


////////////////////////////////////////////////////////////////////////////
//  Loads the calibration coefficients prevously stored in a file.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::LoadCoefficients(string file)
{
  string aux;
  ifstream read_coeffs(file.c_str()); 
  if(read_coeffs.is_open()){
    // Offsets from pulser data.
    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> Threshold[i];
    //    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> q0[i];
    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> oV[i];
    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> mV[i];
    // Position coefficients
    read_coeffs >> aux;
    for (int order=0; order<=PosPolyOrder; order++)
      read_coeffs >> Cd[order];    
    read_coeffs >> aux;
    for (int order=0; order<=PosPolyOrder; order++)
      read_coeffs >> Cu[order];
    // Energy calibration
    read_coeffs >> aux >> Elin;
    read_coeffs >> aux >> Eshift;
    //cout << "> Coefficients from \"" << file << "\" loaded." << endl;
    read_coeffs.close();
    IsOn = 1;  // The detector is ready to be used.
  } else{
    cout << "> Warning: File " << file << " couldn't be opened. Only raw histogrmas will be available." << endl;
    IsOn = 0;  // The detector can't be used.
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
// Prints the current values of the event quantities.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::PrintEvent()
{
  cout << setw(5) << "Type=" << setw(8) << EventType
       << setw(5) << "Vd=" << setw(8) << EventVd 
       << setw(5) << "Vu=" << setw(8) << EventVu
       << setw(5) << "Ef=" << setw(8) << EventFinalE
       << setw(5) << "xf=" << setw(8) << EventFinalX
       << setw(5) << "yf=" << setw(8) << EventFinalY
       << setw(5) << "zf=" << setw(8) << EventFinalZ
       << setw(5) << "t1=" << setw(8) << EventTime1
       << setw(5) << "t2=" << setw(8) << EventTime2 << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Reset the event related quantities.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::Reset()
{
  mult = 0;
  for(int m=0; m<MaxMultiplicity; m++){
    RawE[m]=0;
    Channel[m]=-1;
  }
  EventVd = EventVu = 0;
  EventFinalE = EventE = 0;
  EventFinalX = EventFinalY = 0;
  EventFinalZ = EventZ = -100;
  EventTime1 = EventTime2 = 0;
  EventType = 0;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Simple function that saves the 'voltages' of valid events for this wire.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::SaveEvent(int event)
{
  Vd[event] = EventVd;
  Vu[event] = EventVu;
  FinalE[event] = EventFinalE;
  FinalX[event] = EventFinalX;
  FinalY[event] = EventFinalY;
  FinalZ[event] = EventFinalZ;
  Time1[event] = EventTime1;
  Time2[event] = EventTime2;
  return;
}


////////////////////////////////////////////////////////////////////////////
// This functions helps you find the raw energy threshold, one of the most
// basic parameters of the detectors. When extracting data for this detector
// only events with RawE>Threshold[channel] will be saved in memory.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::SetThresholds(float DisplayLimit)
{
  TCanvas* Display = new TCanvas(Form("%s.Disp.FindThresh",Name.c_str()),
				 "Find thresholds",1000,800);
  cout << ">\tDouble-click on the canvas to set threshold value." << endl;
  for (int ch=0; ch<NumChan; ch++) {
    Display->cd(1);
    HistRawE[ch]->GetXaxis()->SetRangeUser(0, DisplayLimit);
    HistRawE[ch]->Draw();
    Display->Update();
    Display->WaitPrimitive();
    cout << ">\tChannel=" << ch <<". Set raw energy threshold = ";
    cin >> Threshold[ch];
  }
  delete Display;
  return;
}



////////////////////////////////////////////////////////////////////////////
//  Writes the calibration coefficients to a file. Be careful when over
//  writing older calibration files.
////////////////////////////////////////////////////////////////////////////
void PC_Detector::WriteCoefficients(string file){
  string aux;
  ofstream write_coeffs(file.c_str());
  // Raw energy threshold
  for (int i=0; i<NumChan; i++) write_coeffs << "Threshold[" << i << "]= " << Threshold[i] << endl;
  // Offsets from pulser data.
  for (int i=0; i<NumChan; i++) write_coeffs << "oV[" << i << "]= " << oV[i] << endl;
  // Linear coeff. from pulser data.
  for (int i=0; i<NumChan; i++) write_coeffs << "mV[" << i << "]= " << mV[i] << endl;
  // Position caibration
  write_coeffs << "Cd[order]= ";
  for (int order=0; order<=PosPolyOrder; order++)
    write_coeffs << Cd[order] << " ";
  write_coeffs << "\n";
  write_coeffs << "Cu[order]= ";
  for (int order=0; order<=PosPolyOrder; order++)
    write_coeffs << Cu[order] << " ";
  write_coeffs << "\n";
  
  //Energy calibration
  write_coeffs << "Elin= " << Elin << endl;
  write_coeffs << "Eshift= " << Eshift << endl;
  write_coeffs.close();
  cout << "> Calibration coefficients written to \"" << file << "\"." << endl;
  return;
}

									 

