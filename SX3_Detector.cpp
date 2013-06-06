//Methods for the SX3_Detector.

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
#include <TF2.h>
#include <TLine.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>

using namespace std;

//Detector's library
#include "SX3_Detector.hpp"

////////////////////////////////////////////////////////////////////////////
//  Simple constructor just needs the name of the object.
////////////////////////////////////////////////////////////////////////////
SX3_Detector::SX3_Detector(string Name, int MaxMultiplicity)
{
  NumChan = 12;
  int NumBCh = 4; 
  int NumStp = 4;
  this->Name = Name;
  this->MaxMultiplicity = MaxMultiplicity;

  //Arrays to store the event data.
  Channel = new short[MaxMultiplicity];
  RawE = new unsigned short[MaxMultiplicity];
  Time = new unsigned short[MaxMultiplicity];
  //Initialize them.
  mult = 0;
  for (int i=0; i<MaxMultiplicity; i++) {
    Channel[i] = -1;
    RawE[i] = 0;
    Time[i] = 0;
  }

  BackRef = -1;
  EventType = 0;  

  ArraysCreated = 0;
  ForceFullHits = 0;
  ForceBchPosRange = 0;
  IsOn = 1;
  InternalEvents = 0;
  HistDiffQCreated = 0;
  HistPvECreated = 0;
  GetEnergyFrom = 'b';

  rand = new TRandom3();

  // Arrays for calibration coefficients.
  Threshold = new unsigned short[NumChan];
  q0 = new float[NumChan];
  mU = new float[NumStp];
  mD = new float[NumStp];
  oF = new float[NumStp];
  mB = new float[NumBCh];
  oB = new float[NumBCh];
  Cu = new float*[NumStp];
  Cd = new float*[NumStp];
  PosPolyOrder = 1;
  for (int stp=0; stp<NumStp; stp++) {
    Cd[stp] = new float[PosPolyOrder+1];
    Cu[stp] = new float[PosPolyOrder+1];
  }

  //Initialize the calibration coefficients, linears=1 and offsets=0.
  for (int ch=0; ch<NumChan; ch++) {
    Threshold[ch] = 100;
    q0[ch] = 0;
  }

  for(int i=0; i<4; i++){
    mB[i]=1;
    oB[i]=0;
    mU[i]=1;
    mD[i]=1;
    oF[i]=0;    
    for (int j=0; j<PosPolyOrder+1; j++) {
      Cd[i][j]=1;
      Cu[i][j]=1;
    }    
  }
  Elin=1; Eshift=0;
  DetLength = 7.5;        // Front strips length in cm.

  //- Basic histograms -----------------------------------------------------------------------------------
  // Histogram for the channel multiplicity.
  HistMult = new TH1I(Form("%s.HistMult",Name.c_str()),Form("Multiplicity in %s",Name.c_str()),
		      MaxMultiplicity,-0.5,MaxMultiplicity-0.5);
  HistMult->GetXaxis()->SetTitle("Channels fired per event");
  HistMult->GetXaxis()->CenterTitle();
  // Histograms for each channel.
  HistRawE = new TH1I*[NumChan];
  for (int ch=0; ch<NumChan; ch++)
    HistRawE[ch] = new TH1I(Form("%s.HistRawE%d",Name.c_str(),ch),
			    Form("%s | Raw energy (ch=%d)",Name.c_str(),ch),1024,0,16384);
  
  // Summary histograms
  HistCalESummary = new TH2F(Form("%s.HistCalESummary",Name.c_str()), Form("%s CalE summary",Name.c_str()),
			     4,-0.5,4-0.5, 1000,0,10);
  HistCalESummary->GetYaxis()->SetTitle("Energy [MeV]");
  HistCalESummary->GetYaxis()->CenterTitle();
  HistCalESummary->GetXaxis()->SetTitle("Channel");
  HistCalESummary->GetXaxis()->CenterTitle();
  
  HistRawESummary = new TH2F(Form("%s.HistRawESummary",Name.c_str()), Form("%s RawE summary",Name.c_str()),
			     NumChan,-0.5,NumChan-0.5, 4096,0,16384);
  HistRawESummary->GetYaxis()->SetTitle("Raw energy [arb. units]");
  HistRawESummary->GetYaxis()->CenterTitle();
  HistRawESummary->GetXaxis()->SetTitle("Channel");
  HistRawESummary->GetXaxis()->CenterTitle();
  
  
  HistTimeSummary = new TH2F(Form("%s.HistTimeSummary",Name.c_str()), Form("%s Time summary",Name.c_str()),
			     NumChan,-0.5,NumChan-0.5, 4096,0,16384);
  HistTimeSummary->GetYaxis()->SetTitle("Time [arb. units]");
  HistTimeSummary->GetYaxis()->CenterTitle();
  HistTimeSummary->GetXaxis()->SetTitle("Channel");
  HistTimeSummary->GetXaxis()->CenterTitle();

  // Histograms for the back 'charge' (~energy).
  HistQb = new TH1F*[5];
  for(int bch = 0; bch<4; bch++)
    HistQb[bch] = new TH1F(Form("%s.HistQb%d",Name.c_str(),bch),Form("Back charge (bch=%d)",bch),1000,0,10000);
  HistQb[4] = new TH1F(Form("%s.HistQb4",Name.c_str()),"Back charge (all back channels)",1000,0,10000);

  // Histograms for the front 'charge' (~energy).
  HistQf = new TH1F*[5];
  for(int stp = 0; stp<4; stp++)
    HistQf[stp] = new TH1F(Form("%s.HistQf%d",Name.c_str(),stp),Form("Front charge (stp=%d)",stp),1000,0,10000);
  HistQf[4] = new TH1F(Form("%s.HistQf4",Name.c_str()),"Front charge (all back channels)",1000,0,10000);

}

////////////////////////////////////////////////////////////////////////////
//  Destructor
////////////////////////////////////////////////////////////////////////////
SX3_Detector::~SX3_Detector()
{
#if 0
  delete[] Channel;
  delete[] RawE;
  delete[] Time;
  delete rand;
  // Arrays for calibration coefficients.
  delete[] Threshold;
  delete[] q0;
  delete[] mB;
  delete[] oB;
  for (int stp=0; stp<NumStp; stp++) {
    delete[] Cd[stp];
    delete[] Cu[stp];
  }
  delete[] Cd;
  delete[] Cu;
  delete[] mU;
  delete[] mD;
  delete[] oF;  
  
  delete HistMult;
  delete HistCalESummary;
  delete HistRawESummary;
  delete HistTimeSummary;
  for (int ch=0; ch<NumChan; ch++)
    delete HistRawE[ch];
  delete[] HistRawE;

  for (int ch=0; ch<5; ch++)
    delete HistQb[ch];
  delete[] HistQb;

  if (HistPvECreated) {
    DeleteHistPvE();
  }
      
#endif
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::CreateArrays()
{
  Qd = new float[InternalEvents];
  Qu = new float[InternalEvents];
  Qb = new float[InternalEvents];
  BackSeg = new int[InternalEvents];
  FrontStp = new int[InternalEvents];
  TypeOfEvent = new int[InternalEvents];

  for (int e=0; e<InternalEvents; e++) {
    Qd[e] = Qu[e] = Qb[e] = 0;
    BackSeg[e] = FrontStp[e] = -1;
    TypeOfEvent[e] = 0;
  }
  ArraysCreated = 1;
  cout << "> Arrays to store data of internal events created. Dimension = "
       << InternalEvents << "." << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Define histograms used in special cases.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::CreateHistDiffQ(int BinsQ, float MinQ, float MaxQ)
{
  HistDiffQCreated = 1;
  HistDiffQ = new TH2F**[4];
  for (int stp=0; stp<4; stp++) {
    HistDiffQ[stp] = new TH2F*[4];
    for (int bch = 0; bch<4; bch++) {
      HistDiffQ[stp][bch] = new TH2F(Form("%s.HistDiffQ_stp%d_bch%d",Name.c_str(),stp,bch),
				     Form("%s | Stp=%d, Bch=%d",Name.c_str(),stp,bch),BinsQ,MinQ,MaxQ,120,-30,30);
      HistDiffQ[stp][bch]->GetXaxis()->SetTitle(Form("Q_{%c} [arb. units]",GetEnergyFrom));
      HistDiffQ[stp][bch]->GetXaxis()->CenterTitle();
      HistDiffQ[stp][bch]->GetYaxis()->SetTitle("100%(Q_{f}-Q_{b})/Q_{b} [%]");
      HistDiffQ[stp][bch]->GetYaxis()->CenterTitle();
    }
  }  
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Define histograms used in special cases.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::CreateHistPvE(int BinsE, float MinE, float MaxE)
{
  HistPvECreated = 1;
  float ResolutionP = 0.05;    // cm
  int bins_p = (int)((DetLength+2)/ResolutionP);
  HistP = new TH1F**[4];
  for (int stp=0; stp<4; stp++) {
    HistP[stp] = new TH1F*[4];
    for (int bch = 0; bch<4; bch++) {
      HistP[stp][bch] = new TH1F(Form("%s.HistP_stp%d_bch%d",Name.c_str(),stp,bch),
				 Form("Position (stp=%d, bch=%d)",stp,bch),
				 bins_p,-1,DetLength+1);
      HistP[stp][bch]->GetXaxis()->SetTitle("Hit position [cm]");
      HistP[stp][bch]->GetXaxis()->CenterTitle();
    }
  }
  HistPvE = new TH2F**[4];
  for (int stp=0; stp<4; stp++) {
    HistPvE[stp] = new TH2F*[4];
    for (int bch = 0; bch<4; bch++) {
      HistPvE[stp][bch] = new TH2F(Form("%s.HistPvE_stp%d_bch%d",Name.c_str(),stp,bch),
				   Form("Position v. charge (stp=%d, bch=%d)",stp,bch),
				   bins_p,-1,DetLength+1, BinsE,MinE,MaxE);
      HistPvE[stp][bch]->GetXaxis()->SetTitle("Hit position [cm]");
      HistPvE[stp][bch]->GetXaxis()->CenterTitle();
      HistPvE[stp][bch]->GetYaxis()->SetTitle("Energy [MeV]");
      HistPvE[stp][bch]->GetYaxis()->CenterTitle();
    }
  }
  HistPEdge = new TH1F**[4];
  for (int stp=0; stp<4; stp++) {
    HistPEdge[stp] = new TH1F*[4];
    for(int bch = 0; bch<4; bch++){
      HistPEdge[stp][bch] = new TH1F(Form("%s.HistPEdge_stp%d_bch%d",Name.c_str(),stp,bch),
				     Form("Position (stp=%d, bch=%d)",stp,bch), bins_p,-1,DetLength+1);
    }
  }
  HistPvEEdge = new TH2F**[4];
  for (int stp=0; stp<4; stp++) {
    HistPvEEdge[stp] = new TH2F*[4];
    for (int bch = 0; bch<4; bch++) {
      HistPvEEdge[stp][bch] = new TH2F(Form("%s.HistPvEEdge_stp%d_bch%d",Name.c_str(),stp,bch),
				       Form("Position v. charge (stp=%d, bch=%d)",stp,bch),
				       bins_p,-1,DetLength+1, BinsE,MinE,MaxE);
    }
  }  
}


////////////////////////////////////////////////////////////////////////////
//  Define histograms used in special cases.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::CreateHistRFTime(int BinsT, float MinT, float MaxT)
{
  HistRFTimeCreated = 1;
  HistRFTime = new TH1F*[4];
  for (int bch = 0; bch<4; bch++) {
    HistRFTime[bch] = new TH1F(Form("%s.HistRFTime_bch%d",Name.c_str(),bch),
				    Form("%s | Bch=%d",Name.c_str(),bch), BinsT, MinT, MaxT);
    HistRFTime[bch]->GetXaxis()->SetTitle("RF Time [arb. units]");
    HistRFTime[bch]->GetXaxis()->CenterTitle();
  }
  return;
}




////////////////////////////////////////////////////////////////////////////
void SX3_Detector::DeleteArrays()
{
  delete[] Qd;
  delete[] Qu;
  delete[] Qb;
  delete[] BackSeg;
  delete[] FrontStp;
  delete[] TypeOfEvent;
  ArraysCreated = 0;
  return;
}


////////////////////////////////////////////////////////////////////////////
void SX3_Detector::DeleteHistDiffQ()
{
  HistDiffQCreated = 0;
  for (int stp=0; stp<4; stp++) {
    for (int bch = 0; bch<4; bch++) 
      delete HistDiffQ[stp][bch];
    delete[] HistDiffQ[stp];
  }
  delete[] HistDiffQ;
  return;
}


////////////////////////////////////////////////////////////////////////////
void SX3_Detector::DeleteHistPvE()
{
  HistPvECreated = 0;
  for (int stp=0; stp<4; stp++) {
    for (int bch = 0; bch<4; bch++) {
      delete HistP[stp][bch];
      delete HistPvE[stp][bch];
      delete HistPEdge[stp][bch];
      delete HistPvEEdge[stp][bch];	
    }
    delete[] HistP[stp];
    delete[] HistPvE[stp];
    delete[] HistPEdge[stp];
    delete[] HistPvEEdge[stp];
  }
  delete[] HistP;
  delete[] HistPvE;
  delete[] HistPEdge;
  delete[] HistPvEEdge;
  return;
}


////////////////////////////////////////////////////////////////////////////
void SX3_Detector::DeleteHistRFTime()
{
  HistRFTimeCreated = 0;
  for (int bch = 0; bch<4; bch++) 
    delete HistRFTime[bch];
  delete[] HistRFTime;
  return;
}



////////////////////////////////////////////////////////////////////////////
// This method fills the PvE histograms in two different ways and for all
// or just one front strip. When invoked, the method fills the histogrmas
// by looping over the data arrays of this detector. It assumes these arrays
// have been created and filled.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::FillPvE(int strip)
{
  // Reset the histograms.
  for (int stp=0; stp<4; stp++) {
    for (int bch=0; bch<4; bch++) {
      HistP[stp][bch]->Reset();
      HistPEdge[stp][bch]->Reset();
      HistPvE[stp][bch]->Reset();
      HistPvEEdge[stp][bch]->Reset();
    }
  }

  // Loop over all the internal events to fill the HistEvPs.
  for (int n=0; n<InternalEvents; n++) {
    // Get the data from the internal event.
    EventBch = BackSeg[n];
    EventStp = FrontStp[n];
    EventQb = Qb[n];
    EventQd = Qd[n];
    EventQu = Qu[n];
    EventType = TypeOfEvent[n];
    if (EventStp==strip) {
      GetEventEP();
      // Filling position histograms.
      if (EventE>0) {
	if (EventType==1) {
	  HistP[EventStp][EventBch]->Fill(EventP);
	  HistPvE[EventStp][EventBch]->Fill(EventP, EventE);
	} else if (EventType==2 || EventType==3) {
	  HistPEdge[EventStp][EventBch]->Fill(EventP);
	  HistPvEEdge[EventStp][EventBch]->Fill(EventP, EventE);
	}
      }      
    }
  }  

  return;
}


////////////////////////////////////////////////////////////////////////////
// This method fills the DiffQ histograms with data of the current event.
// It is assumed that the GetEventEP method has been invoked and that the
// histograms have been created. 
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::FillHistDiffQ()
{
  int stp = EventStp;
  int bch = EventBch;
  float Qb_cal = GetQbCal();
  float Qf_cal = GetQdCal() + GetQuCal();
  if (EventType==1) 
    HistDiffQ[stp][bch]->Fill((EventE-Eshift)/Elin, 100*(Qf_cal-Qb_cal)/Qb_cal);
  return;
}



////////////////////////////////////////////////////////////////////////////
// This method fills the PvE histograms with the data of the current event.
// It is assumed that the GetEventEP method has been invoked before and 
// that the histograms have been created. 
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::FillHistPvE()
{
  float MinE = HistPvE[0][0]->GetYaxis()->GetXmin();
  float MaxE = HistPvE[0][0]->GetYaxis()->GetXmax();
  int stp = EventStp;
  int bch = EventBch;
  if(EventE>MinE && EventE<MaxE) {
    if (EventType==1) {
      HistP[stp][bch]->Fill(EventP);
      HistPvE[stp][bch]->Fill(EventP, EventE);
    } else if (EventType==2 || EventType==3) {
      HistPEdge[stp][bch]->Fill(EventP);
      HistPvEEdge[stp][bch]->Fill(EventP, EventE);
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Gets the energy calibration coefficients, this is, a linear coeff. and
//  offset in MeV/ch and Mev, respectivelly.
//
//  IMPORTANT: It is assumed that the HistQb[4] histogram has been filled
//             with charge calibrated data, in other words it assumes the
//             relative coefficients have been obtained and that the raw 
//             data was sorted and calibrated with such coefficients.
//
//  Arguments:
//  * PeakSigma.- Related to the peak width.  A relatively small value of 
//      sigma (~2) will give a more accurate position of the peak.
//  * MinPeakAmplitude.- The threshold peak amplitud percent to be used in
//      the Search() method.
//  * Th228Peaks.- array with the emission energy of the alphas from the
//      228Th decay chain in increasing energy.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetEnergyCoeffs(int npeaks, float PeakSigma, float MinPeakAmplitude,
				   float* Th228Peaks, float MinQ, float MaxQ)
{
  cout << "> Getting energy calibration coefficients from a histogram\n"
       << "> with the calibrated back charge data." << endl;
 
  int MaxCounts;                 // Maximum number of counts displayed in the histograms.
  float E_max = 10 ;             // Range in MeV where the data will be fitted.
  string aux;
  int nfound;
  float Temp;
  TCanvas* Display;
  TF1* linear;
  TGraph* FitGraph;

  Display = new TCanvas(Form("D.%s_EneCoeffs",Name.c_str()),"Energy coefficients",0,0,1400,700);
  Display->Divide(2);
  
  // The position of the peaks are found for the rest of the back channels, then 
  // a graph is made with the reference data in the y-axis and the data from the
  // other channels in the x-axis. A straight line is fitted to each of these 
  // graphs from which one gets the relative coefficients.

  MaxCounts = 1.25*(HistQb[4]->GetMaximum());      
  HistQb[4]->SetMaximum(MaxCounts);
  HistQb[4]->SetAxisRange(MinQ, MaxQ);
  
  Display->cd(1);            
  TSpectrum *s = new TSpectrum();
  nfound = s->Search(HistQb[4],PeakSigma," ",MinPeakAmplitude);
  if(nfound != npeaks)
    cout << "> WARNING: number of peaks found (" << nfound 
	 << ") is not the same as the number of expected peaks (" << npeaks << ")." << endl;
  
  cout << "> Showing spectrum from all back channels with peaks found." << endl;
            
  float *xpeaks = s->GetPositionX();
  Temp=0;
  //Reorder the peaks' x-coordinate from lower to higher x-value.
  for(int i=0;i<nfound;i++) {
    for(int j=i;j<nfound;j++) {
      if (xpeaks[j] < xpeaks[i]) {
	Temp = xpeaks[i];
	xpeaks[i] = xpeaks[j];
	xpeaks[j] = Temp;  
      }
    }
  }
   
  FitGraph = new TGraph(nfound,xpeaks,Th228Peaks);
  FitGraph->GetXaxis()->SetTitle("Peak position [arb. units]");
  FitGraph->GetXaxis()->SetLimits(0,MaxQ);
  FitGraph->GetYaxis()->SetTitle("Alpha energy [MeV]");
  FitGraph->GetYaxis()->SetLimits(0,E_max);

  linear = new TF1("linear","pol1",0,MaxQ);
  
  FitGraph->Fit(linear,"qr");
  float chi2 = linear->GetChisquare();
  FitGraph->SetTitle(Form("Detector %s |  chi^2=%f",Name.c_str(), chi2));
  
  Elin = linear->GetParameter(1);
  Eshift = linear->GetParameter(0);
  
  Display->cd(2);
  FitGraph->Draw("AP*");
  Display->Update();
 
  cout << "> Extracted coefficients:" << endl;
  cout << ">\tElin = " << Elin << endl;
  cout << ">\tEshift = " << Eshift << endl;
  
  cout << "> Double-click on the canvas to continue." << endl;
  Display->WaitPrimitive();  // Wait for double click or space bar on the canvas.
  delete Display;

  cout << "> All done. " << endl;
 
}


////////////////////////////////////////////////////////////////////////////
// This method gets the energy (in MeV) and z-position (in cm) of an event
// that took place in this detector (an internal event). Before calling this
// method GetEventQs() must be called. This method must be invoked before
// passing the detector's physical quantities (E and P) to the ANASEN
// physical event.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetEventEP()
{
  float Qb_cal, Qtot, Qd_cal, Qu_cal, PosUS, PosDS;
  int bch,stp;
  float DSEdge, USEdge;
  
  EventE = 0;
  EventP = -100;
  
  if (EventType>0) {
    stp = EventStp;
    bch = EventBch;
    Qb_cal = mB[bch]*EventQb;
    Qd_cal = mD[stp]*(EventQd + oF[stp]/(mD[stp]+mU[stp]));
    Qu_cal = mU[stp]*(EventQu + oF[stp]/(mD[stp]+mU[stp]));

    if (GetEnergyFrom=='b')
      Qtot = Qb_cal;
    else
      Qtot = Qd_cal + Qu_cal;

    // The event energy can be obtained from the back or front sides of the detector.  
    EventE = Elin*Qtot + Eshift; // in MeV.
 
    // Calculate the acutal z-position in cm.
    if (Qtot>0) {
      PosUS = 0;
      PosDS = 0;
      for (int order=0; order<=PosPolyOrder; order++) {
	PosDS += Cd[stp][order]*pow(Qd_cal/Qtot, order);
	PosUS += Cu[stp][order]*pow(Qu_cal/Qtot, order);
      }
      if (EventType==1) 
	EventP = (PosDS + PosUS)/2;
      else if (EventType==2) 
	EventP = PosDS;
      else if (EventType==3) 
	EventP = PosUS;

      if (ForceBchPosRange) {
	DSEdge = bch*DetLength/4;
	USEdge = (bch+1)*DetLength/4;
	if (EventP<DSEdge || EventP>USEdge) {
	  EventE = 0;
	  EventP = -100;
	  EventType = 0;
	}
      }
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetEventQs(string HitType)
{
  Int_t B[12] = {0x001,0x002,0x004,0x008,    //here we are using hex for the hit specific hit key
		 0x010,0x020,0x040,0x080,
		 0x100,0x200,0x400,0x800};
  Int_t HitKey, BackKey, DSKey, USKey;
  int bch, stp, chan, stp_DS, stp_US;
  float Qb,Qd,Qu;

  EventType = 0;
  Qb = Qd = Qu = 0;
  stp = -1;
  bch = -1; 

  //- Step 1 --------------------------------------------------------------------------------
  // Loop over the multiplicity index and sort the RawE into Qb, Qd and Qu. Also fill the 
  // HitKey depending on the combination of channels that fired.
  HitKey = 0x0;
  for (int mi=0; mi<mult; mi++) {
    chan = Channel[mi];
    // If it's a back channel make sure it has a good time.
    if (chan<4) {
      if (Time[mi]>4000 && Time[mi]<5000)
	HitKey = HitKey | B[chan];
    } else {
      if (Time[mi]>1000 && Time[mi]<8000)
	HitKey = HitKey | B[chan];
    }        
    // At this point the channel number is already > 0, so no need to check that again.
    if (chan<4) {
      if (Time[mi]>4000 && Time[mi]<5000) {
	Qb = RawE[mi] - q0[chan];
	EventT = Time[mi];
      }
    } else if (chan<8) 
      Qd = RawE[mi] - q0[chan];
    else if (chan<12) 
      Qu = RawE[mi] - q0[chan];
  }// enf for(mi)


  //- Step 2 --------------------------------------------------------------------------------
  // Read the HitKey to assign 'EventBch' and 'EventStp'.
  BackKey = HitKey & 0x00F;  
  if (HitType=="normal" || HitType=="any") {
    switch(BackKey) {
      // Normal hits
    case 0x1: bch = 0; break;
    case 0x2: bch = 1; break;
    case 0x4: bch = 2; break;
    case 0x8: bch = 3; break;
    }
  } 
  if(HitType=="border" || HitType=="any"){
    switch(BackKey) {
      // Border hits
    case 0x3: bch = 0; break;
    case 0x6: bch = 1; break;
    case 0xC: bch = 2; break;    
    }
  }
  // The 'bch' must have been changed above form its initialization value (-1) in order to
  // continue.
  if (bch>-1) {
    // Get the strip number from the down-stream signal.
    DSKey = HitKey & 0x0F0;
    stp_DS = -1;
    for (int i=0; i<4; i++)
      if (DSKey==B[i+4]) {
	stp_DS = i;
	break;
      }
    // Get the strip number from the up-stream signal.
    USKey = HitKey & 0xF00;
    stp_US = -2;
    for (int i=0; i<4; i++)
      if (USKey==B[i+8]) {
	stp_US = i;
	break;
      }

    //- Step 3 --------------------------------------------------------------------------------
    // Assign the type of event index and save the 'Event' variables. For a full hit we need 1
    // back signal and 2 front ones from the same front strip but events with only one front
    // signal are also taken into account. 

    // This is a full hit (1 back + 2 front).
    if (Qu>0 && Qd>0 && stp_US==stp_DS) {
      EventStp = stp_US;
      EventType = 1;
    }
    // This is a half front hit on the down-stream side.
    else if (Qu==0 && Qd>0 && stp_DS>=0) {
      if (!ForceFullHits) {
	EventStp = stp_DS;
	EventType = 2;
      }
    }
    // This is a half front hit on the up-stream side.
    else if (Qu>0 && Qd==0 && stp_US>=0) {
      if (!ForceFullHits) {
	EventStp = stp_US;
	EventType = 3;
      }
    }
    EventBch = bch;
    EventQd = Qd;
    EventQu = Qu;
    EventQb = Qb;  
  }
  
  return;
}


////////////////////////////////////////////////////////////////////////////
// New method for obtaining the front strip coefficients.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetFrontCoeffs(int strip, float MaxQDifference, int Cuts, float *EnergyCuts, float Uncertainty)
{
  cout << ">\t-< Strip " << strip << " >-" << endl;

  TCanvas *Display;
  TF2 *QfCal;
  TGraph2D *FitGraph;
  TH2F** HistFB_before;
  TH2F** HistFB_after;
  int bch, stp;
  float channels_per_bin = 16;
  float hist_max_Q = 8192;
  int i=0;
  float Qb_cal, Qf;
  int MinPoints = 200; // Minumum number of points used in fits.
 
  FitGraph = new TGraph2D();
  HistFB_before = new TH2F*[4];
  HistFB_after = new TH2F*[4];
  for (bch=0; bch<4; bch++) {
    HistFB_before[bch] = new TH2F(Form("HistFB_b_bch%d",bch),"Before calibration (all back channels)",
				  int(2000/channels_per_bin),-1000,1000, int(hist_max_Q/channels_per_bin),0,hist_max_Q);
    HistFB_before[bch]->GetXaxis()->SetTitle("Q_{f} - Q_{b}^{cal}");
    HistFB_before[bch]->GetYaxis()->SetTitle("Q_{f} [arb. units]");

    HistFB_after[bch] = new TH2F(Form("HistFB_a_bch%d",bch),"After calibration (all back channels)",
				 int(2000/channels_per_bin),-1000,1000, int(hist_max_Q/channels_per_bin),0,hist_max_Q);
    HistFB_after[bch]->GetXaxis()->SetTitle("Q_{f}^{cal} - Q_{b}^{cal} [arb. units]");
    HistFB_after[bch]->GetYaxis()->SetTitle("Q_{f}^{cal} [arb. units]");
  }
  // Variables: x = Qd, y = Qu.
  // Parameters' indices:
  // 0 - mD ~ down-stream gain
  // 1 - mU ~ up-stream gain
  // 2 - oF ~ total front offset
  QfCal = new TF2("QfCal","[0]*x + [1]*y + [2]",200.0,8000.0,200.0,8000.0);
  QfCal->SetParNames("mD","mU");
  QfCal->SetParameters(1.0, 1.0, 0.0);
  QfCal->SetParLimits(0, 0.75, 1.45);
  QfCal->SetParLimits(1, 0.75, 1.45);
   
  // Useful lines that represent the fitting region.
  TLine CenterLine(0,0,0,hist_max_Q);
  TLine MinDiffLine(-MaxQDifference,0, -MaxQDifference,hist_max_Q);
  MinDiffLine.SetLineStyle(2);
  MinDiffLine.SetLineColor(kRed);
  TLine MaxDiffLine(MaxQDifference,0, MaxQDifference,hist_max_Q);
  MaxDiffLine.SetLineStyle(2);
  MaxDiffLine.SetLineColor(kRed);
  
  // Reset the calibration coefficients.
  for (int bch=0; bch<4; bch++) {
    mD[strip] = mU[strip] = 1;
    oF[strip] = 0;
  }   
  cout << ">\tCoefficients mD, mU, and oF reset to default values." << endl;
  cout << ">\tFilling histograms ..." << endl;

  //- Step 1 --------------------------------------------------------------------------------
  // Fill graphs used for fitting and the histograms to show the data before getting the coeffs.
  for(int n=0; n<InternalEvents; n++) {
    EventStp = stp = FrontStp[n];
    EventBch = bch = BackSeg[n];
    EventQb = Qb[n];
    EventQd = Qd[n];
    EventQu = Qu[n];
    EventType = TypeOfEvent[n];   
    Qb_cal = mB[bch]*Qb[n] + oB[bch];
    Qf = Qd[n] + Qu[n];
    if (stp==strip) {
      GetEventEP();
      for (int c=0; c<Cuts; c++) {
	if (fabs(EventE-EnergyCuts[c])<Uncertainty) {
	  HistFB_before[bch]->Fill(Qb_cal - Qf, Qf);
	  if (fabs(Qb_cal - Qf)<MaxQDifference) {
	    FitGraph->SetPoint(i, Qd[n], Qu[n], Qb_cal);
	    i++;
	  }
	  break;
	}
      }      
    }
  }

  cout << ">\tPoints used in fits: "<< i << endl;

  if (i>MinPoints) {
    QfCal->SetParameters(1.0, 1.0, 0.0);
    FitGraph->Fit(QfCal,"RQM");
    mD[strip] = QfCal->GetParameter(0);
    mU[strip] = QfCal->GetParameter(1);
    oF[strip] = QfCal->GetParameter(2);

    //Extract the values of the parameters and store them for future use.
    cout << ">\tmD[" << strip << "]=" << mD[strip] << endl;
    cout << ">\tmU[" << strip << "]=" << mU[strip] << endl;
    cout << ">\toF[" << strip << "]=" << oF[strip] << endl;
  } else {
    cout << ">\tWrining: Not enough points. Skipping fit for bch=" << bch << endl;
  }
  
  
  //- Step 2 --------------------------------------------------------------------------------
  // Fill the histograms with the calibrated data.
  for(int n=0; n<InternalEvents; n++) {
    if (FrontStp[n]==strip) {
      EventStp = stp = FrontStp[n];
      EventBch = bch = BackSeg[n];
      EventQb = Qb[n];
      EventQd = Qd[n];
      EventQu = Qu[n];
      EventType = TypeOfEvent[n];
      if (stp==strip) {
	Qf = QfCal->Eval(Qd[n], Qu[n]);
	Qb_cal = mB[bch]*Qb[n] + oB[bch];
	GetEventEP();
	for (int c=0; c<Cuts; c++) {
	  if (fabs(EventE-EnergyCuts[c])<Uncertainty) {
	    HistFB_after[bch]->Fill(Qb_cal - Qf, Qf);
	    break;
	  }
	}      

      }     
    }
  }
  
  
  //- Step 3 --------------------------------------------------------------------------------
  // Show all histograms
  Display = new TCanvas("Disp_BF","Front-Back calibration",0,0,1500,800);
  Display->Divide(2);

  Display->cd(1);
  HistFB_before[0]->Draw();
  HistFB_before[0]->SetMarkerStyle(7);
  for (bch=1; bch<4; bch++) {
    HistFB_before[bch]->SetMarkerColor(2*bch);
    HistFB_before[bch]->SetMarkerStyle(7);
    HistFB_before[bch]->Draw("same");
  }
  CenterLine.Draw();
  MinDiffLine.Draw();
  MaxDiffLine.Draw();

  Display->cd(2);
  HistFB_after[0]->Draw();
  HistFB_after[0]->SetMarkerStyle(7);
  for (bch=1; bch<4; bch++) {
    HistFB_after[bch]->SetMarkerColor(2*bch);
    HistFB_after[bch]->SetMarkerStyle(7);
    HistFB_after[bch]->Draw("same");
  }
  CenterLine.Draw();
  MinDiffLine.Draw();
  MaxDiffLine.Draw();

  Display->Update();
  Display->WaitPrimitive();

  
  //- Step 5 --------------------------------------------------------------------------------
  // Delete all objects created with 'new'.
  for (bch=0; bch<4; bch++) {
    delete HistFB_before[bch];
    delete HistFB_after[bch];
  }
  delete[] HistFB_after;
  delete[] HistFB_before;
  delete FitGraph;
  delete QfCal;
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////
//   
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetPositionPoly(int strip, float MinE, float MaxE)
{
  TCanvas* Display;
  TF1* PolyDS;
  TF1* PolyUS;
  TGraph* GraphUS;
  TGraph* GraphDS;
  TGraph* FitGraphDS;
  TGraph* FitGraphUS;
  TH1F* HistDSRatio;
  TH1F* HistUSRatio;
  TH2F* BackgroundDS;
  TH2F* BackgroundUS;
  int brd, stp;
  int points = 0, peaks_found;
  float Qtot, pos, Qd_cal, Qu_cal, Temp;
  float p[3] = {DetLength/4, DetLength/2, 3*DetLength/4};
  TSpectrum *s = new TSpectrum();

  BackgroundDS = new TH2F("BackgroundDS","",200,0,1,75,0,DetLength);
  BackgroundDS->GetYaxis()->SetTitle("Hit position [cm]");
  BackgroundDS->SetStats(0);
  BackgroundUS = new TH2F("BackgroundUS","",200,0,1,75,0,DetLength);
  BackgroundUS->GetYaxis()->SetTitle("Hit position [cm]");
  BackgroundUS->SetStats(0);
  PolyDS = new TF1("PolyDS",Form("pol%d",PosPolyOrder),0.2,0.8);
  PolyUS = new TF1("PolyUS",Form("pol%d",PosPolyOrder),0.2,0.8);
  GraphUS = new TGraph(); 
  GraphDS = new TGraph(); 
  HistDSRatio = new TH1F("HistDSRatio","",100,0,1);
  HistDSRatio->GetXaxis()->SetTitle("Q^{cal}_{d}/Q_{tot}");
  HistDSRatio->GetXaxis()->CenterTitle();
  HistUSRatio = new TH1F("HistUSRatio","",100,0,1);
  HistUSRatio->GetXaxis()->SetTitle("Q^{cal}_{u}/Q_{tot}");
  HistUSRatio->GetXaxis()->CenterTitle();


  // Select to get the energy from the front side since the border hits do not have
  // a good back signal.
  GetEnergyFrom = 'f';

  cout << "> Finding position coefficients for strip " << strip << " ..." << endl;

  // Loop over the internal events to fill the graphs and histograms.
  for (int n=0; n<InternalEvents; n++) {
    EventBch = brd = BackSeg[n];
    EventStp = stp = FrontStp[n];
    EventQb = Qb[n];
    EventQd = Qd[n];
    EventQu = Qu[n];
    EventType = TypeOfEvent[n];
    if (stp==strip && TypeOfEvent[n]>0) {
      // Qu[] and Qd[] are the pulser corrected signals, this is, RawE - q0;
      Qd_cal = GetQdCal();
      Qu_cal = GetQuCal();
      Qtot = Qd_cal + Qu_cal;
      // We add to the border position a random number between -0.5 and 0.5 mm to account
      // for the detector uncertainty (~ 1mm).
      pos = p[brd] + rand->Uniform(-0.05,0.05);
      HistDSRatio->Fill(Qd_cal/Qtot);
      HistUSRatio->Fill(Qu_cal/Qtot);
      if (Qd_cal/Qtot>0.15 && Qd_cal/Qtot<0.9 && Qu_cal/Qtot>0.15 && Qu_cal/Qtot<0.9) {
	GraphDS->SetPoint(points, Qd_cal/Qtot, pos);
	GraphUS->SetPoint(points, Qu_cal/Qtot, pos);
	points++;
      }
    }
  }

  // IMPORTANT NOTE:
  // I've seen that to decrease the error in the fit it is better to find the peaks
  // that represent the borders and to fit a 3 point graph to the polynomial. This 
  // method gives better results than just fitting the GraphDS or GraphUS to the 
  // position polynomials.

  // Fitting the down-stream graph.
  HistDSRatio->SetAxisRange(0.2,0.8);
  peaks_found = s->Search(HistDSRatio,PeakSigma," ",MinPeakAmplitude);
  if(peaks_found != NPeaks)
    cout << "> WARNING: number of peaks found (" << peaks_found
	 << ") in down-stream spectrum is not " << NPeaks << "." << endl;
            
  float *peaks_ds = s->GetPositionX();

  //Reorder the peaks' x-coordinate (Qd/Qtot) from higher to lower x-value.
  for(int i=0;i<peaks_found;i++) {
    for(int j=i;j<peaks_found;j++) {
      if (peaks_ds[j] > peaks_ds[i]) {
	Temp = peaks_ds[i];
	peaks_ds[i] = peaks_ds[j];
	peaks_ds[j] = Temp;  
      }
    }
  }
  FitGraphDS = new TGraph(peaks_found, peaks_ds, p);
  FitGraphDS->Fit(PolyDS,"QR");
  FitGraphDS->SetMarkerColor(kBlue);

  if(peaks_found == NPeaks)
    for (int order=0; order<=PosPolyOrder; order++)
      Cd[strip][order] = PolyDS->GetParameter(order);
  else {
    Cd[strip][0] = -100;
    for (int order=1; order<=PosPolyOrder; order++)
      Cd[strip][order] = 0;
  }
  
  // Fitting the up-stream graph.
  HistUSRatio->SetAxisRange(0.2,0.8);
  peaks_found = s->Search(HistUSRatio,PeakSigma," ",MinPeakAmplitude);
  if(peaks_found != NPeaks)
    cout << "> WARNING: number of peaks found (" << peaks_found
	 << ") in up-stream spectrum is not " << NPeaks << "." << endl;
            
  float *peaks_us = s->GetPositionX();

  //Reorder the peaks' x-coordinate (Qu/Qtot) from lower to higher x-value.
  for(int i=0;i<peaks_found;i++) {
    for(int j=i;j<peaks_found;j++) {
      if (peaks_us[j] < peaks_us[i]) {
	Temp = peaks_us[i];
	peaks_us[i] = peaks_us[j];
	peaks_us[j] = Temp;  
      }
    }
  }
  FitGraphUS = new TGraph(peaks_found, peaks_us, p);
  FitGraphUS->Fit(PolyUS,"QR");
  FitGraphUS->SetMarkerColor(kBlue);

  if(peaks_found == NPeaks)
    for (int order=0; order<=PosPolyOrder; order++)
      Cu[strip][order] = PolyUS->GetParameter(order);
  else {
    Cu[strip][0] = -100;
    for (int order=1; order<=PosPolyOrder; order++)
      Cu[strip][order] = 0;
  }


  // Show the graphs and fits.   
  Display = new TCanvas(Form("D.%s.PosCoeffs",Name.c_str()),
			Form("%s strip %d coefficients",Name.c_str(),strip),0,0,1500,700);
  Display->Divide(2,2);

  Display->cd(1);
  BackgroundDS->GetXaxis()->SetTitle("Q^{cal}_{d}/Q_{tot}");
  BackgroundDS->Draw();
  GraphDS->Draw("p");
  FitGraphDS->Draw("* same");
  PolyDS->SetRange(0.0, 1.0);
  PolyDS->Draw("l same");

  Display->cd(2);
  BackgroundUS->GetXaxis()->SetTitle("Q^{cal}_{u}/Q_{tot}");
  BackgroundUS->Draw();
  GraphUS->Draw("p");
  FitGraphUS->Draw("* same");
  PolyUS->SetRange(0.0, 1.0);
  PolyUS->Draw("l same");

  Display->cd(3);
  HistDSRatio->SetAxisRange(0.0,1.0);
  HistDSRatio->Draw();

  Display->cd(4);
  HistUSRatio->SetAxisRange(0.0,1.0);
  HistUSRatio->Draw();

  // Print the obtaied coefficients.
  cout << "> Parameters found:\n>\tOrder" << setw(6) <<"Cd" << setw(13) <<"Cu" << endl;
  for (int order=0; order<=PosPolyOrder; order++) 
    cout << ">\t" << order << setw(13) << Cd[strip][order] << setw(13) << Cu[strip][order] << endl;


  Display->Update();
  Display->WaitPrimitive();

  delete Display;
  delete s;
  delete HistDSRatio;
  delete HistUSRatio;
  delete BackgroundUS;
  delete BackgroundDS;
  delete PolyDS;
  delete PolyUS;
  delete FitGraphUS;
  delete FitGraphDS;
  delete GraphUS;
  delete GraphDS;
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
void SX3_Detector::GetPulserOffsets(int NPeaks, float PeakSigma, float MinPeakAmplitude,
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

  cout << "> Offsets found:\n"
       << "> ch\tq0" << endl;

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
    for(Int_t i=0;i<PeaksFound;i++) {
      for(Int_t j=i;j<PeaksFound;j++) {
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
    FitGraph[ch]->GetXaxis()->SetLimits(0,16000);
    FitGraph[ch]->GetYaxis()->SetTitle("Pulser voltage [V]");
    FitGraph[ch]->GetYaxis()->SetLimits(0,10);
    
    linear = new TF1("linear","pol1",100,16000);
    FitGraph[ch]->Fit(linear,"qrROB");
    chi2[ch] = linear->GetChisquare();
    slope[ch] = linear->GetParameter(1);
    offset[ch] = linear->GetParameter(0);
    //cout << offset[ch] << "\t" << slope[ch] << "\t" << chi2[ch] << "\t" << endl;
    
    FitGraph[ch]->SetTitle(Form("Ch=%d  |  Chi^2=%f", ch, chi2[ch]));
    
    Display->cd(2);
    FitGraph[ch]->Draw("AP*");
    Display->Update();
    if(chi2[ch]<NormalChi2Value){
      q0[ch] = -offset[ch]/slope[ch];
      cout << "> " << ch << "\t" << q0[ch] << endl;
    }

    else 
      cout << "> Warning: ch=" << ch << " has a high chi^2 value " << chi2[ch] << endl;

    Display->WaitPrimitive();

  }// end for(ch)
  delete Display;
  return; 
}


////////////////////////////////////////////////////////////////////////////
// The followin three functions simply return the calibrated back, 
// down-stream and up-stream 'charges', respectively.
////////////////////////////////////////////////////////////////////////////
float SX3_Detector::GetQbCal()
{
  return (mB[EventBch]*EventQb + oB[EventBch]);
}

float SX3_Detector::GetQdCal()
{
  return (mD[EventStp]*(EventQd + oF[EventStp]/(mD[EventStp]+mU[EventStp])));
}

float SX3_Detector::GetQuCal()
{
  return (mU[EventStp]*(EventQu + oF[EventStp]/(mD[EventStp]+mU[EventStp])));
}



////////////////////////////////////////////////////////////////////////////
// Gets the relative coefficients with respect to one reference back
// channel for the rest of the back channels.
// * PeakSimga.- Related to the peak width. A relatively small value of
//   sigma (~2) will give a more accurate position of the peak.
// * MinPeakAmplitude.- The threshold peak amplitud percent to be used in
//   the Search() method.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::GetRelativeCoeffs(int back_ref, float PeakSigma, float MinPeakAmplitude,
				     float MinQ, float MaxQ)
{
  cout << "> Getting relative coefficients for back channels. Reference channel = " 
       << back_ref << "." << endl;
  BackRef = back_ref;
 
  int MaxCounts;         // Maximum number of counts displayed in the histograms.
  string aux;
  TCanvas* Display;
  TSpectrum *ref, *s=0;
  TF1* linear = 0;
  TGraph** FitGraph;
  FitGraph = new TGraph*[4];
  
  Display = new TCanvas("Display","Relative Coefficients",0,0,1400,700);
  Display->Divide(2);

  // The coefficients of the reference channel are by definition 1 for the 
  // linear coefficient and 0 for the offset.
  mB[back_ref]=1;
  oB[back_ref]=0;

  //- Step 1 --------------------------------------------------------------------------------
  // First we need to get the position of the peaks in the histogram of the
  // reference channel.
  ref = new TSpectrum();
  Display->cd(1);
  MaxCounts =  (int)(1.25*(HistQb[back_ref]->GetMaximum()));
  HistQb[back_ref]->SetMaximum(MaxCounts);
  HistQb[back_ref]->SetAxisRange(MinQ, MaxQ);
  int nfound = ref->Search(HistQb[back_ref],PeakSigma," ",MinPeakAmplitude);
  Display->Update();

  cout << "> Showing reference spectrum with " << nfound << " peaks found." << endl;

  float *x_ref = ref->GetPositionX();
  float Temp=0;
  // Reorder the peaks' x-coordinate from lower to higher x-value.
  for(int i=0;i<nfound;i++) {
    for(int j=i;j<nfound;j++) {
      if (x_ref[j] < x_ref[i]) {
	Temp = x_ref[i];
	x_ref[i] = x_ref[j];
	x_ref[j] = Temp;  
      }
    }
  }
  Display->WaitPrimitive();

  //- Step 2 --------------------------------------------------------------------------------
  // The position of the peaks are found for the rest of the back channels, then 
  // a graph is made with the reference data in the y-axis and the data from the
  // other channels in the x-axis. A straight line is fitted to each of these 
  // graphs from which one gets the relative coefficients.
  for(int bch=0; bch<4; bch++){
    if(bch!=back_ref){
      
      HistQb[bch]->SetMaximum(MaxCounts);
      HistQb[bch]->SetAxisRange(MinQ, MaxQ);

      Display->cd(1);            
      s = new TSpectrum();
      nfound = s->Search(HistQb[bch],PeakSigma," ",MinPeakAmplitude);

      cout << "> Showing bch=" << bch << " spectrum with " << nfound << " peaks found." << endl;
            
      float *xpeaks = s->GetPositionX();
      Temp=0;
      //Reorder the peaks' x-coordinate from lower to higher x-value.
      for(int i=0;i<nfound;i++) {
	for(int j=i;j<nfound;j++) {
	  if (xpeaks[j] < xpeaks[i]) {
	    Temp = xpeaks[i];
	    xpeaks[i] = xpeaks[j];
	    xpeaks[j] = Temp;  
	  }
	}
      }
      	  
      FitGraph[bch] = new TGraph(nfound,xpeaks,x_ref);
      FitGraph[bch]->GetXaxis()->SetTitle("Peak position [arb. units]");
      FitGraph[bch]->GetXaxis()->SetLimits(0,MaxQ);
      FitGraph[bch]->GetYaxis()->SetTitle("Reference peak position [arb. units]");
      FitGraph[bch]->GetYaxis()->SetLimits(0,MaxQ);
      
      linear = new TF1("linear","[0]*x",0,MaxQ);
      
      FitGraph[bch]->Fit(linear,"rq");
      float chi2 = linear->GetChisquare();
      FitGraph[bch]->SetTitle(Form("det=%s  ch=%d |  chi^2=%f",Name.c_str(), bch, chi2));
      mB[bch] = linear->GetParameter(0);
      
      Display->cd(2);
      FitGraph[bch]->Draw("AP*");
      Display->Update();
      Display->WaitPrimitive();
      
    }// end if(bch!=back_ref)
  }// end for(bch)

  cout << "> Extracted coefficients:" << endl;
  for (int bch=0; bch<4; bch++)
    cout << ">\tmB[" << bch << "]=" << mB[bch] << endl;
  

  delete ref;
  delete s;
  delete linear;
  for (int bch=0; bch<4; bch++) 
    if (bch!=BackRef)
      delete FitGraph[bch];
  delete[] FitGraph;
  delete Display;
  return;
}



////////////////////////////////////////////////////////////////////////////
//  Loads the calibration coefficients prevously stored in a file.
//  It retunrs 1 if the text file was open succesfully and 0 if not.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::LoadCoefficients(string file){
  string aux;
  ifstream read_coeffs(file.c_str());
  
  if(read_coeffs.is_open()){
    // Offsets from pulser data.
    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> Threshold[i];
    for (int i=0; i<NumChan; i++) read_coeffs >> aux >> q0[i];
    // Energy calibration
    for (int i=0; i<4; i++) read_coeffs >> aux >> mB[i];
    read_coeffs >> aux >> Elin;
    read_coeffs >> aux >> Eshift;
    // Front-back charge calibration    
    for (int i=0; i<4; i++) read_coeffs >> aux >> mD[i];
    for (int i=0; i<4; i++) read_coeffs >> aux >> mU[i];
    for (int i=0; i<4; i++) read_coeffs >> aux >> oF[i];
    // Front position caibration
    read_coeffs >> aux;
    for (int stp=0; stp<4; stp++)
      for (int order=0; order<=PosPolyOrder; order++)
	read_coeffs >> Cd[stp][order];    
    read_coeffs >> aux;
    for (int stp=0; stp<4; stp++)
      for (int order=0; order<=PosPolyOrder; order++)
	read_coeffs >> Cu[stp][order];
          
    read_coeffs.close();

    IsOn = 1;  // The detector is ready to be used.
  }
  else{
    cout << "> Warning: File " << file << " couldn't be opened. Only raw histogrmas will be available." << endl;
    IsOn = 0;  // The detector can't be used.
  }
  
  return;
}


////////////////////////////////////////////////////////////////////////////
// Simple function that resets the quantities used in one event.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::Reset()
{
  mult = 0;
  for (int mi=0; mi<MaxMultiplicity; mi++) {
    Channel[mi] = -1;
    RawE[mi] = 0;
    Time[mi] = 0;
  }
  EventQb = EventQd = EventQu = 0;
  EventBch = EventStp = -1;
  EventP = -100;
  EventE = 0;
  EventT = 0;
  EventType = 0;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Simple function that saves the 'charges' of events of non-zero type.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::SaveEventQs(int event)
{
  Qb[event] = EventQb;
  Qd[event] = EventQd;
  Qu[event] = EventQu;
  BackSeg[event] = EventBch;
  FrontStp[event] = EventStp;
  TypeOfEvent[event] = EventType;
  return;
}



////////////////////////////////////////////////////////////////////////////
// The user provides the parameters related to peak finding. These are used,
// for example, when looking for the pulser offsets and when getting the 
// relative energy coefficients.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::SetPeakFindingParameters(int NPeaks, float PeakSigma, float MinPeakAmplitude, float NormalChi2Value)
{
  this->NPeaks = NPeaks;
  this->PeakSigma = PeakSigma;
  this->MinPeakAmplitude = MinPeakAmplitude;
  this->NormalChi2Value = NormalChi2Value;   // Used in linear fits.
  return;
}



////////////////////////////////////////////////////////////////////////////
// This functions helps you find the raw energy threshold, one of the most
// basic parameters of the detectors. When extracting data for this detector
// only events with RawE>Threshold[channel] will be saved in memory.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::SetThresholds(float DisplayLimit)
{
  TCanvas* Display = new TCanvas(Form("%s.Disp.FindThresh",Name.c_str()),
				 "Find thresholds",0,0,1000,800);
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
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::ShowHistDiffQ(int strip)
{
  TCanvas* Display;
  float MinQ =  HistDiffQ[0][0]->GetXaxis()->GetXmin();
  float MaxQ =  HistDiffQ[0][0]->GetXaxis()->GetXmax();
  TLine* Plus5Percent = new TLine(MinQ, 5, MaxQ, 5);
  TLine* Minus5Percent = new TLine(MinQ, -5, MaxQ, -5);
  TLine* Center = new TLine(MinQ, 0, MaxQ, 0);

  Plus5Percent->SetLineWidth(2);
  Plus5Percent->SetLineStyle(2);
  Plus5Percent->SetLineColor(kGray+2);
  Minus5Percent->SetLineWidth(2);
  Minus5Percent->SetLineStyle(2);
  Minus5Percent->SetLineColor(kGray+2);
  Center->SetLineWidth(2);
  Center->SetLineStyle(2);
  Center->SetLineColor(kBlack);

  cout << "> Showing HistDiffQ gated on back channels for " << Name << ", front strip " << strip << "." << endl;

  Display = new TCanvas(Form("D.%s.DiffQ",Name.c_str()),Form("%s | Qf - Qb",Name.c_str()),0,0,1550,800);
  Display->Divide(2,2);
  
  // Draw the percent difference histograms.
  for(int bch=0; bch<4; bch++) {
    Display->cd(bch+1);  
    HistDiffQ[strip][bch]->Draw("colz");
    Center->DrawClone();     
    Plus5Percent->DrawClone();
    Minus5Percent->DrawClone();
  }

  Display->Update();
  cout << "> Double-click on the canvas to delete it." << endl;
  Display->WaitPrimitive();
  delete Center;
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::ShowHistPvE(int strip, float MinE, float MaxE, string CanvasExtension)
{
  TCanvas* Display;
  string NameExt = "";
  string TitleExt = "";
  if (CanvasExtension!="") {
    NameExt = "_" + CanvasExtension;
    TitleExt = " (" + CanvasExtension + ")";
  }

  TLine* LLeft2  = new TLine(0, MinE, 0, MaxE);
  TLine* LLeft1  = new TLine(DetLength/4, MinE, DetLength/4, MaxE);
  TLine* LCenter = new TLine(DetLength/2, MinE, DetLength/2, MaxE);
  TLine* LRight1 = new TLine(3*DetLength/4, MinE, 3*DetLength/4, MaxE);
  TLine* LRight2 = new TLine(DetLength, MinE, DetLength, MaxE);

  LLeft2->SetLineWidth(2);   LLeft2->SetLineStyle(2);   LLeft2->SetLineColor(2);
  LLeft1->SetLineWidth(2);   LLeft1->SetLineStyle(2);   LLeft1->SetLineColor(2);
  LCenter->SetLineWidth(2);  LCenter->SetLineStyle(2);  LCenter->SetLineColor(2);
  LRight1->SetLineWidth(2);  LRight1->SetLineStyle(2);  LRight1->SetLineColor(2);
  LRight2->SetLineWidth(2);  LRight2->SetLineStyle(2);  LRight2->SetLineColor(2);

  cout << "> Showing PvE gated on back channels for strip " << strip << "." << endl;

  Display = new TCanvas(Form("D.%s.PvE%s",Name.c_str(),NameExt.c_str()),
			Form("%s position%s",Name.c_str(),TitleExt.c_str()),0,0,1550,800);
  Display->Divide(4,2);
  
  for(int bch=0; bch<4; bch++){
    Display->cd(bch+1);
    // Draw the position against charge histogram.
    HistPvE[strip][bch]->SetMarkerColor(1);
    HistPvE[strip][bch]->GetYaxis()->SetRangeUser(MinE, MaxE);
    HistPvE[strip][bch]->DrawClone();
    HistPvEEdge[strip][bch]->SetMarkerColor(4);
    HistPvEEdge[strip][bch]->GetYaxis()->SetRangeUser(MinE, MaxE);
    HistPvEEdge[strip][bch]->DrawClone("same");
    // Draw the lines corresponding to the mean position of the borders.
    LLeft2->DrawClone();
    LLeft1->DrawClone();
    LCenter->DrawClone();
    LRight1->DrawClone();
    LRight2->DrawClone();
    
  }
  for(int bch=0; bch<4; bch++){
    Display->cd(bch+5);
    HistP[strip][bch]->SetLineColor(1);
    HistP[strip][bch]->DrawClone();
    HistPEdge[strip][bch]->SetLineColor(4);
    HistPEdge[strip][bch]->DrawClone("same");
  }
  Display->Update();
  cout << "> Double-click on the canvas to delete it." << endl;
  Display->WaitPrimitive();
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::ShowGraphsQbQf(int Cuts, float *EnergyCuts, float Uncertainty)
{
  TCanvas* Display1;
  TCanvas* Display2;
  TGraph**** DiffQP;
  TGraph**** DiffQQf;
  TH2F** BackGnd1;
  TH2F** BackGnd2;
  int*** points;
  float Qf_cal,Qb_cal;
  int stp, bch;

  BackGnd1 = new TH2F*[4];
  for (stp=0; stp<4; stp++) {
    BackGnd1[stp] = new TH2F(Form("Backgnd1_stp%d",stp),Form("Front strip %d",stp),95,-1,DetLength+1,60,-30,30);
    BackGnd1[stp]->GetYaxis()->SetTitle("Q_{f}^{cal}-Q_{b}^{cal} % difference");
    BackGnd1[stp]->GetYaxis()->CenterTitle();
    BackGnd1[stp]->GetXaxis()->SetTitle("Hit position [cm]");
    BackGnd1[stp]->GetXaxis()->CenterTitle();
    BackGnd1[stp]->SetStats(0);
  }

  BackGnd2 = new TH2F*[4];
  for (stp=0; stp<4; stp++) {
    BackGnd2[stp] = new TH2F(Form("Backgnd2_stp%d",stp),Form("Front strip %d",stp),100,0,9000,60,-30,30);
    BackGnd2[stp]->GetYaxis()->SetTitle("Q_{f}^{cal}-Q_{b}^{cal} % difference");
    BackGnd2[stp]->GetYaxis()->CenterTitle();
    BackGnd2[stp]->GetXaxis()->SetTitle("Q_{f}^{cal} [arb. units]");
    BackGnd2[stp]->GetXaxis()->CenterTitle();
    BackGnd2[stp]->SetStats(0);
  }

  DiffQP = new TGraph***[Cuts];
  DiffQQf = new TGraph***[Cuts];
  points = new int**[Cuts];
  for (int c=0; c<Cuts; c++) {
    DiffQP[c] = new TGraph**[4];
    DiffQQf[c] = new TGraph**[4];
    points[c] = new int*[4];
    for (stp=0; stp<4; stp++) {
      DiffQP[c][stp] = new TGraph*[4];
      DiffQQf[c][stp] = new TGraph*[4];
      points[c][stp] = new int[4];
      for (bch=0; bch<4; bch++) {
	DiffQP[c][stp][bch] = new TGraph();
	DiffQQf[c][stp][bch] = new TGraph();
	points[c][stp][bch] = 0;
      }
    }
  }

  GetEnergyFrom = 'b';
    
  // Loop over all the internal events to fill the HistEvPs.
  for (int n=0; n<InternalEvents; n++) {

    // Get the data from the internal event.
    EventBch = bch = BackSeg[n];
    EventStp = stp = FrontStp[n];
    EventQb = Qb[n];
    EventQd = Qd[n];
    EventQu = Qu[n];
    EventType = TypeOfEvent[n];
    
    if (EventType==1) {
      Qb_cal = mB[bch]*Qb[n] + oB[bch];
      Qf_cal = mU[stp]*Qu[n] + mD[stp]*Qd[n] + oF[stp];
      GetEventEP();      
      // Fill graphs for different energy cuts.
      for (int c=0; c<Cuts; c++) {
	if (fabs(EventE-EnergyCuts[c])<Uncertainty) {
	  DiffQP[c][stp][bch]->SetPoint(points[c][stp][bch], EventP, 100.0*(Qf_cal-Qb_cal)/Qb_cal);
	  DiffQQf[c][stp][bch]->SetPoint(points[c][stp][bch], Qf_cal, 100.0*(Qf_cal-Qb_cal)/Qb_cal);
	  points[c][stp][bch]++;
	  break;
	}
      }      
    }
  }  

  Display1 = new TCanvas(Form("D.%s_DiffQP",Name.c_str()),"DiffQ v. Pos",0,0,1400,700);
  Display1->Divide(2,2);

  for (stp=0; stp<4; stp++) {
    Display1->cd(stp+1);
    BackGnd1[stp]->Draw();
    
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	DiffQP[c][stp][bch]->SetMarkerStyle(7);
	if (c==Cuts-1)
	  DiffQP[c][stp][bch]->SetMarkerColor(kBlack);
	else {
	  if (bch==0)
	    DiffQP[c][stp][bch]->SetMarkerColor(kPink+c);
	  else if (bch==1)
	    DiffQP[c][stp][bch]->SetMarkerColor(kRed+c);
	  else if (bch==2)
	    DiffQP[c][stp][bch]->SetMarkerColor(kGreen+c);
	  else if (bch==3)
	    DiffQP[c][stp][bch]->SetMarkerColor(kBlue+c);
	}
	if (points[c][stp][bch]>0)
	  DiffQP[c][stp][bch]->Draw("p same");
      }
    }
  }
  Display1->Update();
  
  Display2 = new TCanvas(Form("D.%s_DiffQQf",Name.c_str()),"DiffQ v. Qf",0,0,1400,700);
  Display2->Divide(2,2);

  for (stp=0; stp<4; stp++) {
    Display2->cd(stp+1);
    BackGnd2[stp]->Draw();
    
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	DiffQQf[c][stp][bch]->SetMarkerStyle(7);
	if (bch==0)
	  DiffQQf[c][stp][bch]->SetMarkerColor(kPink+1);
	else if (bch==1)
	  DiffQQf[c][stp][bch]->SetMarkerColor(kRed+1);
	else if (bch==2)
	  DiffQQf[c][stp][bch]->SetMarkerColor(kGreen+1);
	else if (bch==3)
	  DiffQQf[c][stp][bch]->SetMarkerColor(kBlue+1);
      
	if (points[c][stp][bch]>0)
	  DiffQQf[c][stp][bch]->Draw("p same");
      }
    }
  }
  Display2->Update();
  Display2->WaitPrimitive();
  

  for (int c=0; c<Cuts; c++) {
    for (stp=0; stp<4; stp++) {
      for (bch=0; bch<4; bch++) {
	delete DiffQP[c][stp][bch];
	delete DiffQQf[c][stp][bch];
      }
      delete[] DiffQP[c][stp];
      delete[] DiffQQf[c][stp];
      delete[] points[c][stp];
    }
    delete[] DiffQP[c];
    delete[] DiffQQf[c];
    delete[] points[c];
  }
  delete[] DiffQP;
  delete[] DiffQQf;
  delete[] points;
  for (stp=0; stp<4; stp++) {
    delete BackGnd1[stp];
    delete BackGnd2[stp];
  }
  delete[] BackGnd1;
  delete[] BackGnd2;
  delete Display1;
  delete Display2;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::ShowGraphsQdQu(int Cuts, float *EnergyCuts, float Uncertainty)
{
  TCanvas* Display;
  TGraph**** CorrQdQu;
  TGraph**** CorrQdQb;
  TGraph**** CorrQuQb;
  TGraph**** CorrQbQuasiP;
  TH2F** Background;
  TH2F** BackgroundEvP;
  int*** points;
  float Qb_cal, QuasiP, Qf_cal;
  int stp, bch;

  Background = new TH2F*[4];
  for (stp=0; stp<4; stp++) {
    Background[stp] = new TH2F(Form("Backgnd%d",stp),Form("Front strip %d",stp),500,0,4500,500,0,4500);
    Background[stp]->GetYaxis()->SetTitle("Q_{d} [arb. units]");
    Background[stp]->GetYaxis()->CenterTitle();
    Background[stp]->GetXaxis()->SetTitle("Q_{u} [arb. units]");
    Background[stp]->GetXaxis()->CenterTitle();
    Background[stp]->SetStats(0);
  }
  BackgroundEvP = new TH2F*[4];
  for (stp=0; stp<4; stp++) {
    BackgroundEvP[stp] = new TH2F(Form("BackgndEvP%d",stp),Form("Front strip %d",stp),500,0,6000,500,-0.5,0.5);
    BackgroundEvP[stp]->GetYaxis()->SetTitle("Q_{b}^{cal} [arb. units]");
    BackgroundEvP[stp]->GetYaxis()->CenterTitle();
    BackgroundEvP[stp]->GetXaxis()->SetTitle("Quasi position [arb. units]");
    BackgroundEvP[stp]->GetXaxis()->CenterTitle();
    BackgroundEvP[stp]->SetStats(0);
  }
  // Correlation between d and u.
  CorrQdQu = new TGraph***[Cuts];
  points = new int**[Cuts];
  for (int c=0; c<Cuts; c++) {
    CorrQdQu[c] = new TGraph**[4];
    points[c] = new int*[4];
    for (stp=0; stp<4; stp++) {
      CorrQdQu[c][stp] = new TGraph*[4];
      points[c][stp] = new int[4];
      for (bch=0; bch<4; bch++) {
	CorrQdQu[c][stp][bch] = new TGraph();
	points[c][stp][bch] = 0;
      }
    }
  }
  // Correlation between d and b.
  CorrQdQb = new TGraph***[Cuts];
  points = new int**[Cuts];
  for (int c=0; c<Cuts; c++) {
    CorrQdQb[c] = new TGraph**[4];
    points[c] = new int*[4];
    for (stp=0; stp<4; stp++) {
      CorrQdQb[c][stp] = new TGraph*[4];
      points[c][stp] = new int[4];
      for (bch=0; bch<4; bch++) {
	CorrQdQb[c][stp][bch] = new TGraph();
	points[c][stp][bch] = 0;
      }
    }
  }
  // Correlation between u and b.
  CorrQuQb = new TGraph***[Cuts];
  points = new int**[Cuts];
  for (int c=0; c<Cuts; c++) {
    CorrQuQb[c] = new TGraph**[4];
    points[c] = new int*[4];
    for (stp=0; stp<4; stp++) {
      CorrQuQb[c][stp] = new TGraph*[4];
      points[c][stp] = new int[4];
      for (bch=0; bch<4; bch++) {
	CorrQuQb[c][stp][bch] = new TGraph();
	points[c][stp][bch] = 0;
      }
    }
  }
  // Correlation between QuasiP and Qb_cal.
  CorrQbQuasiP = new TGraph***[Cuts];
  points = new int**[Cuts];
  for (int c=0; c<Cuts; c++) {
    CorrQbQuasiP[c] = new TGraph**[4];
    points[c] = new int*[4];
    for (stp=0; stp<4; stp++) {
      CorrQbQuasiP[c][stp] = new TGraph*[4];
      points[c][stp] = new int[4];
      for (bch=0; bch<4; bch++) {
	CorrQbQuasiP[c][stp][bch] = new TGraph();
	points[c][stp][bch] = 0;
      }
    }
  }
  
  // Loop over all the internal events to fill the HistEvPs.
  for (int n=0; n<InternalEvents; n++) {

    // Get the data from the internal event.
    EventBch = bch = BackSeg[n];
    EventStp = stp = FrontStp[n];
    EventQb = Qb[n];
    EventQd = Qd[n];
    EventQu = Qu[n];
    EventType = TypeOfEvent[n];
    Qb_cal = mB[bch]*Qb[n] + oB[bch];
    Qf_cal = mD[stp]*Qd[n] + mU[stp]*Qu[n] + oF[stp];
    
    if (EventType==1) {
      QuasiP = (Qf_cal - Qb_cal)/Qb_cal;
      GetEventEP();      
      // Fill graphs for different energy cuts.
      for (int c=0; c<Cuts; c++) {
	if (fabs(EventE-EnergyCuts[c])<Uncertainty) {
	  CorrQdQu[c][stp][bch]->SetPoint(points[c][stp][bch],Qu[n],Qd[n]);
	  CorrQdQb[c][stp][bch]->SetPoint(points[c][stp][bch],Qb_cal,Qd[n]);
	  CorrQuQb[c][stp][bch]->SetPoint(points[c][stp][bch],Qb_cal,Qu[n]);
	  CorrQbQuasiP[c][stp][bch]->SetPoint(points[c][stp][bch],Qd[n],QuasiP);
	  points[c][stp][bch]++;
	  break;
	}
      }      
    }
  }  

  // Show the Qd v. Qu graphs.
  Display = new TCanvas(Form("D.%s_CorrQdQu",Name.c_str()),"QdQu Corr",0,0,1400,700);
  Display->Divide(2,2);
  for (stp=0; stp<4; stp++) {
    Display->cd(stp+1);
    Background[stp]->Draw();
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	CorrQdQu[c][stp][bch]->SetMarkerStyle(7);
	if (bch==0)
	  CorrQdQu[c][stp][bch]->SetMarkerColor(kPink+2);
	else if (bch==1)
	  CorrQdQu[c][stp][bch]->SetMarkerColor(kRed+2);
	else if (bch==2)
	  CorrQdQu[c][stp][bch]->SetMarkerColor(kGreen+2);
	else if (bch==3)
	  CorrQdQu[c][stp][bch]->SetMarkerColor(kBlue+2);
	
	if (points[c][stp][bch]>0)
	  CorrQdQu[c][stp][bch]->Draw("p same");
      }
    }
  }
  // Show the Qd v. Qb graphs.
  Display = new TCanvas(Form("D.%s_CorrQdQb",Name.c_str()),"QdQb Corr",0,0,1400,700);
  Display->Divide(2,2);
  for (stp=0; stp<4; stp++) {
    Display->cd(stp+1);
    Background[stp]->Draw();
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	CorrQdQb[c][stp][bch]->SetMarkerStyle(7);
	if (bch==0)
	  CorrQdQb[c][stp][bch]->SetMarkerColor(kPink+2);
	else if (bch==1)
	  CorrQdQb[c][stp][bch]->SetMarkerColor(kRed+2);
	else if (bch==2)
	  CorrQdQb[c][stp][bch]->SetMarkerColor(kGreen+2);
	else if (bch==3)
	  CorrQdQb[c][stp][bch]->SetMarkerColor(kBlue+2);
	
	if (points[c][stp][bch]>0)
	  CorrQdQb[c][stp][bch]->Draw("p same");
      }
    }
  }
  // Show the Qu v. Qb graphs.
  Display = new TCanvas(Form("D.%s_CorrQuQb",Name.c_str()),"QuQb Corr",0,0,1400,700);
  Display->Divide(2,2);
  for (stp=0; stp<4; stp++) {
    Display->cd(stp+1);
    Background[stp]->Draw();
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	CorrQuQb[c][stp][bch]->SetMarkerStyle(7);
	if (bch==0)
	  CorrQuQb[c][stp][bch]->SetMarkerColor(kPink+2);
	else if (bch==1)
	  CorrQuQb[c][stp][bch]->SetMarkerColor(kRed+2);
	else if (bch==2)
	  CorrQuQb[c][stp][bch]->SetMarkerColor(kGreen+2);
	else if (bch==3)
	  CorrQuQb[c][stp][bch]->SetMarkerColor(kBlue+2);
	
	if (points[c][stp][bch]>0)
	  CorrQuQb[c][stp][bch]->Draw("p same");
      }
    }
  } 
  // Show the Qb_cal v. QuasiP graphs.
  Display = new TCanvas(Form("D.%s_CorrQbQuasiP",Name.c_str()),"QbQuasiP Corr",0,0,1400,700);
  Display->Divide(2,2);
  for (stp=0; stp<4; stp++) {
    Display->cd(stp+1);
    BackgroundEvP[stp]->Draw();
    for (bch=0; bch<4; bch++) {
      for (int c=0; c<Cuts; c++) {
	CorrQbQuasiP[c][stp][bch]->SetMarkerStyle(7);
	if (bch==0)
	  CorrQbQuasiP[c][stp][bch]->SetMarkerColor(kPink+2);
	else if (bch==1)
	  CorrQbQuasiP[c][stp][bch]->SetMarkerColor(kRed+2);
	else if (bch==2)
	  CorrQbQuasiP[c][stp][bch]->SetMarkerColor(kGreen+2);
	else if (bch==3)
	  CorrQbQuasiP[c][stp][bch]->SetMarkerColor(kBlue+2);
	
	if (points[c][stp][bch]>0)
	  CorrQbQuasiP[c][stp][bch]->Draw("p same");
      }
    }
  }
  Display->WaitPrimitive();
  
  delete Background;
  delete BackgroundEvP;
  for (int c=0; c<Cuts; c++) {
    for (stp=0; stp<4; stp++) {
      for (bch=0; bch<4; bch++) {
	delete CorrQdQu[c][stp][bch];
	delete CorrQdQb[c][stp][bch];
	delete CorrQuQb[c][stp][bch];
	delete CorrQbQuasiP[c][stp][bch];
      }
      delete[] CorrQdQu[c][stp];
      delete[] CorrQdQb[c][stp];
      delete[] CorrQuQb[c][stp];
      delete[] CorrQbQuasiP[c][stp];
      delete[] points[c][stp];
    }
    delete[] CorrQdQu[c];
    delete[] CorrQdQb[c];
    delete[] CorrQuQb[c];
    delete[] CorrQbQuasiP[c];
    delete[] points[c];
  }
  delete[] CorrQdQu;
  delete[] CorrQdQb;
  delete[] CorrQuQb;
  delete[] CorrQbQuasiP;
  delete[] points;
  delete Display;
  return;
}




////////////////////////////////////////////////////////////////////////////
//  Writes the calibration coefficients to a file. Be careful when over
//  writing older calibration files.
////////////////////////////////////////////////////////////////////////////
void SX3_Detector::WriteCoefficients(string file)
{
  string aux;
  ofstream write_coeffs(file.c_str());
  // Raw energy threshold
  for (int i=0; i<NumChan; i++) write_coeffs << "Threshold[" << i << "]= " << Threshold[i] << endl;
  // Offsets from pulser data.
  for (int i=0; i<NumChan; i++) write_coeffs << "q0[" << i << "]= " << q0[i] << endl;
  //Energy calibration
  for (int i=0; i<4; i++) write_coeffs << "mB[" << i << "]= " << mB[i] << endl;
  write_coeffs << "Elin= " << Elin << endl;
  write_coeffs << "Eshift= " << Eshift << endl;
  // Front-back charge calibration
  for (int i=0; i<4; i++) write_coeffs << "mD[" << i << "]= " << mD[i] << endl;
  for (int i=0; i<4; i++) write_coeffs << "mU[" << i << "]= " << mU[i] << endl;
  for (int i=0; i<4; i++) write_coeffs << "oF[" << i << "]= " << oF[i] << endl;
  // Front position caibration
  write_coeffs << "Cd[stp][order]= ";
  for (int stp=0; stp<4; stp++) {
    if (stp!=0) write_coeffs << "\t";
    for (int order=0; order<=PosPolyOrder; order++)
      write_coeffs << Cd[stp][order] << " ";
    write_coeffs << "\n";
  }
  write_coeffs << "Cu[stp][order]= ";
  for (int stp=0; stp<4; stp++) {
    if (stp!=0) write_coeffs << "\t";
    for (int order=0; order<=PosPolyOrder; order++)
      write_coeffs << Cu[stp][order] << " ";
    write_coeffs << "\n";
  }
  
  write_coeffs.close();
  
  cout << "> Calibration coefficients written to \"" << file << "\"." << endl;
  return;    
}
