/*******************************************************************
Code: EventProcessor.cpp

Description: Methods for the EventProcessor class. The EventProcessor 
   is meant to provide an interface for the user to calibrate and 
   monitor the detectors that form ANASEN. 
   
Author: Daniel Santiago-Gonzalez
2013-06
Otro comentario
*******************************************************************/

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
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TTree.h>

// Useful libraries
#include "EnergyLoss.cpp"

//Detetectors' libraries
#include "NTupleType.hpp"
#include "ChannelMap.cpp"
#include "SX3_Detector.cpp"
#include "PC_Detector.cpp"

//EventProcessor headers.
#include "EventProcessor.hpp"

using namespace std;


///////////////////////////////////////////////////////////////////////////////////
// Constructor. 
///////////////////////////////////////////////////////////////////////////////////
EventProcessor::EventProcessor()
{
  ASICS = new ASICS_Signal(12*30); 
  // The maximum number of signals is the total number of channels,
  // this is, (30 SX3)*(12 channels) = 360

  CAEN = new CAEN_Signal(500);

  Particle = new Physical_Event(30);
  // I don't expect to have more than 30 (the total number of SX3 detectors) hits
  // to happen in one event.

  Max3DPoints = 0;
  Save3DHits = 0;
  DataExtracted = 0;
  FillFinalEvDist = 0;
  DataFraction = 1.0;
  ApplyTime1Cuts = 0;
  BeamReady = 0;
  MassesReady = 0;
  EnergyLossReady = 0;
  RF_SX3ID = new  TH2F("RF_SX3ID","", 30*4,0,30*4,250,1000,2000);
  RF_SiT = new  TH2F("RF_SiT","", 100,4000,4800,250,1000,2000);
  
}


///////////////////////////////////////////////////////////////////////////////////
// Destructor. 
// Still need to work on the destructor of the SX3 and PC because when calling
// delete it does not release the used memory.
///////////////////////////////////////////////////////////////////////////////////
EventProcessor::~EventProcessor()
{
  cout << "Destructor starts" << endl;
  delete ASICS;
  cout << "ASICS deleted" << endl;
  delete CAEN;
  cout << "CAEN deleted" << endl;
  delete Particle;
  cout << "Particle deleted" << endl;
  delete Chain;
  cout << "Chain deleted " << Chain << endl;
#if 0
  cout << "Deleting SX3s" << endl;
  for (int d=0; d<NumSX3Dets; d++)
    delete SX3[d];  
  delete[] SX3;  
#endif
  cout << "Deleting PC wires" << endl;
  for (int w=0; w<NumPCWires; w++)
    delete PC[w];
  delete[] PC;
  cout << "Destructor ends" << endl;
}  

  
///////////////////////////////////////////////////////////////////////////////////
// File converter, from .evt to .root. Arguments:
// - files_list: used to specify the 
///////////////////////////////////////////////////////////////////////////////////
void EventProcessor::ConvertEVTtoROOT(string files_list, int BufferWords){
  
  TFile* RootFile;
  TTree* DataTree;
  int BufferBytes = BufferWords*2;
  unsigned short buffer[BufferWords];
  unsigned int BufferType = 0;
  unsigned int NBuffers = 0; 
  int BufferPhysics = 0;
  ifstream ListEVT;
  ifstream evtfile;
  string aux;
  string data_dir = "";
  string output_file;
  int run_number;    

  // Global variables
  unsigned int Nevents;
  unsigned int TotEvents;
  unsigned int words;  
  unsigned int EOB_NEvents=0;  
  unsigned int outputCounter=0;
  unsigned int CAENCounter = 0;
  unsigned short *point,*epoint;
  string output_dir;
  string output_filename;
  string evt_filename;
  int BufferScaler;
  int nWords;
  int file_subindex;
  bool go_to_next_run;

  // Open the text file with the list of evt files.  If this file cannot
  // be found end the function.
  ListEVT.open(files_list.c_str());
  if (!ListEVT.is_open()) {
    cout << "*** Error: could not open " << files_list << endl;
    return;
  }
  else{
    cout << "> File " << files_list << " opened." <<endl;
    ListEVT >> aux >> aux >> output_file;
    ListEVT >> aux >> aux >> aux >> data_dir;
  }

  // ROOT output file
  RootFile = new TFile(output_file.c_str(),"RECREATE");  

  // Data Tree
  DataTree = new TTree("DataTree","DataTree");

  // Branches for ASICs' data
  DataTree->Branch("ASICS.Nhits",&ASICS->Nhits,"ASICS_hits/I");
  DataTree->Branch("ASICS.MB",ASICS->MB,"MB[ASICS_hits]/I");
  DataTree->Branch("ASICS.Chip",ASICS->Chip,"Chip[ASICS_hits]/I");
  DataTree->Branch("ASICS.Chan",ASICS->Chan,"Chan[ASICS_hits]/I");
  DataTree->Branch("ASICS.Energy",ASICS->Energy,"Energy[ASICS_hits]/I");
  DataTree->Branch("ASICS.Time",ASICS->Time,"Time[ASICS_hits]/I");

  // Branches for CAEN modules' data
  DataTree->Branch("CAEN.Nhits",&CAEN->Nhits,"CAEN_hits/I");
  DataTree->Branch("CAEN.Module",CAEN->Module,"Module[CAEN_hits]/I");
  DataTree->Branch("CAEN.Channel",CAEN->Channel,"Channel[CAEN_hits]/I");
  DataTree->Branch("CAEN.Value",CAEN->Value,"Value[CAEN_hits]/I");

  cout << "> Loop over evt files " <<endl; //debug

  ListEVT >> run_number;

  //Loop over files in the data file list.
  while(!ListEVT.eof()){

    evt_filename = data_dir + Form("run%d-%d.evt",run_number,BufferWords);
    file_subindex = 0;

    do {
      //Open evt file
      evtfile.open(evt_filename.c_str(),ios::binary);      
      if (!evtfile.is_open()){
	cout << ">\tCould not open " << evt_filename << endl;
	break;
      }
      else 
	cout << ">\tConverting " << evt_filename << " ..." << endl;
    
 
      //Loop over buffers in this file.
      while( !evtfile.eof() ){

	evtfile.read((char*)buffer,BufferBytes);

	if (!evtfile) {
	  //this could be a bad file or the file
	  // is subdivided into parts
	  //fileProblem = 1;
	  //cout << "  ** Warning: evt file problem" << endl;
	  break;
	}
      
	point = buffer;
	// cout << "buffer = " << buffer << endl; //debug
      
	//This is how many elements in the buffer (words) are useful (including this one).
	nWords = buffer[0];
	//From buffer[nWords-1] to buffer[BufferWords-1] is junk data. From buffer[0-15] is the
	//header of the buffer. 

	NBuffers++;
	//      if (NBuffers%50000 == 0) cout << "  Buffers read " << NBuffers << endl;
	epoint = point;
	if (epoint<point+nWords) words = *epoint++;
	//cout<<"Words "<<words<<endl;
	if (epoint<point+nWords) BufferType = *epoint++;
	//cout<<"Type "<<type<<endl;
      
	Nevents = *(point+6);
	TotEvents += Nevents;
     				
	string TitleStr;
      
	switch(BufferType){
	case 11: 
	  // begin run buffer
	  //runNum = *(epoint+1);
	  //cout << "run number="<< runNum << endl;
	  char Run_Title[200];
	  char *pTitle;
	  pTitle = (char*)(point+16);
	  for (Int_t i=0;i<200;i++) {
	    Run_Title[i] = *(pTitle++);
	  }
	  TitleStr.assign(Run_Title);
	  cout << ">\tRun number: " << *(point+3) << endl;
	  cout << ">\tRun title: " << TitleStr << " ..." << endl;
	  printf(">\tStart: %02i/%02i/%02i at %02i:%02i:%02i\n",*(point+58)+1,*(point+59),*(point+60)-100,
		 *(point+61),*(point+62),*(point+63));
	  break;
	
	case 12:
	  // end run buffer
	  break;
	
	case 2:
	  // scaler buffer
	  break;
	
	case 1:
	  // Physics buffer
	  Nevents = *(epoint+4);
	  TotEvents += Nevents;
	  epoint = epoint + 14;

	  for (unsigned int ievent=0; ievent<Nevents; ievent++) {

	    //Set all the elements of ROOT branches to zero before filling them.    
	    ASICS->Reset();    
	    CAEN->Reset();
	  
	    //create pointer inside of each  event
	    unsigned short * fpoint = epoint;
	  
	    words = *fpoint++;

	    int XLMdata1 = *fpoint++;
	  
	    if (XLMdata1==0xaaaa) {

	      fpoint +=2 ;
	      unsigned short Nstrips = *fpoint;
         
	      fpoint+=5;
      
	      for (int istrip=0; istrip<Nstrips; istrip++) {
		unsigned short *gpoint = fpoint;
		unsigned short id = *gpoint;
		unsigned short chipNum = (id&0x1FE0)>>5;
		unsigned short chanNum = id& 0x1F;
		gpoint++;
		int energy = *gpoint;
		gpoint++;
		unsigned short time = *gpoint;
		time = 16384 - time;

		if(chanNum<16){ 
		  if (chipNum == 3 || chipNum == 4 ) energy = 16384 - energy;
		  if (chipNum == 9 || chipNum == 10 ) energy = 16384 - energy;
		  if (chipNum == 13 || chipNum == 14) energy = 16384 - energy;
		  //Chipboard 8.
		  if (chipNum == 15 || chipNum == 16) energy = 16384 - energy;
		  if (chipNum == 17 || chipNum == 18) energy = 16384 - energy;
	  	  
		  if (energy>0 && energy<10000 && ASICS->Nhits<ASICS->Max_Hits) {
		    ASICS->MB[ASICS->Nhits] = 1;
		    ASICS->Chip[ASICS->Nhits] = chipNum;
		    ASICS->Chan[ASICS->Nhits] = chanNum;
		    ASICS->Energy[ASICS->Nhits] = energy;
		    ASICS->Time[ASICS->Nhits] = time;
		    ASICS->Nhits++;
		  }
	  
		  fpoint +=3;
		}
	      }// end for(istrip)       
	    }// end if(XLMdata1)


	    // Second XLM readout if present
	    int XLMdata2 = *fpoint++;
	    int counter;
    
	    counter = 0;
	    while (XLMdata2 != 0xbbbb) {
	      XLMdata2 = *fpoint++;
	      counter++;
	      if (counter>10) break;
	    }
    
	    if (XLMdata2==0xbbbb) {

	      fpoint += 2;
	      unsigned short Nstrips = *fpoint;
      
	      fpoint+=5;
      
	      for (int istrip=0;istrip<Nstrips;istrip++){
		unsigned short *gpoint = fpoint;
		unsigned short id = *gpoint;
		unsigned short chipNum = (id&0x1FE0)>>5;
		unsigned short chanNum = id& 0x1F;
		gpoint++;
		unsigned short energy = *gpoint;
		gpoint++;
		unsigned short time = *gpoint;
		time = 16384 - time;
		if (chanNum<16) {
		  if (chipNum == 3 || chipNum == 4 ) energy = 16384 - energy;
		  if (chipNum == 5 || chipNum == 6 ) energy = 16384 - energy;
		  if (chipNum == 7 || chipNum == 8 ) energy = 16384 - energy;
		  if (chipNum == 11 || chipNum == 12) energy = 16384 - energy;
	  
		  if (energy>0 && energy<10000 && ASICS->Nhits<ASICS->Max_Hits) {
		    ASICS->MB[ASICS->Nhits] = 2;
		    ASICS->Chip[ASICS->Nhits] = chipNum;
		    ASICS->Chan[ASICS->Nhits] = chanNum;
		    ASICS->Energy[ASICS->Nhits] = energy;
		    ASICS->Time[ASICS->Nhits] = time;
		    ASICS->Nhits++;
		  }	
		  fpoint +=3;
		}
	      }// end second for(istrip)
	    }// end if(XMLdata2)

	    //--- CAEN readout section -------------------------------------------
	    int caen = *fpoint++;
	    counter = 0;
	    while (caen != 0xcccc) {
	      caen = *fpoint++;
	      counter++;
	      if (counter>10) break;
	    }
	    //--- CAEN TDC / ADC readout -----------------------------------------------

	    if (caen==0xcccc) CAENCounter++;
	    while (fpoint < epoint + words) {
	      if(*fpoint == 0xffff ){
		fpoint++;
		continue;
	      }
	      unsigned short *gpoint = fpoint;
	      unsigned short chanCount = (*(gpoint++) & 0xff00)>>8;
	      unsigned short GEOaddress = (*(gpoint++) & 0xf800)>>11;
	      unsigned short data[32]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    
	      int i;

	      for (i=0; i<chanCount; i++) {
		if (i>31) continue;
		unsigned short ov  = (*gpoint&0x1000)>>12;
		unsigned short un  = (*gpoint&0x2000)>>13;
		unsigned short dat = (*(gpoint++)&0xfff);
		unsigned short geo = (*gpoint&0xf800)>>11;
		unsigned short chn = (*(gpoint++)&0x1f);
		if (geo == GEOaddress) {
		

		  //		  if (GEOaddress>=3 && GEOaddress<=7 && chn<32 && CAEN->Nhits<CAEN->Max_Hits) {
		  if (CAEN->Nhits<CAEN->Max_Hits) {
		    //cout << CAENCounter << " ov=" << ov << " un=" << un << " d=" << dat << " g=" << geo << " ch=" << chn << " ";

	    
		    CAEN->Module[CAEN->Nhits] = GEOaddress;
		    CAEN->Channel[CAEN->Nhits] = chn;
	    
		    if (ov) {
		      CAEN->Value[CAEN->Nhits] = 5000;
		      CAEN->Nhits++;
		    } else if (un) {
		      //		      CAEN->Value[CAEN->Nhits] = -1000;
		      CAEN->Value[CAEN->Nhits] = dat;
		      CAEN->Nhits++;
		    } else {
		      CAEN->Value[CAEN->Nhits] = dat;
		      CAEN->Nhits++;
		    }   
		  }
	  
		  // TDCs (DSG: there shouldn't be a difference between TDC and ADC).
		  /*
		  if (GEOaddress>=11 && GEOaddress<=12 && chn<32 && CAEN->Nhits<CAEN->Max_Hits) {
		    CAEN->Module[CAEN->Nhits] = GEOaddress;
		    CAEN->Channel[CAEN->Nhits] = chn;
		    CAEN->Value[CAEN->Nhits] = dat;
		    CAEN->Nhits++;
		  }	
		  */
		}
	      }
	    
	      unsigned short EOB_l = *(gpoint++);
	      unsigned short EOB_h = *(gpoint++);
	      unsigned short EOB_bit;
	      unsigned short geo = (EOB_h&0xf800)>>11;
	      EOB_bit = (EOB_h&0x0400)>>10;
	    
	      if (geo == GEOaddress && EOB_bit) {
		EOB_NEvents = EOB_l+(EOB_h&0x00ff)*65536+1;
	      } 
	    
	      while ((gpoint < epoint + words )&&(*gpoint==0xffff)){
		gpoint ++;
	      }
	      // go to next CAEN data
	      fpoint = gpoint;
	    
	    }// end if(fpoint<epoint+words)
	  
	    epoint += words+1; // This skips the rest of the event
	  
	    DataTree->Fill();
	    
	    outputCounter++;
	  
	  }//end for over events
	
	  BufferPhysics++;
	  break; //end of physics buffer
	
	  //      default: break;	
	}//end switch(BufferType)
      
      } //end loop over evtfile
    
      evtfile.close();

      // Check if thre is more than one file for this run number(this is denoted
      // by a '_#' where # is the file index number, starting at 1).
      file_subindex++;
      evt_filename = data_dir + Form("run%d_%d-%d.evt", run_number, file_subindex, BufferWords);
      evtfile.open(evt_filename.c_str(),ios::binary);      
      if (evtfile.is_open()) 
	go_to_next_run = 0;
      else 
	go_to_next_run = 1;
      evtfile.close();
    
    } while (!go_to_next_run);
    
    ListEVT >> run_number;
    
  }// end loop over list of files.
  
  cout << setprecision(3);
  cout << "> Total buffers = " << NBuffers << endl;
  cout << "> Physics buffers = " << BufferPhysics  << " (" 
       << 100.0*BufferPhysics/NBuffers << "\% of total buffers)" << endl;
  cout << "> Output file: " << output_file <<  endl;
  
  DataTree->Write();
  RootFile->Close();
  return;
}





////////////////////////////////////////////////////////////////////////////
void EventProcessor::CreateDetectors(int NumQ3Dets, int NumSX3Dets, int NumPCWires, int NumCsIDets, int MaxMultiplicity,
				     string CombinationsSX3andPC, string CombinationsSX3andCsI)
{
  //- Q3 section ------------------------------------------------------------------------
  this->NumQ3Dets = NumQ3Dets;
  // Q3 dets not yet ready to use.
 

  //- SX3 section -----------------------------------------------------------------------
  this->NumSX3Dets = NumSX3Dets;
  SX3 = new SX3_Detector*[NumSX3Dets];
  for (int d=0; d<NumSX3Dets; d++){
    SX3[d] = new SX3_Detector(Form("SX3_%d",d), MaxMultiplicity);
  }

  // Histogram for correlation between SX3 and CsI hits.
  Corr_SX3_CsI = new TH2F("Corr_SX3_CsI","SX3-CsI hit correlation", 30,-0.5,29.5,36,-0.5,35.5);
  Corr_SX3_CsI->GetXaxis()->SetTitle("SX3 detector ID");
  Corr_SX3_CsI->GetXaxis()->CenterTitle();
  Corr_SX3_CsI->GetYaxis()->SetTitle("CsI detector ID");
  Corr_SX3_CsI->GetYaxis()->CenterTitle();

  // Histogram to summarize the calibrated z-position in each SX3 detector.
  HistSX3PosSummary = new TH2F("HistSX3PosSummary","SX3 position summary", NumSX3Dets*4,0,NumSX3Dets*4, 160,0,8);
  HistSX3PosSummary->GetXaxis()->SetTitle("SX3 front strips");
  HistSX3PosSummary->GetXaxis()->CenterTitle();
  HistSX3PosSummary->GetYaxis()->SetTitle("z-position [cm]");
  HistSX3PosSummary->GetYaxis()->CenterTitle();
  

  //- PC section ------------------------------------------------------------------------
  this->NumPCWires = NumPCWires;
  PC = new PC_Detector*[NumPCWires];
  for (int w=0; w<NumPCWires; w++){
    PC[w] = new PC_Detector(Form("PC_%d",w), MaxMultiplicity);
    PC[w]->CreateHistograms();
  }

  // Histogram to summarize the calibrated z-position in each PC wire.
  HistPCPosSummary = new TH2F("HistPCPosSummary","PC position summary", NumPCWires,-0.5,NumPCWires-0.5, 500,-10,40);
  HistPCPosSummary->GetXaxis()->SetTitle("PC wire");
  HistPCPosSummary->GetXaxis()->CenterTitle();
  HistPCPosSummary->GetYaxis()->SetTitle("Hit z-position [cm]");
  HistPCPosSummary->GetYaxis()->CenterTitle();


  //- SX3 and PC combinations -----------------------------------------------------------
  // In here is where the geometrical constraints are loaded.
  Good_SX3_PC_comb = new bool*[NumSX3Dets];
  for (int d=0;  d<NumSX3Dets; d++){
    Good_SX3_PC_comb[d] = new bool[NumPCWires];
  }

  // Good combinations of detectors obtained from the correlation between Si and PC.
  for(Int_t si=0; si<NumSX3Dets; si++)
    for(Int_t w=0; w<NumPCWires; w++)
      Good_SX3_PC_comb[si][w]=0;

  if (CombinationsSX3andPC!="") {
    string aux;
    Int_t wire_id, sx3_id;
    CombinationsSX3andPC = ParamDirectory + CombinationsSX3andPC;
    ifstream FileDetComb(CombinationsSX3andPC.c_str());
    if(!FileDetComb.is_open()){
      cout << "> ERROR: File with combination \"" << CombinationsSX3andPC
	   << "\" could not be opened." << endl;
    }
    
    FileDetComb >> aux >> aux;
    do{
      FileDetComb >> sx3_id >> wire_id;
      Good_SX3_PC_comb[sx3_id][wire_id] = 1;
    } while(!FileDetComb.eof());
    FileDetComb.close();
    cout << "> Combinations of SX3 detector and PC wires loaded. " << endl;
  }

  //- CsI section -----------------------------------------------------------------------
  this->NumCsIDets = NumCsIDets;
  HistCsIRawESummary = new TH2F("HistCsIRawESummary","CsI RawE summary", 36,-0.5,36-0.5, 1024,0,4096);
  HistCsIRawESummary->GetXaxis()->SetTitle("CsI det");
  HistCsIRawESummary->GetXaxis()->CenterTitle();
  HistCsIRawESummary->GetYaxis()->SetTitle("RawE [arb. units]");
  HistCsIRawESummary->GetYaxis()->CenterTitle();

  //- SX3 and CsI combinations ----------------------------------------------------------
  // In here is where the geometrical constraints are loaded.
  Good_SX3_CsI_comb = new bool*[NumSX3Dets];
  for (int d=0;  d<NumSX3Dets; d++){
    Good_SX3_CsI_comb[d] = new bool[NumCsIDets];
  }
  // Good combinations of detectors obtained from the correlation between Si and CsI.
  for(Int_t si=0; si<NumSX3Dets; si++)
    for(Int_t csi=0; csi<NumCsIDets; csi++)
      Good_SX3_CsI_comb[si][csi]=0;
  if (CombinationsSX3andCsI!="") {
    string aux;
    Int_t csi_id, sx3_id;
    CombinationsSX3andCsI = ParamDirectory + CombinationsSX3andCsI;
    ifstream FileSX3CsIComb(CombinationsSX3andCsI.c_str());
    if(!FileSX3CsIComb.is_open()){
      cout << "> ERROR: File with SX3-CsI combination \"" << CombinationsSX3andCsI
	   << "\" could not be opened." << endl;
    }    
    FileSX3CsIComb >> aux >> aux;
    do{
      FileSX3CsIComb >> sx3_id >> csi_id;
      Good_SX3_CsI_comb[sx3_id][csi_id] = 1;
    } while(!FileSX3CsIComb.eof());
    FileSX3CsIComb.close();
    cout << "> Combinations of SX3 and CsI detectors loaded. " << endl;
  }

  return;
}


////////////////////////////////////////////////////////////////////////////
void EventProcessor::CreatePhysicalEventsFile(string filename)
{
  // ROOT output file
  PE_ROOT_File = new TFile(filename.c_str(),"RECREATE");  
  
  // Data Tree
  PE_Tree = new TTree("PE_Tree","Physical events TTree");

  // Branches for the physical event.
  PE_Tree->Branch("NHits",&Particle->NHits,"NHits/I");
  PE_Tree->Branch("FinalE",Particle->FinalE,"FinalE[NHits]/F");
  PE_Tree->Branch("FinalX",Particle->FinalX,"FinalX[NHits]/F");
  PE_Tree->Branch("FinalY",Particle->FinalY,"FinalY[NHits]/F");
  PE_Tree->Branch("FinalZ",Particle->FinalZ,"FinalZ[NHits]/F");
  PE_Tree->Branch("PC_Vcor",Particle->PC_V,"PC_Vcor[NHits]/F");
  PE_Tree->Branch("PC_Z",Particle->PC_Z,"PC_Z[NHits]/F");
  PE_Tree->Branch("PC_wire",Particle->PC_wire,"PC_wire[NHits]/I");
  PE_Tree->Branch("Time1",Particle->Time1,"Time1[NHits]/I");
  PE_Tree->Branch("Time2",Particle->Time2,"Time2[NHits]/I");


  // Create the PID histogrmas for each PC wire.
  PID_PC = new TH2F*[NumPCWires];
  for (int w=0; w<NumPCWires; w++) {
    PID_PC[w] = new TH2F(Form("PID_wire%d",w), Form("PID wire %d",w), 440,0,22,500,0,0.4);
    PID_PC[w]->GetXaxis()->SetTitle("Final energy [MeV]");
    PID_PC[w]->GetXaxis()->CenterTitle();
    PID_PC[w]->GetYaxis()->SetTitle("Path-corrected PC voltage [V]");
    PID_PC[w]->GetYaxis()->CenterTitle();
  }
  
  /*
  PID_CsI = new TH2F*[32];
  for (int w=0; w<NumPCWires; w++) {
    //      PID_CsI[w] = new TH2F(Form("PID_wire%d",w), Form("PID wire %d",w), 400,0,20,8000,0,8000);
    PID_CsI[w] = new TH2F(Form("PID_CsI%d",w), Form("PID CsI %d",w), 400,0,20,4000,0,2);
    PID_CsI[w]->GetXaxis()->SetTitle("Final energy [MeV]");
    PID_CsI[w]->GetXaxis()->CenterTitle();
    //      PID_CsI[w]->GetYaxis()->SetTitle("PC raw energy [arb. units]");
    PID_CsI[w]->GetYaxis()->SetTitle("PC voltage [V]");
    PID_CsI[w]->GetYaxis()->CenterTitle();
  }
  */
  
  HistTimes = new TH2F("HistTimes","",1024,0,4096, 1024,0,4096);
  HistTimes->GetXaxis()->SetTitle("RF time [arb.units]");
  HistTimes->GetYaxis()->SetTitle("SiPM time [arb.units]");

  WritePEFile = 1;
  cout << "> ROOT file with physical events \"" << filename 
       << "\" is ready to get some data." << endl;
  return;
}






////////////////////////////////////////////////////////////////////////////
void EventProcessor::Enable3DHits(int Max3DPoints)
{
  this->Max3DPoints = Max3DPoints;
  Save3DHits = 1;

  Float_t Lim_X = 12; // cm 
  Float_t Lim_Y = 12; // cm 
  Float_t Lim_Z_US = 40; // cm 
  Float_t Lim_Z_DS = -10; // cm 
  
  // The 'universe' can be delimited by two points at opposite conrnes of a 
  // cuboid (rectangular parallelepiped).
  Universe = new TGraph2D(2);
  Universe->SetPoint(0,-Lim_X,-Lim_Y,Lim_Z_US);
  Universe->SetPoint(1, Lim_X, Lim_Y,Lim_Z_DS);
  // Formatting the axes
  Universe->GetHistogram("empty")->GetXaxis()->SetTitle("x (horizontal) [cm]");
  Universe->GetHistogram("empty")->GetXaxis()->CenterTitle();
  Universe->GetHistogram("empty")->GetXaxis()->SetTitleOffset(1.5);
  Universe->GetHistogram("empty")->GetYaxis()->SetTitle("y (vertical) [cm]");
  Universe->GetHistogram("empty")->GetYaxis()->CenterTitle();
  Universe->GetHistogram("empty")->GetYaxis()->SetTitleOffset(1.5);
  Universe->GetHistogram("empty")->GetZaxis()->SetTitle("z (up stream) [cm]");
  Universe->GetHistogram("empty")->GetZaxis()->CenterTitle();
  Universe->GetHistogram("empty")->GetZaxis()->SetTitleOffset(1.5);
  
  Hits3D = new TGraph2D(Max3DPoints);
  return;
}
  

////////////////////////////////////////////////////////////////////////////
//  Once the ROOT file has been opened and the histograms and data arrays
//  have been created by using CreateHistograms() and CreateArrays(), 
//  respectively, one can use this function to read the data in the ROOT
//  file and store it in the arrays for futher computations and use the
//  histograms for data visualization.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ExtractData(string HitType, bool EmptyChamber){

  cout << "> Extracting data ..." << endl;

  Long64_t ASICS_counter=0, CAEN_counter=0, PE_counter=0;
  Long64_t ASICS_hit_counter=0, SX3_hit_counter=0, Q3_hit_counter=0;
  Long64_t CAEN_hit_counter=0, PC_hit_counter=0, SiPM_hit_counter=0, RF_hit_counter=0;
  Long64_t counter1=0, counter2=0, counter3=0;
  Long64_t TotEntries = Chain->GetEntries();
  Int_t Det_ID, Det_Ch, m, points3D=0;
  Int_t bch, stp, hit, w, d, h;
  Int_t SiPM_time = 0, RF_time=0;
  Float_t x, y, z, dist, Qb_cal;
  string Det_Type;
  Int_t WireWithMaxV;
  Float_t MaxV;
  Float_t CsI_Energy[36];
  bool Veto = 0, GoodTime1 = 0;
  TStopwatch ProcessTime;

  cout << ">\tTotal entries: " << TotEntries;
  if (DataFraction<1.0)
    cout << " (taking only " << DataFraction*100 << "% of total entries)";
  cout << "\n";

  // Reset the number of internal events for all the detectors.
  for (d=0; d<NumSX3Dets; d++) 
    SX3[d]->InternalEvents = 0;
  for (w=0; w<NumPCWires; w++) 
    PC[w]->InternalEvents = 0;

  // Index for the arrays qb and qu. It will increment if a valid hit is found.
  int event=0;
  ProcessTime.Start();
  
  for(Long64_t n=0; n<TotEntries; n++){
      
    Particle->Reset();
    SiPM_time = RF_time = 0;
    for (int det=0; det<36; det++)
      CsI_Energy[det] = 0;

    Chain->GetEntry(n);
    
    if (n == TMath::Nint(0.10*TotEntries))  cout << ">\t10% through the data" << endl;
    if (n == TMath::Nint(0.25*TotEntries))  cout << ">\t25% through the data" << endl;
    if (n == TMath::Nint(0.50*TotEntries))  cout << ">\t50% through the data" << endl;
    if (n == TMath::Nint(0.75*TotEntries))  cout << ">\t75% through the data" << endl;
    if (n == TMath::Nint(0.90*TotEntries))  cout << ">\t90% through the data" << endl;
    if (n>TMath::Nint(DataFraction*TotEntries)) break;
    
    if (ASICS_Nhits>0) 
      ASICS_counter++;
    if (CAEN_Nhits>0)
      CAEN_counter++;

    // CAEN hit section ----------------------------------------------------------
    if (NumPCWires>0) {
      for (w=0; w<NumPCWires; w++) 
	PC[w]->Reset();    
    }

    for (h=0; h<CAEN_Nhits; h++) {
      CAEN_hit_counter++;
      Det_ID = Det_Ch = -1;
      Det_Type = "";
      Map->IdentifyCAEN(CAEN_Module[h], CAEN_Channel[h], Det_Type, Det_ID, Det_Ch);
      
      if (Det_ID>=0 && Det_Ch>=0 && CAEN_Value[h]>0 && CAEN_Value[h]<4096) {
	// Proportional counter
	if (Det_Type=="PC" && NumPCWires>0 && !EmptyChamber) {    
	  PC_hit_counter++;
	  m = PC[Det_ID]->mult;
	  PC[Det_ID]->RawE[m] = CAEN_Value[h];
	  PC[Det_ID]->Channel[m] = Det_Ch;
	  PC[Det_ID]->HistRawE[Det_Ch]->Fill(CAEN_Value[h]);
	  PC[Det_ID]->HistRawESummary->Fill(Det_Ch, CAEN_Value[h]);	    
	  // Increase the multiplicity of this detector.
	  PC[Det_ID]->mult++;	   
	} 
	// RF time
	else if (Det_Type == "RF") {
	  RF_hit_counter++;
	  RF_time = CAEN_Value[h];
	}	
	// Silicon photo-multiplier
	else if (Det_Type == "SiPM") {
	  SiPM_hit_counter++;
	  SiPM_time = CAEN_Value[h];
	}
	// Ionization chamber
	else if (Det_Type == "IC") {
	  
	}
	else if (Det_Type == "CsI") {
	  CsI_Energy[Det_ID] = CAEN_Value[h];
	  HistCsIRawESummary->Fill(Det_ID, CAEN_Value[h]);
	}
      }  
    }// End loop over CAEN_Nhits
  
    // Fill the multiplicity histogram of each PC wire.
    if (NumPCWires>0) {
      for (int w=0; w<NumPCWires; w++) {
	m = PC[w]->mult;
	if (m>0)
	  PC[w]->HistMult->Fill(m);
      }
    }

    // ASICs hit section  ---------------------------------------------------------
    // Befor going into the ASICs hit loop the event related quantities have to
    // be reset.
    for (d=0; d<NumSX3Dets; d++)
      SX3[d]->Reset();
    
    for (h=0; h<ASICS_Nhits; h++) {
      ASICS_hit_counter++;
      Det_ID = Det_Ch = -1;
      Det_Type = "";
      Map->IdentifyASICS(ASICS_MB[h], ASICS_Chip[h], ASICS_Chan[h], Det_Type, Det_ID, Det_Ch);
      
      // Verify that the MB, Chip and Chan combination exists.
      if (Det_ID>=0 && Det_Ch>=0 && ASICS_Energy[h]>0 && ASICS_Energy[h]<15000) {
	if (Det_Type=="Q3") {
	  Q3_hit_counter++;
	  // Do the Q3 stuff.
	} else if (Det_Type=="SX3") {
	  SX3_hit_counter++;
	  m = SX3[Det_ID]->mult;
	 
	  SX3[Det_ID]->Channel[m] = Det_Ch;
	  SX3[Det_ID]->RawE[m] = ASICS_Energy[h];
	  SX3[Det_ID]->Time[m] = ASICS_Time[h];
	  SX3[Det_ID]->HistRawE[Det_Ch]->Fill(ASICS_Energy[h]);
	  SX3[Det_ID]->HistRawESummary->Fill(Det_Ch, ASICS_Energy[h]);
	  SX3[Det_ID]->HistTimeSummary->Fill(Det_Ch, ASICS_Time[h]);	  
	  // Increase the multiplicity of this detector.
	  SX3[Det_ID]->mult++;
	}
      }// end if(Det_ID && Det_Ch)
    }// end for(h)
    

    // After extracting the CAEN and ASICs data -----------------------------------
    // Check if there was a time gate provided by the user.
    if (ApplyTime1Cuts) {
      GoodTime1 = 0;
      for (int c=0; c<NT1Cuts; c++) 
	if (RF_time>=T1CutLowerLimit[c] && RF_time<=T1CutUpperLimit[c]) {
	  GoodTime1 = 1;
	  break;
	}
      if (!GoodTime1)
	continue;
    }
    
    // Actions taken for each SX3 detector once an event is fully read.
    for (d=0; d<NumSX3Dets; d++) {
      m = SX3[d]->mult;
      if (m>0) {
	// Fill the multiplicity histogram.
	SX3[d]->HistMult->Fill(m);
	if (m>1 && HitType!="raw") {
	  counter1++;  // aux counter: SX3 hits with m>1.
	  SX3[d]->GetEventQs(HitType);
	  if (SX3[d]->EventType>0 && SX3[d]->IsOn) {
	    bch = SX3[d]->EventBch;
	    stp = SX3[d]->EventStp;
	    // These histograms are used in the energy calibration.
	    Qb_cal = (SX3[d]->mB[bch])*(SX3[d]->EventQb) + (SX3[d]->oB[bch]);
	    SX3[d]->HistQb[bch]->Fill(Qb_cal);
	    SX3[d]->HistQb[4]->Fill(Qb_cal);
	    SX3[d]->GetEventEP();
	    // Check if the final energy is within a reasonable range.
	    if (SX3[d]->EventE>0.2 && SX3[d]->EventE<30.0) {
	      RF_SX3ID->Fill(4*d + bch, RF_time);
	      RF_SiT->Fill(SX3[d]->EventT, RF_time);
	      SX3[d]->HistCalESummary->Fill(SX3[d]->EventBch, SX3[d]->EventE);
	      // Filling position histograms if necessary (used to check the SX3 calibration).
	      if (SX3[d]->HistPvECreated)
		SX3[d]->FillHistPvE();
	      // Filling charge difference histograms if necessary (used to check the SX3 calibration).
	      if (SX3[d]->HistDiffQCreated)
		SX3[d]->FillHistDiffQ();
	      // Veto SX3 punch-through. A particle is belived to punch-through only when a
	      // silicon detector and it's corresponding CsI detector (located behind the
	      // silicon det.) fire in the same event.
	      Veto = 0;
	      for (int d_CsI=0; d_CsI<NumCsIDets; d_CsI++) 
		if (SX3[d]->EventE>3 && Good_SX3_CsI_comb[d][d_CsI] && CsI_Energy[d_CsI]>100) 
		  Veto = 1;
	      
	      if (Map->WorldCoordinatesLoaded  && !Veto && SX3[d]->IsOn) {
		x = y = 0;
		z = -100;
		// For this particular experiment some SX3s were used as forward detectors (DetID=0-5),
		// thus a special way of getting the world coordinates needed for these detectors.
		if (d<6)
		  Map->GetSX3FwdWorldCoordinates(d, SX3[d]->EventStp, SX3[d]->EventP, x, y, z);
		else 
		  Map->GetSX3WorldCoordinates(d, SX3[d]->EventStp, SX3[d]->EventP, x, y, z);
		if (z>-10 && z<40 && sqrt(x*x+y*y)<15) {
		  counter2++;  // aux counter
		  // This is used in the GlobalGeometryTest function.
		  if (FillFinalEvDist) {
		    dist = sqrt(x*x + y*y + pow(z_228Th-z,2));
		    HistEvD->Fill(dist, SX3[d]->EventE);
		  }
		  // In case we want to show the 3D hits we save them.
		  if (Save3DHits && points3D<Max3DPoints) {
		    Hits3D->SetPoint(points3D, x, y, z);
		    points3D++;
		  }
		  // After the event has proven to have good silicon detector signals we have two
		  // possible scenarios: the chamber is filled with gas or it is empty.
		  if (!EmptyChamber){
		    WireWithMaxV = -1;
		    MaxV = 0;
		    // PC hit analysis when the chamber has gas.
		    for (w=0; w<NumPCWires; w++) {
		      if (PC[w]->mult==2 && Good_SX3_PC_comb[d][w]) {
			counter3++;  // aux counter
			PC[w]->GetEventVs();
			if (PC[w]->EventVd>0 && PC[w]->EventVu>0 && (PC[w]->EventVd + PC[w]->EventVu)>MaxV) {
			  PC[w]->EventFinalE = SX3[d]->EventE;
			  PC[w]->EventFinalX = x;
			  PC[w]->EventFinalY = y;
			  PC[w]->EventFinalZ = z;
			  PC[w]->EventTime1 = RF_time;
			  PC[w]->EventTime2 = SiPM_time;
			  PC[w]->GetEventZ();
			  MaxV = PC[w]->EventVd + PC[w]->EventVu;
			  WireWithMaxV = w;
			}	      		     
		      }
		    } // end for(w)	
		    if (WireWithMaxV>=0) {
		      w = WireWithMaxV;
		      // Save the event to the internal arrays, if they have been created.
		      if (PC[w]->ArraysCreated) 
			PC[w]->SaveEvent(PC[w]->InternalEvents);
		      PC[w]->InternalEvents++;
		      //-----------------------------------------------------------------------------		   
		      // The event is complete. Set the quantities related to the physical event.
		      hit = Particle->NHits;
		      Particle->FinalE[hit] = SX3[d]->EventE;
		      Particle->FinalX[hit] = x;
		      Particle->FinalY[hit] = y;
		      Particle->FinalZ[hit] = z;
		      Particle->PC_wire[hit] = w;
		      Particle->PC_Z[hit] = PC[w]->EventZ;
		      Particle->PC_V[hit] = PC[w]->GetCorrectedV();
		      Particle->Time1[hit] = RF_time;
		      Particle->Time2[hit] = SiPM_time;		    
		      // Check the hit correlation between the SX3 and CsI detectors.
		      for (int d_CsI=0; d_CsI<36; d_CsI++)
			if (CsI_Energy[d_CsI]>290 && SX3[d]->EventE>4)
			  Corr_SX3_CsI->Fill(d, d_CsI);
		      HistPCPosSummary->Fill(w, Particle->PC_Z[hit]);
		      HistSX3PosSummary->Fill(4*d+SX3[d]->EventStp, SX3[d]->EventP);
		      if (WritePEFile) {
			PID_PC[w]->Fill(Particle->FinalE[hit], Particle->PC_V[hit]);
			HistTimes->Fill(RF_time, SiPM_time);
		      }
		      Particle->NHits++; 	    
		    }// end if(WireWithMaxV>=0)
		  } // end if(!EmptyChamber)
		  else {
		    // When the chamber is empty there cannot be PC data, thus we exclude the
		    // PC from the physical event.
		    hit = Particle->NHits;
		    Particle->FinalE[hit] = SX3[d]->EventE;
		    Particle->FinalX[hit] = x;
		    Particle->FinalY[hit] = y;
		    Particle->FinalZ[hit] = z;
		    Particle->Time1[hit] = RF_time;
		    Particle->Time2[hit] = SiPM_time;
		    HistTimes->Fill(RF_time, SiPM_time);
		    Particle->NHits++;		  
		  } // end else if(!EmptyChamber)
		}
	      } // end if(Map->WorldCoordinatesLoaded)	    
	      if (SX3[d]->ArraysCreated)
		SX3[d]->SaveEventQs(SX3[d]->InternalEvents);
	      SX3[d]->InternalEvents++;
	    } // end if(energy within reasonable range)
	  }// end if(EventType>0)
	}
      }
    } // end for(d)
  
    if (Particle->NHits>0 && WritePEFile) {
      PE_Tree->Fill();
      PE_counter++;
    }
    
    // Increase the event number by 1.
    event++;    
  }// end for(n)

  ProcessTime.Stop();
  cout << ">\tTime processing the event loop: " << ProcessTime.RealTime() << " s" << endl;
  
  Entries = event;
  
  cout << ">\tEvents with ASICs data: " << ASICS_counter << " (" 
       << 100*ASICS_counter/TotEntries << "% of total events)"<< endl;
  cout << ">\tEvents with CAEN data: " << CAEN_counter << " (" 
       << 100*CAEN_counter/TotEntries << "% of total events)"<< endl;
  cout << ">\tASICS hits: " << ASICS_hit_counter << endl;
  if (ASICS_hit_counter>0) {
    cout << ">\t\twith SX3 hits: " << SX3_hit_counter << " (" 
	 << 100*SX3_hit_counter/ASICS_hit_counter << "% of ASICs hits)"<< endl;
    cout << ">\t\twith Q3 hits: " << Q3_hit_counter << " (" 
	 << 100*Q3_hit_counter/ASICS_hit_counter << "% of ASICs hits)"<< endl;
  }
  cout << ">\tCAEN hits: " << CAEN_hit_counter << endl; 
  if (CAEN_hit_counter>0) {
    cout << ">\t\twith PC hits: " << PC_hit_counter << " (" 
	 << 100*PC_hit_counter/CAEN_hit_counter << "% of CAEN hits)"<< endl;
    cout << ">\t\twith RF hits: " << RF_hit_counter << " (" 
	 << 100*RF_hit_counter/CAEN_hit_counter << "% of CAEN hits)"<< endl;
    cout << ">\t\twith SiPM hits: " << SiPM_hit_counter << " (" 
	 << 100*SiPM_hit_counter/CAEN_hit_counter << "% of CAEN hits)"<< endl;
  }
    
  if (WritePEFile) {   
    cout << ">\tSX3 hits with mult>1: H=" << counter1 << endl;
    cout << ">\t& SX3 hits with reasonable energy and position: " << counter2 <<" (" << 100*counter2/counter1 
	 << "% of H)"<< endl;
    cout << ">\t& PC hits with Qu>0 and Qd>0 and geometrical constraint: " << endl;
    cout << ">\t=> Physical events = " << PE_counter <<" (" << 100*PE_counter/counter1 
	 << "% of H)"<< endl;
  }

  DataExtracted = 1;
  return;
  
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::GlobalGeometryTest(float SourcePos, int FwdDet, int R1Det, int R2Det, string FileELoss)
{
#if 0
  EnergyLoss* a_in_H;
  float Peaks228Th[6]={5.41, 5.685, 6.051, 6.288, 6.778, 8.785};
  TCanvas* Display;
  TGraph** EvD_graph;

  HistEvD = new TH2F("HistEvD", Form("Showing detectors %d, %d and %d", FwdDet,R1Det,R2Det), 400,9,29,900,0,9);
  HistEvD->GetXaxis()->SetTitle("Distance from source [cm]");
  HistEvD->GetXaxis()->CenterTitle();
  HistEvD->GetYaxis()->SetTitle("SX3 energy [MeV]");
  HistEvD->GetYaxis()->CenterTitle();

  a_in_H = new EnergyLoss(FileELoss);

  EvD_graph = new TGraph*[6];
  for (int i=0; i<6; i++) {
    a_in_H->GetEvDCurve(Peaks228Th[i], 30, 300);
    a_in_H->EvD->SetLineWidth(2);
    a_in_H->EvD->SetLineColor(kBlue);
  }
  
  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++)
    if (d==FwdDet || d==R1Det || d==R2Det) {
      SX3[d]->ForceBchPosRange = 1;
      continue;
    }
    else
      SX3[d]->IsOn = 0;

  // These 2 are used in ExtractData().
  FillFinalEvDist = 1;
  z_228Th = SourcePos;
  
  ExtractData();

  FillFinalEvDist = 0; // Turnning the switch off in case Extract data is used later.
  
  Display = new TCanvas("DispGeomTes","Global geometry test");
  a_in_H->EvD->Draw("colz");
  for (int i=0; i<6; i++) 
    EvD_graph[i]->Draw("l same");
#endif
  return;
}


////////////////////////////////////////////////////////////////////////////
// This fucntion attempts to load all the calibratrion coefficients for all
// the detectors assuming a predefined format for the file names.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::LoadDetectorsParameters()
{
  string filename;
  for (int d=0; d<NumSX3Dets; d++) {
    filename = ParamDirectory + Form("SX3_%d.coeffs",d);
    SX3[d]->LoadCoefficients(filename);
  }
  for (int w=0; w<NumPCWires; w++) {
    filename = ParamDirectory + Form("PCWire_%d.coeffs",w);
    PC[w]->LoadCoefficients(filename);
  }

  return;
}

////////////////////////////////////////////////////////////////////////////
// Similar to the function LoadDetectorsParameters() but for only a type of
// detector and from an initial to a final detector number.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::LoadParamsSelectedDetectors(string DetType, int InitDet, int FinalDet)
{
  string filename;
  if (DetType=="SX3") {
    for (int d=InitDet; d<=FinalDet; d++) {
      filename = ParamDirectory + Form("SX3_%d.coeffs",d);
      SX3[d]->LoadCoefficients(filename);
    }
  } else if (DetType=="PC") {
    for (int w=InitDet; w<=FinalDet; w++) {
      filename = ParamDirectory + Form("PCWire_%d.coeffs",w);
      PC[w]->LoadCoefficients(filename);
    }
  }

  return;
}





////////////////////////////////////////////////////////////////////////////
//  Simple function that opens a ROOT file (read mode).  The method
//  SetTreeName() must be called before using this method.
////////////////////////////////////////////////////////////////////////////
bool EventProcessor::OpenROOTFile(string filename)
{
  ROOTFileName = filename;
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
void EventProcessor::PC_CalibratePosition(float pos_source_DS, float pos_source_US, string AlphaELossFile,
					  int InitWire, int FinalWire)
{
  char GoodParameters;
  cout << "---- PC pos. calibration (begin) -------------------------------" << endl;

  // Before getting to work make sure the energy loss file exists.
  AlphaELossFile = ParamDirectory + AlphaELossFile;
  EnergyLoss ELossTest(AlphaELossFile);
  if (!ELossTest.GoodELossFile) {
    cout << "---- PC pos. calibration (end) ---------------------------------" << endl;
    return;
  }

  // Reset the calibration coefficients to their original value.
  for (int WireID=InitWire; WireID<=FinalWire; WireID++) {
    for (int order=0; order<=PC[WireID]->PosPolyOrder; order++) {
      PC[WireID]->Cd[order] = 1;
      PC[WireID]->Cu[order] = 1;
    }
  }
 
  // Sort the data to count how many internal events are for each PC wire.
  ExtractData();
  
  for (int w=InitWire; w<=FinalWire; w++) {
    cout << "> Working on detector " << PC[w]->Name << endl;
    // Create the internal events' arrays and set some parameters.
    PC[w]->CreateArrays();
 
      // Resort the data and save the internal events.
    cout << "> Saving data in the internal arrays of the detector ..." << endl;
    ExtractData();
    
    GoodParameters = 'n';
    // First we give it a try with the default values of the parameters.
    PC[w]->GetPositionCoeffs(pos_source_DS, pos_source_US, AlphaELossFile);
    cout << "> Are these good values for the parameters (y/n): ";
    cin >> GoodParameters;
    if (GoodParameters=='n') {
      cout << "> Cu and Cd reset." << endl;
      PC[w]->Cd[0] = -50;
      PC[w]->Cu[0] = -50;
      for (int order=1; order<=PC[w]->PosPolyOrder; order++) {
	PC[w]->Cd[order] = 0;
	PC[w]->Cu[order] = 0;
      }        
    }
    else
      PC[w]->WriteCoefficients(ParamDirectory+Form("PCWire_%d.coeffs",w));
    
    // Since the data arrays could be pretty large they are deleted once they have
    // been used to keep the memory usage low.  The same is done for the HistPvE's.
    PC[w]->DeleteArrays();
  }
  
  cout << "---- PC pos. calibration (end) ---------------------------------" << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Finds the pulser offsets and saves them in a text file with the PC wire
// calibration coefficients. It does this from InitWire to FinalWire. The
// first argument is and array with the voltages that were used. The method
// SetPeakFindingParameters() must be called before calling this method.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::PC_GetPulserOffsets(float* Voltages, int InitWire, int FinalWire)
{ 
  cout << "---- PC pulser offsets (begin) ---------------------------------" << endl;

  if (FinalWire>NumPCWires-1)
    FinalWire = NumPCWires-1;
 
  // Sort the data to fill the HistRawE histograms.
  ExtractData("raw");
  
  for (int w=InitWire; w<=FinalWire; w++) {
    cout << "> Working on detector " << PC[w]->Name << endl;
    PC[w]->GetPulserOffsets(NPeaks, PeakSigma, PeakMinAmplitude, NormalChi2Value, Voltages);
    PC[w]->WriteCoefficients(ParamDirectory+Form("PCWire_%d.coeffs",w));
  }
  
  cout << "---- PC pulser offsets (end) -----------------------------------" << endl;
  return; 
}



////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////
void EventProcessor::PC_SetThresholds(float DisplayLimit, int InitWire, int FinalWire)
{ 
  cout << "---- Set PC thresholds (begin) ---------------------------------" << endl;

  if (FinalWire==-1)
    FinalWire = NumPCWires-1;
 
  // Sort the data to fill the HistRawE histograms.
  ExtractData("raw");
  
  for (int w=InitWire; w<=FinalWire; w++) {
    cout << "> Working on detector " << PC[w]->Name << endl;
    PC[w]->SetThresholds(DisplayLimit);
    PC[w]->WriteCoefficients(ParamDirectory+Form("PCWire_%d.coeffs",w));
  }
  
  cout << "---- Set PC thresholds (end) -----------------------------------" << endl;
  return; 
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::PC_ShowPosSummary()
{
  TCanvas* Display;
  Display = new TCanvas("PCPosSumm","PC pos. summary",1400,800);
  HistPCPosSummary->Draw("colz");
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ResetChain()
{
  Chain->Reset();
  delete Chain;
  Chain = new TChain(TreeName);  
  return;
}


////////////////////////////////////////////////////////////////////////////
// This method is meant to be used after all the ROOT files have been opened.  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SetBranches()
{
  // This procedure to is similar to the one used in a TreeSelector. 
  Chain->SetBranchAddress("ASICS.Nhits", &ASICS_Nhits);
  Chain->SetBranchAddress("ASICS.MB", ASICS_MB);
  Chain->SetBranchAddress("ASICS.Chip", ASICS_Chip);
  Chain->SetBranchAddress("ASICS.Chan", ASICS_Chan);
  Chain->SetBranchAddress("ASICS.Energy", ASICS_Energy);
  Chain->SetBranchAddress("ASICS.Time", ASICS_Time);
  Chain->SetBranchAddress("CAEN.Nhits", &CAEN_Nhits);
  Chain->SetBranchAddress("CAEN.Module", CAEN_Module);
  Chain->SetBranchAddress("CAEN.Channel", CAEN_Channel);
  Chain->SetBranchAddress("CAEN.Value", CAEN_Value);
  return;
}



////////////////////////////////////////////////////////////////////////////
//  Set the channel map for ASICs and CAEN channels.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SetChannelMap(string ASICS_map, string CAEN_map)
{
  ASICS_map = ParamDirectory + ASICS_map;
  CAEN_map = ParamDirectory + CAEN_map;
  Map = new ChannelMap(ASICS_map, CAEN_map);
  return;
}

void EventProcessor::SetChannelMap(string ASICS_map, string CAEN_map, string SX3_geometry_map)
{
  ASICS_map = ParamDirectory + ASICS_map;
  CAEN_map = ParamDirectory + CAEN_map;
  SX3_geometry_map = ParamDirectory + SX3_geometry_map;
  Map = new ChannelMap(ASICS_map, CAEN_map);
  Map->InitWorldCoordinates(SX3_geometry_map);
  return;
}


////////////////////////////////////////////////////////////////////////////
//  Set the directory where the parameters such as calibration coefficients
//  can be found.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SetParamDirectory(string ParamDirectory){
  this->ParamDirectory = ParamDirectory;
  return;
}


////////////////////////////////////////////////////////////////////////////
// The user provides the parameters related to peak finding. These are used,
// for example, when looking for the pulser offsets and when getting the 
// relative energy coefficients in the SX3 detectors.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SetPeakFindingParameters(int NPeaks, float PeakSigma, float PeakMinAmplitude, float NormalChi2Value)
{
  this->NPeaks = NPeaks;
  this->PeakSigma = PeakSigma;
  this->PeakMinAmplitude = PeakMinAmplitude;
  this->NormalChi2Value = NormalChi2Value;   // Used in linear fits.
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SetTime1Cuts(int NT1Cuts, float* T1CutLowerLimit, float* T1CutUpperLimit)
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
void EventProcessor::SetTreeName(TString TreeName){
  this->TreeName = TreeName;
  Chain = new TChain(TreeName);  
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::Show3DHits()
{
  TCanvas* Display;
  if (Save3DHits) {
    cout << "> Displaying 3D hits. Double-click on the canvas to delete it." << endl;
    Display = new TCanvas("Disp3DHits","3D hits",1300,800);
    Display->cd(1);
    Display->cd(1)->SetPhi(150);
    Display->cd(1)->SetTheta(10);
    Universe->SetTitle("");
    Universe->Draw("p");
    Hits3D->Draw("same p");
    Display->Update();
    Display->WaitPrimitive();
    delete Display;
  }
  else
    cout << "> ERROR: Didn't call Enable3DHits()." << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ShowCalESummary()
{
  TCanvas* Display;

  Display = new TCanvas("SX3_CalE_R1","SX3 CalE Ring 1",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=0; d<(int)(NumSX3Dets/2); d++) {
    Display->cd(d+1);
    SX3[d]->HistCalESummary->Draw("col");
  }
  
  Display = new TCanvas("SX3_CalE_R2","SX3 CalE Ring 2",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=(int)(NumSX3Dets/2); d<NumSX3Dets; d++) {
    Display->cd(d+1-(int)(NumSX3Dets/2));
    SX3[d]->HistCalESummary->Draw("col");
  }
  return;
}



////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ShowMultiplicity()
{
  TCanvas* Display;

  Display = new TCanvas("D.SX3_mult_1","SX3 Mult Part 1",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=0; d<(int)(NumSX3Dets/2); d++) {
    Display->cd(d+1);
    SX3[d]->HistMult->Draw();
  }

  Display = new TCanvas("D.SX3_mult_2","SX3 Mult Part 2",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=(int)(NumSX3Dets/2); d<NumSX3Dets; d++) {
    Display->cd(d+1-(int)(NumSX3Dets/2));
    SX3[d]->HistMult->Draw();
  }

  Display = new TCanvas("D.PC_mult","PC Mult",1400,800);
  Display->Divide((int)ceil((float)NumPCWires/4),4);
  for (int w=0; w<NumPCWires; w++) {
    Display->cd(w+1);
    PC[w]->HistMult->Draw();
  }

  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ShowSX3CsICorrelation()
{
  TCanvas* Display;
  Display = new TCanvas("SX3_CsI","SX3-CsI correlation",1400,800);
  Corr_SX3_CsI->Draw("colz");
  Display = new TCanvas("CsI_RawE_Summ","CsI RawE",1400,800);
  HistCsIRawESummary->Draw("col");
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ShowRawESummary()
{
  TCanvas* Display;

  Display = new TCanvas("SX3_RawE_1","SX3 RawE Part 1",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=0; d<(int)(NumSX3Dets/2); d++) {
    Display->cd(d+1);
    SX3[d]->HistRawESummary->Draw("col");
  }

  Display = new TCanvas("SX3_RawE_2","SX3 RawE Part 2",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=(int)(NumSX3Dets/2); d<NumSX3Dets; d++) {
    Display->cd(d+1-(int)(NumSX3Dets/2));
    SX3[d]->HistRawESummary->Draw("col");
  }

  Display = new TCanvas("PC_RawE_Summ","PC RawE",1400,800);
  Display->Divide((int)ceil((float)NumPCWires/4),4);
  for (int w=0; w<NumPCWires; w++) {
    Display->cd(w+1);
    PC[w]->HistRawESummary->Draw("col");
  }
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::ShowTimeSummary()
{
  TCanvas* Display;
  Display = new TCanvas("SX3_Time_1","SX3 Time Part 1",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=0; d<(int)(NumSX3Dets/2); d++) {
    Display->cd(d+1);
    SX3[d]->HistTimeSummary->Draw("col");
  }

  Display = new TCanvas("SX3_Time_2","SX3 Time Part 2",1400,800);
  Display->Divide((int)(NumSX3Dets/2/3),3);
  for (int d=(int)(NumSX3Dets/2); d<NumSX3Dets; d++) {
    Display->cd(d+1-(int)(NumSX3Dets/2));
    SX3[d]->HistTimeSummary->Draw("col");
  }
  Display = new TCanvas("RF_Time","RF Time",1400,800);
  RF_SX3ID->Draw("col");
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// This method performs the full energy calibration for the back channels of the
// SX3 detectors. Arguments:
// - *back_ref, is an array with the reference back channels for all the (24)
//    detectors.
// - *SourcePeaks, is an array with the source's alpha energies.
// - InitDet and FinalDet set the range of detactors for which the coeffs. will be 
//   obtained.
// - GetRelCoeffs, its default value is 1 ( 'on'). When set to 0 ('off') it will
//   skip the step in which the relative coefficients are obtained and go straight
//   to the step where the detector's energy coefficientes are obtained.
///////////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_CalibrateEnergy(int* back_ref, float* SourcePeaks, float MinQ, float MaxQ,
					int InitDet, int FinalDet, bool GetRelCoeffs)
{ 
  cout << "---- SX3 energy calibration (begin) ---------------------------" << endl;

  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;

  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d>=InitDet && d<=FinalDet)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }

  // I've seen it's better to get only full hits when working on calibrations and 
  // to reset the calibration coefficients to their original value.
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    SX3[Det_ID]->ForceFullHits = 1;
    if (GetRelCoeffs)
      for (int bch=0; bch<4; bch++) {
	SX3[Det_ID]->mB[bch] = 1;
	SX3[Det_ID]->oB[bch] = 0;
      }   
    SX3[Det_ID]->Elin = 1;
    SX3[Det_ID]->Eshift = 0;
  }
    
  if (GetRelCoeffs) {
    // Sort the data to fill the HistQb histograms.
    ExtractData("normal");
    
    for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
      cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
      SX3[Det_ID]->GetRelativeCoeffs(back_ref[Det_ID], PeakSigma, PeakMinAmplitude, MinQ, MaxQ);
      // It's good to write the coefficients in this intermediate step.  If you
      // wait until the last step some of your work may get lost.
      SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
      
      for (int bch=0; bch<5; bch++)
	SX3[Det_ID]->HistQb[bch]->Reset();
    }
  }
  
  // Resort the data now that the relative coefficients have been obtained.
  ExtractData("normal");
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
    SX3[Det_ID]->GetEnergyCoeffs(NPeaks, PeakSigma, PeakMinAmplitude, SourcePeaks, MinQ, MaxQ);
    SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
  }
  cout << "---- SX3 energy calibration (end) -----------------------------" << endl;

  return; 
}


///////////////////////////////////////////////////////////////////////////////////
// This is the third step in the SX3 calibration procedure.  It fits the total 
// front 'charge' to the calibrated back charge for each front strip  combination.
///////////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_CalibrateFrontBack(float MaxQDiff, int Cuts, float *EnergyCuts, float Uncertainty, int InitDet, int FinalDet)
{  
  cout << "---- SX3 front-back matching (begin) --------------------------" << endl;
  char GoodParameters;
  float InitMaxQDiff = MaxQDiff;

  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;

  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d>=InitDet && d<=FinalDet)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }

  // I've seen it's better to get only full hits when working on calibrations and 
  // to reset the calibration coefficients to their original value.
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    SX3[Det_ID]->ForceFullHits = 1;
    for (int strip=0; strip<4; strip++) {
      SX3[Det_ID]->mD[strip] = 1;
      SX3[Det_ID]->mU[strip] = 1;
      SX3[Det_ID]->oF[strip] = 0;
    }   
  }
  
  cout << "> Getting number of events in each detector ..." << endl;

  ExtractData("normal");

  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
    // Arrays to store qu, qd and qb are created. Their dimension was determined in the first
    // time ExtractData() was used in this method.
    SX3[Det_ID]->CreateArrays();
    
    cout << "> Saving data in the internal arrays of the detector ..." << endl;
    ExtractData("normal");

    for (int strip=0; strip<4; strip++) {
      GoodParameters = 'n';
      // First we give it a try with the default parameters.    
      MaxQDiff = InitMaxQDiff; 
      SX3[Det_ID]->GetFrontCoeffs(strip, MaxQDiff, Cuts, EnergyCuts, Uncertainty);
      cout << "> Are these good values for the parameters (y/n): ";
      cin >> GoodParameters;
      
      while(GoodParameters=='n'){
	cout << "> Type max. value of |Qb-Qf| for strip=" << strip << " (~500): ";
	cin >> MaxQDiff;
	SX3[Det_ID]->GetFrontCoeffs(strip, MaxQDiff, Cuts, EnergyCuts, Uncertainty);
	cout << "> Are these good values for the parameters (y/n): ";
	cin >> GoodParameters;
      }
      SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
    }    
    SX3[Det_ID]->DeleteArrays();
  }

  cout << "---- SX3 front-back matching (end) ----------------------------" << endl;

  return; 
}



///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_CalibratePosition(float MinE, float MaxE, int InitDet, int FinalDet)
{ 
  cout << "---- SX3 position calibration (start) --------------------------" << endl;
  char GoodParameters;
  float InitMaxE = MaxE, InitMinE = MinE;

  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;

  // I've seen it's better to get only full hits when working on calibrations and 
  // to reset the calibration coefficients to their original value.
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    SX3[Det_ID]->ForceFullHits = 1;
    for (int strip=0; strip<4; strip++) {
      for (int order=0; order<=SX3[Det_ID]->PosPolyOrder; order++) {
	SX3[Det_ID]->Cd[strip][order] = 1;
	SX3[Det_ID]->Cu[strip][order] = 1;
      }
    }
  }
 
  // Sort the data to fill the PvE histograms and count how many events were counted
  // by each detector.
  cout << "> Getting number of events in each detector ..." << endl;
  ExtractData("border");
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    cout << "> *** Working on detector " << SX3[Det_ID]->Name << " ***" << endl;
    // Pass the peak finding parameters provided by the user to each active SX3 det.
    SX3[Det_ID]->SetPeakFindingParameters(NPeaks, PeakSigma, PeakMinAmplitude, NormalChi2Value);
    // Arrays to store qu, qd and qb are created. Their dimension was determined in the first
    // time ExtractData() was used in this method.
    SX3[Det_ID]->CreateArrays();
    // Sort the data again now to fill the created arrays.
    ExtractData("border");
    for (int strip=0; strip<4; strip++) {
      GoodParameters = 'n';
      // First we give it a try with the default values of the parameters.
      MinE = InitMinE;  MaxE = InitMaxE;  
      SX3[Det_ID]->GetPositionPoly(strip, MinE, MaxE);
      cout << "> Are these good values for the parameters (y/n): ";
      cin >> GoodParameters;
      if (GoodParameters=='n') {
	SX3[Det_ID]->Cd[strip][0] = -50;
	SX3[Det_ID]->Cu[strip][0] = -50;
	for (int order=1; order<=SX3[Det_ID]->PosPolyOrder; order++) {
	  SX3[Det_ID]->Cd[strip][order] = 0;
	  SX3[Det_ID]->Cu[strip][order] = 0;
	}        
      }
      else
	SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
    }
    // Since the data arrays could be pretty large they are deleted once they have
    // been used to keep the memory usage low.  The same is done for the HistPvE's.
    SX3[Det_ID]->DeleteArrays();
  }
  cout << "---- SX3 position calibration (end) ----------------------------" << endl;
  return; 
}


////////////////////////////////////////////////////////////////////////////
// Simply function that enables for all the SX3 detectors the back-segment
// position restriction.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ForceBchPosRange()
{  
  for (int d=0; d<NumSX3Dets; d++)
    SX3[d]->ForceBchPosRange = 1;
  return;  
}
  

////////////////////////////////////////////////////////////////////////////
// Finds the pulser offsets and saves them in a text file with the SX3
// calibration coefficients. It does this from InitDet to FinalDet. The
// first argument is and array with the voltages that were used. The method
// SetPeakFindingParameters() must be called before calling this method.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_GetPulserOffsets(float* Voltages, int InitDet, int FinalDet)
{ 
  cout << "---- SX3 pulser offsets (begin) --------------------------------" << endl;

  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;
 
  // Sort the data to fill the HistRawE histograms.
  ExtractData("raw");
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
    SX3[Det_ID]->GetPulserOffsets(NPeaks, PeakSigma, PeakMinAmplitude, NormalChi2Value, Voltages);
    SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
  }
  
  cout << "---- SX3 pulser offsets (end) ----------------------------------" << endl;
  return; 
}

// Similar as the function above but only for one channel.
void EventProcessor::SX3_GetPulserOffsetsOneChannel(float* Voltages, int Det_ID, int Chan)
{ 
  cout << "---- SX3 pulser offsets (begin) --------------------------------" << endl;

  // Sort the data to fill the HistRawE histograms.
  ExtractData("raw");
  cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
  SX3[Det_ID]->GetPulserOffsets(NPeaks, PeakSigma, PeakMinAmplitude, NormalChi2Value, Voltages, Chan);
  SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
  
  cout << "---- SX3 pulser offsets (end) ----------------------------------" << endl;
  return; 
}

  
////////////////////////////////////////////////////////////////////////////
// Border hits in SX3 detectors are not very likely to occur, however they
// are the pilars of the position calibration procedure since we know 
// exactly where they hit the detector. To get enough statistics one needs
// to sort a lot of data events making the position calibration slow. To 
// speed it up you can use the function below which filters the data and
// saves in a cloned TTree (with the same structure as 'DataTree') the 
// events that had at least one border hit.
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_MakeBorderHitsFile(string OutputFile, float Min_Energy, float Max_Energy)
{
  Int_t  DetID, DetCh, m, SavedDetID[NumSX3Dets], FiredSX3s;
  string DetType;
  bool DetWithBorderHit[NumSX3Dets], AlreadyFiredDet, FoundOneBorderHit, FillFlag;

  // Two signal objects are created identical to the ones created in the constructor.
  ASICS_Signal* BH_ASICS = new ASICS_Signal(12*30); 
  CAEN_Signal* BH_CAEN = new CAEN_Signal(500);

  Long64_t TotEntries = Chain->GetEntries();

  // Create a root file only for the border hit data.
  TFile* BH_File = new TFile(OutputFile.c_str(), "RECREATE");

  // Identical TTree to the one used in ConvertEVTtoROOT() but addressed to the 'BH' objects.
  //TTree* BH_Tree = new TTree("DataTree","DataTree");
  TTree* BH_Tree = Chain->CloneTree(0);  // Clone the Tree structure but not the entries.  
  BH_Tree->SetBranchAddress("ASICS.Nhits",&BH_ASICS->Nhits);
  BH_Tree->SetBranchAddress("ASICS.MB",BH_ASICS->MB);
  BH_Tree->SetBranchAddress("ASICS.Chip",BH_ASICS->Chip);
  BH_Tree->SetBranchAddress("ASICS.Chan",BH_ASICS->Chan);
  BH_Tree->SetBranchAddress("ASICS.Energy",BH_ASICS->Energy);
  BH_Tree->SetBranchAddress("ASICS.Time",BH_ASICS->Time);
  BH_Tree->SetBranchAddress("CAEN.Nhits",&BH_CAEN->Nhits);
  BH_Tree->SetBranchAddress("CAEN.Module",BH_CAEN->Module);
  BH_Tree->SetBranchAddress("CAEN.Channel",BH_CAEN->Channel);
  BH_Tree->SetBranchAddress("CAEN.Value",BH_CAEN->Value);

  cout << "> Getting border hits ..." << endl;
  cout << ">\tTotal entries: " << TotEntries << endl;

  for(Long64_t n=0; n<TotEntries; n++){
    if (n == TMath::Nint(0.10*TotEntries))  cout << ">\t10% through the data" << endl;
    if (n == TMath::Nint(0.25*TotEntries))  cout << ">\t25% through the data" << endl;
    if (n == TMath::Nint(0.50*TotEntries))  cout << ">\t50% through the data" << endl;
    if (n == TMath::Nint(0.75*TotEntries))  cout << ">\t75% through the data" << endl;
    if (n == TMath::Nint(0.90*TotEntries))  cout << ">\t90% through the data" << endl;

    Chain->GetEntry(n);
    if (ASICS_Nhits<2 || ASICS_Nhits>=ASICS->Max_Hits)
      continue;

    for (int d=0; d<NumSX3Dets; d++) {
      SX3[d]->Reset();
      DetWithBorderHit[d] = 0; // Assume there was no border hit for this detector.
      SavedDetID[d]=-1;
    }
    BH_ASICS->Reset();    
    BH_CAEN->Reset();
    FiredSX3s = 0;  

    // 1. Sort the data.
    for (Int_t h=0; h<ASICS_Nhits; h++) {
      DetID = DetCh = -1;
      DetType = "";
      Map->IdentifyASICS(ASICS_MB[h], ASICS_Chip[h], ASICS_Chan[h], DetType, DetID, DetCh);      
      // Verify that the MB, Chip and Chan combination exists.
      if (DetType=="SX3" && DetID>=0 && DetCh>=0
	  && ASICS_Energy[h]>Min_Energy && ASICS_Energy[h]<Max_Energy) {
	m = SX3[DetID]->mult;
	SX3[DetID]->Channel[m] = DetCh;
	SX3[DetID]->RawE[m] = ASICS_Energy[h];
	SX3[DetID]->Time[m] = ASICS_Time[h];	   
	// Increase the multiplicity of this detector.
	SX3[DetID]->mult++;

	// This few lines optimize the sorting time since it enables us to keep track of
	// which detectors actually fired in this event and then (in step 2) insetead of
	// looping over all the SX3 detectors we only loop over the ones that fired.
	AlreadyFiredDet = 0;
	for (int f=0; f<FiredSX3s; f++) {
	  if (SavedDetID[f]==DetID)
	    AlreadyFiredDet = 1;
	}
	if (!AlreadyFiredDet) {
	  SavedDetID[FiredSX3s] = DetID;
	  FiredSX3s++;
	}
      }
    }// end for(h)
  
    // 2. Once an event is fully read check the SX3 detectors that fired to see if there was a border hit.
    FoundOneBorderHit = 0;
    for (int f=0; f<FiredSX3s; f++) {
      DetID = SavedDetID[f];
      m = SX3[DetID]->mult;
      // We need at least mutiplicity 4 to have a full border hit: 2 neighboring back channels and 
      // two front signals from the same strip.
      if (m>3) {
	SX3[DetID]->GetEventQs("border");
	if (SX3[DetID]->EventType>0) {
	  DetWithBorderHit[DetID] = 1;
	  FoundOneBorderHit = 1;
	}
      }
    }

    // 3. Resort the ASICs data (if there was a border hit) and save it in the BH_Tree.
    FillFlag = 0;
    if (FoundOneBorderHit) {
      for (Int_t h=0; h<ASICS_Nhits; h++) {
	DetID = DetCh = -1;
	DetType = "";
	Map->IdentifyASICS(ASICS_MB[h],ASICS_Chip[h],ASICS_Chan[h], DetType, DetID, DetCh);
	if (DetType=="SX3" && DetID>=0 && DetCh>=0 && DetWithBorderHit[DetID]) {
	  if (ASICS_Energy[h]>SX3[DetID]->Threshold[DetCh] && ASICS_Energy[h]<Max_Energy) {
	    BH_ASICS->MB[BH_ASICS->Nhits] = ASICS_MB[h];
	    BH_ASICS->Chip[BH_ASICS->Nhits] = ASICS_Chip[h];
	    BH_ASICS->Chan[BH_ASICS->Nhits] = ASICS_Chan[h];
	    BH_ASICS->Energy[BH_ASICS->Nhits] = ASICS_Energy[h];
	    BH_ASICS->Time[BH_ASICS->Nhits] = ASICS_Time[h];
	    BH_ASICS->Nhits++;
	    FillFlag = 1;
	  }
	}
      }// end for(h)
    }
    
    if (FillFlag)
      BH_Tree->Fill(); 

  }// end for(n)

  BH_Tree->Write(); 
  //  BH_File->Close();
  cout << ">\tDone!\n>\tUseful entries = " << BH_Tree->GetEntries() << endl;
  cout << ">\tFile \"" << OutputFile  << "\" created." << endl;

  //  delete DetWithBorderHit;
  delete BH_ASICS;
  delete BH_CAEN;
  delete BH_Tree;
  delete BH_File;
  return;
}
  
  
////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_MakeFilteredFile(string FilteredFilesDir, int SelectedDet)
{
  Int_t  DetID, DetCh, m, SavedDetID[NumSX3Dets], FiredSX3s;
  string DetType, OutputFile, OriginalFileName, DetString;
  bool DetWithGoodData[NumSX3Dets], AlreadyFiredDet, FoundOneActiveDet, FillFlag;

  // Two signal objects are created identical to the ones created in the constructor.
  ASICS_Signal* New_ASICS = new ASICS_Signal(12*30); 
  CAEN_Signal* New_CAEN = new CAEN_Signal(500);

  Long64_t TotEntries = Chain->GetEntries();

  // Get the name of the root file with the original data and insert the '_SX3_%d'
  // string before the '.root' extension.
  OriginalFileName = ROOTFileName;
  OriginalFileName = OriginalFileName.substr(OriginalFileName.find_last_of("/")+1);
  DetString = Form("_SX3_%d", SelectedDet);
  OutputFile = FilteredFilesDir + OriginalFileName.insert(OriginalFileName.find(".root"), DetString);
  cout << "> Writing filtered file: " << OutputFile << endl;

  // Create a root file only for the filtered data.
  TFile* FilteredFile = new TFile(OutputFile.c_str(), "RECREATE");

  // Identical TTree to the one used in ConvertEVTtoROOT() but addressed to the 'New' objects.
  TTree* New_Tree = Chain->CloneTree(0);  // Clone the Tree structure but not the entries (hence the '0').  
  New_Tree->SetBranchAddress("ASICS.Nhits",&New_ASICS->Nhits);
  New_Tree->SetBranchAddress("ASICS.MB",New_ASICS->MB);
  New_Tree->SetBranchAddress("ASICS.Chip",New_ASICS->Chip);
  New_Tree->SetBranchAddress("ASICS.Chan",New_ASICS->Chan);
  New_Tree->SetBranchAddress("ASICS.Energy",New_ASICS->Energy);
  New_Tree->SetBranchAddress("ASICS.Time",New_ASICS->Time);
  New_Tree->SetBranchAddress("CAEN.Nhits",&New_CAEN->Nhits);
  New_Tree->SetBranchAddress("CAEN.Module",New_CAEN->Module);
  New_Tree->SetBranchAddress("CAEN.Channel",New_CAEN->Channel);
  New_Tree->SetBranchAddress("CAEN.Value",New_CAEN->Value);
  
  cout << "> Filtering data ..." << endl;
  cout << ">\tTotal entries: " << TotEntries << endl;

  for(Long64_t n=0; n<TotEntries; n++){
    if (n == TMath::Nint(0.10*TotEntries))  cout << ">\t10% through the data" << endl;
    if (n == TMath::Nint(0.25*TotEntries))  cout << ">\t25% through the data" << endl;
    if (n == TMath::Nint(0.50*TotEntries))  cout << ">\t50% through the data" << endl;
    if (n == TMath::Nint(0.75*TotEntries))  cout << ">\t75% through the data" << endl;
    if (n == TMath::Nint(0.90*TotEntries))  cout << ">\t90% through the data" << endl;
    if (n == TotEntries-1)                  cout << ">\tDone." << endl;

    /*if (n>TMath::Nint(0.1*TotEntries))
      continue;
    */
    Chain->GetEntry(n);

    if (ASICS_Nhits<2 || ASICS_Nhits>=ASICS->Max_Hits)
      continue;

    for (int d=0; d<NumSX3Dets; d++) {
      SX3[d]->Reset();
      DetWithGoodData[d] = 0; // Assume there was no hit in this detector.
      SavedDetID[d]=-1;
    }
    New_ASICS->Reset();    
    New_CAEN->Reset();
    FiredSX3s = 0;  

    // 1. Sort the data.
    for (Int_t h=0; h<ASICS_Nhits; h++) {
      DetID = DetCh = -1;
      DetType = "";
      Map->IdentifyASICS(ASICS_MB[h], ASICS_Chip[h], ASICS_Chan[h], DetType, DetID, DetCh);      
      // Verify that the MB, Chip and Chan combination exists.
      if (DetType=="SX3" && DetID>=0 && DetCh>=0)
	if(ASICS_Energy[h]>SX3[DetID]->Threshold[DetCh] && ASICS_Energy[h]<16000) {
	  m = SX3[DetID]->mult;
	  SX3[DetID]->Channel[m] = DetCh;
	  SX3[DetID]->RawE[m] = ASICS_Energy[h];
	  SX3[DetID]->Time[m] = ASICS_Time[h];	   
	  // Increase the multiplicity of this detector.
	  SX3[DetID]->mult++;
	  
	  // These few lines optimize the sorting time since it enables us to keep track of
	  // which detectors actually fired in this event and then (in step 2) insetead of
	  // looping over all the SX3 detectors we only loop over the ones that fired.
	  AlreadyFiredDet = 0;
	  for (int f=0; f<FiredSX3s; f++) {
	    if (SavedDetID[f]==DetID)
	      AlreadyFiredDet = 1;
	  }
	  if (!AlreadyFiredDet) {
	    SavedDetID[FiredSX3s] = DetID;
	    FiredSX3s++;
	  }
	}
    }// end for(h)
  
    // 2. Once an event is fully read loop over the SX3 detectors that fired to see whether
    // the SelecetDet is one of them.
    FoundOneActiveDet = 0;
    for (int f=0; f<FiredSX3s; f++) {
      DetID = SavedDetID[f];
      if (DetID==SelectedDet) {
	m = SX3[DetID]->mult;
	// We need at least mutiplicity 2 to have a good hit: 1 back channel and 1 or 2 front signals
	// from the same strip.
	if (m>1) {
	  SX3[DetID]->GetEventQs();
	  if (SX3[DetID]->EventType>0) {
	    DetWithGoodData[DetID] = 1;
	    FoundOneActiveDet = 1;
	  }
	}
      }
    }
    
    
    // 3. Resort the ASICs data (if there was a an active det. that fired) and save it in the New_Tree.
    FillFlag = 0;
    if (FoundOneActiveDet) {
      for (Int_t h=0; h<ASICS_Nhits; h++) {
	DetID = DetCh = -1;
	DetType = "";
	Map->IdentifyASICS(ASICS_MB[h],ASICS_Chip[h],ASICS_Chan[h], DetType, DetID, DetCh);
	if (DetType=="SX3" && DetID==SelectedDet && DetCh>=0 && DetWithGoodData[DetID]) {
	  if (ASICS_Energy[h]>SX3[DetID]->Threshold[DetCh] && ASICS_Energy[h]<16000) {
	    New_ASICS->MB[New_ASICS->Nhits] = ASICS_MB[h];
	    New_ASICS->Chip[New_ASICS->Nhits] = ASICS_Chip[h];
	    New_ASICS->Chan[New_ASICS->Nhits] = ASICS_Chan[h];
	    New_ASICS->Energy[New_ASICS->Nhits] = ASICS_Energy[h];
	    New_ASICS->Time[New_ASICS->Nhits] = ASICS_Time[h];
	    New_ASICS->Nhits++;
	    FillFlag = 1;
	  }
	}
      }// end for(h)
    }
    
    if (FillFlag)
      New_Tree->Fill(); 

  }// end for(n)

  New_Tree->Write(); 
  cout << "> New TTree written.  Useful entries = " << New_Tree->GetEntries() << endl;
  FilteredFile->Close();
  cout << "> FilteredFile closed " << endl;
  delete FilteredFile;
  cout << "> FilteredFile deleted " << endl;
  delete New_ASICS;
  cout << "> New_ASICS deleted " << endl;
  delete New_CAEN;
  cout << "> New_CAEN deleted " << endl;
  
  //New_Tree->Reset();
  New_Tree = 0;
  cout << "> New_Tree set to  " << New_Tree << endl;
  /*  delete FilteredFile;
  cout << ">\tFilteredFile deleted " << endl;
  */
  //  FilteredFile->Close();  // Close() makes this function crash. I don't know why.

  //  cout << ">\tFile \"" << OutputFile  << "\" created." << endl;

  return;
}
  

////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_SetThresholds(float DisplayLimit, int InitDet, int FinalDet)
{ 
  cout << "---- Set SX3 thresholds (begin) --------------------------------" << endl;

  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;
 
  // Sort the data to fill the HistRawE histograms.
  ExtractData("raw");
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    cout << "> Working on detector " << SX3[Det_ID]->Name << endl;
    SX3[Det_ID]->SetThresholds(DisplayLimit);
    SX3[Det_ID]->WriteCoefficients(ParamDirectory+Form("SX3_%d.coeffs",Det_ID));
  }
  
  cout << "---- Set SX3 thresholds (end) ----------------------------------" << endl;
  return; 
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ShowGraphsQbQf(int DetID, int Cuts, float *EnergyCuts, float Uncertainty)
{
  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d==DetID)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }
  SX3[DetID]->GetEnergyFrom = 'b';
  // I've seen it's better to get only full hits when working on calibrations.
  // SX3[DetID]->ForceFullHits = 1;

  // Sort the data to fill the PvE histograms and count how many events were counted
  // by each detector.
  cout << "> Getting number of events in detector " << DetID << " ..." << endl;
  ExtractData();
  
  // Arrays to store qu, qd and qb are created. Their dimension was determined in the first
  // time ExtractData() was used in this method.
  SX3[DetID]->CreateArrays();

  ExtractData();

  SX3[DetID]->ShowGraphsQbQf(Cuts,EnergyCuts,Uncertainty);

  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ShowGraphsQdQu(int DetID, int Cuts, float *EnergyCuts, float Uncertainty)
{
  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d==DetID)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }
  SX3[DetID]->GetEnergyFrom = 'b';
  // I've seen it's better to get only full hits when working on calibrations.
  // SX3[DetID]->ForceFullHits = 1;

  // Sort the data to fill the PvE histograms and count how many events were counted
  // by each detector.
  cout << "> Getting number of events in detector " << DetID << " ..." << endl;
  ExtractData();
  
  // Arrays to store qu, qd and qb are created. Their dimension was determined in the first
  // time ExtractData() was used in this method.
  SX3[DetID]->CreateArrays();

  ExtractData();

  SX3[DetID]->ShowGraphsQdQu(Cuts,EnergyCuts,Uncertainty);

  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ShowHistQb(float MinQ, float MaxQ, int InitDet, int FinalDet)
{
  TCanvas* Display;

  // We turn off all the SX3 detectors except the ones selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d>=InitDet && d<=FinalDet)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }

  ExtractData();

  Display = new TCanvas("SX3HistQb","SX3 Qb histograms",0,0,1400,800);
  Display->Divide(2);
  for (int d=InitDet; d<=FinalDet; d++) {
    Display->cd(1);
    for (int bch=3; bch>=0; bch--) {
      SX3[d]->HistQb[bch]->SetLineColor(bch+1);
      SX3[d]->HistQb[bch]->SetAxisRange(MinQ, MaxQ);
      if (bch==3)
	SX3[d]->HistQb[bch]->Draw();
      else 
	SX3[d]->HistQb[bch]->Draw("same");
    }
    Display->cd(2);
    SX3[d]->HistQb[4]->SetAxisRange(MinQ, MaxQ);
    SX3[d]->HistQb[4]->Draw();
    Display->Update();
    cout << "> Showing HistQb of detector " << SX3[d]->Name << "." << endl;
    cout << "> Double-click on the canvas to continue." << endl;
    Display->WaitPrimitive();
  }
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ShowPosSummary()
{
  TCanvas* Display;
  TLine** DetDivision;
  
  Display = new TCanvas("SX3PosSumm","SX3 pos. summary",1400,800);
  DetDivision = new TLine*[NumSX3Dets];
  for (int d=0; d<NumSX3Dets; d++) {
    DetDivision[d] = new TLine(4*d,0,4*d,8);
    DetDivision[d]->SetLineStyle(1);
  }
  
  HistSX3PosSummary->Draw("colz");
  for (int d=1; d<NumSX3Dets-1; d++)
    DetDivision[d]->Draw();
  return;
}



////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_ShowPvE(int Det, char GetEnergyFrom, int HistEBins, float HistMinE, float HistMaxE)
{
  TCanvas* Display;
  SX3[Det]->CreateHistPvE(HistEBins, HistMinE, HistMaxE);
  // Here you select from which side the energy of the SX3 detector will be 
  // extracted (front-'f' or back-'b')
  SX3[Det]->GetEnergyFrom = GetEnergyFrom;
  
  // We turn off all the SX3 detectors except the one selected by the user.
  for (int d=0; d<NumSX3Dets; d++) {
    if (d==Det)
      SX3[d]->IsOn = 1;
    else
      SX3[d]->IsOn = 0;
  }
    
  ExtractData();
  
  Display = new TCanvas("SX3HistPvE",Form("%s HistPvE",SX3[Det]->Name.c_str()),0,0,1400,800);
  Display->Divide(4);
  for (int stp=0; stp<4; stp++) {
    Display->cd(stp+1);
    SX3[Det]->HistPvE[stp][0]->SetTitle(Form("%s | Front strip %d", SX3[Det]->Name.c_str(),stp));
    SX3[Det]->HistPvE[stp][0]->SetMarkerColor(1);
    SX3[Det]->HistPvE[stp][0]->SetMarkerStyle(7);
    SX3[Det]->HistPvE[stp][0]->GetYaxis()->SetTitle(Form("Q_{%c} [arb. units]",GetEnergyFrom));
    SX3[Det]->HistPvE[stp][0]->DrawClone();
    SX3[Det]->HistPvEEdge[stp][0]->SetMarkerColor(1);
    SX3[Det]->HistPvEEdge[stp][0]->SetMarkerStyle(5);
    SX3[Det]->HistPvEEdge[stp][0]->DrawClone("same");
    for (int bch=1; bch<4; bch++) {
      SX3[Det]->HistPvE[stp][bch]->SetMarkerColor(bch+1);
      SX3[Det]->HistPvE[stp][bch]->SetMarkerStyle(7);
      SX3[Det]->HistPvE[stp][bch]->DrawClone("same");
      SX3[Det]->HistPvEEdge[stp][bch]->SetMarkerColor(bch+1);
      SX3[Det]->HistPvEEdge[stp][bch]->SetMarkerStyle(5);
      SX3[Det]->HistPvEEdge[stp][bch]->DrawClone("same");
    }
  }
  Display->Update();
  cout << "> Showing overlayed HistPvE from detector " << SX3[Det]->Name << "." << endl;
  cout << "> Double-click on the canvas to delete it." << endl;
  Display->WaitPrimitive();

  SX3[Det]->DeleteHistPvE();
  delete Display;
  return;
}


////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
void EventProcessor::SX3_TestFrontBackMatch(char GetEnergyFrom, int BinsQ, float MinQ, float MaxQ, int InitDet, int FinalDet)
{ 
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    if (SX3[Det_ID]->HistDiffQCreated)
      SX3[Det_ID]->DeleteHistDiffQ();
    SX3[Det_ID]->CreateHistDiffQ(BinsQ, MinQ, MaxQ);
    SX3[Det_ID]->GetEnergyFrom = GetEnergyFrom;
  }
  
  ExtractData();
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    for (int stp=0; stp<4; stp++) {
      SX3[Det_ID]->ShowHistDiffQ(stp);
    }
  }
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++)
    SX3[Det_ID]->DeleteHistDiffQ();

  return;
}
    

////////////////////////////////////////////////////////////////////////////
//  
////////////////////////////////////////////////////////////////////////////
  void EventProcessor::SX3_TestPosCal(char GetEnergyFrom, float MinE, float MaxE, int InitDet, int FinalDet)
{ 
  if (FinalDet==-1)
    FinalDet = NumSX3Dets-1;
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    if (SX3[Det_ID]->HistPvECreated)
      SX3[Det_ID]->DeleteHistPvE();
    SX3[Det_ID]->CreateHistPvE(500,MinE,MaxE);
    SX3[Det_ID]->GetEnergyFrom = GetEnergyFrom;
  }
  
  ExtractData();
  
  for (int Det_ID=InitDet; Det_ID<=FinalDet; Det_ID++) {
    for (int stp=0; stp<4; stp++) {
      SX3[Det_ID]->ShowHistPvE(stp, MinE, MaxE, Form("stp%d",stp));
    }
  }

  return;
}
    

////////////////////////////////////////////////////////////////////////////
// Simple function that writes the data in the physical event TTree to a
// ROOT file (prevoiusly created in CreatePhysicalEventsFile()).
////////////////////////////////////////////////////////////////////////////
void EventProcessor::WritePhysicalEventsFile()
{
  PE_ROOT_File->Write();
  PE_ROOT_File->Close();
  cout << "> File with physical events \"" << PE_ROOT_File->GetName() << "\" has been written." << endl;
  return;
}



