/*******************************************************************
Code: EnergyLoss.cpp

Description: Simple class that calculates the energy loss of an ion
  in a gas target. The input is a three-column text file with the
  energy of the ion and the electrical and nuclear stopping powers
  (dE/dx) of the target for that ion energy. The units are assumed
  to be MeV and MeV/mm for the ion's energy and the stopping powers,
  respectively. This information (E and dE/dx) can be obtained from
  SRIM. Also, it is assumed that the first line are three strings
  describing the columns.

Author: Daniel Santiago-Gonzalez
2012-Sep
*******************************************************************/

#include <iostream>
#include <fstream>
#include <string.h>
#include <TGraph.h>

#include "EnergyLoss.hpp"

using namespace std;

// Constructor.
EnergyLoss::EnergyLoss(string Eloss_file, float IonMass/*MeV/c^2*/)
{
  double IonEnergy, dEdx_e, dEdx_n;
  string aux;
  ifstream Read(Eloss_file.c_str());
  last_point = 0;
  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: File " << Eloss_file << " was not found." << endl;
    GoodELossFile = 0;
  } 
  else {
    GoodELossFile = 1;        
    // The first line has three strings (columns' description).
    Read >> aux >> aux >> aux;
    // Cout the number of points.
    points = 0;
    do{
      Read >> IonEnergy >> dEdx_e >> dEdx_n;
      points++;
    }while(!Read.eof());
    Read.close();
    // Create the arrays depending on the number rows in the file.
    this->IonEnergy = new double[points];
    this->dEdx_e = new double[points];
    this->dEdx_n = new double[points];    
    // Go to the begining of the file and read it again to now save the info in the
    // newly created arrays.
    Read.open(Eloss_file.c_str());
    Read >> aux >> aux >> aux;
    for(int p=0; p<points; p++){
      Read >> IonEnergy >> dEdx_e >> dEdx_n;
      this->IonEnergy[p] = IonEnergy;
      this->dEdx_e[p] = dEdx_e;
      this->dEdx_n[p] = dEdx_n;
    }    
    Energy_in_range = 1;
    this->IonMass = IonMass;  // In MeV/c^2
    c = 29.9792458;           // Speed of light in cm/ns.
    EvD = new TGraph();
  }
}

// Get the energy loss of the gas for an initial ion energy and a certain
// distance through the target.
double EnergyLoss::GetEnergyLoss(float energy /*MeV*/, float distance /*cm*/)
{
  int i = -1;
  // Look for two points for which the initial energy lays in between.
  // This for-loop should find the points unless there was a big jump from
  // the energy used in the last point and the energy used now.
  for(int p=last_point-1; p<points-1; p++){
    if (last_point>=0)
      if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	i = p+1;
	last_point = p;
	break;
      }
  }
  // It is probable that if the point wasn't found could have been because of
  // a big jump in the energy (see above), so we need to look in the remaining
  // points.
  if (i==-1) {
    for(int p=0; p<last_point-1; p++){
      if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	i = p+1;
	last_point = p;
	break;
      }
    }
  }
  // If after those two loops i is still -1 it means the energy was out of range.
  if(i==-1){
    cout << "*** EnergyLoss Error: energy not within range: " << energy << endl;
    Energy_in_range = 0;
    return 0;
  }
  // If the initial energy is within the range of the function, get the stopping power
  // for the initial energy.  
  double E1 = IonEnergy[i-1];        double E2 = IonEnergy[i];
  double dEdx_e1 = dEdx_e[i-1];      double dEdx_e2 = dEdx_e[i];  
  double dEdx_n1 = dEdx_n[i-1];      double dEdx_n2 = dEdx_n[i];  
  // Interpolating the electric stopping power (from point 1 to 'e').
  double dEdx_e1e = dEdx_e1 + (energy - E1)*(dEdx_e2 - dEdx_e1)/(E2 - E1);
  // Interpolating the nuclear stopping power (usually negligable).
  double dEdx_n1e = dEdx_n1 + (energy - E1)*(dEdx_n2 - dEdx_n1)/(E2 - E1);
  // The stopping power units are in MeV/mm so we multiply by 10 to convert to MeV/cm.
  return ( (dEdx_e1e+dEdx_n1e)*10*distance );  
}


void EnergyLoss::GetEvDCurve(float InitEne/*MeV*/, float FinalDepth/*cm*/, int steps)
{
  Double_t current_ene, current_depth;
  current_ene = InitEne;
  current_depth = 0;
  for(int s=0; s<steps; s++){
    current_ene -= GetEnergyLoss(current_ene,FinalDepth/steps);
    current_depth += FinalDepth/steps;
    EvD->SetPoint(s, current_ene, current_depth);
  }
}



double EnergyLoss::GetInitialEnergy(float FinalEnergy /*MeV*/, float PathLength /*cm*/, float StepSize/*cm*/)
{
  double Energy = FinalEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  last_point = 0;
  // The function starts by assuming FinalEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {
    Energy = Energy + GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range)
      break;
  } 
  Energy = Energy + GetEnergyLoss(Energy,PathLength-Steps*StepSize);
  if (!Energy_in_range)
    Energy = -1000; // Return an unrealistic value.
  
  //  cout << "d: K_lf=" << FinalEnergy << "  K_lr=" << Energy << "  l=" << PathLength <<endl;

  return Energy;
}


double EnergyLoss::GetFinalEnergy(float InitialEnergy /*MeV*/, float PathLength /*cm*/, float StepSize/*cm*/)
{
  double Energy = InitialEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  // The function starts by assuming InitialEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;
  for (int s=0; s<Steps; s++) {
    Energy = Energy - GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range)
      break;
  }  
  Energy = Energy - GetEnergyLoss(Energy,PathLength-Steps*StepSize);
  if (!Energy_in_range) 
    Energy = -1000;
  //  cout << "O: K_bw=" << InitialEnergy << "  K_br=" << Energy << "  l=" << PathLength <<endl;
  return Energy;
}


////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's time of flight in ns.
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetTimeOfFlight(float InitialEnergy, float PathLength, float StepSize)
{
  double TOF = 0;
  double Kn = InitialEnergy;
  int Steps = PathLength/StepSize;
  if (IonMass==0)
    cout << "*** EnergyLoss Error: Time of flight cannot be calculated for IonMass = 0." << endl;
  else {
    // The TOF is proportional to 1/sqrt(Kn). After the sum, TOF will be multiplied by
    // the proportionality factor.
    for (int n=0; n<Steps; n++) {
      TOF += 1/sqrt(Kn);               // DeltaT going from point n to n+1.
      Kn -= GetEnergyLoss(Kn, StepSize); // After the TOF is added the kinetic energy at point n+1 is calc.
    }
    TOF *= sqrt(IonMass/2)*StepSize/c;
  }
  return TOF;
}
