/*******************************************************************
Header file: EnergyLoss.hpp

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

class EnergyLoss{
public:
  EnergyLoss(string Eloss_file, float IonMass=0);
  double GetEnergyLoss(float initial_energy, float distance);
  void GetEvDCurve(float InitEne, float FinalDepth, int steps);
  double GetInitialEnergy(float FinalEnergy, float PathLength, float StepSize);
  double GetFinalEnergy(float InitialEnergy, float PathLength, float StepSize);  
  double GetTimeOfFlight(float InitialEnergy, float PathLength, float StepSize);
  bool GoodELossFile;
  TGraph* EvD;

private:
  double c;
  double* IonEnergy;
  double IonMass;
  double* dEdx_e;
  double* dEdx_n;
  int points;
  int last_point;
  bool Energy_in_range;
};


