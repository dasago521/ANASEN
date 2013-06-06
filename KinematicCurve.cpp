/*******************************************************************
Code: KinematicCurve.cpp

Description: This class gives the kinematic curves as a TGraph
  object. The masses for the reaction and the energy of excitation
  must be set before getting the curve.

Author: Daniel Santiago-Gonzalez
2013-March
*******************************************************************/

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TGraph.h>
#include "KinematicCurve.hpp"

using namespace std;

KinematicCurve::KinematicCurve()
{

}

void KinematicCurve::SetReaction(float Kb, float mb, float mt, float mh, float ml, float Eexc)
{
  // All masses must be in MeV/c^2 and energies in MeV.
  this->Kb = Kb;     // beam kinetic energy.
  this->mb = mb;     // beam mass
  this->mt = mt;     // target mass
  this->mh = mh;     // heavy ion mass
  this->ml = ml;     // light ion mass
  this->Eexc = Eexc; // excitation energy of heavy ion.
  return;
}


void KinematicCurve::GetCurve(int Points)
{
  Float_t A=0, B=0, C=0;
  Float_t Kl=0, sqrtKl_plus=0, sqrtKl_minus=0, Q0=0, theta_deg=0;
  Float_t pi = TMath::Pi();
  //  TGraph* KC = new TGraph(Points);

  // Initialize the curve.
  for(Int_t p=0; p<Points; p++){
    theta_deg = 180.0*p/Points;
    this->SetPoint(p, theta_deg, 0.0);
  }
  
  // Fill the points of the kinematic curve
  for(Int_t p=0; p<Points; p++){
    Q0 = mb + mt - mh - ml;
    theta_deg = 180.0*p/Points;
      
    A = 1 + ml/mh;
    B = -2*sqrt(mb*ml*Kb)*cos(theta_deg*pi/180.0)/mh;
    C = -Kb*(1-mb/mh) - Q0 + Eexc;
    
    if(B*B-4*A*C>=0){
      sqrtKl_plus  = (-B+sqrt(B*B-4*A*C))/(2*A);
      sqrtKl_minus = (-B-sqrt(B*B-4*A*C))/(2*A);
      Kl = pow(sqrtKl_plus,2);
      this->SetPoint(p, theta_deg, Kl);
    }
    else
      break;
  }// end for(p)    

  return;
}
