/*******************************************************************
Code: KinematicCurve.hpp

Description: This class gives the kinematic curves as a TGraph
  object. The masses for the reaction and the energy of excitation
  must be set before getting the curve.

Author: Daniel Santiago-Gonzalez
2013-March
*******************************************************************/

class KinematicCurve: public TGraph{
public:
  KinematicCurve();
  void SetReaction(float Kb, float mb, float mt, float mh, float ml, float Eexc=0);
  void GetCurve(int Points);

private:
  float mb, mt, mh, ml, Eexc, Kb;

  // Need to uncomment this when using this class directly in CINT.
  //  ClassDef(KinematicCurve,1);
  // I've seen that when this class is included in another class or code
  // and ACLiC is used to compile such code you don't need to uncomment the
  // line above.
};
