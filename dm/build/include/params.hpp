#ifndef __PARAMS__
#define __PARAMS__

#include <vector>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#include "spline.hpp"

// Struct of parameters
struct PARAMETERS {
  int Nk;
  int Nc;
  int Nb;
  int Nl;
  int Ik;
  int Ic;
  int Ib;
  int Il;
  int NEQ;
  int NEQ2;
  int numOutputSteps;
  bool bridge_on;
  bool scale_bubr;
  bool scale_brqd;
  bool scale_buqd;
  bool parabolicCoupling;
  bool torsion;
  int torsionSite;

  realtype kBandEdge;
  realtype kBandTop;
  realtype kBandWidth;
  realtype tout;

  realtype gamma1;
  realtype gamma2;

  std::vector<realtype> energies;
  std::vector<realtype> Vbridge;
  std::vector<realtype> Vnobridge;
  std::vector<realtype> H;
  std::vector<realtype> times;

  Spline * torsionV;
};

#endif
