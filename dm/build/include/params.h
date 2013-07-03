#ifndef __PARAMS__
#define __PARAMS__

#include <vector>

// Struct of parameters
struct PARAMETERS {
 int Nk;
 int Nb;
 int Nc;
 int Ik;
 int Ib;
 int Ic;
 int NEQ;
 int NEQ2;
 int numOutputSteps;
 int bridge_on;
 int scale_bubr;
 int scale_brqd;
 int scale_buqd;
 int parabolicCoupling;

 realtype kBandEdge;
 realtype kBandTop;
 realtype tout;

 std::vector<realtype> energies;
 std::vector<realtype> Vbridge;
 std::vector<realtype> Vnobridge;
 std::vector<realtype> times;
};

#endif
