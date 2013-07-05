#ifndef __PARAMS__
#define __PARAMS__

#include <vector>

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
 int bridge_on;
 int scale_bubr;
 int scale_brqd;
 int scale_buqd;
 int parabolicCoupling;

 realtype kBandEdge;
 realtype kBandTop;
 realtype tout;

 realtype gamma1;
 realtype gamma2;

 std::vector<realtype> energies;
 std::vector<realtype> Vbridge;
 std::vector<realtype> Vnobridge;
 std::vector<realtype> H;
 std::vector<realtype> times;
};

#endif
