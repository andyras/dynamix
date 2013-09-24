#ifndef __RHS__
#define __RHS__


#include <iostream>
#include <vector>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#include "params.hpp"
#include "numerical.hpp"

void updateHamiltonian(PARAMETERS * p, realtype t);

int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data);

void buildFDD(struct PARAMETERS * p, N_Vector y, std::vector<double> & fdd);

double b13(double bm, double ekin, double ne, double K1, double K2, double K3, double X);

int RHS_DM_RTA(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data);

#endif
