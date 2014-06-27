#ifndef __DYNAMIX__
#define __DYNAMIX__

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <numeric>
#include <complex>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>
#include <map>
#include <omp.h>

#include "libdynamix_input_parser.hpp"
#include "libdynamix_outputs.hpp"
#include "output.hpp"
#include "numerical.hpp"
#include "params.hpp"
#include "userdata.hpp"
#include "rhs.hpp"
#include "plots.hpp"
#include "constants.hpp"
#include "conversions.hpp"
#include "analytic.hpp"

void updateDM(N_Vector dm, int timeStep, Params * p);

void updateWfn(N_Vector wfn, int timeStep, Params * p);

int bandStartIdx(int bandFlag, Params * p);

int bandEndIdx(int bandFlag, Params * p);

int bandNumStates(int bandFlag, Params * p);

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, Params * p);

void updateHamiltonian(Params * p, realtype t);

void initialize(Params * p);

void initHamiltonian(Params * p);

void initWavefunction(Params * p);

#endif