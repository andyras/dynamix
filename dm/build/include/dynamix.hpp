#pragma once

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <vector>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <nvector/nvector_serial.h>

#ifdef __USE_MKL__
  #include <mkl.h>
#endif

#include "analytic.hpp"
#include "constants.hpp"
#include "conversions.hpp"
#include "numerical.hpp"
#include "output.hpp"
#include "params.hpp"
#include "parser.hpp"
#include "plots.hpp"
#include "rhs.hpp"

void updateDM(N_Vector dm, int timeStep, Params * p);

void updateWfn(N_Vector wfn, int timeStep, Params * p);

int bandStartIdx(int bandFlag, Params * p);

int bandEndIdx(int bandFlag, Params * p);

int bandNumStates(int bandFlag, Params * p);

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, Params * p);

void updateHamiltonian(Params * p, realtype t);

void initialize(Params * p, const bool readFiles = true);

void initHamiltonian(Params * p, const bool readFiles = true);

void initWavefunction(Params * p);