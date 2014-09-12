#pragma once

#include <iostream>
#include <vector>
#include <numeric>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>

#include "params.hpp"
#include "numerical.hpp"
#include "dynamix.hpp"
#include "omp.h"

int RHS_WFN(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_WFN_SPARSE(realtype t, N_Vector y, N_Vector ydot, void * data);

void RELAX_KINETIC(int bandFlag, realtype * yp, realtype * ydotp, Params * p);

int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data);

double findDynamicMu(double pop, double T, int bandFlag, Params * p);

double FDDBinarySearch(double lower, double upper, double T, double n,
    int bandFlag, Params * p);

double FDDSum(double mu, double T, int bandFlag, Params * p);

void FDD(double mu, double T, double * fdd, double * E, int N, double P);