#pragma once

#include <iostream>
#include <numeric>
#include <vector>

#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#ifdef __USE_MKL__
#include <mkl.h>
#endif
extern "C" void sgemm_( char *, char *, int *, int *, int *, float *, float *,
  int *, float *, int *, float *, float *, int *);
extern "C" void dsymv_(const char *, int *, double *, double *, int *, double *, int *,
  double *, double *, int *);
extern "C" void dsymm_(char *, char *, int *, int *, double *, double *, int *,
  double *, int *, double *, double *, int *);
extern "C" void dgemv_(char *, int *, int *, double *, double *, int *, double *,
  int *, double *, double *, int *);
extern "C" double ddot_(int *, double *, int *, double *, int *);

#include "params.hpp"
#include "numerical.hpp"
#include "dynamix.hpp"
#include "omp.h"

int RHS_WFN(realtype t, N_Vector y, N_Vector ydot, void * data);

#ifdef __USE_MKL__
int RHS_WFN_SPARSE(realtype t, N_Vector y, N_Vector ydot, void * data);
#endif

void RELAX_KINETIC(int bandFlag, realtype * yp, realtype * ydotp, Params * p);

int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data);

double findDynamicMu(double pop, double T, int bandFlag, Params * p);

double FDDBinarySearch(double lower, double upper, double T, double n,
    int bandFlag, Params * p);

double FDDSum(double mu, double T, int bandFlag, Params * p);

void FDD(double mu, double T, double * fdd, double * E, int N, double P);