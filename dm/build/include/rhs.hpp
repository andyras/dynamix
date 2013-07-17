#ifndef __RHS__
#define __RHS__


#include <iostream>
#include <vector>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

#include "params.hpp"

int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data);

int RHS_DM_RTA(realtype t, N_Vector y, N_Vector ydot, void * data);

#endif