#ifndef __PLOTS__
#define __PLOTS

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>

#include "output.hpp"
#include "params.hpp"

void makePlots(Params * p);

void plotPopulations(char * fileName, Params * p);

void plotKProbs(char * fileName, Params * p);

void plotCProbs(Params * p);

void plotMuFromPops(char * fileName, Params * p);

void plotDMt_z(char * fileName, Params * p);

void plotKProbsMovie(char * fileName, Params * p);

void plotCProbsMovie(char * fileName, Params * p);

#endif
