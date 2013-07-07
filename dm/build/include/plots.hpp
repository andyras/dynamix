#ifndef __PLOTS__
#define __PLOTS

#define DEBUG_PLOT

#include <iostream>
#include <fstream>
#include <map>

#include "params.hpp"

void makePlots(std::map<std::string, bool> &outs, struct PARAMETERS * p);

void plotPopulations(char * fileName, struct PARAMETERS * p);

void plotKProbs(char * fileName, struct PARAMETERS * p);

void plotCProbs(struct PARAMETERS * p);

#endif
