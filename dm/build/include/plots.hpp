#ifndef __PLOTS__
#define __PLOTS

#include <iostream>
#include <fstream>
#include <map>

#include "params.hpp"

void makePlots(std::map<const std::string, bool> &outs, struct PARAMETERS * p);

void plotPopulations(char * fileName, struct PARAMETERS * p);

void plotKProbs(char * fileName, struct PARAMETERS * p);

void plotCProbs(struct PARAMETERS * p);

void plotDMt_z(char * fileName, struct PARAMETERS * p);

void plotKProbsMovie(char * fileName, struct PARAMETERS * p);

#endif
