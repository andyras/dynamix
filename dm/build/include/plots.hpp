#ifndef __PLOTS__
#define __PLOTS

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>

#include "output.hpp"
#include "params.hpp"

void makePlots(std::map<const std::string, bool> &outs, struct Params * p);

void plotPopulations(char * fileName, struct Params * p);

void plotKProbs(char * fileName, struct Params * p);

void plotCProbs(struct Params * p);

void plotMuFromPops(char * fileName, struct Params * p);

void plotDMt_z(char * fileName, struct Params * p);

void plotKProbsMovie(char * fileName, struct Params * p);

void plotCProbsMovie(char * fileName, struct Params * p);

#endif
