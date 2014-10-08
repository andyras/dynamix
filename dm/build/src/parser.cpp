#include "parser.hpp"

// #define DEBUG
// #define DEBUG_INPUT_PARSER

void assignOutputs(const char * inputFile, std::map<const std::string, bool> &outs,
    Params * p) {
  std::string line;
  std::ifstream input;

  input.open(inputFile, std::ios::in);

#ifdef DEBUG_INPUT_PARSER
  std::cout << "Parsing file " << inputFile << "\n";
#endif

  if (input.good() == false) {
    std::cerr << "[assignOutputs] ERROR: bad input file " << inputFile << "\n";
  }

  // get first line
  getline(input, line);

  // skip lines until header specifying start of outputs
  while (getline(input, line)) {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "line is " << line << "\n";
#endif
    if (line != "[[Output Files]]") {
      continue;
    }
    else {
      break;
    }
  }

  // read inputs until [[End]] line or EOF
  while (getline(input, line)) {
    if (line != "[[End]]") {

      // skip comments
      if (line.substr(0,1) == "#") {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Skipping comment line: " << line << "\n";
#endif
	continue;
      }

      // skip whitespace (space/tab) and blank lines
      else if ((line.find_first_not_of(' ') == std::string::npos)
	  || (line.find_first_not_of('\t') == std::string::npos)) {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Skipping blank/whitespace line\n";
#endif
	continue;
      }

      // turn on outputs which are in this section of the input
      else {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Creating output file:  " << line << "\n";
#endif
	outs.insert(std::pair<const std::string,bool>(line,true));
	if ((line.substr(line.length()-4, line.length()) != ".out")
	    && (line.substr(line.length()-4, line.length()) != ".plt")) {
	  std::cerr << "WARNING [" << __FUNCTION__ << "]: output file extension is not '.out' or '.plt'\n";
	}
      }
    }
  }

  // if creating plot files, make sure that the appropriate outputs are being created.
  if (isOutput(outs, "projections.plt")) {
    outs.insert(std::pair<const std::string,bool>("tkprob.out", true));
    outs.insert(std::pair<const std::string,bool>("tcprob.out", true));

    // create bridge output if bridge is on
    if (p->bridge_on) {
      outs.insert(std::pair<const std::string,bool>("tbprob.out", true));
    }
  }

  if (isOutput(outs, "cprobs.plt")) {
    outs.insert(std::pair<const std::string,bool>("cprobs.out", true));
  }

  if ((isOutput(outs, "kprobs.plt")) || (isOutput(outs, "kprobs_movie.plt"))) {
    outs.insert(std::pair<const std::string,bool>("kprobs.out", true));
  }

  if ((isOutput(outs, "dmt_z.plt")) && (! p->wavefunction)) {
    outs.insert(std::pair<const std::string,bool>("dmt_z.out", true));
  }

  // create output directory if it does not exist //////////////////////////////

  if (!fileExists(p->outputDir)) {
    mkdir(p->outputDir.c_str(), 0755);
  }
}

/* assigns params to the Params struct from the input file */
void assignParams(std::string inputFile, Params * p) {
  std::string line;
  std::string input_param;
  std::string param_val;
  size_t equals_pos;
  size_t space_pos;

#ifdef DEBUG_INPUT_PARSER
  std::cout << "\nParsing file: " << inputFile << "\n";
#endif

  std::ifstream bash_in;  // declare input file stream

  bash_in.open(inputFile.c_str(), std::ios::in);	// open file as input stream
  if (bash_in.good() == false) {
    fprintf(stderr, "ERROR [Inputs]: file '%s' not available for reading\n", inputFile.c_str());
    exit(-1);
  }

  // read first line of input file
  getline (bash_in,line);

  // skip non-parameter lines
  while ( line != "## START INPUT PARAMETERS ##") {
    getline (bash_in,line);
  }

  while ( line != "## END INPUT PARAMETERS ##") {
    // skip comment lines
    if ( line.substr(0,1) == "#" ) {
      getline (bash_in,line);
      continue;
    }
    // find first equals sign
    equals_pos=line.find("=");
    // find first whitespace
    space_pos=(line.find(" ") > line.find("\t") ? line.find("\t") : line.find(" "));
    // parameter name is before equals sign
    input_param = line.substr(0,int(equals_pos));
    // parameter is after equals sign, before space
    param_val = line.substr(int(equals_pos)+1,int(space_pos)-int(equals_pos));
    // extract parameters
#ifdef DEBUG
    std::cout << "Parameter: " << input_param << std::endl << "New value: " << atof(param_val.c_str()) << std::endl;
#endif
    if (input_param == "timedepH") { p->timedepH = atoi(param_val.c_str()); }
    else if (input_param == "justPlots") { p->justPlots = atoi(param_val.c_str()); }
    else if (input_param == "nproc") { p->nproc = atoi(param_val.c_str()); }
    else if (input_param == "wavefunction") { p->wavefunction = atoi(param_val.c_str()); }
    else if (input_param == "coherent") { p->coherent = atoi(param_val.c_str()); }
    else if (input_param == "analytical") { p->analytical = atoi(param_val.c_str()); }
    else if (input_param == "kinetic") { p->kinetic = atoi(param_val.c_str()); }
    else if (input_param == "kineticQD") { p->kineticQD = atoi(param_val.c_str()); }
    else if (input_param == "dynamicMu") { p->dynamicMu = atoi(param_val.c_str()); }
    else if (input_param == "dynamicMuQD") { p->dynamicMuQD = atoi(param_val.c_str()); }
    else if (input_param == "progressFile") { p->progressFile = atoi(param_val.c_str()); }
    else if (input_param == "progressStdout") { p->progressStdout = atoi(param_val.c_str()); }
    else if (input_param == "abstol") { p->abstol = atof(param_val.c_str()); }
    else if (input_param == "reltol" ) { p->reltol = atof(param_val.c_str()); }
    else if (input_param == "tout" ) { p->tout = atof(param_val.c_str()); }
    else if (input_param == "numsteps" ) { p->numsteps = atoi(param_val.c_str()); }
    else if (input_param == "numOutputSteps" ) { p->numOutputSteps = atoi(param_val.c_str()); }
    else if (input_param == "kBandEdge" ) { p->kBandEdge = atof(param_val.c_str()); }
    else if (input_param == "kBandTop" ) { p->kBandTop = atof(param_val.c_str()); }
    else if (input_param == "lBandTop" ) { p->lBandTop = atof(param_val.c_str()); }
    else if (input_param == "bulk_gap" ) { p->bulk_gap = atof(param_val.c_str()); }
    else if (input_param == "Nk" ) { p->Nk = atoi(param_val.c_str()); }
    else if (input_param == "Nk_first" ) { p->Nk_first = atoi(param_val.c_str()); }
    else if (input_param == "Nk_final" ) { p->Nk_final = atoi(param_val.c_str()); }
    else if (input_param == "Nc_first" ) { p->Nc_first = atoi(param_val.c_str()); }
    else if (input_param == "Nc_final" ) { p->Nc_final = atoi(param_val.c_str()); }
    else if (input_param == "valenceBand" ) { p->valenceBand = atof(param_val.c_str()); }
    else if (input_param == "Nl" ) { p->Nl = atoi(param_val.c_str()); }
    else if (input_param == "bulkGaussSigma" ) { p->bulkGaussSigma = atof(param_val.c_str()); }
    else if (input_param == "bulkGaussMu" ) { p->bulkGaussMu = atof(param_val.c_str()); }
    else if (input_param == "me" ) { p->me = atof(param_val.c_str()); }
    else if (input_param == "mh" ) { p->mh = atof(param_val.c_str()); }
    else if (input_param == "me_c" ) { p->me_c = atof(param_val.c_str()); }
    else if (input_param == "mh_c" ) { p->mh_c = atof(param_val.c_str()); }
    else if (input_param == "X2" ) { p->X2 = atof(param_val.c_str()); }
    else if (input_param == "temperature" ) { p->temperature = atof(param_val.c_str()); }
    else if (input_param == "EF" ) { p->EF = atof(param_val.c_str()); }
    else if (input_param == "gamma1" ) { p->gamma1 = atof(param_val.c_str()); }
    else if (input_param == "gamma1_c" ) { p->gamma1_c = atof(param_val.c_str()); }
    else if (input_param == "gamma2" ) { p->gamma2 = atof(param_val.c_str()); }
    else if (input_param == "muLK" ) { p->muLK = atof(param_val.c_str()); }
    else if (input_param == "pumpFWHM" ) { p->pumpFWHM = atof(param_val.c_str()); }
    else if (input_param == "pumpPeak" ) { p->pumpPeak = atof(param_val.c_str()); }
    else if (input_param == "pumpFreq" ) { p->pumpFreq = atof(param_val.c_str()); }
    else if (input_param == "pumpAmpl" ) { p->pumpAmpl = atof(param_val.c_str()); }
    else if (input_param == "pumpPhase" ) { p->pumpPhase = atof(param_val.c_str()); }
    else if (input_param == "CBPopFlag" ) { p->CBPopFlag = atoi(param_val.c_str()); }
    else if (input_param == "VBPopFlag" ) { p->VBPopFlag = atoi(param_val.c_str()); }
    else if (input_param == "QDPopFlag" ) { p->QDPopFlag = atoi(param_val.c_str()); }
    else if (input_param == "laser_on" ) { p->laser_on = atoi(param_val.c_str()); }
    else if (input_param == "parabolicCoupling" ) { p->parabolicCoupling = atoi(param_val.c_str()); }
    else if (input_param == "scale_bubr" ) { p->scale_bubr = atoi(param_val.c_str()); }
    else if (input_param == "scale_brqd" ) { p->scale_brqd = atoi(param_val.c_str()); }
    else if (input_param == "scale_buqd" ) { p->scale_buqd = atoi(param_val.c_str()); }
    else if (input_param == "scale_laser" ) { p->scale_laser = atoi(param_val.c_str()); }
    else if (input_param == "bridge_on" ) { p->bridge_on = atoi(param_val.c_str()); }
    else if (input_param == "random_phase" ) { p->random_phase = atoi(param_val.c_str()); }
    else if (input_param == "random_seed" ) { p->random_seed = atoi(param_val.c_str()); }
    else if (input_param == "torsion" ) { p->torsion = atoi(param_val.c_str()); }
    else if (input_param == "torsionFile" ) { p->torsionFile = param_val; }
    else if (input_param == "torsionSite" ) { p->torsionSite = atoi(param_val.c_str()); }
    else if (input_param == "torsionSin2" ) { p->torsionSin2 = atoi(param_val.c_str()); }
    else if (input_param == "torsionCos2Pulse" ) { p->torsionCos2Pulse = atoi(param_val.c_str()); }
    else if (input_param == "torsionGaussianPulse" ) { p->torsionGaussianPulse = atoi(param_val.c_str()); }
    else if (input_param == "torsionCouplingV0" ) { p->torsionCouplingV0 = atof(param_val.c_str()); }
    else if (input_param == "torsionCouplingV1" ) { p->torsionCouplingV1 = atof(param_val.c_str()); }
    else if (input_param == "torsionCouplingOmega" ) { p->torsionCouplingOmega = atof(param_val.c_str()); }
    else if (input_param == "torsionCouplingPhi" ) { p->torsionCouplingPhi = atof(param_val.c_str()); }
    else {
      std::cerr << "ERROR [" << __FUNCTION__ << "]: input parameter " << input_param << " not recognized." << std::endl;
      _exit(-1);
    }
    getline (bash_in,line);
  }

  // close input file
  bash_in.close();

#ifdef DEBUG
  std::cout << std::endl;
  std::cout << "justPlots is " << p->justPlots << std::endl;
  std::cout << "timedepH is " << p->timedepH << std::endl;
  std::cout << "nproc is " << p->nproc << std::endl;
  std::cout << "wavefunction is " << p->wavefunction << std::endl;
  std::cout << "coherent is " << p->coherent << std::endl;
  std::cout << "analytical is " << p->analytical << std::endl;
  std::cout << "kinetic is " << p->kinetic << std::endl;
  std::cout << "kineticQD is " << p->kineticQD << std::endl;
  std::cout << "dynamicMu is " << p->dynamicMu << std::endl;
  std::cout << "dynamicMuQD is " << p->dynamicMuQD << std::endl;
  std::cout << "progressFile is " << p->progressFile << std::endl;
  std::cout << "progressStdout is " << p->progressStdout << std::endl;
  std::cout << "abstol is " << p->abstol << std::endl;
  std::cout << "reltol is " << p->reltol << std::endl;
  std::cout << "tout is " << p->tout << std::endl;
  std::cout << "numsteps is " << p->numsteps << std::endl;
  std::cout << "numOutputSteps is " << p->numOutputSteps << std::endl;
  std::cout << "kBandEdge is " << p->kBandEdge << std::endl;
  std::cout << "kBandTop is " << p->kBandTop << std::endl;
  std::cout << "lBandTop is " << p->lBandTop << std::endl;
  std::cout << "bulk_gap is " << p->bulk_gap << std::endl;
  std::cout << "Nk is " << p->Nk << std::endl;
  std::cout << "valenceBand is " << p->valenceBand << std::endl;
  std::cout << "Nl is " << p->Nl << std::endl;
  std::cout << "bulkGaussSigma is " << p->bulkGaussSigma << std::endl;
  std::cout << "bulkGaussMu is " << p->bulkGaussMu << std::endl;
  std::cout << "me is " << p->me << std::endl;
  std::cout << "mh is " << p->mh << std::endl;
  std::cout << "me_c is " << p->me_c << std::endl;
  std::cout << "mh_c is " << p->mh_c << std::endl;
  std::cout << "X2 is " << p->X2 << std::endl;
  std::cout << "temperature is " << p->temperature << std::endl;
  std::cout << "EF is " << p->EF << std::endl;
  std::cout << "gamma1 is " << p->gamma1 << std::endl;
  std::cout << "gamma1_c is " << p->gamma1_c << std::endl;
  std::cout << "gamma2 is " << p->gamma2 << std::endl;
  std::cout << "muLK is " << p->muLK << std::endl;
  std::cout << "pumpFWHM is " << p->pumpFWHM << std::endl;
  std::cout << "pumpPeak is " << p->pumpPeak << std::endl;
  std::cout << "pumpFreq is " << p->pumpFreq << std::endl;
  std::cout << "pumpAmpl is " << p->pumpAmpl << std::endl;
  std::cout << "pumpPhase is " << p->pumpPhase << std::endl;
  std::cout << "CBPopFlag is " << p->CBPopFlag << std::endl;
  std::cout << "VBPopFlag is " << p->VBPopFlag << std::endl;
  std::cout << "QDPopFlag is " << p->QDPopFlag << std::endl;
  std::cout << "Nk_first is " << p->Nk_first << std::endl;
  std::cout << "Nk_final is " << p->Nk_final << std::endl;
  std::cout << "Nc_first is " << p->Nc_first << std::endl;
  std::cout << "Nc_final is " << p->Nc_final << std::endl;
  std::cout << "laser_on is " << p->laser_on << std::endl;
  std::cout << "parabolicCoupling is " << p->parabolicCoupling << std::endl;
  std::cout << "scale_bubr is " << p->scale_bubr << std::endl;
  std::cout << "scale_brqd is " << p->scale_brqd << std::endl;
  std::cout << "scale_buqd is " << p->scale_buqd << std::endl;
  std::cout << "scale_laser is " << p->scale_laser << std::endl;
  std::cout << "bridge_on is " << p->bridge_on << std::endl;
  std::cout << "random_phase is " << p->random_phase << std::endl;
  std::cout << "random_seed is " << p->random_seed << std::endl;
  std::cout << "torsion is " << p->torsion << std::endl;
  std::cout << "torsionFile is " << p->torsionFile << std::endl;
  std::cout << "torsionSite is " << p->torsionSite << std::endl;
  std::cout << "torsionSin2 is " << p->torsionSin2 << std::endl;
  std::cout << "torsionCos2Pulse is " << p->torsionCos2Pulse << std::endl;
  std::cout << "torsionGaussianPulse is " << p->torsionGaussianPulse << std::endl;
  std::cout << "torsionCouplingV0 is " << p->torsionCouplingV0 << std::endl;
  std::cout << "torsionCouplingV1 is " << p->torsionCouplingV1 << std::endl;
  std::cout << "torsionCouplingOmega is " << p->torsionCouplingOmega << std::endl;
  std::cout << "torsionCouplingPhi is " << p->torsionCouplingPhi << std::endl;
#endif

  // Error checking
  if ((p->VBPopFlag && p->CBPopFlag) || (p->VBPopFlag && p->QDPopFlag) || (p->CBPopFlag && p->QDPopFlag)) {
    std::cerr << "\nWARNING: population starting in multiple locations.\n";
  }

  if (p->CBPopFlag == POP_CONSTANT) {
    if (p->Nk_final > p->Nk || p->Nk_final < 1) {
      fprintf(stderr, "ERROR [Inputs]: Nk_final greater than Nk or less than 1.\n");
      exit(-1);
    }
    if (p->Nk_final < p->Nk_first) {
      fprintf(stderr, "ERROR [Inputs]: Nk_final is less than Nk_first.\n");
      exit(-1);
    }
  }

  if (p->QDPopFlag == POP_CONSTANT) {
    if (p->Nc_final > p->Nc || p->Nc_final < 1) {
      fprintf(stderr, "ERROR [Inputs]: Nc_final greater than Nc or less than 1.\n");
      exit(-1);
    }
    if (p->Nc_final < p->Nc_first) {
      fprintf(stderr, "ERROR [Inputs]: Nc_final is less than Nc_first.\n");
      exit(-1);
    }
  }

  if (p->Nl < 0) {
    fprintf(stderr, "ERROR [Inputs]: Nl less than 0.\n");
    exit(-1);
  }

  if (p->random_seed < -1) {
    std::cerr << "\nERROR: random_phase must be -1 or greater.\n";
    exit(-1);
  }

  // Decide which output files to make /////////////////////////////////////////

#ifdef DEBUG
  std::cout << "Assigning outputs as specified in " << p->inputFile << "\n";
#endif
  assignOutputs(p->inputFile.c_str(), p->outs, p);

  // error checking on various parameters //////////////////////////////////////

  // check torsion parameters, set up torsion spline ///////////////////////////

  if (p->torsion) {
#ifdef DEBUG
    std::cout << "Torsion is on." << std::endl;
#endif

    // error checking
    if (p->torsionSite < 0) {
      std::cerr << "ERROR: torsion site is less than zero." << std::endl;
      exit(-1);
    }

    if (!p->torsionSin2 && !p->torsionCos2Pulse && !p->torsionGaussianPulse) {
      if (!fileExists(p->torsionFile)) {
      std::cerr << "ERROR: torsion file " << p->torsionFile << " does not exist." << std::endl;
      exit(-1);
      }

      // create spline
      p->torsionV.readFile(p->torsionFile.c_str());
      if (p->torsionV.getFirstX() != 0.0) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should start at 0.0." << std::endl;
        exit(-1);
      }

      if (p->torsionV.getLastX() < p->tout) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should be >= tout." << std::endl;
        exit(-1);
      }
    }

    if ((p->torsionSin2 + p->torsionCos2Pulse + p->torsionGaussianPulse) > 1) {
      std::cerr << "ERROR: more than one torsional coupling turned on." << std::endl;
      exit(-1);
    }
  }

  return;
}

// checks if a file exists (can be opened)
bool fileExists(std::string fileName) {
  std::ifstream ifile(fileName.c_str());
  return ifile.good();
}
