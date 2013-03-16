#include <fstream>
#include <iostream>
#include <string>
#include <map>

#include "libdynamix_input_parser.h"

//#define DEBUG_INPUT_PARSER

// Methods for the outputFile class

// Construct an outputFile object
outputFile::outputFile(char * fileName) {
 name = fileName;
}

// Get the name of an output file.
char * outputFile::getName() {
 return name;
}

// Toggle the creation of an output file
void outputFile::create() {
 createMe = true;
}

//void parseInput(const char * inputFile)

void assignOutputs(const char * inputFile, std::map<std::string, bool> outputs) {
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

 // read inputs until [[End]] header or EOF
 while (getline(input, line)) {
  if (line != "[[End]]") {
   // skip comments
   if (line.substr(0,1) == "#") {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "Skipping comment line: " << line << "\n";
#endif
    getline(input, line);
    continue;
   }
#ifdef DEBUG_INPUT_PARSER
   else {
    std::cout << "This line is not a comment.\n";
    std::cout << "Output toggled: " << line << "\n";
   }
#endif
   // turn on outputs which are in this section of the input
   outputs[line] = true;
  }
 }
}
