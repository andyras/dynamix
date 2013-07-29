#ifndef __LIBDYNAMIX_INPUT_PARSER_H__
#define __LIBDYNAMIX_INPUT_PARSER_H__

#include <fstream>
#include <iostream>
#include <string>
#include <map>

class outputFile {
 private:
  // variables
  bool createMe;
  char * name;
  //constructor
  outputFile() {};

 public:
  // public constructor
  outputFile(char * fileName);
  // methods
  void create();
  char * getName();
};

// initiator for a string->bool map
std::map<std::string, bool> initOutputMap();

// turns on output files in 'outputs' map
void assignOutputs(const char * inputFile, std::map<std::string, bool> &outputs);

bool fileExists(std::string fileName);

#endif
