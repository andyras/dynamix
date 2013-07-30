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
std::map<const std::string, bool> initOutputMap();

const std::string makeConstString(std::string inputString);

// turns on output files in 'outputs' map
void assignOutputs(const char * inputFile, std::map<const std::string, bool> &outputs);

bool fileExists(std::string fileName);

#endif
