#include <iostream>
#include <string>
#include <stdio.h>

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
