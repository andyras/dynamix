#include <iostream>
#include <map>
#include <string>

#include "libdynamix_input_parser.h"

int main () {
 std::map<std::string, bool> outputs;

 std::cout << "\n";
 std::cout << "testing getting outputs to create from file ins/parameters.in\n";

 assignOutputs("ins/outputs.in", outputs);

 return 0;
}
