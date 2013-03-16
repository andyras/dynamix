
#ifndef __LIBDYNAMIX_INPUT_PARSER_H__
#define __LIBDYNAMIX_INPUT_PARSER_H__

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

#endif
