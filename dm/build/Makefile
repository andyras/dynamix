# name of main executable
### This line should be just about all you have to change ###
BIN = dynamix
#############################################################

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# flags for compiler and linker
CXX = icpc
CPP = $(CXX)
# optimized CPPFLAGS
CPPFLAGS = -O3 -Wall -xHOST -ipo -no-prec-div
# non-optimized CPPFLAGS
CPPFLAGS = -Wall -g -O0 -openmp
INCLUDES = -I./$(INCDIR) -I$(MKLROOT)/include/fftw
# addition for radon
INCLUDES += -I/home/andyras/bin/include -I/home/andyras/bin/include/cvode -L/home/andyras/bin/lib
# addition for quest
INCLUDES += -I/home/ara132/bin/include -I/home/ara132/bin/include/cvode -L/home/ara132/bin/lib
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CPP) -o $@ $^ $(LDFLAGS) $(INCLUDES)

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# make object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean cleantest build install uninstall $(BINDIR)
clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)

cleantest:
	rm -f *.out
	rm -f *.plt
	rm -f *test

build:
	make clean
	make

install: $(BIN) | $(BINDIR)
	cp dynamix ../bin

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(BINDIR):
	mkdir -p ../bin

# TODO update testing section of this file
tests: dynamix libdynamix_outputs_test.o libdynamix_input_parser_test.o
	$(CPP) -o libdynamix_outputs_test \
	 libdynamix_outputs_test.o libdynamix_outputs.o $(LDFLAGS)
	$(CPP) -o libdynamix_input_parser_test \
	 libdynamix_input_parser_test.o libdynamix_input_parser.o $(LDFLAGS)

libdynamix_outputs_test.o: libdynamix_outputs_test.cpp \
        libdynamix_outputs.cpp libdynamix_outputs.h
	$(CPP) $(CPPFLAGS) -c libdynamix_outputs_test.cpp

libdynamix_input_parser_test.o: libdynamix_input_parser_test.cpp \
        libdynamix_input_parser.cpp libdynamix_input_parser.h
	$(CPP) $(CPPFLAGS) -c libdynamix_input_parser_test.cpp
