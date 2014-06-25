# name of main executable
### This line should be just about all you have to change ###
BIN = dynamix
#############################################################

#DEBUG_MAKE = "foo"

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# compiler
CPP = icpc

CPPFLAGS = -O3 -Wall -std=c++11 -fopenmp
LDFLAGS = -lsundials_cvode -lsundials_nvecserial
INCLUDES = -I$(INCDIR)

# basic compiler-dependent flags
ifeq ($(CPP),icpc)
  CPPFLAGS += -fast -xHOST -no-prec-div -mkl -no-multibyte-chars
  LDFLAGS += -mkl
else # g++*
  LDFLAGS += -lgomp
  LDFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
  ifeq ($(shell hostname),tim.selfip.org)
    #INCLUDES += -L/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64/
    LDFLAGS += -lpthread -lm
  endif
endif

# includes for SUNDIALS/CVODE, useful on Arch
#INCLUDES += -I/usr/include/sundials -I/usr/include/cvode

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
  LDFLAGS += -pg -g -debug
endif


#INCLUDES += -I/usr/local/Cellar/gcc48/4.8.3/lib/gcc/x86_64-apple-darwin13.2.0/4.8.3/include/
#INCLUDES += -L/opt/intel/composer_xe_2013_sp1.1.103/mkl/lib/
#INCLUDES += -L/usr/local/Cellar/gcc48/4.8.3/lib/gcc/x86_64-apple-darwin13.2.0/4.8.3/
#LDFLAGS = -lsundials_cvode -lsundials_nvecserial -lgomp -lfftw3

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CPP) -o $@ $^ $(LDFLAGS)

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# make object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean build install uninstall $(BINDIR)

clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)

build:
	make clean
	make

install: $(BIN) | $(BINDIR)
	cp dynamix $(BINDIR)

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(BINDIR):
	mkdir -p ../bin
