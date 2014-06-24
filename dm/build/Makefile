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

# flags for compiler
# optimized CPPFLAGS
# Intel
#CPPFLAGS = -O0 -Wall --std=c++11 -openmp
CPPFLAGS = -O3 -fast -xHOST -no-prec-div -Wall -std=c++11 -fopenmp -mkl -no-multibyte-chars
CPPFLAGS = -O3 -Wall -std=c++11 -fopenmp

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
endif

# non-optimized CPPFLAGS
#CPPFLAGS = -Wall -g -O0 -openmp --std=c++11

# flags for linker
INCLUDES = -I$(INCDIR) -I$(MKLROOT)/include/fftw
# includes for CVODE
INCLUDES += -I$(HOME)/bin/include -I$(HOME)/bin/include/cvode -L$(HOME)/bin/lib -I/usr/include/sundials -I/usr/include/cvode -L/usr/lib
INCLUDES += -I/usr/local/Cellar/gcc48/4.8.3/lib/gcc/x86_64-apple-darwin13.2.0/4.8.3/include/
INCLUDES += -L/opt/intel/composer_xe_2013_sp1.1.103/mkl/lib/
INCLUDES += -L/usr/local/Cellar/gcc48/4.8.3/lib/gcc/x86_64-apple-darwin13.2.0/4.8.3/
#LDFLAGS = -lsundials_cvode -lsundials_nvecserial -lgomp -lfftw3
# Intel flags for linker
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -lmkl_intel_ilp64 -lgomp
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl
ifdef DEBUG_MAKE
  LDFLAGS += -pg -g -debug
endif

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
