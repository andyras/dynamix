### This line should be just about all you have to change ######################
# name of main executable
BIN = dynamix
################################################################################

#DEBUG_MAKE = "foo"
# optional #####################################################################
# USE_MKL = yes
# USE_BOOST = yes
# DEBUG_MAKE = yes
################################################################################

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
LIBDIR = ./lib
BINDIR = ./bin

# compiler and archiver
CPP = g++
AR = ar

CPPFLAGS = -O3 -Wall -std=c++11 -fopenmp
LDFLAGS = -lsundials_cvode -lsundials_nvecserial
INCLUDES = -I$(INCDIR)

ifdef USE_BOOST
  CPPFLAGS += -D__BOOST_SERIALIZE__
endif

ifdef USE_MKL
  CPPFLAGS += -D__USE_MKL__
  LDFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
endif

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
  LDFLAGS += -pg -g -debug
endif

# basic compiler-dependent flags
ifeq ($(CPP),icpc)
  CPPFLAGS += -fast -xHOST -no-prec-div -mkl -no-multibyte-chars
  LDFLAGS += -mkl
else # g++*
  LDFLAGS += -lblas -lgfortran
  ifeq ($(shell hostname),tim.selfip.org)
    LDFLAGS += -lpthread -lm
  endif
endif

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CPP) -o $@ $^ $(LDFLAGS) -fopenmp

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

# add all obj files except main to library
library: $(OBJECTS) $(LIBDIR)
	find $(OBJDIR) -name *.o ! -name main.o | xargs $(AR) vrs $(LIBDIR)/lib$(BIN).a

.PHONY: clean install uninstall $(BINDIR) $(OBJDIR) $(LIBDIR	)

clean:
	rm -f $(BIN)
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

install: $(BIN) | $(BINDIR)
	cp dynamix $(BINDIR)

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(LIBDIR):
	mkdir -p $(LIBDIR)
