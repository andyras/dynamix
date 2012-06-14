CPP = g++
CFLAGS = -O2 -Wall -I/home/andyras/bin/include -L/home/andyras/bin/lib \
	 -lm -lsundials_cvode -lsundials_nvecserial

# Test if new g++ is installed as a different version
#EXISTS := $(shell which g++-4.7.0)
#ifneq ($(EXISTS),)
 #CPP = g++-4.7.0
#endif
#EXISTS := $(shell which g++-4.7.1)
#ifneq ($(EXISTS),)
 #CPP = g++-4.7.1
#endif

dynamix: dynamix.cpp
	$(CPP) dynamix.cpp -o dynamix $(CFLAGS)

.PHONY : clean
clean:
	rm -f dynamix
