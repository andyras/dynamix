CPP = g++
CFLAGS = -O2 -Wall -I/home/andyras/bin/include -L/home/andyras/bin/lib \
	 -lm -lsundials_cvode -lsundials_nvecserial

dynamix: dynamix.cpp
	$(CPP) dynamix.cpp -o dynamix $(CFLAGS)

.PHONY : clean
clean:
	rm -f dynamix
