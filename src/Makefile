# original line
# #CPPFLAGS = -O3 -Wall -I/home/andyras/bin/include -L/home/andyras/bin/lib \
	 -lsundials_cvode -lsundials_nvecserial -xHOST -ipo -no-prec-div \
	 -mkl

CPP = icpc
CPPFLAGS = -O3 -Wall -xHOST -ipo -no-prec-div
LDFLAGS = -lsundials_cvode -lsundials_nvecserial -mkl

#CPP = g++
#CPPFLAGS = -O3 -Wall -lsundials_cvode -lsundials_nvecserial -llapack

dynamix: dynamix.o libdynamix_outputs.o
	$(CPP) -o dynamix dynamix.o libdynamix_outputs.o $(LDFLAGS)

dynamix.o: dynamix.cpp libdynamix_outputs.cpp libdynamix_outputs.h
	$(CPP) $(CPPFLAGS) -c dynamix.cpp

libdynamix_outputs.o: libdynamix_outputs.cpp libdynamix_outputs.h
	$(CPP) $(CPPFLAGS) -c libdynamix_outputs.cpp

clean:
	rm -f *.o

cleanall:
	make clean
	rm -f dynamix

build:
	make cleanall
	make
