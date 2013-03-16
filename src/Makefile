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

libdynamix_input_parser.o: libdynamix_input_parser.cpp libdynamix_input_parser.h
	$(CPP) $(CPPFLAGS) -c libdynamix_input_parser.cpp

clean:
	rm -f *.o
	rm -f dynamix
	rm -f *.out
	rm -f tests

build:
	make clean
	make

tests: dynamix libdynamix_outputs_test.o
	$(CPP) -o tests libdynamix_outputs_test.o libdynamix_outputs.o $(LDFLAGS)

libdynamix_outputs_test.o: libdynamix_outputs_test.cpp libdynamix_outputs.cpp libdynamix_outputs.h
	$(CPP) $(CPPFLAGS) -c libdynamix_outputs_test.cpp
