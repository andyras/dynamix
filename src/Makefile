CPP = icpc
CFLAGS = -O3 -Wall -I/home/andyras/bin/include -L/home/andyras/bin/lib \
	 -lm -lsundials_cvode -lsundials_nvecserial -xHOST -ipo -no-prec-div \
	 -mkl

dynamix: dynamix.cpp
	$(CPP) dynamix.cpp -o dynamix $(CFLAGS)

.PHONY : clean
clean:
	rm -f dynamix
