dynamix: dynamix.cpp
	g++ -Wall dynamix.cpp -O3 -o dynamix -I/home/andyras/bin/include -L/home/andyras/bin/lib -lm -lsundials_cvode -lsundials_nvecserial

clean:
	rm -f dynamix
