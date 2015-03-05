CXX=g++
CXXFLAGS=-O3 -Wall
.PHONY: pred
pred:
	rm -f pred
	$(CXX) $(CXXFLAGS) pred.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o pred 

hw_counters:
	rm -f hw
	$(CXX) hw_counters.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o hw

small:
	rm -f small
	$(CXX) -O3 small_test.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o small

matrix:
	rm -f matrix
	$(CXX) -03 matrix.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -i matrix