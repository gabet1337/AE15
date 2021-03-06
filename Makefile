CXX=g++
CXXFLAGS=-Wall -O2 -std=c++11
.PHONY: pred
pred:	clean
	$(CXX) $(CXXFLAGS) pred.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o pred 

hw_counters: clean
	$(CXX) hw_counters.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o hw

small:	clean
	$(CXX) -O3 small_test.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o small

matrix:	clean
	$(CXX) -03 -fopenmp -Wall -std=c++11 matrix.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o matrix

count:	clean
	$(CXX) $(CXXFLAGS) counting_sort.cpp -o count

radix:	clean
	$(CXX) -O3 -fopenmp -Wall -std=c++11 radix_sort.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o radix

clean:
	rm -f count
	rm -f pred
	rm -f hw
	rm -f small
	rm -f matrix
	rm -f radix
