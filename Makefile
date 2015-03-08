CXX=g++
CXXFLAGS=-Wall -O2
.PHONY: pred
pred:
	rm -f pred
	$(CXX) $(CXXFLAGS) pred.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o pred 

hw_counters:
	rm -f hw
	$(CXX) hw_counters.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o hw

matrix: clean
	rm -f matrix
	$(CXX) $(CXXFLAGS) matrix.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o matrix

small:
	rm -f small
	$(CXX) -O3 small_test.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o small

clean:
	rm -f matrix
	rm -f small
	rm -f pred
	rm -f hw
	$(CXX) -03 matrix.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -i matrix
