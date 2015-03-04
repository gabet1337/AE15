.PHONY: pred
pred:
	rm -f pred
	g++ -O3 pred.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o pred 

hw_counters:
	rm -f hw
	g++ hw_counters.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o hw
