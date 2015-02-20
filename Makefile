.PHONY: pred
pred:
	rm -f pred
	g++ pred.cpp /usr/lib/x86_64-linux-gnu/libpapi.so -o pred
