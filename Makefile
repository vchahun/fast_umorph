prefsuf: prefsuf.cc prob.h vocabulary.h
	clang++ -std=c++11 -O3 -I/usr/include/x86_64-linux-gnu/c++/4.7/ $< -o $@
