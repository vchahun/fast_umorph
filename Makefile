segment: segment.cc vocabulary.h corpus.h prob.h trie.h banana.h pss_model.h
	clang++ -std=c++11 -O3 -I/usr/include/x86_64-linux-gnu/c++/4.7/ $< -o $@ -lfst -ldl

prefsuf: prefsuf.cc prob.h vocabulary.h corpus.h
	clang++ -std=c++11 -O3 -I/usr/include/x86_64-linux-gnu/c++/4.7/ $< -o $@
