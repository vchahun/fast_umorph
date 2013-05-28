segment: segment.cc vocabulary.h corpus.h prob.h trie.h banana.h pss_model.h thread_pool.h
	g++-4.7 -std=c++11 -O3 -pthread $< -o $@ -lfst -ldl

prefsuf: prefsuf.cc prob.h vocabulary.h corpus.h
	g++-4.7 -std=c++11 -O3 $< -o $@
