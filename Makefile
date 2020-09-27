.PHONY: clean, tsk1
tsk1: tsk1_graph_prepare.cpp tsk1_graph_prepare.h tsk1_real.cpp tsk1_real.h \
    tsk1_utils.h
	g++ -O3 -o tsk1 -fopenmp --std=c++11 tsk1_graph_prepare.cpp tsk1_real.cpp
clean: 
	rm tsk1