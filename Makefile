.PHONY: clean, tsk1
CFLAGS=-O3 -fopenmp --std=c++11
LLIB = 
ifeq ($(OS),Windows_NT)
LLIB += -lpsapi
TESTS_DIR := tests\
else
TESTS_DIR := tests/
endif
tsk1: tsk1_graph_prepare.cpp tsk1_real.cpp tsk1_vector.cpp\
$(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp
	g++ $(CFLAGS) -o tsk1 tsk1_graph_prepare.cpp tsk1_vector.cpp\
    $(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp tsk1_real.cpp $(LLIB)
tsk1_Measure: tsk1_graph_prepare.cpp tsk1_real.cpp tsk1_vector.cpp\
$(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp
	g++ $(CFLAGS) -DMEASURE_GENERATE -DMEASURE_FILL -DMEASURE_SOLVER -o tsk1_msr\
    tsk1_graph_prepare.cpp tsk1_vector.cpp $(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp tsk1_real.cpp $(LLIB)
tsk1_Measure_Solver: tsk1_graph_prepare.cpp tsk1_real.cpp tsk1_vector.cpp\
$(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp
	g++ $(CFLAGS) -DMEASURE_VECTOR_OPS -DMEASURE_SOLVER -o tsk1_msr_slv\
    tsk1_graph_prepare.cpp tsk1_vector.cpp $(TESTS_DIR)test_Vector.cpp tsk1_solver.cpp tsk1_real.cpp $(LLIB)
clean: 
	rm tsk1
