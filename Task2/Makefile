.PHONY: clean, tsk2
CFLAGS=-O3 -lmsmpi --std=c++11
LLIB = 
ifeq ($(OS),Windows_NT)
LLIB += -lpsapi
TESTS_DIR := tests\\
else
TESTS_DIR := tests\/
endif
tsk2:
	g++ $(CFLAGS) -o tsk2 tsk2_graph_prepare.cpp tsk2_vector.cpp\
    $(TESTS_DIR)test_Vector.cpp tsk2_solver.cpp tsk2_real.cpp $(LLIB)
tsk2_Measure:
	g++ $(CFLAGS) -DMEASURE_GENERATE -DMEASURE_FILL -DMEASURE_SOLVER -o tsk2_msr\
    tsk2_graph_prepare.cpp tsk2_vector.cpp $(TESTS_DIR)test_Vector.cpp tsk2_solver.cpp tsk2_real.cpp $(LLIB)
tsk2_Measure_Solver:
	g++ $(CFLAGS) -DMEASURE_VECTOR_OPS -DMEASURE_SOLVER -o tsk2_msr_slv\
    tsk2_graph_prepare.cpp tsk2_vector.cpp $(TESTS_DIR)test_Vector.cpp tsk2_solver.cpp tsk2_real.cpp $(LLIB)
clean: 
	rm tsk2
