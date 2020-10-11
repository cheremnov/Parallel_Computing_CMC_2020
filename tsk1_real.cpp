/**
 * A main module for the task1.
 */
#include <cstring>
#include <stdint.h>
#include "omp.h"
#include "tsk1_utils.h"
#include "tsk1_vector.h"
#include "tsk1_solver.h"
#include "tests/test_Vector.h"
 
const double CONVERGENCE_EPS = 0.00001;
/**
 * A wrapper over a graph generation
 * Warning: 
 * 	   Manually delete the arrays afterwards
 */
void generateWrapper( size_t row_len, size_t column_len, size_t undivided, size_t divided,
	int* N, int** IA, int** JA){
	MatrixParameters matrix_param( row_len, column_len, undivided, divided);
	NetGraph graph( &matrix_param);
	graph.generate( &matrix_param, 1);
	*N = graph.getNodesCount();
	// Allocate memory
	size_t edges_max = *N * NETGRAPH_MAX_EDGES_NODE;
	*IA = new int[edges_max];
	*JA = new int[edges_max];
	memcpy( *IA, graph.getIA(), edges_max);
	memcpy( *JA, graph.getJA(), edges_max);
}
int main( int argc, char **argv){
    if( argc == 1 ){
        printHelp();
        return 0;
    } 
    MatrixParameters matrix_param;
    ProgramEnv program_env;
    int parse_env = parseCMDArguments( argc, argv, &matrix_param, &program_env);
    if( parse_env == -1 ){
        return -1;
    }
    omp_set_num_threads( program_env.getThreadsNum());
    // Run the tests
    launchTests();
    double start = omp_get_wtime();
    // Measure the phases time if the parameter is set
#ifdef MEASURE_GENERATE
    #ifdef MEASURE_MEMORY
    uint64_t generate_before_mem = getMemoryUsage();
    #endif
    double generate_start = omp_get_wtime();
#endif
    NetGraph graph( &matrix_param);
    graph.generate( &matrix_param, program_env.getThreadsNum());
#ifdef MEASURE_GENERATE
    #ifdef MEASURE_MEMORY
    uint64_t generate_after_mem = getMemoryUsage();
    std::cout << "Generation memory usage: " << generate_after_mem -
    generate_before_mem << std::endl;
    #endif
    double generate_end = omp_get_wtime();
    std::cout << "Generation time: " << generate_end - generate_start << std::endl;
#endif
#ifdef MEASURE_FILL
    #ifdef MEASURE_MEMORY
    uint64_t fill_before_mem = getMemoryUsage();
    #endif
    double fill_start = omp_get_wtime();
#endif
    graph.fillMatrix( program_env.getThreadsNum());
    MathVector b_vec( graph.getNodesCount());
    b_vec.fillVector( program_env.getThreadsNum());
#ifdef MEASURE_FILL
    #ifdef MEASURE_MEMORY
    uint64_t fill_after_mem = getMemoryUsage();
    std::cout << "Fill memory usage: " << fill_after_mem -
    fill_before_mem << std::endl;
    #endif
    double fill_end = omp_get_wtime();
    std::cout << "Fill time: " << fill_end - fill_start << std::endl;
#endif
#ifdef MEASURE_SOLVER
    #ifdef MEASURE_MEMORY
    uint64_t solver_before_mem = getMemoryUsage();
    #endif
    double solver_start = omp_get_wtime();
#endif
    solverCG( graph, b_vec, program_env.isDebugPrint(), CONVERGENCE_EPS);
#ifdef MEASURE_SOLVER
    #ifdef MEASURE_MEMORY
    uint64_t solver_after_mem = getMemoryUsage();
    std::cout << "Solver memory usage: " << solver_after_mem -
    solver_before_mem << std::endl;
    #endif
    double solver_end = omp_get_wtime();
    std::cout << "Solver time: " << solver_end - solver_start << std::endl;
#endif
    double end = omp_get_wtime();
#ifdef MEASURE_MEMORY
    std::cout << "Memory usage: " << getMemoryUsage() << std::endl;
#endif
    std::cout << "Time: " << end - start << std::endl;
    if( program_env.isDebugPrint() ){
        graph.printGraph();
        b_vec.printVector();
    }
    return 0;
 }
