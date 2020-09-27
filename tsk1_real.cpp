/**
 * A main module for the task1.
 */
#include <cstring>
#include "omp.h"
#include "tsk1_utils.h"
 
/**
 * A wrapper over a graph generation
 * Warning: 
 * 	   Delete 
 */
void generateWrapper( size_t row_len, size_t column_len, size_t undivided, size_t divided,
	int* N, int** IA, int** JA){
	MatrixParameters matrix_param( row_len, column_len, undivided, divided);
	NetGraph graph( &matrix_param);
	graph.generate( &matrix_param);
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
    NetGraph graph( &matrix_param);
    double start = omp_get_wtime();
    graph.generate( &matrix_param);
    double end = omp_get_wtime();
    std::cout << "Time: " << end - start << std::endl;
    if( program_env.isDebugPrint() ){
        graph.printGraph();
    }
    return 0;
 }
