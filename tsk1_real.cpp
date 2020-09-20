/**
 * A main module for the task1.
 */
 #include "omp.h"
 #include "tsk1_utils.h"

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
