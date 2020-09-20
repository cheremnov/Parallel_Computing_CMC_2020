/**
 * Different utilities for the parsing
 */
#include <cstring>
#include <fstream>
#include "tsk1_real.h"
#include "tsk1_graph_prepare.h"
enum { 
    // A file argument number
    FILE_ARG_NUM = 1 
};
/**
 * Print a program help
 */
void printHelp(){
    std::cout << "tsk1 FILE_NAME [-d]" << std::endl;
    std::cout << "File must be put at the same directory" << std::endl;
    std::cout << "-d enables a debug print" << std::endl;
}
/**
 * Read the parameters from the file
 */
int readMatrixParametersFile( const char* file_name, 
                              MatrixParameters *matrix_params_p){
    std::ifstream matrix_param_file( file_name);
    // Was the file opened
    if( !matrix_param_file.is_open() ){
        std::cout << "Can't open the file" << std::endl;
        return -1;
    }
    size_t row_len, column_len, divided, not_divided;
    // Was the input correct
    if( !( matrix_param_file >> row_len >> column_len >> not_divided >> divided) ){
        std::cout << "Can't parse the input" << std::endl;
        return -1;
    }
    size_t max_matrix_dimension = MatrixParameters::MatrixConstraints_t::MAX_MATRIX_DIMENSION;
    size_t max_cells_cnt = MatrixParameters::MatrixConstraints_t::MAX_CELLS;
    if( row_len > max_matrix_dimension || column_len > max_matrix_dimension ){
        std::cout << "The matrix parameters exceeds limitations" << std::endl;
        return -1;
    }
    if( divided > max_cells_cnt || not_divided > max_cells_cnt ){
        std::cout << "A number of cells exceeds boundaries" << std::endl;
        return -1;
    }
    matrix_params_p->setParameters( row_len, column_len, not_divided, divided);
    return 0;
}
/**
 * Parse the cmd arguments
 * Results:
 *     -1, if failure. 0 otherwise
 */
int parseCMDArguments( int argc, char** argv, 
                       MatrixParameters* matrix_params_p, // Matrix parameters
                       ProgramEnv* program_env_p){        // Environment state
    if( argc <= FILE_ARG_NUM ){
        std::cout << "No filename specified" << std::endl;
        return -1;
    }
    int file_read = readMatrixParametersFile( argv[FILE_ARG_NUM], matrix_params_p);
    if( file_read == -1 ){
        std::cout << "Can't parse command-line arguments" << std::endl;
        return -1;
    }
    // Find the other options
    for( int arg_idx = 1; arg_idx < argc; ++arg_idx ){
        if( !strcmp( "-d", argv[arg_idx]) ){
            program_env_p->setDebugPrint( true);
        }
    }
    return 0;
}
