/**
 * A main module for the task1.
 */
#include <cstring>
#include <stdint.h>
#include <chrono>
#include <mpi.h>
#include "tsk2_utils.h"
#include "tsk2_vector.h"
#include "tsk2_solver.h"
#include "tests/test_Vector.h"
 

typedef std::chrono::milliseconds ms;
const double CONVERGENCE_EPS = 0.00001;
/**
 * Divide a matrix onto the blocks
 * Get the parameters of the current block,
 * based on the processor number
 */
void divideMatrix( MatrixParameters* params_p, ProgramEnv* env_p){
    size_t row_len = params_p->getRowLen();
    size_t column_len = params_p->getColumnLen();
    // How many rows and columns are there in block representation of the matrix
    size_t block_rows = params_p->getBlockRows();
    size_t block_columns = params_p->getBlockColumns();
    // A current processor rank
    size_t processor_rank = env_p->getProcessRank();
    /**
     * A block numeration example:
     * [0 1 2 3]
     * [4 5 6 7]
     * [8 9 10 11]
     */
    // Current block indexes
    size_t block_row_idx = processor_rank / block_columns;
    size_t block_column_idx = processor_rank % block_columns;
    /**
     * Divide the net onto blocks,
     * so that they has almost the same number of rows and columns
     * If the rows or columns can't be splitted evenly, 
     * then the modula is distributed to the first blocks
     * The blocks have rows [row_start_idx; row_end_idx)
     * and columns [column_start_idx; column_end_idx)
     */
    // How many rows and columns are in the block
    size_t rows_in_block = (row_len + 1) / block_rows;
    size_t columns_in_block = (column_len + 1) / block_columns;
    size_t row_start_idx = block_row_idx * rows_in_block 
        + std::min( block_row_idx, (row_len + 1) % block_rows);
    size_t row_end_idx = (block_row_idx + 1) * rows_in_block
        + std::min( block_row_idx + 1, (row_len + 1) % block_rows);
    size_t column_start_idx = block_column_idx * columns_in_block 
        + std::min( block_column_idx, (column_len + 1) % block_columns);
    size_t column_end_idx = (block_column_idx + 1) * columns_in_block
        + std::min( block_column_idx + 1, (column_len + 1) % block_columns);
    env_p->setBlockParams( row_start_idx, row_end_idx,
        column_start_idx, column_end_idx);
    std::map<int,int>& local_to_global = env_p->getLocalToGlobal();
    std::map<int,int>& global_to_local = env_p->getGlobalToLocal();
    std::vector<int>& parts = env_p->getParts(); 
    // Calculate local nodes
    size_t local_nodes_num = 0;
    for( size_t row_idx = row_start_idx; row_idx < row_end_idx; ++row_idx ){
        for( size_t column_idx = column_start_idx; column_idx < column_end_idx;
            ++column_idx){
            global_to_local[row_idx * (column_len + 1) + column_idx] 
                = local_nodes_num;
            local_to_global[local_nodes_num] 
                = row_idx * (column_len + 1) + column_idx;
            // The local vertexes belong to the current processor
            parts.push_back( processor_rank); 
            ++local_nodes_num;
        }
    }
    // Create halo
    size_t halo_nodes_num = 0;
    // The neighbor block on the previous row
    if( block_row_idx > 0 ){
        size_t row_idx = row_start_idx - 1;
        for( size_t column_idx = column_start_idx; column_idx < column_end_idx;
            ++column_idx){
            global_to_local[row_idx * (column_len + 1) + column_idx] 
                = local_nodes_num + halo_nodes_num;
            local_to_global[local_nodes_num + halo_nodes_num] 
                = row_idx * (column_len + 1) + column_idx;
            // Get the processor on the previous row
            parts.push_back( processor_rank - block_columns);
            ++halo_nodes_num;
        }
    }
    // The neighbor block on the next row
    if( block_row_idx < block_rows - 1 ){
        size_t row_idx = row_end_idx;
        for( size_t column_idx = column_start_idx; column_idx < column_end_idx;
            ++column_idx){
            global_to_local[row_idx * (column_len + 1) + column_idx] 
                = local_nodes_num + halo_nodes_num;
            local_to_global[local_nodes_num + halo_nodes_num] 
                = row_idx * (column_len + 1) + column_idx;
            parts.push_back( processor_rank + block_columns);
            ++halo_nodes_num;
        }
    }
    // The neighbor block on the previous column
    if( block_column_idx > 0 ){
        size_t column_idx = column_start_idx - 1;
        for( size_t row_idx = row_start_idx; row_idx < row_end_idx;
            ++row_idx){
            global_to_local[row_idx * (column_len + 1) + column_idx] 
                = local_nodes_num + halo_nodes_num;
            local_to_global[local_nodes_num + halo_nodes_num] 
                = row_idx * (column_len + 1) + column_idx;
            parts.push_back( processor_rank - 1);
            ++halo_nodes_num;
        }
    }
    // The neighbor block on the next column
    if( block_column_idx < block_columns - 1 ){
        size_t column_idx = column_end_idx;
        for( size_t row_idx = row_start_idx; row_idx < row_end_idx; ++row_idx )
        {
            global_to_local[row_idx * (column_len + 1) + column_idx] 
                = local_nodes_num + halo_nodes_num;
            local_to_global[local_nodes_num + halo_nodes_num] 
                = row_idx * (column_len + 1) + column_idx;
            parts.push_back( processor_rank + 1);
            ++halo_nodes_num;
        }
    }
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
    MPI_Init( nullptr, nullptr);
    int process_rank, process_num;
    MPI_Comm_rank( MPI_COMM_WORLD, &process_rank);
    MPI_Comm_size( MPI_COMM_WORLD, &process_num);
    program_env.setProcessNum( process_num);
    program_env.setProcessRank( process_rank);
    divideMatrix( &matrix_param, &program_env);
    // Run the tests
    launchTests( &program_env);
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    // Measure the phases time if the parameter is set
#ifdef MEASURE_GENERATE
    #ifdef MEASURE_MEMORY
    uint64_t generate_before_mem = getMemoryUsage();
    #endif
    std::chrono::high_resolution_clock::time_point generate_start =
    std::chrono::high_resolution_clock::now();
#endif
    NetGraph graph( &matrix_param, &program_env);
    graph.generate( &matrix_param);
#ifdef MEASURE_GENERATE
    #ifdef MEASURE_MEMORY
    uint64_t generate_after_mem = getMemoryUsage();
    std::cout << "Generation memory usage: " << generate_after_mem -
    generate_before_mem << std::endl;
    #endif
    std::chrono::high_resolution_clock::time_point generate_end =
    std::chrono::high_resolution_clock::now();
    if( program_env.getProcessRank() == 0 ){
        std::cout << "Generation time: " <<
        (std::chrono::duration_cast<ms>(generate_end - generate_start)).count() << std::endl;
    }
#endif
#ifdef MEASURE_FILL
    #ifdef MEASURE_MEMORY
    uint64_t fill_before_mem = getMemoryUsage();
    #endif
    std::chrono::high_resolution_clock::time_point fill_start =
    std::chrono::high_resolution_clock::now();
#endif
    graph.fillMatrix();
    MathVector b_vec( graph.getNodesCount(), &program_env);
    b_vec.fillVector();
    graph.createComScheme();
    ComScheme* com_scheme_p = graph.getComScheme();
#ifdef MEASURE_FILL
    #ifdef MEASURE_MEMORY
    uint64_t fill_after_mem = getMemoryUsage();
    if( program_env.getProcessRank() == 0 ){
        std::cout << "Fill memory usage: " << fill_after_mem -
    fill_before_mem << std::endl;
    }
    #endif
    std::chrono::high_resolution_clock::time_point fill_end =
    std::chrono::high_resolution_clock::now();
    if( program_env.getProcessRank() == 0 ){
        std::cout << "Fill time: " << (std::chrono::duration_cast<ms>(fill_end
        - fill_start)).count() << std::endl;
    }
#endif
#ifdef MEASURE_SOLVER
    #ifdef MEASURE_MEMORY
    uint64_t solver_before_mem = getMemoryUsage();
    #endif
    std::chrono::high_resolution_clock::time_point solver_start =
    std::chrono::high_resolution_clock::now();
#endif
    solverCG( graph, b_vec, program_env.isDebugPrint(), CONVERGENCE_EPS,
    &program_env);
#ifdef MEASURE_SOLVER
    #ifdef MEASURE_MEMORY
    uint64_t solver_after_mem = getMemoryUsage();
    std::cout << "Solver memory usage: " << solver_after_mem -
    solver_before_mem << std::endl;
    #endif
    std::chrono::high_resolution_clock::time_point solver_end =
    std::chrono::high_resolution_clock::now();
    if( program_env.getProcessRank() == 0 ){
        std::cout << "Solver time: " << (std::chrono::duration_cast<ms>(solver_end -
    solver_start)).count() << std::endl;
    }
#endif
    std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
#ifdef MEASURE_MEMORY
    std::cout << "Memory usage: " << getMemoryUsage() << std::endl;
#endif
    if( program_env.getProcessRank() == 0 ){
        std::cout << "Time: " << (std::chrono::duration_cast<ms>(end - start)).count() << std::endl;
    }
    if( program_env.isDebugPrint() ){
        graph.printGraph();
        com_scheme_p->print( &program_env);
        //b_vec.printVector();
    }
    MPI_Finalize();
    return 0;
 }
