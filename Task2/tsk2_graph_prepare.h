#ifndef GRAPH_PREPARE_H
    #define GRAPH_PREPARE_H
#include <iostream>
#include <cstddef>
#include <mpi.h>
#include <cmath>
#include "tsk2_real.h"
#include "tsk2_com.h"
enum { 
    NETGRAPH_MAX_EDGES_NODE = 6,
    NETGRAPH_NOT_DIVIDED_EDGES = 2,
    NETGRAPH_DIVIDED_EDGES = 3
};
class MatrixParameters{
public:
    typedef enum{ 
        // A number of columns and row can't exceed this
        MAX_MATRIX_DIMENSION = 10000,
        MAX_CELLS = MAX_MATRIX_DIMENSION * MAX_MATRIX_DIMENSION
    } MatrixConstraints_t;
    size_t getRowLen(){
        return row_len_;
    }
    size_t getColumnLen(){
        return column_len_;
    }
    size_t getNotDivided(){
        return not_divided_;
    }
    size_t getDivided(){
        return divided_;
    }
    size_t getBlockRows(){
        return block_rows_;
    }
    size_t getBlockColumns(){
        return block_columns_;
    }
    /**
     * Set the parameters in one batch
     */
    void setParameters( size_t row_len, size_t column_len, size_t not_divided,
    size_t divided, size_t block_rows, size_t block_columns ){
        row_len_ = row_len;
        column_len_ = column_len;
        not_divided_ = not_divided;
        divided_ = divided;
        block_rows_ = block_rows;
        block_columns_ = block_columns;
    }
    MatrixParameters( size_t row_len, size_t column_len, size_t not_divided,
    size_t divided, size_t block_rows, size_t block_columns): 
        row_len_( row_len), column_len_( column_len),
        not_divided_( not_divided), divided_( divided), 
        block_rows_ (block_rows), block_columns_ (block_columns) {}
    MatrixParameters(): 
        row_len_( 0), column_len_( 0),
        not_divided_( 0), divided_( 0) {}
private:
    // How many block rows and block columns exist 
    size_t block_rows_;
    size_t block_columns_;

    // A number of rows
    size_t row_len_;
    // A number of columns
    size_t column_len_;
    /**
     * Every matrix cell can be cut in half or not
     * It depends on the pattern: first the cells are not divided,
     * then divided, then again not divided, and so on.
     * These parameteres dependent how many divided and not divided cells
     * are in one rotation
     */
    // A number of not cut in half cells
    size_t not_divided_;
    // A number of cut in half cells
    size_t divided_;
};
class NetGraph{
/** 
 * A matrix describes the graph.
 * Matrix rows are the graph nodes.
 */
public:
NetGraph( MatrixParameters *params_p, 
          ProgramEnv* env_p ){       // program environment
    env_p_ = env_p;
    // Get the block boundaries
    size_t start_row_idx = env_p_->getStartRow();
    size_t end_row_idx = env_p->getEndRow();
    size_t start_column_idx = env_p_->getStartColumn();
    size_t end_column_idx = env_p->getEndColumn();

    nodes_count_ = (end_row_idx - start_row_idx) * (end_column_idx -
    start_column_idx);
    size_t edges_count_max = nodes_count_ * NETGRAPH_MAX_EDGES_NODE;
    IA = new int[nodes_count_ + 1];
    JA = new int[2 * edges_count_max];
    A = new double[2 * edges_count_max];
    // The communication scheme is not prepared
    com_scheme_p_ = nullptr;
}
/** 
 * Not generate a graph(a matrix) conventional way
 * Instead, manually insert already prepared
 */
NetGraph( size_t nodes_count, size_t edges_count, int* IA, int *JA, double *A,
          ProgramEnv* env_p, ComScheme* scheme_p){        // program environment
    nodes_count_ = nodes_count;
    edges_count_ = edges_count;
    this->IA = IA;
    this->JA = JA;
    this->A = A;
    env_p_ = env_p;
    com_scheme_p_ = scheme_p;
}
~NetGraph(){
    delete[] IA;
    // Not-null matrix columns
    delete[] JA;
    // An array that stores matrix coefficients
    delete[] A;
    if( com_scheme_p_ ){
        delete com_scheme_p_;
    }
}
void printGraph(){
    std::map<int, int>& local_to_global = env_p_->getLocalToGlobal();
    for( int processor_rank = 0; processor_rank < env_p_->getProcessNum();
            ++processor_rank ){
        // Print the block, located on the processor
        if( env_p_->getProcessRank() == processor_rank ){
            for( int node_idx = 0; node_idx < nodes_count_; ++node_idx ){
                std::cout << "Node index: " << local_to_global[node_idx] << std::endl << "Edges to: ";        
                size_t end_idx = node_idx + 1 < nodes_count_ ? IA[node_idx + 1] :
                    edges_count_;
                for( int edge_idx = IA[node_idx]; edge_idx < end_idx; ++edge_idx){
                    std::cout << local_to_global[JA[edge_idx]] << " ";
                }
                std::cout << std::endl << "Matrix elements: ";
                for( int edge_idx = IA[node_idx]; edge_idx < end_idx; ++edge_idx){
                    std::cout << A[edge_idx] << " ";
                }
                std::cout << std::endl;
            }
        }
        MPI_Barrier( MPI_COMM_WORLD);
    }
}
size_t getNodesCount(){
	return nodes_count_;
}
size_t getEdgesCount(){
    return edges_count_;
}
int* getIA(){
	return IA;
}
int* getJA(){
	return JA;
}
double* getA(){
    return A;
}
ComScheme* getComScheme(){
    return com_scheme_p_;
}
void generate( MatrixParameters *params_p);
void fillMatrix();
NetGraph makeDiagonalMatrix( bool is_reverse);
void createComScheme();
private:
    /** 
     * JA and A stores information about all rows.
     * To get information about a single row,
     * you have to know the index, where
     * JA and A store it.
     */
    int* IA;
    int* JA;
    double* A;
    size_t nodes_count_;
    // A number of the edges( not-null cells) in the graph
    size_t edges_count_;
    ProgramEnv* env_p_;
    ComScheme* com_scheme_p_;
};
#endif
