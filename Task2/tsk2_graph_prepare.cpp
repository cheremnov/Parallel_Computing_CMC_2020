/**
 * Generate a graph 
 */ 
#include <cassert>
#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <algorithm>
#include "tsk2_graph_prepare.h"
#include "tsk2_com.h"
/**
 * Generate a graph of a network from a matrix
 * of a pre-set parameters.
 */
void
NetGraph::generate( MatrixParameters *params_p){ 
    size_t row_len = params_p->getRowLen();
    size_t column_len = params_p->getColumnLen();
    size_t not_divided = params_p->getNotDivided();
    size_t divided = params_p->getDivided();
    /**
     * Load a parameters of the net block
     * from the program environment
     */
    size_t start_row_idx = env_p_->getStartRow();
    size_t end_row_idx = env_p_->getEndRow();
    size_t start_column_idx = env_p_->getStartColumn();
    size_t end_column_idx = env_p_->getEndColumn();
    std::map<int, int>& global_to_local = env_p_->getGlobalToLocal();
    /**
     * An edge index in a simple case is a data dependency
     * Can't sum all the edges from the previous rows — a formula is required.
     * In the block the edges from the previous column may be inavailable,
     * so create formula for each cell.
     * Calculate the indexes for each cell independently
     */
    // A number of edges in the block
    size_t block_edge_idx = 0;
    for( size_t row_idx = start_row_idx; row_idx < end_row_idx; ++row_idx )
    {            
        /**
         * N------N-------N---------N    row_idx
         * |      |       |         |
         * N      N       N         N    row_idx + 1
         * Calculate a number of edges for a matrix row.
         * Include edges, for which any node is at the row_idx row
         * For each divided cell there are 4 edges, for each non-divided — 3. 
         * But one of the edges was already included in the left cell
         * So, calculate as if it had 3 and 2 edges accordingly.
         */
        for( size_t column_idx = start_column_idx; column_idx < end_column_idx; ++column_idx )
        {
			size_t node_idx = row_idx * (column_len + 1) + column_idx;
            assert( global_to_local[node_idx] < nodes_count_);
            IA[global_to_local[node_idx]] = block_edge_idx;
			
            if( row_idx > 0 ){
                // Edge from (row_idx; column_idx) to (row_idx - 1; column_idx)
                JA[block_edge_idx] = 
                    global_to_local[node_idx - (column_len + 1)]; 
                A[block_edge_idx] = 1;
                ++block_edge_idx;
            }
            if( row_idx > 0 && column_idx < column_len ){
                // Include an edge from a cell division
                size_t upper_cell_idx = (row_idx - 1) * column_len + column_idx;
				// Cell is divided
                if( upper_cell_idx % (divided + not_divided) >= not_divided ){
                    // An edge from (row_idx; column_idx) to (row_idx - 1; column_idx + 1)
                    JA[block_edge_idx] = global_to_local[node_idx - column_len];
                    A[block_edge_idx] = 1;
                    ++block_edge_idx;
                }
            }
            // Edge from (row_idx; column_idx) to (row_idx; column_idx - 1)
            if( column_idx > 0 ){
                JA[block_edge_idx] = global_to_local[node_idx - 1];
                A[block_edge_idx] = 1;
                ++block_edge_idx;
            }
			// Edge from a node to itself
			JA[block_edge_idx] = global_to_local[node_idx];
			A[block_edge_idx] = 1;
            ++block_edge_idx;
            // Edge from (row_idx; column_idx) to (row_idx; column_idx + 1)
            if( column_idx < column_len ){
                JA[block_edge_idx] = global_to_local[node_idx + 1];
                A[block_edge_idx] = 1;
                ++block_edge_idx;
            }
            // Look at the cell below
            if( row_idx < row_len && column_idx > 0 ){
                size_t below_cell_idx = row_idx * column_len + column_idx - 1;
                if( below_cell_idx % (divided + not_divided) >= not_divided ){
                    // An edge from (row_idx; column_idx) to (row_idx + 1; column_idx - 1)
                    JA[block_edge_idx] = global_to_local[node_idx + column_len];
                    A[block_edge_idx] = 1;
                    ++block_edge_idx;
                }
            }
            if( row_idx < row_len ){
                // Edge from (row_idx; column_idx) to (row_idx + 1; column_idx)
                JA[block_edge_idx] = global_to_local[node_idx + (column_len + 1)];
                A[block_edge_idx] = 1;
                ++block_edge_idx;
            }
			
        }
    }
    edges_count_ = block_edge_idx;
	IA[nodes_count_] = edges_count_;
}

/** 
 * Fill the matrix
 * Make it diagonally dominant
 */
void NetGraph::fillMatrix(){
    std::map<int, int>& local_to_global = env_p_->getLocalToGlobal();
    const double DOMINANCE_COEFF = 2;
    for( size_t node_idx = 0; node_idx < nodes_count_; ++node_idx){
        // A sum of all matrix cells on the row, except for the diagonal
        double row_sum = 0;
        /**
         * For the rightmost node IA doesn't specify the edges index.
         * Instead, get it from the graph
         */
        size_t end_idx = node_idx + 1 < nodes_count_ ? IA[node_idx + 1] :
            edges_count_;
        /** 
         * A position of a diagonal element:
         * of a node that has an edge to itself
         */
        size_t diagonal_idx = 0;
        for( size_t edge_idx = IA[node_idx]; edge_idx < end_idx; ++edge_idx){
            // An other node of the edge
            size_t neighbor_idx = JA[edge_idx];
            if( neighbor_idx == node_idx ){
                diagonal_idx = edge_idx;
            } else{
                size_t global_node_idx = local_to_global[node_idx];
                size_t global_neighbor_idx = local_to_global[neighbor_idx];
                A[edge_idx] = cos( global_node_idx + global_neighbor_idx
                    + global_node_idx * global_neighbor_idx );
                row_sum += fabs( A[edge_idx]);
            }
        }
        A[diagonal_idx] = DOMINANCE_COEFF * row_sum;
    }
}
/**
 * Make a diagonal matrix from the current graph
 * If is_reverse = true, then make it reversed
 */
NetGraph NetGraph::makeDiagonalMatrix( bool is_reverse){
    int* diagonal_IA = new int[nodes_count_];
    int* diagonal_JA = new int[nodes_count_];
    double* diagonal_A = new double[nodes_count_];
    for( size_t node_idx = 0; node_idx < nodes_count_; ++node_idx ){
        diagonal_IA[node_idx] = node_idx;
        /**
         * For the rightmost node IA doesn't specify the edges index.
         * Instead, get it from the graph
         */
        size_t end_idx = node_idx + 1 < nodes_count_ ? IA[node_idx + 1] :
            edges_count_;
        for( size_t edge_idx = IA[node_idx]; edge_idx < end_idx;
        ++edge_idx){
            if(JA[edge_idx] == node_idx ){
                diagonal_JA[node_idx] = node_idx;
                if( is_reverse ){
                    diagonal_A[node_idx] = 1.0 / A[edge_idx];
                } else{
                    diagonal_A[node_idx] = A[edge_idx];
                }
            }
        }
    }
    ComScheme* com_scheme_p = new ComScheme;
    return NetGraph( nodes_count_, nodes_count_, diagonal_IA, diagonal_JA,
    diagonal_A, env_p_, com_scheme_p);
}
/**
 * Create a communication scheme for the graph
 * Warning:
 *     Free the communication sc
 */
void
NetGraph::createComScheme(){
    com_scheme_p_ = new ComScheme;  
    std::map<int, int>& global_to_local = env_p_->getGlobalToLocal();
    std::map<int, int>& local_to_global = env_p_->getLocalToGlobal();
    std::vector<int>& parts = env_p_->getParts();
    int processor_rank = env_p_->getProcessRank();
    ComArray_t& send_to_process = com_scheme_p_->getSendToProcess();
    ComArray_t& recv_from_process = com_scheme_p_->getRecvFromProcess();
    std::vector<int>& neighbors = com_scheme_p_->getNeighbors();
    /**
     * Add the halo to recv_from_process array.
     * Add the local nodes, that neighbor halo, to the send_to_process array.
     */
    for( size_t node_idx = 0; node_idx < nodes_count_; ++node_idx ){
        size_t end_idx = node_idx + 1 < nodes_count_ ? IA[node_idx + 1] :
            edges_count_;
        for( size_t edge_idx = IA[node_idx]; edge_idx < end_idx;
        ++edge_idx){
            size_t neighbor_node = JA[edge_idx];
            size_t neighbor_node_process = parts[neighbor_node];
            if( neighbor_node_process != processor_rank ){
                // The node belongs to another process
                send_to_process[neighbor_node_process].insert(
                local_to_global[node_idx]);
                recv_from_process[neighbor_node_process].insert( 
                local_to_global[neighbor_node]);
            }
        }
    }
    // Fill neighbors array
    for( auto send_iter = send_to_process.begin(); send_iter != 
        send_to_process.end(); ++send_iter){
        neighbors.push_back( send_iter->first);
    }
    std::sort( neighbors.begin(), neighbors.end());
    std::vector<int>& send = com_scheme_p_->getSend();
    std::vector<int>& recv = com_scheme_p_->getRecv();
    std::vector<int>& send_offset = com_scheme_p_->getSendOffset();
    std::vector<int>& recv_offset = com_scheme_p_->getRecvOffset();
    send_offset.push_back( 0);
    recv_offset.push_back( 0);
    for( int neighbor_idx = 0; neighbor_idx < neighbors.size(); ++neighbor_idx
    ){
        int neighbor_process = neighbors[neighbor_idx];
        for( int global_node_idx: send_to_process[neighbor_process] ){
            send.push_back( global_to_local[global_node_idx]);
        }
        for( int global_node_idx: recv_from_process[neighbor_process] ){
            recv.push_back( global_to_local[global_node_idx]);
        }
        send_offset.push_back( send.size());
        recv_offset.push_back( recv.size());
    }
}
