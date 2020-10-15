/**
 * Generate a graph 
 */ 
#include <iostream>
#include <omp.h>
#include "tsk1_graph_prepare.h"
/**
 * Generate a graph of a network from a matrix
 * of a pre-set parameters.
 */
void
NetGraph::generate( MatrixParameters *params_p, 
                    int threads_num = 1){  // A number of threads
    size_t row_len = params_p->getRowLen();
    size_t column_len = params_p->getColumnLen();
    size_t not_divided = params_p->getNotDivided();
    size_t divided = params_p->getDivided();
    /**
     * An edge index in a simple case is a data dependency
     * Can't sum all the edges from the previous rows — a formula is required.
     * Calculate the indexes for the rows independently
     */
    //omp_set_num_threads( threads_num);
    #pragma omp parallel for
    for( size_t row_idx = 0; row_idx <= row_len; ++row_idx )
    {            
        std::pair<int, int> cells = countDividedCells( row_idx, params_p);
        size_t row_not_divided_nodes = cells.first;
        size_t row_divided_nodes = cells.second;
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
        size_t edge_idx = 2 * (row_not_divided_nodes * NETGRAPH_NOT_DIVIDED_EDGES + 
            row_divided_nodes * NETGRAPH_DIVIDED_EDGES);
        if( row_idx > 1 ){
            // The leftmost edge wasn't included, fix it
            edge_idx += 2 * (row_idx - 1);
        }
        if( row_idx > 0 ){
            /**
             * Include all edges between the nodes
             * (row_idx - 1, column_idx) and (row_idx, column_idx)
             */
            edge_idx += column_len;
            // Include other edges from the row_idx only once
            cells = countDividedCells( row_idx + 1, params_p);
            size_t new_row_not_divided_nodes = cells.first - row_not_divided_nodes;
            size_t new_row_divided_nodes = cells.second - row_divided_nodes;
            edge_idx += new_row_not_divided_nodes * NETGRAPH_NOT_DIVIDED_EDGES + 
                new_row_divided_nodes * NETGRAPH_DIVIDED_EDGES + 1;
        }
		// Add edges from the node to itself
		edge_idx += row_idx * (column_len + 1);
        for( size_t column_idx = 0; column_idx <= column_len; ++column_idx )
        { 
            size_t node_idx = row_idx * (column_len + 1) + column_idx;
            IA[node_idx] = edge_idx;
            if( row_idx > 0 ){
                // Edge from (row_idx; column_idx) to (row_idx - 1; column_idx)
                JA[edge_idx] = node_idx - (column_len + 1); 
                A[edge_idx] = 1;
                ++edge_idx;
            }
            if( row_idx > 0 && column_idx < column_len ){
                // Include an edge from a cell division
                size_t upper_cell_idx = (row_idx - 1) * column_len + column_idx;
                if( upper_cell_idx % (divided + not_divided) >= not_divided ){
                    // An edge from (row_idx; column_idx) to (row_idx - 1; column_idx + 1)
                    JA[edge_idx] = node_idx - column_len;
                    A[edge_idx] = 1;
                    ++edge_idx;
                }
            }
            // Edge from (row_idx; column_idx) to (row_idx; column_idx - 1)
            if( column_idx > 0 ){
                JA[edge_idx] = node_idx - 1;
                A[edge_idx] = 1;
                ++edge_idx;
            }
			// Edge from a node to itself
			JA[edge_idx] = node_idx;
			A[edge_idx] = 1;
			++edge_idx;
            // Edge from (row_idx; column_idx) to (row_idx; column_idx + 1)
            if( column_idx < column_len ){
                JA[edge_idx] = node_idx + 1;
                A[edge_idx] = 1;
                ++edge_idx;
            }
            // Look at the cell below
            if( row_idx < row_len && column_idx > 0 ){
                size_t below_cell_idx = row_idx * column_len + column_idx - 1;
                if( below_cell_idx % (divided + not_divided) >= not_divided ){
                    // An edge from (row_idx; column_idx) to (row_idx + 1; column_idx - 1)
                    JA[edge_idx] = node_idx + column_len;
                    A[edge_idx] = 1;
                    ++edge_idx;
                }
            }
            if( row_idx < row_len ){
                // Edge from (row_idx; column_idx) to (row_idx + 1; column_idx)
                JA[edge_idx] = node_idx + (column_len + 1);
                A[edge_idx] = 1;
                ++edge_idx;
            }
        }
        // On the last row write the total number of edges
        if( row_idx == row_len ){
            edges_count_ = edge_idx;
        }
    }
	IA[nodes_count_] = edges_count_;
}
/**
 * Calculate a number of the not-divided and divided cells
 * on the rows [0; row_idx-2)
 * Results:
 *      A pair (not-divided, divided)
 */
std::pair<int,int> 
NetGraph::countDividedCells( size_t row_idx, MatrixParameters* params_p){
    size_t row_len = params_p->getRowLen();
    size_t column_len = params_p->getColumnLen();
    size_t not_divided = params_p->getNotDivided();
    size_t divided = params_p->getDivided();
    /**
     * For the [0; row_idx-2) rows the number of edges can be calculated
     * independent of the column_idx
     */
    size_t row_cells = row_idx >= 1 ? (row_idx - 1) * column_len : 0;
    // How many times the not-divided and divided cells rotated
    size_t row_rotate_cnt = row_cells / (not_divided + divided);
    // Index in a sequence of a current rotation
    size_t rotseq_idx = row_cells - row_rotate_cnt * (not_divided + divided);
    size_t row_not_divided_nodes = row_rotate_cnt * not_divided;
    if( rotseq_idx >= not_divided ){
        // The current rotation already has all not-divided cells
        row_not_divided_nodes += not_divided;
    } else {
        row_not_divided_nodes += rotseq_idx;
    } 
    size_t row_divided_nodes = row_cells - row_not_divided_nodes;
    return std::make_pair( row_not_divided_nodes, row_divided_nodes);
}
/** 
 * Fill the matrix
 * Make it diagonally dominant
 */
void NetGraph::fillMatrix( int threads_num ){ // A number of threads
    const double DOMINANCE_COEFF = 2;
    //omp_set_num_threads( threads_num);
    //
    #pragma omp parallel for
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
                A[edge_idx] = cos( node_idx + neighbor_idx 
                    + node_idx * neighbor_idx );
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
    #pragma omp parallel for
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
    return NetGraph( nodes_count_, nodes_count_, diagonal_IA, diagonal_JA,
    diagonal_A);
}
