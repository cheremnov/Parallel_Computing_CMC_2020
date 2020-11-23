#include "tsk2_vector.h"
#include "tests/test_Vector.h"
#include <mpi.h>
#include <chrono>
#include <cassert>
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;
/** 
 * Calculate a dot product of the two vectors
 * They must be the same size
 */
double dotProduct( MathVector& vec_a, MathVector& vec_b){
    ProgramEnv* env_p = vec_a.getEnv();
    assert( vec_a.getVecLen() == vec_b.getVecLen());
    if( vec_a.getVecLen() == vec_b.getVecLen() ){
        size_t vec_len = vec_a.getVecLen();
        float sum = 0;
        for( size_t vec_idx = 0; vec_idx < vec_len; ++vec_idx){
            sum += vec_a[vec_idx] * vec_b[vec_idx];
        }
        float global_sum = 0;
        // Accumulate the sum on the first processor
        MPI_Reduce( &sum, &global_sum, 1, MPI_FLOAT, MPI_SUM, 0,
        MPI_COMM_WORLD);
        MPI_Bcast( &global_sum, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        return global_sum;
    } else{
        // TODO: Change the error handling
        std::cout << "Two vectors doesn't have the same size" << std::endl;
        return 0;
    }
}
/** 
 * Calculate a linear combination of the two vectors
 * They must be the same size
 */
MathVector linearCombination( MathVector& vec_a, MathVector& vec_b, 
                   double alpha_coeff, double beta_coeff){ // Linear coefficients
    assert( vec_a.getVecLen() == vec_b.getVecLen());
    size_t vec_len = vec_a.getVecLen();
    MathVector new_vec( vec_len, vec_a.getEnv());
    for( size_t vec_idx = 0; vec_idx < vec_len; ++vec_idx){
        new_vec[vec_idx] = alpha_coeff * vec_a[vec_idx] + 
            beta_coeff * vec_b[vec_idx];
    }
    return new_vec;
}
/**
 * Multiply a graph matrix to the vector
 * A graph matrix is in the sparse form
 */
MathVector sparseMV( NetGraph& graph, MathVector& vec){
    /**
     * Use communication scheme
     * to send neighbor nodes and receive halo
     */
    ComScheme* com_scheme_p = graph.getComScheme();
    assert( com_scheme_p);
    //std::map<int, int>& global_to_local = graph.getEnv()->getGlobalToLocal();
    // Values of the received halo
    std::map<int, float> received_vec;
    std::vector<int>& send = com_scheme_p->getSend();
    std::vector<int>& send_offset = com_scheme_p->getSendOffset();
    std::vector<int>& recv = com_scheme_p->getRecv();
    std::vector<int>& recv_offset = com_scheme_p->getRecvOffset();
    std::vector<int>& neighbors = com_scheme_p->getNeighbors();
    float* sent_vec = new float[graph.getNodesCount()];
    MPI_Request* request = new MPI_Request[recv.size() + send.size()];
    MPI_Status* status = new MPI_Status[recv.size() + send.size()];
    for( int neighbor_idx = 0; neighbor_idx < neighbors.size();
        ++neighbor_idx ){
        // Send the nodes to the neighbors
        for( int node_idx = send_offset[neighbor_idx]; node_idx <
        send_offset[neighbor_idx + 1]; ++node_idx ){
            int local_node_idx = send[node_idx];
            sent_vec[local_node_idx] = vec[local_node_idx];
            MPI_Isend( &sent_vec[local_node_idx], 1, MPI_FLOAT, neighbors[neighbor_idx], 0,
            MPI_COMM_WORLD, &request[node_idx]);
        }
        // Get the halo from the neighbors
        for( int node_idx = recv_offset[neighbor_idx]; node_idx <
        recv_offset[neighbor_idx + 1]; ++node_idx ){
            int local_node_idx = recv[node_idx];
            received_vec[local_node_idx] = 0;
            // Assume the pointer to the map element is stable
            MPI_Irecv( &received_vec[local_node_idx], 1, MPI_FLOAT,
            neighbors[neighbor_idx], 0, MPI_COMM_WORLD,
            &request[send.size()+node_idx]);
        }
    }
    MPI_Waitall( recv.size() + send.size(), request, status);
    delete request;
    delete status;
    delete sent_vec;
    int* IA = graph.getIA(), *JA = graph.getJA();
    double* A = graph.getA();
    size_t nodes_count = graph.getNodesCount();
    size_t edges_count = graph.getEdgesCount();
    size_t vec_len = vec.getVecLen();
    assert( vec_len == nodes_count );
    MathVector new_vec( nodes_count, vec.getEnv());
    for( size_t node_idx = 0; node_idx < nodes_count; ++node_idx){
        new_vec[node_idx] = 0;
        /**
         * For the rightmost node IA doesn't specify the edges index.
         * Instead, get it from the graph
         */
        size_t end_idx = node_idx + 1 < nodes_count ? IA[node_idx + 1] :
            edges_count;
        for( size_t edge_idx = IA[node_idx]; edge_idx < end_idx; ++edge_idx){
            size_t column_idx = JA[edge_idx];
            double vec_elem = 0;
            if( column_idx < vec.getVecLen() ){
                vec_elem = vec[column_idx];
            } else{
                vec_elem = received_vec[column_idx];
            }
            new_vec[node_idx] += vec_elem * A[edge_idx]; 
        }
    }
    return new_vec;
}
/**
 * The wrappers over the basic operations.
 * Calculate the time of the basic operations and append it to the time.
 * The time is calculated only if MEASURE_VECTOR_OPS define is set
 */
double dotProductWithMeasure( MathVector& vec_a, MathVector& vec_b, double&
time){
	#ifdef MEASURE_VECTOR_OPS
    std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	#endif
    double dot_product = dotProduct( vec_a, vec_b);
	#ifdef MEASURE_VECTOR_OPS
	std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds new_time = std::chrono::duration_cast<ms>(t1 - t0);
    std::cout << new_time.count() << std::endl;
    time += new_time.count();
	#endif
    return dot_product;
}
MathVector linearCombinationWithMeasure( MathVector& vec_a, MathVector& vec_b, 
                   double alpha_coeff, double beta_coeff, double& time){
	#ifdef MEASURE_VECTOR_OPS				
    std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	#endif
    MathVector linear_combination = linearCombination( vec_a, vec_b,
    alpha_coeff, beta_coeff);
	#ifdef MEASURE_VECTOR_OPS
    std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds new_time = std::chrono::duration_cast<ms>(t1 - t0);
    std::cout << new_time.count() << std::endl;
    time += new_time.count();
	#endif
    return linear_combination;
}
MathVector sparseMVWithMeasure( NetGraph& graph, MathVector& vec, double& time){
	#ifdef MEASURE_VECTOR_OPS
    std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	#endif
    MathVector sparse_mv = sparseMV( graph, vec);
	#ifdef MEASURE_VECTOR_OPS
    std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds new_time = std::chrono::duration_cast<ms>(t1 - t0);
    std::cout << new_time.count() << std::endl;
    time += new_time.count();
	#endif
    return sparse_mv;
}
