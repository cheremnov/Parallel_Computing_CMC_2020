#include "tsk1_vector.h"
#include "tests/test_Vector.h"
#include <cassert>
/** 
 * Calculate a dot product of the two vectors
 * They must be the same size
 */
double dotProduct( MathVector& vec_a, MathVector& vec_b){
    assert( vec_a.getVecLen() == vec_b.getVecLen());
    if( vec_a.getVecLen() == vec_b.getVecLen() ){
        size_t vec_len = vec_a.getVecLen();
        double sum = 0;
        #pragma omp parallel for reduction( +:sum)
        for( size_t vec_idx = 0; vec_idx < vec_len; ++vec_idx){
            sum += vec_a[vec_idx] * vec_b[vec_idx];
        }
        return sum;
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
    MathVector new_vec( vec_len);
    #pragma omp parallel for
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
    int* IA = graph.getIA(), *JA = graph.getJA();
    double* A = graph.getA();
    size_t nodes_count = graph.getNodesCount();
    size_t edges_count = graph.getEdgesCount();
    size_t vec_len = vec.getVecLen();
    assert( vec_len == nodes_count );
    MathVector new_vec( nodes_count);
    #pragma omp parallel for
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
            new_vec[node_idx] += vec[column_idx] * A[edge_idx]; 
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
    double start_time = omp_get_wtime();
	#endif
    double dot_product = dotProduct( vec_a, vec_b);
	#ifdef MEASURE_VECTOR_OPS
    double end_time = omp_get_wtime();
    time += (end_time - start_time);
	#endif
    return dot_product;
}
MathVector linearCombinationWithMeasure( MathVector& vec_a, MathVector& vec_b, 
                   double alpha_coeff, double beta_coeff, double& time){
	#ifdef MEASURE_VECTOR_OPS				
    double start_time = omp_get_wtime();
	#endif
    MathVector linear_combination = linearCombination( vec_a, vec_b,
    alpha_coeff, beta_coeff);
	#ifdef MEASURE_VECTOR_OPS
    double end_time = omp_get_wtime();
    time += end_time - start_time;
	#endif
    return linear_combination;
}
MathVector sparseMVWithMeasure( NetGraph& graph, MathVector& vec, double& time){
	#ifdef MEASURE_VECTOR_OPS
    double start_time = omp_get_wtime();
	#endif
    MathVector sparse_mv = sparseMV( graph, vec);
	#ifdef MEASURE_VECTOR_OPS
    double end_time = omp_get_wtime();
    time += end_time - start_time;
	#endif
    return sparse_mv;
}