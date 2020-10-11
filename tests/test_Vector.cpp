#include "../tsk1_vector.h"
/**
 * A module for testing the vector
 */
const double DOUBLE_COMPARISON_ACCURACY = 0.00001;
/**
 * A basic test of a dot product
 * Results:
 *      A control value
 */
static double testDotProduct(){
    MathVector vec_first(2), vec_second(2);
    double res_product = 11;
    vec_first[0] = 1;
    vec_first[1] = 2;
    vec_second[0] = 3;
    vec_second[1] = 4;
    if ( fabs( dotProduct( vec_first, vec_second) - res_product) >=
    DOUBLE_COMPARISON_ACCURACY ){
        std::cout << "A dot product test failed" << std::endl;
    }
    return dotProduct( vec_first, vec_second);
}
/** 
 * A basic test of a linear combination
 * Results:
 *       A control value( a vector sum)
 */
static double testLinearCombination(){    
    MathVector vec_first(2), vec_second(2), res_vec(2);
    double alpha_coeff = 2, beta_coeff = 3;
    vec_first[0] = 1;
    vec_first[1] = 2;
    vec_second[0] = 3;
    vec_second[1] = 4;
    res_vec[0] = 11;
    res_vec[1] = 16;
    MathVector comb_res = linearCombination( vec_first, vec_second,
                                             alpha_coeff, beta_coeff);
    for( size_t vec_idx = 0; vec_idx < 2; ++vec_idx ){
        if ( fabs( res_vec[vec_idx] - comb_res[vec_idx]) >=
    DOUBLE_COMPARISON_ACCURACY ){
        std::cout << "A dot product test failed" << std::endl;
        }
    }
    return comb_res[0] + comb_res[1];
}
/**
 * A basic test of a sparse multiplication
 * Results:
 *      A control value( a vector sum)        
 */
static double testSparseMV(){
    /**
     * Matrix: [0 3] * [5]
     *         [2 0]   [6]
     */
    int* IA = new int[3];
    int* JA = new int[3];
    double* A = new double[3];
    IA[0] = 0;
    IA[1] = 1;
    JA[0] = 1;
    JA[1] = 0;
    A[0] = 3.0;
    A[1] = 2.0;
    NetGraph graph( 2, 2, IA, JA, A);
    MathVector vec_mult(2);
    vec_mult[0] = 5;
    vec_mult[1] = 6;
    MathVector res_vec(2);
    res_vec[0] = 18;
    res_vec[1] = 10;
    MathVector sparse_mult = sparseMV( graph, vec_mult);
    for( size_t vec_idx = 0; vec_idx < 2; ++vec_idx ){
        if ( fabs( res_vec[vec_idx] - sparse_mult[vec_idx]) >=
    DOUBLE_COMPARISON_ACCURACY ){
        std::cout << sparse_mult[vec_idx] << std::endl;
        std::cout << "A sparse MV test failed" << std::endl;
        }
    }
    return sparse_mult[0] + sparse_mult[1];
}
/**
 * Launch all tests
 */
void launchTests(){
    testDotProduct();
    testLinearCombination();
    testSparseMV();
}
