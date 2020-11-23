#include "../tsk2_vector.h"
/**
 * A module for testing the vector
 */
const double DOUBLE_COMPARISON_ACCURACY = 0.00001;
/**
 * A basic test of a dot product
 * Results:
 *      A control value
 */
static double testDotProduct( ProgramEnv* env_p ){
    double res_product = 11;
    if( env_p->getProcessNum() < 2 ){
        MathVector vec_first(2, env_p), vec_second(2, env_p);
        vec_first[0] = 1;
        vec_first[1] = 2;
        vec_second[0] = 3;
        vec_second[1] = 4;
        if ( fabs( dotProduct( vec_first, vec_second) - res_product) >=
        DOUBLE_COMPARISON_ACCURACY ){
            std::cout << "A dot product test failed" << std::endl;
        }
        return dotProduct( vec_first, vec_second);
    } else{
        MathVector vec_first(1, env_p), vec_second(1, env_p);
        if( env_p->getProcessRank() == 0){
            vec_first[0] = 1;
            vec_second[0] = 3;
        } else if( env_p->getProcessRank() == 1){
            vec_first[0] = 2;
            vec_second[0] = 4;
        }
        MPI_Barrier( MPI_COMM_WORLD);
        double dot_product = dotProduct( vec_first, vec_second);
        if ( fabs( dot_product - res_product) >=
        DOUBLE_COMPARISON_ACCURACY ){
            std::cout << "A dot product test failed" << std::endl;
        }
        return dot_product;
    }
}
/** 
 * A basic test of a linear combination
 * Results:
 *       A control value( a vector sum)
 */
static double testLinearCombination( ProgramEnv* env_p){    
    double res_l2 = sqrt( 11 * 11 + 16 * 16);
    double res_comb_l2 = 0;
    double alpha_coeff = 2, beta_coeff = 3;
    if( env_p->getProcessNum() < 2 ){
        MathVector vec_first(2, env_p), vec_second(2, env_p);
        vec_first[0] = 1;
        vec_first[1] = 2;
        vec_second[0] = 3;
        vec_second[1] = 4;
        MathVector comb_res = linearCombination( vec_first, vec_second,
                                                 alpha_coeff, beta_coeff);
        res_comb_l2 = comb_res.calculateL2();
    } else{
        if( env_p->getProcessRank() < 2 ){
            MathVector vec_first(1, env_p), vec_second(1, env_p);
            if( env_p->getProcessRank() == 0 ){
                vec_first[0] = 1;
                vec_second[0] = 3;
            } else{
                vec_first[0] = 2;
                vec_second[0] = 4;
            }
            MathVector comb_res = linearCombination( vec_first, vec_second,
                                                 alpha_coeff, beta_coeff);
            res_comb_l2 = comb_res.calculateL2();
        } else{
            MathVector vec_first(0, env_p), vec_second(0, env_p);
            MathVector comb_res = linearCombination( vec_first, vec_second,
                                                 alpha_coeff, beta_coeff);
            res_comb_l2 = comb_res.calculateL2();
        }
    }
    if( env_p->getProcessRank() == 0 ){
        if ( fabs( res_l2 - res_comb_l2) >=
        DOUBLE_COMPARISON_ACCURACY || fabs( res_comb_l2) <
        DOUBLE_COMPARISON_ACCURACY ){
            std::cout << "A linear combination test failed" << std::endl;
        }
    }
    return fabs( res_l2 - res_comb_l2);
}
/**
 * A basic test of a sparse multiplication
 * Results:
 *      A control value( a vector sum)        
 */
static double testSparseMV( ProgramEnv* env_p ){
    /**
     * Matrix: [0 3] * [5]
     *         [2 0]   [6]
     */
    int* IA = new int[3];
    int* JA = new int[3];
    double* A = new double[3];
    if( env_p->getProcessNum() < 2 ){
        IA[0] = 0;
        IA[1] = 1;
        JA[0] = 1;
        JA[1] = 0;
        A[0] = 3.0;
        A[1] = 2.0;
        NetGraph graph( 2, 2, IA, JA, A, env_p, new ComScheme);
        MathVector vec_mult(2, env_p);
        vec_mult[0] = 5;
        vec_mult[1] = 6;
        double res_l2 = sqrt( 18*18 + 10*10);
        MathVector sparse_mult = sparseMV( graph, vec_mult);
        if ( fabs( sparse_mult.calculateL2() - res_l2) >=
    DOUBLE_COMPARISON_ACCURACY ){
            std::cout << "A sparse MV test failed" << std::endl;
        }
        return sparse_mult[0] + sparse_mult[1];
    } else{
        return 0;
    }
}
/**
 * Launch all tests
 */
void launchTests( ProgramEnv* env_p){
    testDotProduct( env_p);
    testLinearCombination( env_p);
    testSparseMV( env_p);
}
