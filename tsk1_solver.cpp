#include <iostream>
#include "tsk1_graph_prepare.h"
#include "tsk1_solver.h"
enum { MAX_ITERATIONS = 10000 };
/**
 * A CG solver for a matrix
 */
SolverSolution solverCG( NetGraph& matrix, MathVector& right_part, bool print_debug,
               double convergence_accuracy){ // The convergence accuracy
#ifdef MEASURE_VECTOR_OPS
    double dotproduct_time = 0;
    double linearcombination_time = 0;
    double sparsemv_time = 0;
    double start_measurement = 0;
    double end_measurement = 0;
#endif
    size_t row_count = matrix.getNodesCount();
    size_t not_null_cells = matrix.getEdgesCount();
    // Generate an initial guess vector
    MathVector initial_guess( row_count);
    for( size_t vec_idx = 0; vec_idx < row_count; ++vec_idx){
        initial_guess[vec_idx] = 0;
    }
    bool has_converged = false;
    size_t iteration_num = 1;
    /** 
     * Create a preconditioner matrix from a matrix,
     * borrowing only a diagonal
     */
    NetGraph reverse_preconditioner = matrix.makeDiagonalMatrix( true);
#ifdef MEASURE_VECTOR_OPS
    MathVector current_approximation = sparseMVWithMeasure( matrix, initial_guess,
    sparsemv_time);
#else
    MathVector current_approximation = sparseMV( matrix, initial_guess);
#endif
#ifdef MEASURE_VECTOR_OPS
    MathVector r_iter = linearCombinationWithMeasure( right_part, current_approximation, 1,
    -1, linearcombination_time);
#else
    MathVector r_iter = linearCombination( right_part, current_approximation, 1,
    -1);
#endif
    double rho_prev = 0, rho_iter = 0;
    MathVector p_iter( r_iter.getVecLen());
    // A conjugate gradient algorithm
    while( !has_converged ){
#ifdef MEASURE_VECTOR_OPS
        MathVector z_iter = sparseMVWithMeasure( reverse_preconditioner, r_iter,
        sparsemv_time);
#else
        MathVector z_iter = sparseMV( reverse_preconditioner, r_iter);
#endif
        rho_prev = rho_iter;
#ifdef MEASURE_VECTOR_OPS
        rho_iter = dotProductWithMeasure( r_iter, z_iter, dotproduct_time);
#else
        rho_iter = dotProduct( r_iter, z_iter);
#endif
        if( iteration_num == 1 ){
            p_iter.copyValues( z_iter);
        } else{
            if( !rho_prev ){
                std::cout << "Zero dot product" << std::endl;
                break;
            }
            double b_iter = rho_iter / rho_prev;
            p_iter.copyValues( linearCombination( z_iter, p_iter, 1, b_iter));
        }
#ifdef MEASURE_VECTOR_OPS 
        MathVector q_iter = sparseMVWithMeasure( matrix, p_iter, sparsemv_time);
#else
        MathVector q_iter = sparseMV( matrix, p_iter);
#endif
#ifdef MEASURE_VECTOR_OPS
        double pq_product = dotProductWithMeasure( p_iter, q_iter,
        dotproduct_time);
#else
        double pq_product = dotProduct( p_iter, q_iter);
#endif
        if( !pq_product ){
            std::cout << "Product of p_{k} and q_{k} is zero" << std::endl;
            break;
        }
        double alpha_iter = rho_iter / pq_product;
#ifdef MEASURE_VECTOR_OPS
        initial_guess.copyValues( linearCombinationWithMeasure( initial_guess, p_iter, 1,
        alpha_iter, linearcombination_time));   
#else
        initial_guess.copyValues( linearCombination( initial_guess, p_iter, 1,
        alpha_iter));
#endif
#ifdef MEASURE_VECTOR_OPS
        r_iter.copyValues( linearCombinationWithMeasure( r_iter, q_iter, 
        1, -alpha_iter, linearcombination_time));
#else
        r_iter.copyValues( linearCombination( r_iter, q_iter, 1, -alpha_iter));
#endif
        if( print_debug ){
            std::cout << "Iterations:" << iteration_num << " " << rho_iter << std::endl;
        }
        if( rho_iter < convergence_accuracy || iteration_num >= MAX_ITERATIONS ){
            has_converged = true;
        } else{
            iteration_num++;
        }
    }
    if( print_debug ){
        initial_guess.printVector();
        std::cout << "Number of iterations: " << iteration_num << std::endl;
        std::cout << "L2 norm: " << r_iter.calculateL2() << std::endl;
    }
#ifdef MEASURE_VECTOR_OPS
    std::cout << "Dot product time: " << dotproduct_time << std::endl;
    std::cout << "Linear combination time: " << linearcombination_time <<
    std::endl;
    std::cout << "Sparse multiplication time: " << sparsemv_time << std::endl;
#endif
    return SolverSolution( initial_guess, iteration_num, r_iter.calculateL2());
}
