#include <iostream>
#include "tsk2_graph_prepare.h"
#include "tsk2_solver.h"
enum { MAX_ITERATIONS = 10000 };
/**
 * A CG solver for a matrix
 */
SolverSolution solverCG( NetGraph& matrix, MathVector& right_part, bool print_debug,
               double convergence_accuracy,  // The convergence accuracy
               ProgramEnv* env_p){
	/**
	 * The variables for time measurement
	 * Matter only if the time of the basic operations is measure.
	 */
    double dotproduct_time = 0;
    double linearcombination_time = 0;
    double sparsemv_time = 0;
    double start_measurement = 0;
    double end_measurement = 0;
	// Matrix information
    size_t row_count = matrix.getNodesCount();
    size_t not_null_cells = matrix.getEdgesCount();
    // Generate an initial guess vector
    MathVector initial_guess( row_count, env_p);
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
    MathVector current_approximation = sparseMVWithMeasure( matrix, initial_guess,
    sparsemv_time);
    MathVector r_iter = linearCombinationWithMeasure( right_part, current_approximation, 1,
    -1, linearcombination_time);
    double rho_prev = 0, rho_iter = 0;
    MathVector p_iter( r_iter.getVecLen(), env_p);
    // A conjugate gradient algorithm
    while( !has_converged ){
        MathVector z_iter = sparseMVWithMeasure( reverse_preconditioner, r_iter,
        sparsemv_time);
        rho_prev = rho_iter;
        rho_iter = dotProductWithMeasure( r_iter, z_iter, dotproduct_time);
        if( iteration_num == 1 ){
            p_iter.copyValues( z_iter);
        } else{
            if( !rho_prev ){
                std::cout << "Zero dot product" << std::endl;
                break;
            }
            double b_iter = rho_iter / rho_prev;
            p_iter.copyValues( linearCombinationWithMeasure( z_iter, p_iter, 1, b_iter, linearcombination_time));
        }
        MathVector q_iter = sparseMVWithMeasure( matrix, p_iter, sparsemv_time);
        double pq_product = dotProductWithMeasure( p_iter, q_iter,
        dotproduct_time);
        if( !pq_product ){
            std::cout << "Product of p_{k} and q_{k} is zero" << std::endl;
            break;
        }
        double alpha_iter = rho_iter / pq_product;
        initial_guess.copyValues( linearCombinationWithMeasure( initial_guess, p_iter, 1,
        alpha_iter, linearcombination_time));   
        r_iter.copyValues( linearCombinationWithMeasure( r_iter, q_iter, 
        1, -alpha_iter, linearcombination_time));
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
