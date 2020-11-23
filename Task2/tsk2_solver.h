#include "tsk2_vector.h"
// A solver result
class SolverSolution{
public:
    SolverSolution( MathVector approximate_solution, int iterations_number,
        double solution_l2): approximate_solution_( approximate_solution),
        iterations_number_(iterations_number), solution_l2_(solution_l2) {}
    MathVector getApproximateSolution(){
        return approximate_solution_;
    }
    int getIterationsNumber(){
        return iterations_number_;
    }
    double getSolutionL2(){
        return solution_l2_;
    }
private:
    // The approximate vector solution
    MathVector approximate_solution_;
    // A number of iterations for solving the system
    int iterations_number_;
    // An l2 norm of the solution
    double solution_l2_;
};

SolverSolution solverCG( NetGraph& matrix, MathVector& right_part, bool print_debug,
               double convergence_accuracy, ProgramEnv* env_p);
