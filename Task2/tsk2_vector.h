#pragma once
#ifndef MATHVECTOR_H
    #define MATHVECTOR_H
#include <cmath>
#include <cassert>
#include <iostream>
#include "omp.h"
#include "tsk2_graph_prepare.h"
#include "tsk2_real.h"
/**
 * A mathematical vector
 */
#include <cassert>
class MathVector{
public:
MathVector( size_t vec_len, ProgramEnv* env_p){
    values_ = new double[vec_len];
    vec_len_ = vec_len;
    env_p_ = env_p;
}
~MathVector(){
    delete values_;
}
MathVector( const MathVector& source){
    vec_len_ = source.getVecLen();
    values_ = new double[vec_len_];
    double* source_values = source.getValues();
    for( size_t vec_idx = 0; vec_idx < vec_len_; ++vec_idx ){
        values_[vec_idx] = source_values[vec_idx];
    }
    env_p_ = source.getEnv();
}
/**
 * Access the values array through the subscript operator.
 */
double& operator[]( size_t vec_idx ){
    assert( vec_idx < vec_len_);
    return values_[vec_idx];
}
const double& operator[]( size_t vec_idx) const{
    assert( vec_idx < vec_len_);
    return values_[vec_idx];
}
double* getValues() const{
    return values_;
}
size_t getVecLen() const{
    return vec_len_;
}
void setVecLen( int vec_len){
    vec_len_ = vec_len;
}
ProgramEnv* getEnv() const{
    return env_p_;
}
/**
 * Fill a vector
 */
void fillVector(){
    for( size_t vec_idx = 0; vec_idx < vec_len_; ++vec_idx){
        values_[vec_idx] = sin( vec_idx);
    }
}
/**
 * Print a vector
 */
void printVector(){
    std::cout << "Vector elements: "; 
    for( size_t vec_idx = 0; vec_idx < vec_len_; ++vec_idx){
        std::cout << values_[vec_idx] << " ";
    }
    std::cout << std::endl;
}
/** 
 * Copy values of another vector
 */
void copyValues( const MathVector& source){
    assert( source.getVecLen() == vec_len_);
    for( size_t vec_idx = 0; vec_idx < vec_len_; ++vec_idx ){
        values_[vec_idx] = source[vec_idx];
    }
}
/** 
 * Get a L2-norm of a vector
 */
double calculateL2(){
    float l2_norm = 0;
    for( size_t vec_idx = 0; vec_idx < vec_len_; ++vec_idx ){
        l2_norm += values_[vec_idx] * values_[vec_idx];
    }
    /**
     * For process with the rank 0, calculate the l2-norm of the whole vector
     * For other processes calculate the l2-norm of their part of the vector
     */
    float global_l2 = 0;
    MPI_Reduce( &l2_norm, &global_l2, 1, MPI_FLOAT, MPI_SUM, 0,
    MPI_COMM_WORLD);
    MPI_Bcast( &global_l2, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    return sqrt( global_l2);
}
private:
    double* values_;
    size_t vec_len_;
    ProgramEnv* env_p_;
};
double dotProduct( MathVector& vec_a, MathVector& vec_b);
MathVector linearCombination( MathVector& vec_a, MathVector& vec_b, 
                   double alpha_coeff, double beta_coeff);
MathVector sparseMV( NetGraph& graph, MathVector& vec);
/**
 * The wrappers over the basic operations.
 * Calculate the time of the basic operations and append it to the time.
 */
double dotProductWithMeasure( MathVector& vec_a, MathVector& vec_b, double&
time);
MathVector linearCombinationWithMeasure( MathVector& vec_a, MathVector& vec_b, 
                   double alpha_coeff, double beta_coeff, double& time);
MathVector sparseMVWithMeasure( NetGraph& graph, MathVector& vec, double& time);
#endif
