#include <cmath>
/**
 * A mathematical vector
 */
class MathVector{
public:
MathVector( size_t vec_len){
    values_ = new double[vec_len];
    vec_len_ = vec_len;
}
~MathVector(){
    delete values_;
}
int getVecLen(){
    return vec_len_;
}
void setVecLen( int vec_len){
    vec_len_ = vec_len;
}
/**
 * Fill a vector
 */
void fillVector( int threads_num = 1){ // A number of threads
    omp_set_num_threads( threads_num);
    #pragma omp parallel for
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
private:
    double* values_;
    size_t vec_len_;
};
