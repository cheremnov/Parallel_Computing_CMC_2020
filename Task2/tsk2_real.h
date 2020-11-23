#ifndef REAL_H
    #define REAL_H
#include <map>
#include <vector>
/**
 * A class that stores information about the program environment
 */
class ProgramEnv{
private:
    // Is a debug print enabled
    bool debug_print_;
    // A number of processes
    int process_num_;
    // A current process
    int process_rank_;
    // Convertation of node array from local to global and vice versa
    std::map<int,int> local_to_global_;
    std::map<int,int> global_to_local_;
    // Which processor holds the vertex
    std::vector<int> parts_;
    // The block parameters
    size_t start_row_idx_;
    size_t end_row_idx_;
    size_t start_column_idx_;
    size_t end_column_idx_;
public:
    void setDebugPrint( bool debug_print){
        debug_print_ = debug_print;
    }
    bool isDebugPrint(){
        return debug_print_;
    }
    void setProcessNum( int process_num){
        process_num_ = process_num;
    }
    int getProcessNum(){
        return process_num_;
    }
    void setProcessRank( int process_rank){
        process_rank_ = process_rank;
    }
    int getProcessRank(){
        return process_rank_;
    }
    void setBlockParams( size_t start_row_idx, size_t end_row_idx, 
        size_t start_column_idx, size_t end_column_idx){
        start_row_idx_ = start_row_idx;
        start_column_idx_ = start_column_idx;
        end_row_idx_ = end_row_idx;
        end_column_idx_ = end_column_idx;
    }
    size_t getStartRow(){
        return start_row_idx_;
    }
    size_t getEndRow(){
        return end_row_idx_;
    }
    size_t getStartColumn(){
        return start_column_idx_;
    }
    size_t getEndColumn(){
        return end_column_idx_;
    }
    std::map<int, int>& getGlobalToLocal(){
        return global_to_local_;
    }
    std::map<int, int>& getLocalToGlobal(){
        return local_to_global_;
    }
    std::vector<int>& getParts(){
        return parts_;
    }
    ProgramEnv(): debug_print_( false), process_num_(1), process_rank_(0) {}
};
#endif
