/**
 * A class that stores information about the program environment
 */
class ProgramEnv{
private:
    // Is a debug print enabled
    bool debug_print_;
    // Number of threads
    int threads_num_;
public:
    void setDebugPrint( bool debug_print){
        debug_print_ = debug_print;
    }
    bool isDebugPrint(){
        return debug_print_;
    }
    void setThreadsNum( int threads_num){
        threads_num_ = threads_num;
    }
    int getThreadsNum(){
        return threads_num_;
    }
    ProgramEnv(): debug_print_( false), threads_num_( 1){}
};
