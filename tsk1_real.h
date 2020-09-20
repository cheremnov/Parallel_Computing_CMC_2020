// An environemnt
class ProgramEnv{
private:
    // Is a debug print enabled
    bool debug_print_;
public:
    void setDebugPrint( bool debug_print){
        debug_print_ = debug_print;
    }
    bool isDebugPrint(){
        return debug_print_;
    }
    ProgramEnv(): debug_print_( false) {}
};
