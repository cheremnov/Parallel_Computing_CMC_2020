#include <map>
#include <set>
#include "tsk2_real.h"
/**
 * A scheme of communication between two processes
 */
#ifndef COM_H
// Container for the sent and received nodes
typedef std::map<int, std::set<int> > ComArray_t;
    #define COM_H
class ComScheme{
public:
    ComArray_t& getSendToProcess(){
        return send_to_process_;
    }
    ComArray_t& getRecvFromProcess(){
        return recv_from_process_;
    }
    std::vector<int>& getNeighbors(){
        return neighbors_;
    }
    std::vector<int>& getSend(){
        return send_;
    }
    std::vector<int>& getRecv(){
        return recv_;
    }
    std::vector<int>& getSendOffset(){
        return send_offset_;
    }
    std::vector<int>& getRecvOffset(){  
        return recv_offset_;
    }
    void print( ProgramEnv* env_p){
        std::map<int, int>& local_to_global = env_p->getLocalToGlobal();
        for( int processor_rank = 0; processor_rank < env_p->getProcessNum();
                ++processor_rank ){
            // Print the block, located on the processor
            if( env_p->getProcessRank() == processor_rank ){
                std::cout << "Current process: " << processor_rank << std::endl;
                for( int neighbor_idx = 0; neighbor_idx < neighbors_.size();
                ++neighbor_idx ){
                    std::cout << "Neighbor: " << neighbors_[neighbor_idx] <<
                    std::endl << "Sent nodes: ";
                    for( int node_idx = send_offset_[neighbor_idx]; node_idx <
                    send_offset_[neighbor_idx + 1]; ++node_idx ){
                        std::cout << local_to_global[send_[node_idx]] << " ";
                    }
                    std::cout << std::endl << "Received nodes: ";
                    for( int node_idx = recv_offset_[neighbor_idx]; node_idx <
                    recv_offset_[neighbor_idx + 1]; ++node_idx ){
                        std::cout << local_to_global[recv_[node_idx]] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            MPI_Barrier( MPI_COMM_WORLD);
        }
    }
private:
    /**
     * The communication scheme:
     * Which nodes receive from neigbor processes,
     * Which nodes send to neighbor process processes
     */
    ComArray_t send_to_process_;
    ComArray_t recv_from_process_;
    // Neighbors of the current process
    std::vector<int> neighbors_;
    /**
     * Store the received nodes in the recv array
     * Stort the send nodes in the send array
     * The nodes sent to the neighbor process are located
     * from the SendOffset[neigbor_idx] to SendOffset[neighbor_idx+1] - 1
     */
     std::vector<int> send_;
     std::vector<int> recv_;
     std::vector<int> send_offset_;
     std::vector<int> recv_offset_;
};
#endif
