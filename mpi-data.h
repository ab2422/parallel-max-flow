#ifndef MPI_STRUCTS_CODE_H
#define MPI_STRUCTS_CODE_H

#include <mpi.h>
#include <vector>
#include <queue>

struct network {
    int src;  // index of src
    int sink; // index of sink
    int n;    // total num verts
    int m;    // num edges
    int npp;  // nodes on THIS processor
    int std_npp; // ''std'' # nodes/proc
    int *odeg; // array of out degs
    int *ideg; // array of in degs
    
    
    // adj[0][v] = {w_1, j_1, w_2, j_2 ...} 
    //          for w_i adj to v, where j_i is ind of v in w's adj[1] list
    // adj[1][v] = {v_1, i_1, v_2, i_2, ...}
    //          for w adj to v_i, where i_i is ind of w in v's adj[0] list
    std::vector<std::vector<std::vector<int>>> adj; 
    
    //holds capacities: cap[v][w] (local coords) is capacity of v->w
    std::vector<std::vector<int>> cap;
};

struct resgraph : public network {
    int s_proc;
    int t_proc;
    int *ex;  // array of excess
    int *dex; // array of delta excess
    int *hght;// array of heights (d)
    int *hght_p; // array of new heights (d')

    // flow[0][v] = {f(v,w_1), f(v,w_2), ...} for same order as adj[0][v]
    // flow[1][w] = {f(v_1,w), f(v_2,w), ...} for same order as adj[1][w]
    // note: flow[0] all >=0, flow[1] all <=0
    std::vector<std::vector<std::vector<int>>> flow;

    // adj_d[0][v] = list of heights of adj verts, same order as adj[0]
    // adj_d[1][v] = list of hghts of bwd adj verts, same order as adj[1]
    std::vector<std::vector<std::vector<int>>> adj_d;

    int n_act;
    // active is in LOCAL coordinates
    std::queue<int> active;
    std::queue<int> active_p;
};

struct comm_data {
    // for send requests: one for every edge
    MPI_Request ***edge_req;
    std::vector<std::vector<std::vector<int>>> edge_bi;
    // OUT_FLAG[V][I] ==  nothing WHEN (V,NET->ADJ[3*I]) HAS NO PENDING OPERATIONS
    // == 0 WHEN PENDING A SENT QUERY, == 1 FOR REC QUERY
    // == 2 FOR SENT RESPONSE, == 3 FOR REC RESPONSE
    std::vector<std::vector<std::vector<unsigned char>>> edge_flag;

    // fordistance comm: one for every edge
    MPI_Request ***dist_req;
    std::vector<std::vector<std::vector<int>>> dist_bi;

    // for completion checking comm: one for each proc
    MPI_Request *fin_req;
    bool *proc_done;
    bool ring_flag; // true if send in progress
    unsigned char ring_stat;
    int num_times_done;
    int all_done;
    int prev_total;
    int nxt_total;
    MPI_Request ring_req;


    // for all
    int buff_size;
    int **buff;
    std::queue<int> avail;
    int max_deg;
    int *arr_of_inds;
};


#endif
