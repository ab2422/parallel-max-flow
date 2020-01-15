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
    // adj[v] = {w_1, c(v,w_1), j_1, w_2, c(v,w_2), j_2 ...} 
    //          for w_i adj to v, where j_i is ind of v in w's badj list
    std::vector<std::vector<int>> adj; 
    // badj[w] = {v_1, j_1, v_2, j_2, ...} for w adj to v_i, j_i ind of v
    std::vector<std::vector<int>> badj;
    // each proc stores ceil( n/size ) vertices. So, if n = 100, size = 10, rk p = 1, then it stores
    // verts 10-19, and to access vert 12, use adj[2], badj[2]. The adj list holds actual vert #s.
};

struct resgraph : public network {
    int s_proc;
    int t_proc;
    int *ex;  // array of excess
    int *dex; // array of delta excess
    int *hght;// array of heights (d)
    int *hght_p; // array of new heights (d')

    //std::vector<std::vector<int>> badj;
    // aflow[v] = {f(v,w_1), f(v,w_2), ...} for same order as adj
    std::vector<std::vector<int>> aflow;
    // bflow[w] = {f(v_1,w), f(v_2,w), ...} for same order as badj
    // "backward flows" are negative
    std::vector<std::vector<int>> bflow;

    // adj_d[v] = list of heights of adj verts, same order as adj
    std::vector<std::vector<int>> adj_d;
    // badj_d[v] = list of heights of back-adj verts, same order as badj
    std::vector<std::vector<int>> badj_d;

    int n_act;
    std::queue<int> active;
    std::queue<int> active_p;
};

struct comm_data {
    // for send requests: one for every edge
    MPI_Request **out_req;
    MPI_Request **in_req;
    std::vector<std::vector<int>> out_bi;
    std::vector<std::vector<int>> in_bi;
    std::vector<std::vector<int>> dist_adj_bi;
    std::vector<std::vector<int>> dist_badj_bi;
    // OUT_FLAG[V][I] ==  nothing WHEN (V,NET->ADJ[3*I]) HAS NO PENDING OPERATIONS
    // == 0 WHEN PENDING A SENT QUERY, == 1 FOR REC QUERY
    // == 2 FOR SENT RESPONSE, == 3 FOR REC RESPONSE
    std::vector<std::vector<unsigned char>> out_flag;
    std::vector<std::vector<unsigned char>> in_flag;

    // for query receives & distance comm: one for every proc
    MPI_Request *dist_req;

    // for both
    int buff_size;
    int **buff;
    std::queue<int> avail;
    int max_deg;
    int *arr_of_inds;
};


#endif
