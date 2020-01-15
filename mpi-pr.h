#ifndef MPI_PUSH_RELABEL_CODE_H
#define MPI_PUSH_RELABEL_CODE_H

#include <vector>
#include <queue>
#include <iostream> 
#include <mpi.h>

// in binary order, our tags
#define FWD_QUERY 123
#define BWD_QUERY 125
#define FWD_RESPONSE 124
#define BWD_RESPONSE 126
#define DIST_UPDATE 300
// REQUIRE: NOTHING >= 16
#define NOTHING 32

// For binary 4-digit d_2d_1d_0: d_3 = distance / affirmation ,
//                               d_2 = send/receive,  
//                               d_1 = response/query,
//                               d_0 = bwd/fwd (where 1/0).

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

network parse(std::string filename, int rank, int size);

resgraph setup(network *net, int rank, int size);

void pulse(resgraph *net_ptr);

void handle_comm(resgraph *net, int v, int w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi, std::queue<int> *avail, int rank, int size);

void check_comm(resgraph *net, int v, std::vector<int> *flowv, std::vector<int> *adj_dv, MPI_Request *reqv, std::vector<int> *arr_biv, std::vector<unsigned char> *arr_flagv, int buff[][4], std::queue<int> *avail, int arr_of_inds[], MPI_Status arr_of_stat[], int rank, int size);

void async_pr(resgraph *net,int rank, int size);

void output(resgraph net, std::string filename, double time);

bool check(std::string correct, std::string test);

void cleanup(resgraph *net, network *orig);

#endif
