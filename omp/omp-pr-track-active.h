#ifndef OMP_PUSH_RELABEL_CODE_H
#define OMP_PUSH_RELABEL_CODE_H

#include <vector>
#include <iostream> 

struct network {
    int src;  // index of src
    int sink; // index of sink
    int n;    // num verts
    int m;    // num edges
    int *odeg; // array of out degs
    int *ideg; // array of in degs
    // adj[v] = {w_1, c(v,w_1), j_1, w_2, c(v,w_2), j_2 ...} 
    //          for w_i adj to v, where j_i is ind of v in w's list
    std::vector<std::vector<int>> adj; 
    // badj[w] = {v_1, j_1, v_2, j_2, ...} for w adj to v_i, j_i ind of v
    std::vector<std::vector<int>> badj;
};

struct resgraph {
    int src;  // index of src
    int sink; // index of sink
    int n;    // num verts
    int m;    // num edges
    int *odeg; // array of out degs
    int *ideg; // array of in degs
    int *ex;  // array of excess
    int *dex; // array of delta excess
    int *hght;// array of heights (d)
    int *hght_p; // array of new heights (d')

    // adj[v] = {w_1, c(v,w_1), j_1, w_2, c(v,w_2), j_2, ...} 
    //            for w_i adj to v
    std::vector<std::vector<int>> adj; 
    // badj[w] = {v_1, j_1, v_2, j_2, ...}
    //              for v_i adj to w
    std::vector<std::vector<int>> badj;
    // aflow[v] = {f(v,w_1), f(v,w_2), ...} for same order as adj
    std::vector<std::vector<int>> aflow;
    // bflow[w] = {f(v_1,w), f(v_2,w), ...} for same order as badj
    std::vector<std::vector<int>> bflow;

    int n_act;
};

network parse(std::string filename);

void setup(network *net_ptr, resgraph *res_ptr);

void pulse(resgraph *net_ptr);

void output(resgraph net, std::string filename, double time);

bool check(std::string correct, std::string test);

void cleanup(resgraph net, network orig);

#endif
