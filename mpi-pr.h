#ifndef MPI_PUSH_RELABEL_CODE_H
#define MPI_PUSH_RELABEL_CODE_H

#include <vector>
#include <queue>
#include <iostream> 
#include <mpi.h>
#include "mpi-data.h"

// in binary order, our tags
#define FWD_QUERY 100
#define BWD_QUERY 105
#define FWD_RESPONSE 110
#define BWD_RESPONSE 115
#define DIST_UPDATE 300
#define FINISH 250
// REQUIRE: NOTHING >= 16
#define NOTHING 32

// For binary 4-digit d_2d_1d_0: d_3 = distance / affirmation ,
//                               d_2 = send/receive,  
//                               d_1 = response/query,
//                               d_0 = bwd/fwd (where 1/0).
#define AFF_RCV_QRY_FWD 0




bool is_src_loc(resgraph *net, int loc_v, int rank, int size);

bool is_sink_loc(resgraph *net, int loc_v, int rank, int size);

void print_flow(resgraph *net, int rank);

network parse(std::string filename, int rank, int size);

resgraph setup(network *net, int rank, int size);

comm_data setup_cd(resgraph *net, int rank, int size);

//void pulse(resgraph *net_ptr);

//void handle_comm(resgraph *net, int v, int w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi, std::queue<int> *avail, int rank, int size);

//void check_comm(resgraph *net, int v, std::vector<int> *flowv, std::vector<int> *adj_dv, MPI_Request *reqv, std::vector<int> *arr_biv, std::vector<unsigned char> *arr_flagv, int buff[][4], std::queue<int> *avail, int arr_of_inds[], MPI_Status arr_of_stat[], int rank, int size);

void async_pr(resgraph *net, comm_data *cd, int rank, int size);


void output(resgraph net, std::string filename, double time);

bool check(std::string correct, std::string test);

void cleanup(resgraph *net, network *orig);

#endif
