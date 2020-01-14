#ifndef MPI_COMMUNICATION_CODE_H
#define MPI_COMMUNICATION_CODE_H

#include <vector>
#include <queue>
#include <iostream> 


void handle_comm(resgraph *net, int v, int w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi, std::queue<int> *avail, int rank, int size);

void check_comm(resgraph *net, int v, std::vector<int> *flowv, std::vector<int> *adj_dv, MPI_Request *reqv, std::vector<int> *arr_biv, std::vector<unsigned char> *arr_flagv, int buff[][4], std::queue<int> *avail, int arr_of_inds[], MPI_Status arr_of_stat[], int rank, int size);

#endif
