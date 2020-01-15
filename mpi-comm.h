#ifndef MPI_COMMUNICATION_CODE_H
#define MPI_COMMUNICATION_CODE_H

#include <vector>
#include <queue>
#include <iostream> 
#include "mpi-data.h"


// Queries, from v to w: (v, ch, d(v),w,i), where ch = change in flow
//          and i is the index of v in w's adj or badj list (as appropriate)
// Responses, frow v to w: (acc/rej, v, ch, d(v),i), where ch = change in flow
//          and i is index of v in w's adj or badj list (as appropriate)
// for both of these, v is measured in GLOBAL index
// Distance updates, frow v to w: (bwd/fwd, v,dv,w,i), where fwd from v's persp,
//          i is ind of v in w's adj or badj list (as appropriate)
 
/*
 * Checks whether node v needs to be relabelled. If so, relabels & updates nbhrs
 */
void check_dist(resgraph *net, int v, comm_data *cd, int rank, int size);

void handle_comm(resgraph *net, int v, int w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi[], std::queue<int> *avail, comm_data *cd, int rank, int size);

void check_comm_helper(resgraph *net, int v, std::vector<int> *flowv, std::vector<int> *adj_dv, MPI_Request *reqv, std::vector<int> *arr_biv, std::vector<unsigned char> *arr_flagv, int **buff, std::queue<int> *avail, int arr_of_inds[], comm_data *cd, int rank, int size);

/*
 * Checks status of all pending communications, and if completed, processes them
 */ 
void check_comm(resgraph *net, int v, comm_data *cd,int rank, int size);

void listen_helper(resgraph *net, comm_data *cd, int bi, std::vector<std::vector<int>> *flow, std::vector<std::vector<int>> *adj, std::vector<std::vector<unsigned char>> *flagarr, int sc, MPI_Request **req, unsigned char **flagv_pt, int rank, int size);

/*
 * Listens for incoming distance updates & processes them.
 */
void listen_distance(resgraph *net, comm_data *cd, int rank, int size);


/*
 * Listens for incoming queries, processes them, and posts response
 */
void listen(resgraph *net, comm_data *cd, int rank, int size);

#endif
