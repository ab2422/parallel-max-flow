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

/*
 * Handles "planned" communication from local v to global gl_w, 
 * which corresponds to edge in direction dir (1=bwd, 0=fwd)
 * s.t. j is w's location in adj[dir][v]/2 (i.e. the index used
 *      for indexing into the edge, flow, etc, arrays).
 */
void handle_comm(resgraph *net, int v, int gl_w, int dir, int j, comm_data *cd, int rank, int size); 

/*
 * Checks status of all planned communications in direction dir involving v.
 * Incount is the number of possible such communications.
 * arr_of_inds is purely local, but since (possibly) large is stored in cd
 */
void check_comm_helper(resgraph *net, int v, int dir, int incount, int arr_of_inds[], comm_data *cd, int rank, int size);

/*
 * Checks status of all pending "planned" communications, 
 * and if completed, processes them.
 */ 
void check_comm(resgraph *net, int v, comm_data *cd,int rank, int size);

/*
 * Helps the listener function. Dir is the direction that the v w edge goes, 
 * from w's (the receiver's) perspective
 */
void listen_helper(resgraph *net, comm_data *cd, int bi, int dir, MPI_Request **req, unsigned char **flagv_pt, int rank, int size);

/*
 * Listens to see if any processors have finished.
 */
void handle_finish(resgraph *net, comm_data *cd, int rank, int size);

/*
 * Listens for incoming distance updates & processes them.
 */
void listen_distance(resgraph *net, comm_data *cd, int rank, int size);


/*
 * Listens for incoming queries, processes them, and posts response
 */
void listen(resgraph *net, comm_data *cd, int rank, int size);

#endif
