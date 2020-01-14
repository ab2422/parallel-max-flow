#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <vector>
#include <queue>
#include "mpi-comm.h"

void handle_outgoing_comm(resgraph *net, int v, int w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi[], queue<int> *avail, int rank, int size){
    if ((*flagv)/2 == 0){
        if ((*flagv)%2 == 0) {
            // done sending query,wait for response
            MPI_Irecv(buffi,4,MPI_INT,w/net->std_npp,RESPONSE,MPI_COMM_WORLD, req);
            (*flagv) = 3; // 11
        } else {
            //received query, send response
            if (buffi[2] == 1 + net->hght[v]){
                // accept
                net->ex[v] += buffi[1];
                (*flowvj) -= buffi[1]; 
                buffi[0] = 1;
            } else{
                // reject
                buffi[0]=0;
            }
            buffi[2] = buffi[1];
            buffi[1]= v+rank*net->std_npp;
            buffi[3] = net->hght[v];
            MPI_Isend(buffi,4,MPI_INT,w/net->std_npp,RESPONSE,MPI_COMM_WORLD, req);
            (*flagv) = 2; // 10
        }
    } else if ((*flagv)/2 == 1){
        if (((*flagv)%2 == 0)&&(!(buffi[0]))){
            // done receiving response: rejected
            net->ex[v] += buffi[2];
            (*flowvj) -= buffi[2]; 
            (*adj_dvj) = buffi[3];
            // TODO: when to update dv?????
        }
        // cleanup
        (*flagv) = NOTHING;
        (*avail).push(*bi);
        (*bi) = -1;
    } else {
        printf("Invalid flag! Flag was %d\n", (*flagv));
    }
}


void handle_incoming_comm(resgraph *net, int v, MPI_Request ireq[], int buff[][4], queue<int> *avail){
    int flag = 1;
    int bi;
    MPI_Status stat;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &stat); 
    while (flag){
        bi = avail.front();
        avail.pop();
        if (stat.MPI_TAG == QUERY){
            MPI_IRecv(buff[bi], 3, MPI_INT, stat.MPI_SOURCE, QUERY, MPI_COMM_WORLD, 
        } else if (stat.MPI_TAG == RESPONSE){

        } else { printf("Invalid tag! Tag was %d\n" stat.MPI_TAG); }
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &stat); 
    }
}


void check_comm(resgraph *net, int v, vector<int> *flowv, vector<int> *adj_dv, MPI_Request *reqv, vector<int> *arr_biv, vector<unsigned char> *arr_flagv, int buff[][4], queue<int> *avail, int arr_of_inds[], int rank, int size){
        int outcount,w,j,bi;
        MPI_Testsome(net->odeg[v],reqv, &outcount, arr_of_inds,MPI_STATUSES_IGNORE);
        // remember: aflows are positive & bflows are neg, hence the signs working out!
        if (outcount != MPI_UNDEFINED){
            for (int i=0; i<outcount; i++){
                j = arr_of_inds[i];
                bi = (*arr_biv)[j]; //what buffer are we using?
                w = net->adj[v][j];
                handle_comm(net, v, w, &((*flowv)[j]), &((*adj_dv)[j]), &(reqv[j]), &bi, &((*arr_flagv)[j]), buff[bi], avail, rank,size);
            }
        }
}
