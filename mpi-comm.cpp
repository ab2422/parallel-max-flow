#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <vector>
#include <queue>
#include "mpi-comm.h"
#include "mpi-pr.h"
using namespace std;

#define w_in(w,npp) ( (rank*(npp) <= (w)) && ((w) < ((rank+1)*(npp))) )

void check_dist(resgraph *net, int v, comm_data *cd, int rank, int size){
    int min_dist = 2*net->n;
    bool update = 1;
    int w,dw,z;
    z = 2*net->npp/3;
    int dv = net->hght[v];
    if ( (net->ex[v]>0) && (v!=net->src) && (v!=net->sink) ){
        for (int i=0; i<net->odeg[v]; i++){
            if (net->adj[v][3*i+1] - net->aflow[v][i] >0) {
                dw = net->adj_d[v][i];
                update = update && (dv <= dw);
                min_dist = min(min_dist,dw);
            }
        }
        for (int i=0; i<net->ideg[v]; i++){
            if (0 - net->bflow[v][i] >0) {
                dw = net->badj_d[v][i];
                update = update && (dv <= dw);
                min_dist = min(min_dist,dw);
            }
        }
    } else {
        update = 0;
    }
    
    // maybe do update
    bool win;
    int loc_w;
    int j;
    int bi;
    if (update) {
        net->hght[v] = min_dist+1;

        for (int i=0; i<net->odeg[v]; i++){
            w = net->adj[v][3*i];
            loc_w = w - rank*net->std_npp;
            j = net->adj[v][3*i+2];
            win = w_in(w,net->std_npp);
            if (win) {
                net->badj_d[loc_w][j/2];
            } else {
                while (cd->avail.empty()){
                    check_comm(net,z,cd,rank,size);
                    z = (z+1)%net->npp;
                }
                bi = cd->avail.front();
                cd->avail.pop();
                cd->buff[bi][0] = 0; //fwd
                cd->buff[bi][1] = v+rank*net->std_npp;
                cd->buff[bi][2] = net->hght[v];
                cd->buff[bi][3] = w;
                cd->buff[bi][4] = j;
                printf("Sending [%d, %d, %d, %d, %d] from %d to %d\n", cd->buff[bi][0], cd->buff[bi][1], cd->buff[bi][2], cd->buff[bi][3], cd->buff[bi][4], rank, i);
                if (cd->out_dist_bi[v][i] != -1){
                    MPI_Wait(&(cd->out_dist_req[v][i]), MPI_STATUS_IGNORE);
                    cd->avail.push(cd->out_dist_bi[v][i]);
                    cd->out_dist_bi[v][i]=-1;
                }
                MPI_Isend(cd->buff[bi],5,MPI_INT,w/net->std_npp, DIST_UPDATE, MPI_COMM_WORLD, &(cd->out_dist_req[v][i]));
            }
        }

        for (int i=0; i<net->ideg[v]; i++){
            w = net->badj[v][2*i];
            loc_w = w - rank*net->std_npp;
            j = net->badj[v][2*i+1];
            win = w_in(w,net->std_npp);
            if (win) {
                net->adj_d[loc_w][j/3];
            } else {
                while (cd->avail.empty()){
                    check_comm(net,z,cd,rank,size);
                    z = (z+1)%net->npp;
                }
                bi = cd->avail.front();
                cd->avail.pop();
                cd->buff[bi][0] = 1; //bwd
                cd->buff[bi][1] = v+rank*net->std_npp;
                cd->buff[bi][2] = net->hght[v];
                cd->buff[bi][3] = w;
                cd->buff[bi][4] = j;
                if (cd->in_dist_bi[v][i] != -1){
                    MPI_Wait(&(cd->in_dist_req[v][i]), MPI_STATUS_IGNORE);
                    cd->avail.push(cd->in_dist_bi[v][i]);
                    cd->in_dist_bi[v][i]=-1;
                }
                MPI_Isend(cd->buff[bi],5,MPI_INT,w/net->std_npp, DIST_UPDATE, MPI_COMM_WORLD, &(cd->in_dist_req[v][i]));
            }
        }
    }
}

void handle_comm(resgraph *net, int v, int gl_w, int *flowvj, int *adj_dvj, MPI_Request *req,  int *bi, unsigned char *flagv, int buffi[], queue<int> *avail, comm_data *cd, int rank, int size){
    if ((*flagv) == NOTHING){
        // do nothing!
    } else if ( (((*flagv)/2)%2 == 0) && ( (*flagv)/8 == 0) ){
        // just finished query
        if ( (*flagv)/4 ==0) {
            //just finished receiving query 
            printf("Error! Shouldn't be waiting to receive query\n");
        } else{
            // just finished sending query, wait for response
            if ((*flagv)%2 == 0) {
                MPI_Irecv(buffi,5,MPI_INT,gl_w/net->std_npp,FWD_RESPONSE,MPI_COMM_WORLD, req);
                (*flagv) = 2; // 010 (receive response fwd)
                // use same buff to receive response
            } else {
                MPI_Irecv(buffi,5,MPI_INT,gl_w/net->std_npp,BWD_RESPONSE,MPI_COMM_WORLD, req);
                (*flagv) = 3; // 011 (receive response bwd)
            }
        }
    } else if ( (((*flagv)/2)%2 == 1) && ( (*flagv)/8 == 0)){
        // just finished response
        if (((*flagv)/4 == 0)&&(!(buffi[0]))){
            // done receiving response: rejected
            if ((net->ex[v] == 0) && (buffi[2]>0) && (!is_src_loc(net,v,rank,size)) && (!is_sink_loc(net, v,rank,size))){
                net->active.push(v);
            } 
            net->ex[v] += buffi[2];
            (*flowvj) -= buffi[2]; 
            (*adj_dvj) = buffi[3];
            check_dist(net,v,cd,rank,size);
        }
        // cleanup
        (*flagv) = NOTHING;
        (*avail).push(*bi);
        (*bi) = -1;
    } else if ( ((*flagv)/8 == 1) && ( (*flagv)/4 == 1) ){
        // finished sending distance update
        (*flagv) = NOTHING;
        (*avail).push(*bi);
        (*bi)=-1;
    } else {
        printf("Invalid flag! Flag was %d\n", (*flagv));
    }
}

int listen_helper(resgraph *net, comm_data *cd, int bi, vector<vector<int>> *flow, vector<vector<int>> *adj, vector<vector<unsigned char>> *flagarr, int sc, MPI_Request ***reqarr, MPI_Request **req_ptr, unsigned char** flagv_ptr, int rank, int size){
    //received query, send response
    int v = cd->buff[bi][0]; // global!
    int ch = cd->buff[bi][1];
    int dv = cd->buff[bi][2];
    int w = cd->buff[bi][3];
    int i = cd->buff[bi][4];
    int loc_w = w - rank*net->std_npp;

    if (dv == 1 + net->hght[loc_w]){
        // accept
        if ((net->ex[loc_w]==0) && (ch>0) && (w!=net->src) && (w!=net->sink)){
            net->active.push(loc_w);
        }
        net->ex[loc_w] += ch;
        (*flow)[loc_w][i/sc] -= ch; 
        cd->buff[bi][0] = 1;
    } else{
        // reject
        cd->buff[bi][0]=0;
    }

    cd->buff[bi][1] = w;
    cd->buff[bi][2] = ch;
    cd->buff[bi][3] = net->hght[loc_w];
    cd->buff[bi][4] = (*adj)[loc_w][i+sc-1]; // +2 or +1, as approp.
    *req_ptr = &((*reqarr)[loc_w][i/sc]);
    *flagv_ptr = &((*flagarr)[loc_w][i/sc]);
    // bi updated in listen

    return v;
}

void listen_finish(resgraph *net, comm_data *cd, int rank, int size){
    int flag=1;
    int src;
    MPI_Status stat;
    MPI_Iprobe(MPI_ANY_SOURCE, FINISH, MPI_COMM_WORLD, &flag, &stat);
    while (flag){
        src = stat.MPI_SOURCE;
        MPI_Recv(&(cd->proc_done[src]), 1, MPI_CXX_BOOL, src, FINISH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Iprobe(MPI_ANY_SOURCE, FINISH, MPI_COMM_WORLD, &flag, &stat);
    }
}

void listen_distance(resgraph *net, comm_data *cd, int rank, int size){
    int flag=1;
    int src;
    int i,bi;
    int gl_v,w,loc_w,dv;
    int z = (net->std_npp/2)%net->npp;
    MPI_Status stat;
    MPI_Iprobe(MPI_ANY_SOURCE, DIST_UPDATE, MPI_COMM_WORLD, &flag, &stat);
    while (flag) {
        while (cd->avail.empty()) {
            check_comm(net,z,cd,rank,size);
            z = (z+1)%net->npp;
        }
        src = stat.MPI_SOURCE;
        bi = cd->avail.front();
        cd->avail.pop();
        MPI_Recv(cd->buff[bi], 5, MPI_INT, src, DIST_UPDATE, MPI_COMM_WORLD, &stat);
        gl_v = cd->buff[bi][1];
        dv= cd->buff[bi][2];
        w = cd->buff[bi][3];
        i = cd->buff[bi][4];
        loc_w = w - rank*net->std_npp;
        if (cd->buff[bi][0] == 0){
            // forward edge for sender
            net->badj_d[loc_w][i/2] = dv;
        } else {
            // bwd edge for sender
            net->adj_d[loc_w][i/3]= dv;
        }
        cd->avail.push(bi);
        printf("done dist upd: v=%d, w=%d, new dv=%d\n", gl_v, w,dv);

        MPI_Iprobe(MPI_ANY_SOURCE, DIST_UPDATE, MPI_COMM_WORLD, &flag, &stat);
    }
}

void listen(resgraph *net, comm_data *cd, int rank, int size){
    int flag = 1;
    unsigned char *flagv; // pts to appropriate flag in cd
    int bi;
    int src;
    int tag;
    int v,loc_v,dv,sc,ch,i;
    int loc_w,j;
    int z=0;
    int orig_proc;
    vector<vector<int>>* flow, adj;
    MPI_Status stat;
    MPI_Request *req_pt;
    
    // look for messages COMING FROM v s.t. v-> here is fwd
    MPI_Iprobe(MPI_ANY_SOURCE, FWD_QUERY, MPI_COMM_WORLD, &flag, &stat); 
    while (flag){
        while (cd->avail.empty()){
            check_comm(net, z, cd ,rank,size);
            z = (z+1)%net->npp;
        }
        bi = cd->avail.front();
        cd->avail.pop();
        src = stat.MPI_SOURCE;
        
        MPI_Recv(cd->buff[bi], 5, MPI_INT, src, FWD_QUERY, MPI_COMM_WORLD, &stat);
        loc_w = cd->buff[bi][3]-rank*net->std_npp;
        j = cd->buff[bi][4];
        orig_proc = listen_helper(net,cd,bi,&(net->bflow), &(net->badj), &(cd->in_flag), 2, &(cd->in_req), &req_pt, &flagv, rank, size);
        MPI_Isend(cd->buff[bi], 5, MPI_INT, src , FWD_RESPONSE, MPI_COMM_WORLD, req_pt);
        cd->in_bi[loc_w][j/2] = bi;
        (*flagv) = 6; // 110 (send response fwd)
        MPI_Iprobe(MPI_ANY_SOURCE, FWD_QUERY, MPI_COMM_WORLD, &flag, &stat);

    }
    
    // look for messages COMING from v s.t. v->here is bwd
    MPI_Iprobe(MPI_ANY_SOURCE, BWD_QUERY, MPI_COMM_WORLD, &flag, &stat);
    while (flag){
        while (cd->avail.empty()){
            check_comm(net, z, cd ,rank,size);
            z++;
        }
        bi = cd->avail.front();
        cd->avail.pop();
        src = stat.MPI_SOURCE;
    
        MPI_Recv(cd->buff[bi], 5, MPI_INT, src, BWD_QUERY, MPI_COMM_WORLD, &stat);
        loc_w = cd->buff[bi][3]-rank*net->std_npp;
        j = cd->buff[bi][4];
        listen_helper(net,cd,bi,&(net->aflow), &(net->adj), &(cd->out_flag), 2, &(cd->out_req), &req_pt, &flagv, rank, size);
        MPI_Isend(cd->buff[bi], 5, MPI_INT, src, BWD_RESPONSE, MPI_COMM_WORLD, req_pt);
        (*flagv) = 7; // 111 (send response bwd) 
        cd->out_bi[loc_w][j/3] = bi;
        MPI_Iprobe(MPI_ANY_SOURCE, BWD_QUERY, MPI_COMM_WORLD, &flag, &stat); 
    }
}



void check_comm_helper(resgraph *net, int v, int incount, std::vector<int> *flowv, vector<int> *adjv, std::vector<int> *adj_dv, MPI_Request *reqv, std::vector<int> *arr_biv, std::vector<unsigned char> *arr_flagv, int **buff, std::queue<int> *avail, int arr_of_inds[], comm_data *cd, int rank, int size){
    int outcount = 0;
    int j, bi, w;
    MPI_Testsome(incount,reqv, &outcount, arr_of_inds,MPI_STATUSES_IGNORE);
    // remember: aflows are positive & bflows are neg, hence the signs working out!
    if (outcount != MPI_UNDEFINED){
        for (int i=0; i<outcount; i++){
            j = arr_of_inds[i];
            bi = (*arr_biv)[j]; //what buffer are we using?
            w = (*adjv)[j];
            handle_comm(net, v, w, &((*flowv)[j]), &((*adj_dv)[j]), &(reqv[j]), &bi, &((*arr_flagv)[j]), buff[bi], avail, cd, rank,size);
        }
    }


}

void check_comm(resgraph *net, int v, comm_data *cd, int rank, int size){
    // for requests we point out to
    check_comm_helper(net, v, net->odeg[v], &(net->aflow[v]), &(net->adj[v]), &(net->adj_d[v]),cd->out_req[v], 
                      &(cd->out_bi[v]), &(cd->out_flag[v]), cd->buff, &(cd->avail), 
                      cd->arr_of_inds,cd, rank,size);

    // for requests pointing in to us
    check_comm_helper(net, v, net->ideg[v], &(net->bflow[v]), &(net->badj[v]), &(net->badj_d[v]), cd->in_req[v], 
                      &(cd->in_bi[v]), &(cd->in_flag[v]), cd->buff, &(cd->avail), 
                      cd->arr_of_inds, cd, rank,size);
        
}
