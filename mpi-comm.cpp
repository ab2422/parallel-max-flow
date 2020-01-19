#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <vector>
#include <queue>
#include "mpi-comm.h"
#include "mpi-pr.h"
//#define NDEBUG
#include <assert.h>
using namespace std;

#define w_in(w,npp) ( (rank*(npp) <= (w)) && ((w) < ((rank+1)*(npp))) )

void check_dist(resgraph *net, int v, comm_data *cd, int rank, int size){
    int min_dist = 2*net->n;
    bool update = 1;
    int w,dw,z;
    z = 2*net->npp/3;
    int dv = net->hght[v];
    if ( (net->ex[v]>0) && (!is_src_loc(net,v,rank,size)) && (!is_sink_loc(net,v,rank,size)) ){
        for (int i=0; i<net->odeg[v]; i++){
            if (net->cap[v][i] - net->flow[0][v][i] >0) {
                dw = net->adj_d[0][v][i];
                update = update && (dv <= dw);
                min_dist = min(min_dist,dw);
            }
        }
        for (int i=0; i<net->ideg[v]; i++){
            if (0 - net->flow[1][v][i] >0) {
                dw = net->adj_d[1][v][i];
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
            w = net->adj[0][v][2*i];
            loc_w = w - rank*net->std_npp;
            j = net->adj[0][v][2*i+1];
            win = w_in(w,net->std_npp);
            if (win) {
                net->adj_d[1][loc_w][j/2];
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
                if (cd->dist_bi[0][v][i] != -1){
                    MPI_Wait(&(cd->dist_req[0][v][i]), MPI_STATUS_IGNORE);
                    cd->avail.push(cd->dist_bi[0][v][i]);
                    cd->dist_bi[0][v][i]=-1;
                }
                MPI_Isend(cd->buff[bi],5,MPI_INT,w/net->std_npp, DIST_UPDATE, MPI_COMM_WORLD, &(cd->dist_req[0][v][i]));
            }
        }

        for (int i=0; i<net->ideg[v]; i++){
            w = net->adj[1][v][2*i];
            loc_w = w - rank*net->std_npp;
            j = net->adj[1][v][2*i+1];
            win = w_in(w,net->std_npp);
            if (win) {
                net->adj_d[0][loc_w][j/2];
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
                if (cd->dist_bi[1][v][i] != -1){
                    MPI_Wait(&(cd->dist_req[1][v][i]), MPI_STATUS_IGNORE);
                    cd->avail.push(cd->dist_bi[1][v][i]);
                    cd->dist_bi[1][v][i]=-1;
                }
                MPI_Isend(cd->buff[bi],5,MPI_INT,w/net->std_npp, DIST_UPDATE, MPI_COMM_WORLD, &(cd->dist_req[1][v][i]));
            }
        }
    }
}

void handle_comm(resgraph *net, int v, int gl_w, int dir, int j, comm_data *cd, int rank, int size){ 
    int bi = cd->edge_bi[dir][v][j];
    if (cd->edge_flag[dir][v][j] == NOTHING){
        // do nothing!
    } else if ( (((cd->edge_flag[dir][v][j])/2)%2 == 0) && ( (cd->edge_flag[dir][v][j])/8 == 0) ){
        // just finished query
        if ( (cd->edge_flag[dir][v][j])/4 ==0) {
            //just finished receiving query 
            printf("Error! Shouldn't be waiting to receive query\n");
        } else{
            // just finished sending query, wait for response
            if ((cd->edge_flag[dir][v][j])%2 == 0) {
                assert (dir==0);
                MPI_Irecv(cd->buff[bi],5,MPI_INT,gl_w/net->std_npp,FWD_RESPONSE,MPI_COMM_WORLD, &(cd->edge_req[dir][v][j]));
                (cd->edge_flag[dir][v][j]) = 2; // 010 (receive response fwd)
                // use same buff to receive response
            } else {
                assert (dir==1);
                MPI_Irecv(cd->buff[bi],5,MPI_INT,gl_w/net->std_npp,BWD_RESPONSE,MPI_COMM_WORLD, &(cd->edge_req[dir][v][j]));
                (cd->edge_flag[dir][v][j]) = 3; // 011 (receive response bwd)
            }
        }
    } else if ( (((cd->edge_flag[dir][v][j])/2)%2 == 1) && ( (cd->edge_flag[dir][v][j])/8 == 0)){
        // just finished response
        if (((cd->edge_flag[dir][v][j])/4 == 0)&&(cd->buff[bi][0]==0)){
            // done receiving response: rejected
            if ((net->ex[v] == 0) && (cd->buff[bi][2]>0) && (!is_src_loc(net,v,rank,size)) && (!is_sink_loc(net, v,rank,size))){
                net->active.push(v);
            } 
            net->ex[v] += cd->buff[bi][2];
            net->flow[dir][v][j] -= cd->buff[bi][2]; 
            net->adj_d[dir][v][j] = cd->buff[bi][3];
            check_dist(net,v,cd,rank,size);
        }
        // cleanup
        (cd->edge_flag[dir][v][j]) = NOTHING;
        cd->avail.push(bi);
        cd->edge_bi[dir][v][j] = -1;
    } else if ( ((cd->edge_flag[dir][v][j])/8 == 1) && ( (cd->edge_flag[dir][v][j])/4 == 1) ){
        // finished sending distance update
        (cd->edge_flag[dir][v][j]) = NOTHING;
        cd->avail.push(bi);
        cd->edge_bi[dir][v][j]=-1;
    } else {
        printf("Invalid flag! Flag was %d\n", (cd->edge_flag[dir][v][j]));
    }
}

void listen_helper(resgraph *net, comm_data *cd, int bi, int dir, MPI_Request **req_ptr, unsigned char**flagv_ptr, int rank, int size){ 
    //received query, send response
    int v = cd->buff[bi][0]; // global!
    int ch = cd->buff[bi][1];
    int dv = cd->buff[bi][2];
    int gl_w = cd->buff[bi][3];
    int i = cd->buff[bi][4];
    int loc_w = gl_w - rank*net->std_npp;

    if (dv == 1 + net->hght[loc_w]){
        // accept
        if ((net->ex[loc_w]==0) && (ch>0) && (gl_w!=net->src) && (gl_w!=net->sink)){
            net->active.push(loc_w);
        }
        net->ex[loc_w] += ch;
        net->flow[dir][loc_w][i/2] -= ch; 
        cd->buff[bi][0] = 1;
    } else{
        // reject
        cd->buff[bi][0]=0;
    }

    cd->buff[bi][1] = gl_w;
    cd->buff[bi][2] = ch;
    cd->buff[bi][3] = net->hght[loc_w];
    cd->buff[bi][4] = net->adj[dir][loc_w][i+1]; 
    *req_ptr = &(cd->edge_req[dir][loc_w][i/2]);
    *flagv_ptr = &(cd->edge_flag[dir][loc_w][i/2]);
    cd->edge_bi[dir][loc_w][i/2] = bi;
}

/*
 * Handle the ring for completion, both initiating & passing.
 */
void handle_finish(resgraph *net, comm_data *cd, int rank, int size){
    bool f_done;
    int f_tst;
    int f_prev = (rank-1+size)%size;
    int f_nxt = (rank+1)%size;
    int f_buf;
    if (cd->proc_done[rank]){
        if (!net->active.empty()){
            // we're undone
            cd->proc_done[rank]=0;
        }
        
    } else {
        cd->proc_done[rank]=1;
        cd->num_times_done++;
    }
    if ((rank==FIN_PROC)&&(cd->proc_done[rank])){
        // handle initition of ring, & the handoffs
        if (cd->ring_stat==0){

            cd->prev_total=cd->num_times_done;
            MPI_Isend(&cd->prev_total, 1, MPI_INT, f_nxt, F_RING, MPI_COMM_WORLD, &(cd->ring_req));
            cd->ring_flag=1;
            cd->ring_stat=1;

        } else if (cd->ring_stat == 1){

            MPI_Iprobe(f_prev, F_RING, MPI_COMM_WORLD, &f_tst, MPI_STATUS_IGNORE);
            if (f_tst){
                // finished 1st ring, start 2nd
                MPI_Wait(&(cd->ring_req), MPI_STATUS_IGNORE);
                cd->ring_flag=0;
                MPI_Recv(&(cd->prev_total), 1, MPI_INT, f_prev, F_RING, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                cd->nxt_total=cd->num_times_done;
                MPI_Isend(&(cd->nxt_total), 1, MPI_INT, f_nxt, F_RING, MPI_COMM_WORLD, &(cd->ring_req));
                cd->ring_stat=2;
                cd->ring_flag=1;
            }

        } else if (cd->ring_stat ==2){

            MPI_Iprobe(f_prev, F_RING, MPI_COMM_WORLD, &f_tst, MPI_STATUS_IGNORE);
            if (f_tst){
                MPI_Wait(&(cd->ring_req), MPI_STATUS_IGNORE);
                cd->ring_flag=0;
                MPI_Recv(&(cd->nxt_total), 1, MPI_INT, f_prev, F_RING, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (cd->prev_total==cd->nxt_total){
                    // we're done! no one has restarted & finished again since our first ring.
                    cd->all_done=1;
                    for (int j=0; j<size; j++){
                        if (j!=rank){
                            MPI_Isend(&(cd->all_done), 1, MPI_CXX_BOOL, j, FINISH, MPI_COMM_WORLD, &(cd->fin_req[j]));
                        }
                    }
                    MPI_Waitall(size, cd->fin_req, MPI_STATUSES_IGNORE);
                    
                } else {
                    // looks like some stuff has happened, let's restart
                    cd->ring_stat=0;
                }
            }

        } 
    } else {
        MPI_Iprobe(f_prev, F_RING, MPI_COMM_WORLD, &f_tst, MPI_STATUS_IGNORE);
        if (f_tst){
            MPI_Recv(&f_buf, 1, MPI_INT, f_prev, F_RING, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            f_buf += cd->num_times_done;
            MPI_Wait(&(cd->ring_req), MPI_STATUS_IGNORE);
            MPI_Isend(&f_buf, 1, MPI_INT, f_nxt, F_RING, MPI_COMM_WORLD, &(cd->ring_req));
            cd->ring_flag=1;
        }
        MPI_Iprobe(FIN_PROC, FINISH, MPI_COMM_WORLD, &f_tst, MPI_STATUS_IGNORE);
        if (f_tst){
            MPI_Recv(&(cd->all_done), 1, MPI_CXX_BOOL, FIN_PROC, FINISH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

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
            net->adj_d[1][loc_w][i/2] = dv;
        } else {
            // bwd edge for sender
            net->adj_d[0][loc_w][i/2]= dv;
        }
        cd->avail.push(bi);

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
        // since for us, the dir is bwd
        listen_helper(net,cd,bi,1, &req_pt, &flagv, rank, size);
        MPI_Isend(cd->buff[bi], 5, MPI_INT, src , FWD_RESPONSE, MPI_COMM_WORLD, req_pt);
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
        listen_helper(net,cd,bi, 0, &req_pt, &flagv, rank, size);
        MPI_Isend(cd->buff[bi], 5, MPI_INT, src, BWD_RESPONSE, MPI_COMM_WORLD, req_pt);
        (*flagv) = 7; // 111 (send response bwd) 
        MPI_Iprobe(MPI_ANY_SOURCE, BWD_QUERY, MPI_COMM_WORLD, &flag, &stat); 
    }
}



void check_comm_helper(resgraph *net, int v, int dir, int incount, comm_data *cd, int rank, int size){
    int outcount = 0;
    int j, bi, w;
    // TODO find error here
    MPI_Testsome(incount,cd->edge_req[dir][v], &outcount, cd->arr_of_inds,MPI_STATUSES_IGNORE);
    // remember: aflows are positive & bflows are neg, hence the signs working out!
    if (outcount != MPI_UNDEFINED){
        for (int i=0; i<outcount; i++){
            j = cd->arr_of_inds[i];
            bi = cd->edge_bi[dir][v][j]; //what buffer are we using?
            w = net->adj[dir][v][2*j];
            handle_comm(net, v, w, dir, j, cd, rank,size);
        }
    }


}

void check_comm(resgraph *net, int v, comm_data *cd, int rank, int size){
    // for requests we point out to
    check_comm_helper(net, v, 0, net->odeg[v], cd, rank,size);

    // for requests pointing in to us
    check_comm_helper(net, v, 1, net->ideg[v], cd, rank,size);
        
}
