#include "mpi-pr.h"
#include "mpi-data.h"
#include "mpi-comm.h"
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <queue>
//#define NDEBUG
#include <assert.h>
using namespace std;

#define max(a,b) ( ((a)>(b)) ? (a) : (b) )
#define min(a,b) ( ((a)<(b)) ? (a) : (b) )
#define w_in(w,npp) ( (rank*(npp) <= (w)) && ((w) < ((rank+1)*(npp))) )

bool is_src_loc(resgraph *net, int loc_v, int rank, int size){
    return  (net->s_proc == rank) && ( (loc_v +rank*net->std_npp) == net->src);
}

bool is_sink_loc(resgraph *net, int loc_v, int rank, int size){
    return (net->sink == (loc_v + rank*net->std_npp));
}

static void wait_for_debugger(){
    volatile int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i){}
}

/*
Parses a file filename provided in DIMACS netflow format. 
* Stores #verts in *n, #edges in *m.
* nbhr[v] = {w, c(v,w), w', c(v,w'), ... } for nbhrs w
* net.ideg[v] = in degree of v
* net.odeg[v] = out degree of v
* Stores src ind in *src, sink ind in *sink
* Requires: only 1 p max, only 1 src, only 1 sink
* Requires: barring comments, order is p max link, src/sink lines, arc lines
* To ensure zero indexing, shifts indices in file all down by one
*/
network parse(string filename, int rank, int size){
    printf("starting parse\n");
    network net;
    string line;
    ifstream file;
    file.open(filename, ios::in);
    int i=0;
    int npp=0;
    bool win=0;
    int v,w,cap,f,jv,jw;
    int rbuf = 0; // v,w,c(v,w), pos
    int sbuf = 0; // v,w,c(v,w), pos
    MPI_Request sreq = MPI_REQUEST_NULL;
    MPI_Request rreq = MPI_REQUEST_NULL;
    MPI_Status sstat, rstat;
    bool ready = 1;
    if (file.is_open()) {
        while (getline(file, line)){
            if (line[0]=='p'){
                if (!((line[1]==' ')&&(line[2]=='m')&&(line[3]=='a')&&(line[4]=='x')&&(line[5]==' '))) {
                    //return 0;
                } else {
                    i=6;
                    net.n = atoi( &(line[6]));
                    npp = (net.n + size-1  )/size;
                    net.std_npp = npp;
                    net.npp = min(npp, max(0,net.n-rank*npp));
                    while (isspace(line[i])){
                        i++;
                    }
                    while ((isalnum(line[i]))){
                        i++;
                    }
                    net.m = atoi(&(line[i]));
                    net.odeg = (int*) malloc((net.npp)*sizeof(int));
                    net.ideg = (int*) malloc((net.npp)*sizeof(int));
                    net.adj = vector<vector<int>>(net.npp);
                    net.badj = vector<vector<int>>(net.npp);
                    for (int j=0; j<net.npp; j++){
                        net.adj[j].reserve(4*net.m/(size*net.n)+1);
                        net.badj[j].reserve(4*net.m/(size*net.n)+1);
                        net.odeg[j]=0;
                        net.ideg[j]=0;
                    }
                }

            } else if (line[0]=='n'){

                //dealing with nodes
                i=2;
                while (isspace(line[i])){
                    i++;
                }
                while ((isalnum(line[i]))){
                    i++;
                }
                if (line[i+1]=='s'){
                    net.src = atoi(&(line[2])) - 1 ;
                } else if (line[i+1]=='t'){
                    net.sink = atoi(&(line[2])) - 1 ;
                } else {
                    //return 0;
                }

            } else if (line[0]=='a'){

                // dealing with arcs
                i = 2;
                v = atoi(&(line[i]))-1;
                while (isspace(line[i])){
                    i++;
                }
                while ( isalnum(line[i]) ) {
                    i++;
                }
                w = atoi(&(line[i]))-1;
                i++; // now i points to (ws prior to) start of next number
                while (isspace(line[i])){
                    i++;
                }
                while ( isalnum(line[i]) ) {
                    i++;
                }
                cap = atoi(&(line[i]));
                win = (( (rank*npp<=w) && (w < (rank+1)*npp) ));
                // if vin: 
                if (((rank*npp)<=v) && (v < (rank+1)*npp)){
                    net.odeg[v-rank*npp]++;
                    jv = net.adj[v-rank*npp].size();
                    if (win) {
                        net.ideg[w-rank*npp]++;
                        jw = net.badj[w-rank*npp].size();
                        net.adj[v-rank*npp].push_back(w); // adj vert
                        net.adj[v-rank*npp].push_back(cap); // capacity
                        net.adj[v-rank*npp].push_back(jw); //jw = v pos in w list
        
                        net.badj[w-rank*npp].push_back(v); // bacwards edge
                        //net.badj[w].push_back(cap); //cap for back
                        net.badj[w-rank*npp].push_back(jv); //jv = w pos in v list
                    }else{
                        MPI_Irecv(&rbuf,1,MPI_INT, w/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                        // wait for prev send to finish, so know sbuf avail
                        MPI_Wait(&sreq, &sstat);
                        sbuf = jv;
                        MPI_Isend(&sbuf, 1, MPI_INT, w/npp, v*net.n + w, MPI_COMM_WORLD,&sreq);
                        net.adj[v-rank*npp].push_back(w);
                        net.adj[v-rank*npp].push_back(cap);
                        // wait for curr rec to finish, so can use data
                        MPI_Wait(&rreq, &rstat);
                        net.adj[v-rank*npp].push_back(rbuf); 
                    }
                } else if (win) {
                    MPI_Irecv(&rbuf,1,MPI_INT,v/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                    net.ideg[w-rank*npp]++;
                    jw = net.badj[w-rank*npp].size();
                    net.badj[w-rank*npp].push_back(v);
                    MPI_Wait(&sreq, &sstat);
                    sbuf = jw;
                    MPI_Isend(&sbuf,1,MPI_INT, v/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                    MPI_Wait(&rreq, &rstat);
                    jw = rbuf;
                    net.badj[w-rank*npp].push_back(rbuf);

                }
            } else if (line[0]!='c'){
                cout << "Invalid file. Bad line is: " << line << endl;
            }
        }
        file.close();
    }

    for (int j=0; j< (net.npp); j++){ 
        (net.adj[j]).shrink_to_fit();
        (net.badj[j]).shrink_to_fit();
    }

    // deal w/ pending MPI_Sends...
    MPI_Wait(&sreq,&sstat);
    
    return net;
}


/*
* Creates residual graph from network.
* Does NOT copy any arrays, just points to same place.
* Also performs "pulse 1", i.e. put excess through src.
*/
resgraph setup(network *inet, int rank, int size){
    MPI_Status sstat,rstat;
    resgraph onet;
    onet.src = inet->src;
    onet.sink=inet->sink;
    onet.s_proc = onet.src / inet->std_npp;
    onet.t_proc = onet.sink/inet->std_npp;
    onet.n = inet->n;
    onet.m = inet->m;
    onet.npp = inet->npp;
    onet.std_npp = inet->std_npp;
    onet.odeg = inet->odeg;
    onet.ideg = inet->ideg;

    onet.adj = inet->adj;
    onet.badj = inet->badj;
    onet.ex = (int*) malloc(onet.npp*sizeof(int));
    onet.dex = (int*) malloc(onet.npp * sizeof(int));
    onet.hght = (int*) malloc(onet.npp * sizeof(int));
    onet.hght_p = (int*) malloc(onet.npp * sizeof(int));
    
    onet.aflow = vector<vector<int>>(onet.npp);
    onet.bflow = vector<vector<int>>(onet.npp);
    //onet.active = vector<int>();
    //onet.active_p = vector<int>();

    onet.adj_d = vector<vector<int>>(onet.npp);
    onet.badj_d = vector<vector<int>>(onet.npp);

    //printf("about to loop, proc %d\n", rank);
    for (int v=0; v<onet.npp; v++){
        onet.ex[v]=0;
        onet.dex[v]=0;
        onet.hght[v]=0;
        onet.hght_p[v]=0;

        onet.aflow[v] = vector<int>(onet.odeg[v], 0);
        onet.bflow[v] = vector<int>(onet.ideg[v], 0);

        onet.adj_d[v] = vector<int>(onet.odeg[v],0);
        onet.badj_d[v] = vector<int>(onet.ideg[v],0);
    }

    int num_sent=0;    
    bool win;
    int c,w;
    int s = onet.src - rank*onet.std_npp;
    if (rank == onet.s_proc){
        onet.n_act =0;
        for (int i=0; i<onet.odeg[s]; i++){
            c = onet.adj[s][3*i+1];
            w = onet.adj[s][3*i];
            win = ( (rank*onet.std_npp) <= w ) && (w <((rank+1)*onet.std_npp));
            onet.aflow[s][i] = c;
            if (win){
                onet.bflow[w][ onet.adj[s][3*i+2] / 3 ] = -c;
                onet.ex[w] = c;
                onet.badj_d[w][onet.adj[s][3*i+2]/3] = onet.n;
                if (w != onet.sink) {
                    onet.n_act += 1;
                    onet.active.push(w);
                }
            } else {
                num_sent++;
            }
        }
        onet.hght[s] = onet.n;

        for (int i=1; i< size; i++){
            //printf("sending to %d, from proc %d\n", i, rank);
            MPI_Send(&num_sent,1,MPI_INT,i,i,MPI_COMM_WORLD);
        }
    } else {
        onet.n_act = 0;
        //printf("recv on proc %d\n",rank);
        MPI_Recv(&num_sent,1,MPI_INT,onet.s_proc,rank,MPI_COMM_WORLD, &rstat);
    }

    //MPI_Bcast(&num_sent, 1, MPI_INT, onet.s_proc, MPI_COMM_WORLD);
    
    int buffer[num_sent*3];
    int count=0;
    if (rank == onet.s_proc){
        for (int i=0; i<onet.odeg[s]; i++){
            c = onet.adj[s][3*i+1];
            w = onet.adj[s][3*i];
            win = ( (rank*onet.std_npp) <= w ) && (w <((rank+1)*onet.std_npp));
            if (!win) {
                buffer[3*count] = w;
                buffer[3*count + 1] = c;
                buffer[3*count + 2] = onet.adj[s][3*i+2];
                count++;
            }
        }
    }

    printf("about to bcast, proc %d\n", rank);
    MPI_Bcast(buffer, num_sent, MPI_INT, onet.s_proc, MPI_COMM_WORLD);
    printf("finished bcast, proc %d\n", rank);
    if (rank != onet.s_proc){
        for (int i=0; i< num_sent; i++){
            w = buffer[3*i];
            c = buffer[3*i+1];
            onet.bflow[w][buffer[3*i+2]/3] = -c;
            if (!(w_in(w,onet.std_npp))){
                onet.bflow[w-rank*onet.std_npp][buffer[3*i+2]/3]=-c;
                onet.ex[w-rank*onet.std_npp]=c;
                if (w != onet.sink){
                    onet.n_act +=1;
                    onet.active.push(w);
                }
            }
        }
    }
        
    return onet;

}

/*
* Assumes net & orig are associated
*/
void cleanup(resgraph *net, network *orig){
    free(net->hght_p);
    free(net->hght);
    free(net->dex);
    free(net->ex);

    free(orig->ideg);
    free(orig->odeg);    
}

/*
 * A single pulse, from algo as described in Goldberg & Tarjan
 */
void pulse(resgraph *net_ptr){
    int w,r,ch,e,d_p,dw;
    {
    for (int v=0; v<net_ptr->n; v++){
        net_ptr->dex[v] = 0;
    }
    for (int v=0; v<net_ptr->n; v++){
        net_ptr->hght_p[v] = 4*(net_ptr->n)*(net_ptr->n)*(net_ptr->n); // shouldn't ever be bigger!!
        e = net_ptr->ex[v];
        if ( (e==0) || (v==net_ptr->src) || (v==net_ptr->sink) ){
            net_ptr->hght_p[v] = net_ptr->hght[v];
        } else {
            for (int i=0; i<net_ptr->odeg[v]; i++){
                w = net_ptr->adj[v][3*i];
                dw = net_ptr->hght[w];
                r = net_ptr->adj[v][3*i+1]-net_ptr->aflow[v][i];
                if ( (dw == (net_ptr->hght[v]-1)) && (r>0) ) {
                    ch = (e>=r) ? r : e; //update by min(excess, residue)
                    net_ptr->ex[v] = net_ptr->ex[v]-ch;
                    net_ptr->aflow[v][i] += ch; //update flow!
                    net_ptr->bflow[w][ (net_ptr->adj[v][3*i+2])/2 ] -= ch;
    
                    net_ptr->dex[w]+=ch; // don't update e(w) yet!
                    
                }
                if (net_ptr->ex[v]==0){
                    net_ptr->hght_p[v] = net_ptr->hght[v]; //don't update label
                    break;
                } else if (net_ptr->adj[v][3*i+1] > net_ptr->aflow[v][i]){
                    d_p = net_ptr->hght_p[v];
                    net_ptr->hght_p[v] = (d_p > (dw+1)) ? (dw+1) : d_p;
                }
            }
            if (net_ptr->ex[v]!=0){ 
                for (int i=0; i<net_ptr->ideg[v]; i++){
                    w = net_ptr->badj[v][2*i];
                    e = net_ptr->ex[v];
                    dw = net_ptr->hght[w];
                    if (dw == (net_ptr->hght[v]-1)) {
                        r = 0 - net_ptr->bflow[v][i]; // a backwards edge!
                        ch = (e>=r) ? r : e; //update by min(excess, residue)
                        net_ptr->ex[v] -= ch;
                        net_ptr->bflow[v][i] += ch;
                        net_ptr->aflow[w][ (net_ptr->badj[v][2*i+1])/3 ] -= ch; 
                        
                        net_ptr->dex[w] += ch; // don't update e(w) yet!
                        
                    }
                    if (net_ptr->ex[v] == 0) {
                        net_ptr->hght_p[v] = net_ptr->hght[v]; //don't update label
                        break;
                    } else if (net_ptr->bflow[v][i] < 0) {
                        d_p = net_ptr->hght_p[v];
                        net_ptr->hght_p[v] = (d_p > (dw+1)) ? (dw+1) : d_p;
                    }
                }
            }
        }
    }
    
    net_ptr -> n_act= 0;

    for (int v=0; v<net_ptr->n; v++){
        net_ptr->hght[v] = net_ptr->hght_p[v];
        net_ptr->ex[v] += net_ptr->dex[v];
        if ( (net_ptr->ex[v] > 0) && (v!=net_ptr->src) && (v!=net_ptr->sink)){
            net_ptr -> n_act += 1;
        }
    }
    }
    net_ptr->ex[0] = 0;
}

comm_data setup_cd(resgraph *net, int rank, int size){
    comm_data cds;
    cds.out_req = (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));
    cds.in_req = (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));
    cds.out_dist_req = (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));
    cds.in_dist_req = (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));

    cds.max_deg=0;

    cds.out_bi = vector<vector<int>>(net->npp);
    cds.in_bi = vector<vector<int>>(net->npp);  
    cds.out_flag = vector<vector<unsigned char>>(net->npp);
    cds.in_flag = vector<vector<unsigned char>>(net->npp);
    cds.out_dist_bi = vector<vector<int>>(net->npp);
    cds.in_dist_bi = vector<vector<int>>(net->npp);

    // buffers for send & response to update queries
    // query buff[i] = [v,change,d(v), index]
    // response buff[i] = [acc/rej, w, flow change, d(w), index]
    cds.buff_size = 3*net->npp;
    cds.buff = (int**) malloc((cds.buff_size)*sizeof(int*));
    for (int i=0; i<cds.buff_size; i++){
        cds.buff[i] = (int*) malloc((5*sizeof(int)));
        for (int k=0; k<5; k++){
            cds.buff[i][k]=0;
        }
    }


    for (int v=0; v<net->npp; v++){
        cds.max_deg = max(max(cds.max_deg, net->odeg[v]),net->ideg[v]);
        cds.out_req[v] = (MPI_Request*) malloc((net->odeg[v])*sizeof(MPI_Request));
        cds.in_req[v] = (MPI_Request*) malloc((net->ideg[v])*sizeof(MPI_Request));

        cds.out_dist_req[v] = (MPI_Request*) malloc((net->odeg[v])*sizeof(MPI_Request));
        cds.in_dist_req[v] = (MPI_Request*) malloc((net->ideg[v])*sizeof(MPI_Request));

        for (int i=0; i<net->odeg[v]; i++){
            cds.out_req[v][i] = (MPI_REQUEST_NULL);
            cds.out_flag[v].push_back(NOTHING);
            cds.out_bi[v].push_back(-1);
            cds.out_dist_req[v][i] = MPI_REQUEST_NULL;
            cds.out_dist_bi[v].push_back(-1);
        }

        for (int i=0; i<net->ideg[v]; i++){
            cds.in_req[v][i] = MPI_REQUEST_NULL;
            cds.in_flag[v].push_back(NOTHING);
            cds.in_bi[v].push_back(-1);
            cds.in_dist_req[v][i] = MPI_REQUEST_NULL;
            cds.in_dist_bi[v].push_back(-1);
        }
        cds.buff[v][0] = 0;
        cds.buff[v][1] = 0;
        cds.buff[v][2] = 0;
    }

    cds.fin_req = (MPI_Request*) malloc(size*sizeof(MPI_Request));
    cds.proc_done = (bool*) malloc(size*sizeof(bool));
    for (int i=0; i<size; i++){
        cds.fin_req[i] = MPI_REQUEST_NULL;
        cds.proc_done[i] = 0;
    }

    //cds.avail = queue<int>();
    for (int i =0; i<cds.buff_size; i++){
        cds.avail.push(i);
    }
    cds.arr_of_inds = (int*) malloc(cds.max_deg*sizeof(int));
    return cds;
}

void print_queue(queue<int> *q, int rank, int size){
    int tmp = q->size();
    int v;
    int count =0;
    for (int i=0; i<tmp; i++){
        v = q->front();
        printf("%dnd elem is %d, proc %d\n", count, v, rank);
        count++;
        q->pop();
        q->push(v);
    }
}

/*
 * Asynchronous implementation, as described in Goldberg & Tarjan
 */
void async_pr(resgraph *net, comm_data *cd, int rank,int size){
    wait_for_debugger();
    int gl_w,v,z,loc_w; // verts
    z=0;
    int tally=0;
    int r,c,ch,e,d_p,dw; // p-r vars
    int j,bi,tst; // comm vars

    bool done = 0;
    bool all_done=0;
    MPI_Request done_req=MPI_REQUEST_NULL;
    int progress = 1;
    int dv;

    // start actual algo
    bool win;
    
    while (!all_done){
        while (  (!(net->active.empty()) && (progress != 0))   ) {
            progress=0;
            v = net->active.front();
            net->active.pop();

            // check up on all pending communication to forward nodes
            check_comm(net,v, cd, rank,size);
            listen_finish(net,cd,rank,size);
            listen_distance(net,cd,rank,size);
            listen(net,cd,rank,size);

            for (int i=0; i<net->odeg[v]; i++){
                assert( ( (cd->out_flag[v][i]==NOTHING)&&(cd->out_bi[v][i]==-1) )  ||  ( (cd->out_flag[v][i]!=NOTHING)&&(cd->out_bi[v][i]!=-1)));
                listen(net,cd,rank,size);
                listen_distance(net,cd,rank,size);
                gl_w = net->adj[v][3*i];
                dw = net->adj_d[v][i];
                c = net->adj[v][3*i+1];
                r = c - net->aflow[v][i];
                e = net->ex[v];
                ch = min(e,r);
                win = w_in(gl_w,net->npp);
                printf("gl_v: %d, gl_w: %d, dv: %d, dw: %d, res: %d, ex: %d, ch: %d\n", v+rank*net->std_npp, gl_w, net->hght[v], dw,r, e, ch);
                if ((r>0) && (net->hght[v] == dw+1)) {
                    if (win){
                        progress++;
                        loc_w = gl_w-rank*(net->std_npp);
                        net->aflow[v][i] += ch;
                        net->bflow[loc_w][net->adj[v][3*i+2]/3] -= ch;
                        net->ex[v] -= ch;
                        net->ex[loc_w] += ch;
                    } else {
                        if ( (cd->out_flag[v][i] == NOTHING)){
                            while ( (cd->avail.empty())){
                                MPI_Test(&(cd->out_req[v][i]), &tst, MPI_STATUS_IGNORE);
                                if (tst){
                                    handle_comm(net,v,gl_w, &(net->aflow[v][i]), &(net->adj_d[v][i]), &(cd->out_req[v][i]), &(cd->out_bi[v][i]), &(cd->out_flag[v][i]), cd->buff[bi], &(cd->avail),  cd , rank,size);
                                } else {
                                    // check comm backlog while we wait!
                                    check_comm(net, z,  cd , rank,size);
                                    listen(net, cd ,rank,size);
                                    z = (z+1)%(net->npp);
                                }
                            }
                            progress++;
                            bi = cd->avail.front();
                            cd->avail.pop();
                            cd->buff[bi][0] = v+rank*net->std_npp;
                            cd->buff[bi][1] = ch;
                            cd->buff[bi][2] = net->hght[v];
                            cd->buff[bi][3] = gl_w;
                            cd->buff[bi][4] = net->adj[v][3*i+2]; 
                            MPI_Isend(cd->buff[bi], 5, MPI_INT, gl_w/net->std_npp, FWD_QUERY, MPI_COMM_WORLD, &(cd->out_req[v][i]));
                            cd->out_bi[v][i] = bi;
                            cd->out_flag[v][i] = 4; // 100 (send query fwd)
                        }
                    }
                }
            }
            
            for (int i=0; i<net->ideg[v]; i++){
                assert(  ( (cd->in_bi[v][i]==-1)&&(cd->in_flag[v][i]==NOTHING) )   ||   ((cd->in_bi[v][i]!=-1)&&(cd->in_flag[v][i]!=NOTHING)));
                listen(net, cd ,rank,size);
                listen_distance(net, cd ,rank,size);
                gl_w = net->badj[v][2*i];
                dw = net->badj_d[v][i];
                r = 0 - net->bflow[v][i];
                e = net->ex[v];
                ch = min(e,r);
                win = w_in(gl_w,net->npp);
                printf("gl_v: %d, gl_w: %d, dv: %d, dw: %d, res: %d, ex: %d, ch: %d\n", v+rank*net->std_npp, gl_w, net->hght[v], dw,r, e, ch);
                if ((r>0) && (net->hght[v] == dw+1)) {
                    if (win){
                        progress++;
                        loc_w = gl_w-rank*(net->std_npp);
                        net->aflow[v][i] += ch;
                        net->bflow[loc_w][net->badj[v][2*i+1]/2] -= ch;
                        net->ex[v] -= ch;
                        net->ex[loc_w] += ch;
                    } else {
                        if (cd->in_flag[v][i] ==NOTHING){
                            while ((cd->avail.empty())){
                                MPI_Test(&(cd->in_req[v][i]), &tst, MPI_STATUS_IGNORE);
                                if (tst){
                                    handle_comm(net,v,gl_w, &(net->bflow[v][i]), &(net->badj_d[v][i]), &(cd->in_req[v][i]), &(cd->in_bi[v][i]), &(cd->in_flag[v][i]), cd->buff[bi], &(cd->avail),  cd , rank,size);
                                } else {
                                    // check comm backlog while we wait!
                                    check_comm(net, z,  cd , rank,size);
                                    listen(net, cd ,rank,size);
                                    z = (z+1)%(net->npp);
                                }
                            }
                            progress++;
                            bi = cd->avail.front();
                            cd->avail.pop();
                            cd->buff[bi][0] = v+rank*net->std_npp;
                            cd->buff[bi][1] = ch;
                            cd->buff[bi][2] = net->hght[v];
                            cd->buff[bi][3] = gl_w;
                            cd->buff[bi][4] = net->badj[v][2*i+1]; 
                            MPI_Isend(cd->buff[bi], 5, MPI_INT, gl_w/net->std_npp, BWD_QUERY, MPI_COMM_WORLD, &(cd->in_req[v][i]));
                            cd->in_bi[v][i] = bi;
                            cd->in_flag[v][i] = 5; // 101 (send query bwd)
                        }
                    }
                }
            }
            if (net->ex[v] >0){
                net->active.push(v);
            }
            dv = net->hght[v];
            check_dist(net,v,cd,rank,size);
            if (dv>net->hght[v]){
                progress++;
            }
        }
        listen(net,cd,rank,size);
        listen_distance(net,cd,rank,size);
        listen_finish(net, cd, rank,size);
        if ( cd->proc_done[rank==1]){
            printf("I'm done! proc %d\n", rank);
            if (!net->active.empty()){
                cd->proc_done[rank]=0;
                for (int i=0; i<size; i++){
                    if (i!=rank){
                        MPI_Wait(&(cd->fin_req[i]), MPI_STATUS_IGNORE);
                        MPI_Isend(&(cd->proc_done[rank]), 1, MPI_CXX_BOOL, i, FINISH, MPI_COMM_WORLD, &(cd->fin_req[i]));
                    }
                }
            }
            
        } else {
            cd->proc_done[rank]=1;
            for (int i=0; i<size; i++){
                if (i!= rank){
                    MPI_Wait(&(cd->fin_req[i]), MPI_STATUS_IGNORE);
                    MPI_Isend(&(cd->proc_done[rank]), 1, MPI_CXX_BOOL, i,  FINISH, MPI_COMM_WORLD, &(cd->fin_req[i]));
                }
            }
        }
        done = 1;
        for (int i=0; i<size; i++){
            done = (done && cd->proc_done[i]); 
        }

        progress = 1;
    }
}


/*
* Writes output to file <file>.
* Assumes that <net> has finished computation
* i.e. e[v] == 0 for all v
*/
void output(resgraph net, string filename, double time){
    ofstream outf;
    outf.open(filename);
    outf << "c " << filename << endl;
    outf << "c" << endl;
    outf << "c Time: " << time << endl;
    outf << "c" << endl;
    outf << "c Solution" << endl;
    outf << "s " << net.ex[net.sink] << endl;
    outf << "c" << endl;
    
    outf << "c SOURCE DEST FLOW" << endl;

    int w, f;
    for (int v=0; v< net.n; v++){
        for (int i=0; i< net.odeg[v]; i++){
            w = net.adj[v][3*i];
            f = net.aflow[v][i];
            outf << "f " << v << " " << w << " " << f << endl;
        }
    }
}

/*
* Checks if the two provided solution files give the same flow
*/
bool check(string correct, string test){
    int c_flow=-1;
    int t_flow=-1;
    
    string line;
    ifstream file;
    
    file.open(correct, ios::in);
    if (file.is_open()) {
        while (getline(file, line)){
            if (line[0]=='s'){
                c_flow = atoi(&(line[2]));
            } 
        }
        file.close();
    }
    file.open(test, ios::in);
    if (file.is_open()) {
        while (getline(file, line)){
            if (line[0]=='s'){
                t_flow = atoi(&(line[2]));
            } 
        }
        file.close();
    }

    cout << "Correct: " << c_flow << endl;
    cout << "Test: " << t_flow << endl;

    return (c_flow == t_flow);
}

