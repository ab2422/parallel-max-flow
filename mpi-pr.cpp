#include "mpi-pr.h"
#include "mpi-data.h"
#include "mpi-comm.h"
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <string>
#include <cstring>
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

#define MAX_TAG 32767
#define INIT_ADJ 52
// prev MAX_TAG 32767

bool is_src_loc(resgraph *net, int loc_v, int rank, int size){
    return  (net->s_proc == rank) && ( (loc_v +rank*net->std_npp) == net->src);
}

bool is_sink_loc(resgraph *net, int loc_v, int rank, int size){
    return (net->sink == (loc_v + rank*net->std_npp));
}

void wait_for_debugger(){
    volatile int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i){}
}

void print_vector(vector<int> vec){
    string line = "Size: "+ to_string(vec.size())+ ", Elems: ";
    for (int i=0; i<vec.size(); i++){
        line = line + to_string(vec[i]) + " ";
    }
    line  = line + "\n";
    cout << (line);
}

void print_flow(resgraph *net, int rank){
    string line = "";
    for (int v=0; v<net->npp; v++){
        line = "Vert: "+ to_string(v+rank*net->std_npp) +" Out flow: ";
        for (int j=0; j<net->odeg[v]; j++){
            line = line + to_string(net->flow[0][v][j]) + " ";
        }
        line = line + "\n";
        cout << line;
    }
    for (int v=0; v<net->npp; v++){
        line = "Vert: "+ to_string(v+rank*net->std_npp) +" IN flow: ";
        for (int j=0; j<net->ideg[v]; j++){
            line = line + to_string(net->flow[1][v][j]) + " ";
        }
        line = line + "\n";
        cout << line;
    }
}

void print_total(resgraph *net, int rank){
    string line = "";
    int flow=0;
    int loc_s;
    if (rank==net->s_proc){
        loc_s = net->src - rank*net->std_npp;
        for (int i=0; i<net->odeg[loc_s]; i++){
            flow += net->flow[0][loc_s][i];
        }
        line = "Total flow: "+to_string(flow) +"\n";
        cout << line;
    }
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
    network net;
    string line;
    ifstream file;
    file.open(filename, ios::in);
    int i=0;
    int npp=0;
    bool win=0;
    int v,w,cap,f,jv,jw,tag;
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
                    cout << "Invalid file format! Bad line is: " << line << endl;
                } else {
                    i=6;
                    net.n = atoi( &(line[6]));
                    npp = (net.n + size-1  )/size;
                    net.std_npp = npp;
                    net.npp = min(npp, max(1,net.n-rank*npp));
                    while (isspace(line[i])){
                        i++;
                    }
                    while ((isalnum(line[i]))){
                        i++;
                    }
                    net.m = atoi(&(line[i]));
                    net.odeg = (int*) malloc((net.npp)*sizeof(int));
                    net.ideg = (int*) malloc((net.npp)*sizeof(int));
                    net.cap = vector<vector<int>>(net.npp);
                    net.adj = vector<vector<vector<int>>>(2);
                    net.adj[0].resize(net.npp);
                    net.adj[1].resize(net.npp);
                    for (int j=0; j<net.npp; j++){
                        net.adj[0][j].reserve(4*net.m/(size*net.n)+1);
                        net.adj[1][j].reserve(4*net.m/(size*net.n)+1);
                        net.cap[j].reserve(4*net.m/(size*net.n)+1);
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
                tag = INIT_ADJ;//((v%MAX_TAG)*(w%MAX_TAG)+w)%MAX_TAG;
                // if vin: 
                if (((rank*npp)<=v) && (v < (rank+1)*npp)){
                    net.odeg[v-rank*npp]++;
                    net.cap[v-rank*npp].push_back(cap);
                    jv = net.adj[0][v-rank*net.std_npp].size();//pos w in v's
                    net.adj[0][v-rank*net.std_npp].push_back(w);
                    if (win) {
                        net.ideg[w-rank*npp]++;
                        jw = net.adj[1][w-rank*npp].size();
                        
                        net.adj[0][v-rank*npp].push_back(jw); //jw = v pos in w list
        
                        net.adj[1][w-rank*npp].push_back(v); // bacwards edge
                        net.adj[1][w-rank*npp].push_back(jv); //jv = w pos in v list
                    }else{
                        sbuf = jv;
			MPI_Sendrecv(&sbuf,1,MPI_INT, w/net.std_npp,tag,&rbuf, 1, MPI_INT, w/net.std_npp,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			/*
                        MPI_Irecv(&rbuf,1,MPI_INT, w/net.std_npp, tag, MPI_COMM_WORLD, &rreq);
                        // wait for prev send to finish, so know sbuf avail
                        MPI_Wait(&sreq, MPI_STATUS_IGNORE);
                        sbuf = jv;
                        MPI_Isend(&sbuf, 1, MPI_INT, w/net.std_npp, tag, MPI_COMM_WORLD,&sreq);
                        // wait for curr rec to finish, so can use data
                        MPI_Wait(&rreq,MPI_STATUS_IGNORE);
			*/
                        net.adj[0][v-rank*npp].push_back(rbuf); 
                    }
                } else if (win) {
                    net.ideg[w-rank*npp]++;
		    jw = net.adj[1][w-rank*net.std_npp].size();
		    sbuf = jw;
                    MPI_Sendrecv(&sbuf,1,MPI_INT, v/net.std_npp, tag, &rbuf,1,MPI_INT,v/net.std_npp, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    /* 
                    MPI_Irecv(&rbuf,1,MPI_INT,v/npp, tag, MPI_COMM_WORLD, &rreq);
                    jw = net.adj[1][w-rank*npp].size(); // pos of v in w's
                    net.adj[1][w-rank*npp].push_back(v);
                    if ( (w< 10+net.std_npp)&&(v==0)){
			printf("Curr %d, to src, Deg: %d, jw: %d\n", w,net.ideg[w-rank*npp], jw);
		    }
                    MPI_Wait(&sreq, MPI_STATUS_IGNORE);
                    sbuf = jw;
                    MPI_Isend(&sbuf,1,MPI_INT, v/npp, tag, MPI_COMM_WORLD, &rreq);
                    MPI_Wait(&rreq, MPI_STATUS_IGNORE);
                    jv = rbuf;
		    */
                    net.adj[1][w-rank*npp].push_back(rbuf);

                }
            } else if (line[0]!='c'){
                cout << "Invalid file. Bad line is: " << line << endl;
            }
        }
        file.close();
    } else {
        cerr << "File open error: " << strerror(errno) << endl;
    }

    printf("adj size: %d, on proc %d\n", net.adj.size(),rank);
    printf("Adj[0] size: %d, Adj[1] size: %d, on proc %d\n", net.adj[0].size(), net.adj[1].size(), rank);

    for (int j=0; j< (net.npp); j++){ 
        (net.adj[0][j]).shrink_to_fit();
        (net.adj[1][j]).shrink_to_fit();
        net.cap[j].shrink_to_fit();
    }

    //print_vector(net.adj[1][60000-rank*net.std_npp]);
 
    
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
    onet.cap = inet->cap;
    onet.ex = (int*) malloc(onet.npp*sizeof(int));
    onet.dex = (int*) malloc(onet.npp * sizeof(int));
    onet.hght = (int*) malloc(onet.npp * sizeof(int));
    onet.hght_p = (int*) malloc(onet.npp * sizeof(int));
    
    onet.flow = vector<vector<vector<int>>>(2);
    onet.flow[0].resize(onet.npp);
    onet.flow[1].resize(onet.npp);
    onet.adj_d = vector<vector<vector<int>>>(2);
    onet.adj_d[0].resize(onet.npp);
    onet.adj_d[1].resize(onet.npp);

    //printf("about to loop, proc %d\n", rank);
    for (int v=0; v<onet.npp; v++){
        onet.ex[v]=0;
        onet.dex[v]=0;
        onet.hght[v]=0;
        onet.hght_p[v]=0;

        onet.flow[0][v].resize(onet.odeg[v],0);
        onet.flow[1][v].resize(onet.ideg[v], 0);

        onet.adj_d[0][v].resize(onet.odeg[v],0);
        onet.adj_d[1][v].resize(onet.ideg[v],0);
    }

    int num_sent=0;    
    bool win;
    int c,gl_w, loc_w,jw;
    int s = onet.src - rank*onet.std_npp;
    if (rank == onet.s_proc){
        onet.n_act =0;
        for (int i=0; i<onet.odeg[s]; i++){
            c = onet.cap[s][i];
            gl_w = onet.adj[0][s][2*i];
            win = ( (rank*onet.std_npp) <= gl_w ) && (gl_w <((rank+1)*onet.std_npp));
            onet.flow[0][s][i] = c;
            if (win){
                jw = onet.adj[0][s][2*i+1];
                loc_w = gl_w - rank*onet.std_npp;
                onet.flow[1][loc_w][jw/2 ] = -c;
                onet.ex[loc_w] = c;
                onet.adj_d[1][loc_w][jw/2] = onet.n;
                if (gl_w != onet.sink) {
                    onet.n_act += 1;
                    onet.active.push(gl_w);
                }
            } else {
                num_sent++;
            }
        }
        onet.hght[s] = onet.n;

        for (int i=1; i< size; i++){
            //printf("sending to %d, from proc %d\n", i, rank);
            //MPI_Send(&num_sent,1,MPI_INT,i,i,MPI_COMM_WORLD);
        }
    } else {
        onet.n_act = 0;
        //printf("recv on proc %d\n",rank);
        //MPI_Recv(&num_sent,1,MPI_INT,onet.s_proc,rank,MPI_COMM_WORLD, &rstat);
    }

    MPI_Bcast(&num_sent, 1, MPI_INT, onet.s_proc, MPI_COMM_WORLD);
    
    int buffer[num_sent*3];
    int count=0;
    if (rank == onet.s_proc){
        for (int i=0; i<onet.odeg[s]; i++){
            c = onet.cap[s][i];
            gl_w = onet.adj[0][s][2*i];
            win = ( (rank*onet.std_npp) <= gl_w ) && (gl_w <((rank+1)*onet.std_npp));
            if (!win) {
                buffer[3*count] = gl_w;
                buffer[3*count + 1] = c;
                buffer[3*count + 2] = onet.adj[0][s][2*i+1];
                count++;
            }
        }
    }

    //printf("about to bcast, proc %d\n", rank);
    MPI_Bcast(buffer, num_sent, MPI_INT, onet.s_proc, MPI_COMM_WORLD);
    //printf("finished bcast, proc %d\n", rank);
    if (rank != onet.s_proc){
        for (int i=0; i< num_sent; i++){
            gl_w = buffer[3*i];
            loc_w = gl_w - rank*onet.std_npp;
            c = buffer[3*i+1];
            jw = buffer[3*i+2];
            if ((w_in(gl_w,onet.std_npp))){
                if ( (jw/2>onet.ideg[loc_w]) || (jw<0)){
                    printf("w: %d, jw: %d, indeg: %d, proc: %d\n", gl_w, jw, onet.ideg[loc_w], rank);
            }    
                onet.ex[loc_w]=c;
                onet.flow[1][loc_w][jw/2]=-c;
                if (gl_w != onet.sink){
                    onet.n_act +=1;
                    onet.active.push(gl_w);
                }
            }
        }
    }
        
    return onet;

}

/*
* Assumes net & orig are associated
*/
void cleanup(resgraph *net, network *orig, comm_data *cd){
    for (int v=0; v<net->npp; v++){
        free(cd->edge_req[0][v]);
        free(cd->dist_req[0][v]);
    }
    for (int v=0; v<net->npp; v++){
        free(cd->edge_req[1][v]);
        free(cd->dist_req[1][v]);
    }
    free(cd->edge_req[0]);
    free(cd->edge_req[1]);
    free(cd->edge_req);
    free(cd->dist_req[0]);
    free(cd->dist_req[1]);
    free(cd->dist_req);

    free(cd->fin_req);
    free(cd->proc_done);
    free(cd->arr_of_inds);
    for (int bi=0; bi<cd->buff_size; bi++){
        free(cd->buff[bi]);
    }
    free(cd->buff);

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
/*
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
*/

comm_data setup_cd(resgraph *net, int rank, int size){
    comm_data cds;

    cds.num_times_done = 0;

    cds.edge_bi =   vector<vector<vector<int>>>(2);
    cds.edge_bi[0].resize(net->npp);
    cds.edge_bi[1].resize(net->npp);
    cds.edge_flag = vector<vector<vector<unsigned char>>>(2);
    cds.edge_flag[0].resize(net->npp);
    cds.edge_flag[1].resize(net->npp);
    cds.dist_bi =   vector<vector<vector<int>>>(2);
    cds.dist_bi[0].resize(net->npp);
    cds.dist_bi[1].resize(net->npp);

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

    cds.edge_req = (MPI_Request***) malloc(2*sizeof(MPI_Request**));
    cds.edge_req[0]= (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));
    cds.edge_req[1]= (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));

    cds.dist_req = (MPI_Request***) malloc((2)*sizeof(MPI_Request*));
    cds.dist_req[0]= (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));
    cds.dist_req[1]= (MPI_Request**) malloc((net->npp)*sizeof(MPI_Request*));

    cds.max_deg=0;

    for (int v=0; v<net->npp; v++){
        cds.max_deg = max(max(cds.max_deg, net->odeg[v]),net->ideg[v]);
        cds.edge_req[0][v] = (MPI_Request*) malloc((net->odeg[v])*sizeof(MPI_Request));
        cds.dist_req[0][v] = (MPI_Request*) malloc((net->odeg[v])*sizeof(MPI_Request));
        cds.edge_req[1][v] = (MPI_Request*) malloc((net->ideg[v])*sizeof(MPI_Request));
        cds.dist_req[1][v] = (MPI_Request*) malloc((net->ideg[v])*sizeof(MPI_Request));

        for (int i=0; i<net->odeg[v]; i++){
            cds.edge_req[0][v][i] = (MPI_REQUEST_NULL);
            cds.edge_flag[0][v].push_back(NOTHING);
            cds.edge_bi[0][v].push_back(-1);
            cds.dist_req[0][v][i] = MPI_REQUEST_NULL;
            cds.dist_bi[0][v].push_back(-1);
        }
        

        for (int i=0; i<net->ideg[v]; i++){
            cds.edge_req[1][v][i] = MPI_REQUEST_NULL;
            cds.edge_flag[1][v].push_back(NOTHING);
            cds.edge_bi[1][v].push_back(-1);
            cds.dist_req[1][v][i] = MPI_REQUEST_NULL;
            cds.dist_bi[1][v].push_back(-1);
        }
        cds.buff[v][0] = 0;
        cds.buff[v][1] = 0;
        cds.buff[v][2] = 0;
    }

    cds.fin_req = (MPI_Request*) malloc(size*sizeof(MPI_Request));
    cds.proc_done = (bool*) malloc(size*sizeof(bool));
    cds.num_times_done=0;
    cds.all_done=0;
    cds.ring_stat=0;
    cds.prev_total=0;
    cds.nxt_total=0;
    cds.ring_req = MPI_REQUEST_NULL;
    cds.ring_flag=0;
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
    int gl_w,v,z,loc_w; // verts
    z=0;
    int tally=0;
    int r,c,ch,e,d_p,dw; // p-r vars
    int dv;
    int j,bi,tst; // comm vars

    // start actual algo
    bool win;
    
    while (!(cd->all_done)){
        while (  (!(net->active.empty()))   ) {
            v = net->active.front();
            net->active.pop();

            // check up on all pending communication to forward nodes
            check_comm(net,v, cd, rank,size);
            listen_distance(net,cd,rank,size);
            listen(net,cd,rank,size);

            for (int i=0; i<net->odeg[v]; i++){
                assert( ( (cd->edge_flag[0][v][i]==NOTHING)&&(cd->edge_bi[0][v][i]==-1) )  ||  ( (cd->edge_flag[1][v][i]!=NOTHING)&&(cd->edge_bi[1][v][i]!=-1)));
                listen(net,cd,rank,size);
                listen_distance(net,cd,rank,size);
                gl_w = net->adj[0][v][2*i];
                dw = net->adj_d[0][v][i];
                c = net->cap[v][i];
                r = c - net->flow[0][v][i];
                e = net->ex[v];
                ch = min(e,r);
                win = w_in(gl_w,net->npp);
                //printf("gl_v: %d, gl_w: %d, dv: %d, dw: %d, res: %d, ex: %d, ch: %d\n", v+rank*net->std_npp, gl_w, net->hght[v], dw,r, e, ch);
                if ((r>0) && (net->hght[v] == dw+1)) {
                    if (win){
                        loc_w = gl_w-rank*(net->std_npp);
                        net->flow[0][v][i] += ch;
                        net->flow[1][loc_w][ net->adj[0][v][2*i+1]/2] -= ch;
                        net->ex[v] -= ch;
                        net->ex[loc_w] += ch;
                    } else {
                        if ( (cd->edge_flag[0][v][i] == NOTHING)){
                            while ( (cd->avail.empty())){
                                MPI_Test(&(cd->edge_req[0][v][i]), &tst, MPI_STATUS_IGNORE);
                                if (tst){
                                    handle_comm(net,v,gl_w, 0, i, cd,rank,size);
                                } else {
                                    // check comm backlog while we wait!
                                    check_comm(net, z,  cd , rank,size);
                                    listen(net, cd ,rank,size);
                                    z = (z+1)%(net->npp);
                                }
                            }
                            loc_w = gl_w-rank*(net->std_npp);
                            net->flow[0][v][i] += ch;
                            net->ex[v] -= ch;

                            bi = cd->avail.front();
                            cd->avail.pop();
                            cd->buff[bi][0] = v+rank*net->std_npp;
                            cd->buff[bi][1] = ch;
                            cd->buff[bi][2] = net->hght[v];
                            cd->buff[bi][3] = gl_w;
                            cd->buff[bi][4] = net->adj[0][v][2*i+1]; 
                            MPI_Isend(cd->buff[bi], 5, MPI_INT, gl_w/net->std_npp, FWD_QUERY, MPI_COMM_WORLD, &(cd->edge_req[0][v][i]));
                            cd->edge_bi[0][v][i] = bi;
                            cd->edge_flag[0][v][i] = 4; // 100 (send query fwd)
                        }
                    }
                }
            }
            
            for (int i=0; i<net->ideg[v]; i++){
                assert(  ( (cd->edge_bi[1][v][i]==-1)&&(cd->edge_flag[1][v][i]==NOTHING) )   ||   ((cd->edge_bi[1][v][i]!=-1)&&(cd->edge_flag[1][v][i]!=NOTHING)));
                listen(net, cd ,rank,size);
                listen_distance(net, cd ,rank,size);
                gl_w = net->adj[1][v][2*i];
                dw = net->adj_d[1][v][i];
                r = 0 - net->flow[1][v][i];
                e = net->ex[v];
                ch = min(e,r);
                win = w_in(gl_w,net->npp);
                //printf("gl_v: %d, gl_w: %d, dv: %d, dw: %d, res: %d, ex: %d, ch: %d\n", v+rank*net->std_npp, gl_w, net->hght[v], dw,r, e, ch);
                if ((r>0) && (net->hght[v] == dw+1)) {
                    if (win){
                        loc_w = gl_w-rank*(net->std_npp);
                        net->flow[1][v][i] += ch;
                        net->flow[0][loc_w][net->adj[1][v][2*i+1]/2] -= ch;
                        net->ex[v] -= ch;
                        net->ex[loc_w] += ch;
                    } else {
                        if (cd->edge_flag[1][v][i] ==NOTHING){
                            while ((cd->avail.empty())){
                                MPI_Test(&(cd->edge_req[1][v][i]), &tst, MPI_STATUS_IGNORE);
                                if (tst){
                                    handle_comm(net,v,gl_w, 1, i, cd,rank,size);
                                } else {
                                    // check comm backlog while we wait!
                                    check_comm(net, z,  cd , rank,size);
                                    listen(net, cd ,rank,size);
                                    z = (z+1)%(net->npp);
                                }
                            }
                            net->flow[1][v][i] += ch;
                            net->ex[v] -= ch;

                            bi = cd->avail.front();
                            cd->avail.pop();
                            cd->buff[bi][0] = v+rank*net->std_npp;
                            cd->buff[bi][1] = ch;
                            cd->buff[bi][2] = net->hght[v];
                            cd->buff[bi][3] = gl_w;
                            cd->buff[bi][4] = net->adj[1][v][2*i+1]; 
                            MPI_Isend(cd->buff[bi], 5, MPI_INT, gl_w/net->std_npp, BWD_QUERY, MPI_COMM_WORLD, &(cd->edge_req[1][v][i]));
                            cd->edge_bi[1][v][i] = bi;
                            cd->edge_flag[1][v][i] = 5; // 101 (send query bwd)
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
            }
        }
        check_comm(net,z,cd,rank,size);
        z=(z+1)%net->npp;
        listen(net,cd,rank,size);
        listen_distance(net,cd,rank,size);
        handle_finish(net, cd, rank,size);
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
            w = net.adj[0][v][2*i];
            f = net.flow[0][v][i];
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

