#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <vector>
#include "mpi-pr.h"
using namespace std;

#define max(a,b) ( ((a)>(b)) ? (a) : (b) )
#define min(a,b) ( ((a)<(b)) ? (a) : (b) )
#define w_in(w,npp) ( (rank*(npp) <= (w)) && ((w) < ((rank+1)*(npp))) )

#define QUERY 123
#define RESPONSE 124


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
    int rbuf[1] = {0}; // v,w,c(v,w), pos
    int sbuf[1] = {0}; // v,w,c(v,w), pos
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
                        MPI_Irecv(rbuf,1,MPI_INT, w/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                        // wait for prev send to finish, so know sbuf avail
                        MPI_Wait(&sreq, &sstat);
                        sbuf[0] = jv;
                        MPI_Isend(sbuf, 1, MPI_INT, w/npp, v*net.n + w, MPI_COMM_WORLD,&sreq);
                        net.adj[v-rank*npp].push_back(w);
                        net.adj[v-rank*npp].push_back(cap);
                        // wait for curr rec to finish, so can use data
                        MPI_Wait(&rreq, &rstat);
                        net.adj[v-rank*npp].push_back(rbuf[0]); 
                    }
                } else if (win) {
                    MPI_Irecv(rbuf,1,MPI_INT,v/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                    net.ideg[w-rank*npp]++;
                    jw = net.badj[w-rank*npp].size();
                    net.badj[w-rank*npp].push_back(v);
                    MPI_Wait(&sreq, &sstat);
                    sbuf[0] = jw;
                    MPI_Isend(sbuf,1,MPI_INT, v/npp, v*net.n+w, MPI_COMM_WORLD, &rreq);
                    MPI_Wait(&rreq, &rstat);
                    jw = rbuf[0];
                    net.badj[w-rank*npp].push_back(rbuf[0]);

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
    onet.active = vector<int>();
    onet.active_p = vector<int>();

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
                    onet.active.push_back(w);
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
                    onet.active.push_back(w);
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


/*
 * Asynchronous implementation, as described in Goldberg & Tarjan
 */
void async_pr(resgraph *net,int rank,int size){
    int w,loc_w,r,c,ch,e,d_p,dw;
    MPI_Status qstat, rstat;
    MPI_Request *out_req[net->npp];
    MPI_Request *in_req[net->npp];
    // out_flag[v][i] >= 0 when (v,net->adj[3*i]) has a pending query
    // and if >=0, this int is ind of buffer in buff
    vector<vector<int>> out_flag = vector<vector<int>>(net->npp);
    vector<vector<int>> in_flag = vector<vector<int>>(net->npp);
    // buffers for send & response to update queries
    int buff[net->npp][3];
    // indices of available buffers in buff
    vector<int> avail = vector<int>(0); 
    for (int v=0; v<net->npp; v++){
        out_req[v] = (MPI_Request*) malloc((net->odeg[v])*sizeof(MPI_Request));
        in_req[v] = (MPI_Request*) malloc((net->ideg[v])*sizeof(MPI_Request));
        for (int i=0; i<net->odeg[v]; i++){
            out_req[v][i] = (MPI_REQUEST_NULL);
            out_flag[v].push_back(-1);
            buff[v][0] = buff[v][1] = buff[v][2]=0;
        }
        for (int i=0; i<net->ideg[v]; i++){
            in_req[v][i] = MPI_REQUEST_NULL;
            in_flag[v].push_back(-1);
        }
    }
    bool win;
    for (int v=0; v<net->npp; v++){
        e = net->ex[v];
        // if v is active:
        if ( ! ( (e==0) || (v==net->src) || (v==net->sink) ) ) {
            // push flow on outgoing
            for (int i=0; i<net->odeg[v]; i++){
                w = net->adj[v][3*i];
                dw = net->adj_d[v][i];
                c = net->adj[v][3*i+1];
                r = c - net->aflow[v][i];
                win = w_in(w,net->npp);
                if ((r>0) && (net->hght[v] == dw+1)) {
                    if (win){
                        loc_w = w-rank*(net->std_npp);
                        ch = min(e,r);
                        net->aflow[v][i] += ch;
                        net->bflow[loc_w][net->adj[v][3*i+2]/3] -= ch;
                        net->ex[v] -= ch;
                        net->ex[loc_w] += ch;
                    } else {
                        if (out_req[v][i] != MPI_REQUEST_NULL){
                            // since send completed, buff available
                            // deal w/ prev update query
                            MPI_Recv(buff[out_flag[v][i]], 3, MPI_INT, w/(net->std_npp), RESPONSE, MPI_COMM_WORLD, &rstat);
                            
                        }
                    }
                }
            }
        }

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

