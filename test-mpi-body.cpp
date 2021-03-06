#include "catch2/catch.hpp"
#include "mpi-pr.h"
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

//#define DPREFIX "parallel-max-flow/data/"
#define DPREFIX "data/"

TEST_CASE("Parser: serial basic net test", "[1proc]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    SECTION("check node setup"){
        REQUIRE(net.n == 6);
        REQUIRE(net.m == 8);
        REQUIRE(net.src == 0);
        REQUIRE(net.sink == 5);
        REQUIRE(net.npp == 6);
    }
    SECTION("check out degree setup"){
        REQUIRE(net.odeg[0]==2);
        REQUIRE(net.odeg[1]==2);
        REQUIRE(net.odeg[2]==2);
        REQUIRE(net.odeg[3]==1);
        REQUIRE(net.odeg[4]==1);
        REQUIRE(net.odeg[5]==0);
    }
    SECTION("check in degree setup"){
        REQUIRE(net.ideg[0]==0);
        REQUIRE(net.ideg[1]==1);
        REQUIRE(net.ideg[2]==1);
        REQUIRE(net.ideg[3]==2);
        REQUIRE(net.ideg[4]==2);
        REQUIRE(net.ideg[5]==2);
    }
    SECTION("check adjacency matrix"){
        for (int v=0; v<net.n; v++){
            REQUIRE(3*net.odeg[v]==(net.adj[0][v]).size());
        }
        REQUIRE(((net.adj[0][0]==vector<int>({1,0, 2,0}) ) || (net.adj[0][0] == vector<int>({2,0, 1,0}))));
        //REQUIRE(((net.adj[0][1]==vector<int>({3,5,4, 5}) ) || (net.adj[0][0] == vector<int>({4, 5,3,5}))));
        //REQUIRE(((net.adj[0][2]==vector<int>({3,5,4, 5}) ) || (net.adj[0][0] == vector<int>({4, 5,3,5}))));
        REQUIRE(((net.adj[0][3][0]==5) && (net.adj[0][3][1]==15)));
        REQUIRE(((net.adj[0][4][0]==5) && (net.adj[0][4][1]==5)));
        REQUIRE(net.adj[0][5]==vector<int>({}));
    }
    SECTION("check backwards adj matrix"){
        for (int w=0; w<net.n; w++){
            REQUIRE(2*net.ideg[w] == net.adj[1][w].size() );
        }
        REQUIRE(net.adj[1][0]==vector<int>({}));
        REQUIRE(((net.adj[1][1][0]==0)));
        REQUIRE(((net.adj[1][2][0]==0)));
        //REQUIRE(((net.adj[1][3]==vector<int>({1,5,2,5})) || (net.adj[1][3]==vector<int>({2,5,1,5}))));
        //REQUIRE(((net.adj[1][4]==vector<int>({1,5,2,5})) || (net.adj[1][3]==vector<int>({2,5,1,5}))));
        REQUIRE(((net.adj[1][5]==vector<int>({3,0, 4,0}))     ));//|| (net.adj[1][3]==vector<int>({4,5,1, 3,15,1}))));
    }
    SECTION("check compatibility"){
        for (int v=0; v<net.n; v++){
            for (int i=0; i<net.odeg[v]; i++){
                REQUIRE(v == net.adj[1][ net.adj[0][v][2*i] ][ net.adj[0][v][2*i+1] ]);
            }
            for (int i=0; i<net.ideg[v]; i++){
                REQUIRE(v == net.adj[0][ net.adj[1][v][2*i] ][ net.adj[1][v][2*i+1] ]);
            }
        }
    }
}


TEST_CASE("Parser: parallel basic net test", "[2proc]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    SECTION("check node setup"){
        REQUIRE(net.n == 6);
        REQUIRE(net.m == 8);
        REQUIRE(net.src == 0);
        REQUIRE(net.sink == 5);
        REQUIRE(net.npp == 3);
    }
    SECTION("check out degree setup"){
        if (rank==0){
            REQUIRE(net.odeg[0]==2);
            REQUIRE(net.odeg[1]==2);
            REQUIRE(net.odeg[2]==2);
        } else{
            REQUIRE(net.odeg[0]==1);
            REQUIRE(net.odeg[1]==1);
            REQUIRE(net.odeg[2]==0);
        }
    }
    SECTION("check in degree setup"){
        if (rank==0){
            REQUIRE(net.ideg[0]==0);
            REQUIRE(net.ideg[1]==1);
            REQUIRE(net.ideg[2]==1);
        } else {
            REQUIRE(net.ideg[0]==2);
            REQUIRE(net.ideg[1]==2);
            REQUIRE(net.ideg[2]==2);
        }
    }

    SECTION("check capacities"){
        if (rank==0) {
            REQUIRE(((net.cap[0]==vector<int>({5,15}))||(net.cap[0]==vector<int>({15,5})))); 
            REQUIRE(net.cap[1]==net.cap[2]);
            REQUIRE(net.cap[1] == vector<int>({5,5}));
        }
    }
    SECTION("check adjacency matrix"){
        for (int v=0; v<net.npp; v++){
            REQUIRE(2*net.odeg[v]==(net.adj[0][v]).size());
        }
        if (rank==0){
            REQUIRE(((net.adj[0][0]==vector<int>({1,0, 2,0}) ) || (net.adj[0][0] == vector<int>({2,0, 1,0}))));
            //REQUIRE(((net.adj[0][1]==vector<int>({3,5,4, 5}) ) || (net.adj[0][0] == vector<int>({4, 5,3,5}))));
            //REQUIRE(((net.adj[0][2]==vector<int>({3,5,4, 5}) ) || (net.adj[0][0] == vector<int>({4, 5,3,5}))));
        } else {
            REQUIRE(((net.adj[0][0][0]==5) && (net.adj[0][0][1]==0)));
            REQUIRE(((net.adj[0][1][0]==5) && (net.adj[0][1][1]==2)));
            REQUIRE(net.adj[0][2]==vector<int>({}));
        }
    }
    SECTION("check backwards adj matrix"){
        for (int w=0; w<net.npp; w++){
            REQUIRE(2*net.ideg[w] == net.adj[1][w].size() );
        }
        if (rank==0){
            REQUIRE(net.adj[1][0]==vector<int>({}));
            REQUIRE(net.adj[1][1][0]==0);
            REQUIRE(net.adj[1][2][0]==0);
        } else {
    //        REQUIRE(((net.adj[1][0][0]==vector<int>({1,5,2,5})) || (net.adj[1][0][0]==vector<int>({2,5,1,5}))));
    //        REQUIRE(((net.adj[1][0][1]==vector<int>({1,5,2,5})) || (net.adj[1][0][1]==vector<int>({2,5,1,5}))));
            REQUIRE(  net.adj[1][2]==vector<int>({3,0, 4,0})  );
        }
    }
    //SECTION("check compatibility"){
    //    for (int v=0; v<net.npp; v++){
    //        for (int i=0; i<net.odeg[v]; i++){
    //            REQUIRE(v == net.adj[1][0][ net.adj[0][v][3*i] ][ net.adj[0][v][3*i+2] ]);
    //        }
    //        for (int i=0; i<net.ideg[v]; i++){
    //            REQUIRE(v == net.adj[0][ net.adj[1][0][v][2*i] ][ net.adj[1][0][v][2*i+1] ]);
    //        }
    //    }
    //}

}


TEST_CASE("SETUP: serial basic net test", "[1proc]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    resgraph graph = setup(&net,rank,size);
    SECTION("shared attributes"){
        REQUIRE(net.n==graph.n);
        REQUIRE(net.m==graph.m);
        REQUIRE(net.src==graph.src);
        REQUIRE(net.sink==graph.sink);
        REQUIRE(net.odeg == graph.odeg);
        REQUIRE(net.ideg == graph.ideg);
        for (int v=0; v<graph.n; v++){
            REQUIRE(net.adj[0][v] == graph.adj[0][v]);
            REQUIRE(net.adj[1][v] == graph.adj[1][v]);
        }
    }
    
    SECTION("excess"){
        REQUIRE(graph.n_act == 2);
        REQUIRE(graph.ex[0] == 0);
        REQUIRE(graph.ex[1] == 5);
        REQUIRE(graph.ex[2] == 15);
        REQUIRE(graph.ex[3] == 0);
        REQUIRE(graph.ex[4] == 0);
        REQUIRE(graph.ex[5] == 0);
    }

    SECTION("labels"){
        REQUIRE(graph.hght[0]==6);
        for (int v=1; v<graph.n; v++){
            REQUIRE(graph.hght[v]==0);
        }
    }

    SECTION("flows"){
        REQUIRE(graph.flow[0][0][0] == 5);
        REQUIRE(graph.flow[1][1][0] == -5);

        REQUIRE(graph.flow[0][0][1] == 15);
        REQUIRE(graph.flow[1][2][0] == -15);
        
        int w, j;
        for (int v=1; v< graph.n; v++) {
            for (int i=0; i<graph.odeg[v]; i++){
                w = graph.adj[0][v][2*i];
                j = graph.adj[0][v][2*i+1]; // pos in badj
                REQUIRE(graph.flow[0][v][i] == 0);
                REQUIRE(graph.flow[1][w ][ (j-1)/2 ] == 0);
            }
        }
    }

    //cleanup(&graph,&net);

}


TEST_CASE("SETUP: parallel basic net test", "[2proc]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net;
    net=parse(file,rank,size);
    resgraph graph = setup(&net,rank,size);
    SECTION("shared attributes"){
        REQUIRE(net.n==graph.n);
        REQUIRE(net.m==graph.m);
        REQUIRE(net.src==graph.src);
        REQUIRE(net.sink==graph.sink);
        REQUIRE(net.odeg == graph.odeg);
        REQUIRE(net.ideg == graph.ideg);
        for (int v=0; v<graph.npp; v++){
            REQUIRE(net.adj[0][v] == graph.adj[0][v]);
            REQUIRE(net.adj[1][v] == graph.adj[1][v]);
        }
    }
    
    SECTION("excess"){
        int v,w;
        if (rank==0){
            REQUIRE(graph.n_act == 2);
            REQUIRE(graph.ex[0] == 0);
            REQUIRE(graph.ex[1] == 5);
            REQUIRE(graph.ex[2] == 15);

            REQUIRE(graph.active.size() == 2);
            v = graph.active.front();
            graph.active.pop();
            w = graph.active.front();
            REQUIRE( (( (v==1)&&(w==2) )  ||  ( (v==2)&&(w==1) )) );
            
        } else {
            REQUIRE(graph.n_act == 0);
            REQUIRE(graph.ex[0] == 0);
            REQUIRE(graph.ex[1] == 0);
            REQUIRE(graph.ex[2] == 0);
        }
    }

    SECTION("labels"){
        if (rank==0){
            REQUIRE(graph.hght[0]==6);
            REQUIRE(graph.hght[1]==0);
            REQUIRE(graph.hght[2]==0);
        } else {
            for (int v=0; v<graph.npp; v++){
                REQUIRE(graph.hght[v]==0);
            }
        }
    }

    SECTION("flows"){
        if (rank==0){
            REQUIRE(graph.flow[0][0].size() == 2);
            REQUIRE(graph.flow[0][1].size() == 2);
            REQUIRE(graph.flow[0][2].size() == 2);
    
            REQUIRE(graph.flow[1][0].size() == 0);
            REQUIRE(graph.flow[1][1].size() == 1);
            REQUIRE(graph.flow[1][2].size() == 1);
    
            REQUIRE(graph.flow[0][0][0] == 5);
            REQUIRE(graph.flow[1][1][0] == -5);
    
            REQUIRE(graph.flow[0][0][1] == 15);
            REQUIRE(graph.flow[1][2][0] == -15);
    
            REQUIRE(graph.flow[0][1][0]==0);
            REQUIRE(graph.flow[0][1][1]==0);
            REQUIRE(graph.flow[0][1]==graph.flow[0][2]);
    
        } else {
            REQUIRE(graph.flow[0][0].size()==1);
            REQUIRE(graph.flow[0][1].size()==1);
            REQUIRE(graph.flow[0][2].size()==0);
            for (int v=0; v< graph.npp; v++) {
                REQUIRE(graph.flow[1][v].size()==2);
                for (int i=0; i<graph.odeg[v]; i++){
                    REQUIRE(graph.flow[0][v][i] == 0);
                } 
                for (int i=0; i<graph.ideg[v]; i++){
                    REQUIRE(graph.flow[1][v][ i ] == 0);
                }
            }
        }
    }
    //cleanup(&graph,&net);
}

/*
TEST_CASE("1 PULSE: serial basic net test", "[sync]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    resgraph graph = setup(&net,rank,size);
    pulse(&graph);

    SECTION("labels"){
        REQUIRE(graph.hght[0] == 6);
        REQUIRE(graph.hght[1] == 1);
        REQUIRE(graph.hght[2] == 1);
        REQUIRE(graph.hght[3] == 0);
        REQUIRE(graph.hght[4] == 0);
        REQUIRE(graph.hght[5] == 0);
    }

    SECTION("excess"){
        REQUIRE(graph.n_act == 2);
        REQUIRE(graph.ex[0] == 0);
        REQUIRE(graph.ex[1] == 5);
        REQUIRE(graph.ex[2] == 15);
        REQUIRE(graph.ex[3] == 0);
        REQUIRE(graph.ex[4] == 0);
        REQUIRE(graph.ex[5] == 0);
    }

    SECTION("flow[0]s"){
        REQUIRE(graph.flow[0][0] == vector<int>({5,15}) );
        REQUIRE(graph.flow[0][1] == vector<int>({0, 0} ));
        REQUIRE(graph.flow[0][2] == vector<int>({0, 0} ));
        REQUIRE(graph.flow[0][3] == vector<int>({0}) );
        REQUIRE(graph.flow[0][4] == vector<int>({0} ));
        REQUIRE(graph.flow[0][5] == vector<int>({} ));
    }

    SECTION("flow[1]s"){
        REQUIRE(graph.flow[1][0] == vector<int>({}));
        REQUIRE(graph.flow[1][1] == vector<int>({-5}));
        REQUIRE(graph.flow[1][2] == vector<int>({-15}));
        REQUIRE(graph.flow[1][3] == vector<int>({0, 0}));
        REQUIRE(graph.flow[1][4] == vector<int>({0, 0}));
        REQUIRE(graph.flow[1][5] == vector<int>({0, 0}));
    }

    cleanup(&graph,&net);
}



TEST_CASE("2 PULSE: basic net test", "[1proc],[sync]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    resgraph graph = setup(&net,rank,size);
    pulse(&graph);
    pulse(&graph);

    SECTION("labels"){
        REQUIRE(graph.hght[0]==6);
        REQUIRE(graph.hght[1]==1);
        REQUIRE(graph.hght[2]==7);
        REQUIRE(graph.hght[3]==0);
        REQUIRE(graph.hght[4]==0);
        REQUIRE(graph.hght[5]==0);
    }

    SECTION("excess"){
        REQUIRE(graph.ex[0]==0);
        REQUIRE(graph.ex[1]==0);
        REQUIRE(graph.ex[2]==5);
        REQUIRE(graph.ex[3]==10);
        REQUIRE(graph.ex[4]==5);
        REQUIRE(graph.ex[5]==0);
    }
    SECTION("n_act"){
        REQUIRE(graph.n_act == 3);
    }
    SECTION("flow[0]s"){
        REQUIRE(graph.flow[0][0] == vector<int>({5,15}));
        REQUIRE(graph.flow[0][1] == vector<int>({ 5,0 })); 
        REQUIRE(graph.flow[0][2] == vector<int>({ 5,5 })); 
        REQUIRE(graph.flow[0][3] == vector<int>({ 0 })); 
        REQUIRE(graph.flow[0][4] == vector<int>({ 0 })); 
        REQUIRE(graph.flow[0][5] == vector<int>({  })); 
    }
    SECTION("flow[1]s"){
        REQUIRE(graph.flow[1][0] == vector<int>({   }));
        REQUIRE(graph.flow[1][1] == vector<int>({ -5 })); 
        REQUIRE(graph.flow[1][2] == vector<int>({ -15})); 
        REQUIRE(graph.flow[1][3] == vector<int>({ -5,-5 })); 
        REQUIRE(graph.flow[1][4] == vector<int>({ 0,-5 })); 
        REQUIRE(graph.flow[1][5] == vector<int>({ 0,0 })); 
    }

}

TEST_CASE("3 PULSE: basic net test", "[1proc],[sync]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network net = parse(file,rank,size);
    resgraph graph = setup(&net,rank,size);
    pulse(&graph);
    pulse(&graph);
    pulse(&graph);

    SECTION("labels"){
        REQUIRE(graph.hght[0]==6);
        REQUIRE(graph.hght[1]==1);
        REQUIRE(graph.hght[2]==7);
        REQUIRE(graph.hght[3]==1);
        REQUIRE(graph.hght[4]==1);
        REQUIRE(graph.hght[5]==0);
    }

    SECTION("excess"){
        REQUIRE(graph.ex[0]==0);
        REQUIRE(graph.ex[1]==0);
        REQUIRE(graph.ex[2]==0);
        REQUIRE(graph.ex[3]==10);
        REQUIRE(graph.ex[4]==5);
        REQUIRE(graph.ex[5]==0);
    }
    SECTION("n_act"){
        REQUIRE(graph.n_act == 2);
    }
    SECTION("aflows"){
        REQUIRE(graph.aflow[0] == vector<int>({5,10}));
        REQUIRE(graph.aflow[1] == vector<int>({ 5,0 })); 
        REQUIRE(graph.aflow[2] == vector<int>({ 5,5 })); 
        REQUIRE(graph.aflow[3] == vector<int>({ 0 })); 
        REQUIRE(graph.aflow[4] == vector<int>({ 0 })); 
        REQUIRE(graph.aflow[5] == vector<int>({  })); 
    }
    SECTION("bflows"){
        REQUIRE(graph.bflow[0] == vector<int>({   }));
        REQUIRE(graph.bflow[1] == vector<int>({ -5 })); 
        REQUIRE(graph.bflow[2] == vector<int>({ -10})); 
        REQUIRE(graph.bflow[3] == vector<int>({ -5,-5 })); 
        REQUIRE(graph.bflow[4] == vector<int>({ 0,-5 })); 
        REQUIRE(graph.bflow[5] == vector<int>({ 0,0 })); 
    }

}
*/

TEST_CASE("COMPLETE: basic net test", "[2proc]"){
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    string file = ((string) DPREFIX) + ((string) "basic-net.max");
    network ntwk = parse(file,rank,size);
    resgraph net = setup(&ntwk,rank,size);
    comm_data cds = setup_cd(&net,rank,size);
    async_pr(&net, &cds, rank, size);

    SECTION("check forward flows"){
        if (rank==0){
            REQUIRE(net.flow[0][0][0] == 5);
            REQUIRE(net.flow[0][0][1] == 10);
            REQUIRE(net.flow[0][1][0] == 5);
            REQUIRE(net.flow[0][1][1] == 0);
            REQUIRE(net.flow[0][2][0] == 5);
            REQUIRE(net.flow[0][2][1] == 5);
        } else {
            REQUIRE(net.flow[0][0][0] == 10);
            REQUIRE(net.flow[0][1][0] == 5);
            REQUIRE(net.flow[0][2].size() == 0);
        }
    }
    
    SECTION("check backward flows"){
        if (rank==0){
            REQUIRE(net.flow[1][0].size()==0);
            REQUIRE(net.flow[1][1][0] == -5);
            REQUIRE(net.flow[1][2][0] == -10);
        } else {
            REQUIRE(net.flow[1][0][0] == -5);
            REQUIRE(net.flow[1][0][1] == -5);
            REQUIRE(net.flow[1][1][0] == 0);
            REQUIRE(net.flow[1][1][1] == -5);
            REQUIRE(net.flow[1][2][0] == -10);
            REQUIRE(net.flow[1][2][1] == -5);
        }
    }

    SECTION("check anything else"){
        if (rank == 0){
            //REQUIRE(net->ex[0] == 15);
            REQUIRE(net.ex[1] == 0);
            REQUIRE(net.ex[2] == 0);
        } else {
            REQUIRE(net.ex[0] == 0);
            REQUIRE(net.ex[1] == 0);
            REQUIRE(net.ex[2] == 15);
        }
        //REQUIRE(net->n_act == 0);
        REQUIRE(net.active.size()==0);
    }
    cleanup(&net, &ntwk, &cds);
}
