#include "catch2/catch.hpp"
#include "omp-pr-track-active.h"
#include <omp.h>
#include <iostream>
#include <vector>
using namespace std;

TEST_CASE("Parser: basic net test"){
    network net = parse("data/basic-net.max");
    SECTION("check node setup"){
        REQUIRE(net.n == 6);
        REQUIRE(net.m == 8);
        REQUIRE(net.src == 0);
        REQUIRE(net.sink == 5);
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
            REQUIRE(3*net.odeg[v]==(net.adj[v]).size());
        }
        REQUIRE(((net.adj[0]==vector<int>({1,5,0, 2,15,0}) ) || (net.adj[0] == vector<int>({2,15,0, 1,5,0}))));
//        REQUIRE(((net.adj[1]==vector<int>({3,5,4, 5}) ) || (net.adj[0] == vector<int>({4, 5,3,5}))));
//        REQUIRE(((net.adj[2]==vector<int>({3,5,4, 5}) ) || (net.adj[0] == vector<int>({4, 5,3,5}))));
        REQUIRE(((net.adj[3][0]==5) && (net.adj[3][1]==15)));
        REQUIRE(((net.adj[4][0]==5) && (net.adj[4][1]==5)));
        REQUIRE(net.adj[5]==vector<int>({}));
    }
    SECTION("check backwards adj matrix"){
        for (int w=0; w<net.n; w++){
            REQUIRE(2*net.ideg[w] == net.badj[w].size() );
        }
        REQUIRE(net.badj[0]==vector<int>({}));
        REQUIRE(((net.badj[1][0]==0)));
        REQUIRE(((net.badj[2][0]==0)));
//        REQUIRE(((net.badj[3]==vector<int>({1,5,2,5})) || (net.badj[3]==vector<int>({2,5,1,5}))));
//        REQUIRE(((net.badj[4]==vector<int>({1,5,2,5})) || (net.badj[3]==vector<int>({2,5,1,5}))));
        REQUIRE(((net.badj[5]==vector<int>({3,0, 4,0}))     ));//|| (net.badj[3]==vector<int>({4,5,1, 3,15,1}))));
    }
    SECTION("check compatibility"){
        for (int v=0; v<net.n; v++){
            for (int i=0; i<net.odeg[v]; i++){
                REQUIRE(v == net.badj[ net.adj[v][3*i] ][ net.adj[v][3*i+2] ]);
            }
            for (int i=0; i<net.ideg[v]; i++){
                REQUIRE(v == net.adj[ net.badj[v][2*i] ][ net.badj[v][2*i+1] ]);
            }
        }
    }
}


TEST_CASE("SETUP: basic net test"){
    network net = parse("data/basic-net.max");
    resgraph graph = setup(&net);
    SECTION("shared attributes"){
        REQUIRE(net.n==graph.n);
        REQUIRE(net.m==graph.m);
        REQUIRE(net.src==graph.src);
        REQUIRE(net.sink==graph.sink);
        REQUIRE(net.odeg == graph.odeg);
        REQUIRE(net.ideg == graph.ideg);
        for (int v=0; v<graph.n; v++){
            REQUIRE(net.adj[v] == graph.adj[v]);
            REQUIRE(net.badj[v] == graph.badj[v]);
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
        REQUIRE(graph.aflow[0][0] == 5);
        REQUIRE(graph.bflow[1][0] == -5);

        REQUIRE(graph.aflow[0][1] == 15);
        REQUIRE(graph.bflow[2][0] == -15);
        
        int w, j;
        for (int v=1; v< graph.n; v++) {
            for (int i=0; i<graph.odeg[v]; i++){
                w = graph.adj[v][3*i];
                j = graph.adj[v][3*i+2]; // pos in badj
                REQUIRE(graph.aflow[v][i] == 0);
                REQUIRE(graph.bflow[w ][ (j-1)/2 ] == 0);
            }
        }
    }

    cleanup(graph,net);

}



TEST_CASE("1 PULSE: basic net test"){
    network net = parse("data/basic-net.max");
    omp_set_num_threads(4);
    resgraph graph = setup(&net);
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

    SECTION("aflows"){
        REQUIRE(graph.aflow[0] == vector<int>({5,15}) );
        REQUIRE(graph.aflow[1] == vector<int>({0, 0} ));
        REQUIRE(graph.aflow[2] == vector<int>({0, 0} ));
        REQUIRE(graph.aflow[3] == vector<int>({0}) );
        REQUIRE(graph.aflow[4] == vector<int>({0} ));
        REQUIRE(graph.aflow[5] == vector<int>({} ));
    }

    SECTION("bflows"){
        REQUIRE(graph.bflow[0] == vector<int>({}));
        REQUIRE(graph.bflow[1] == vector<int>({-5}));
        REQUIRE(graph.bflow[2] == vector<int>({-15}));
        REQUIRE(graph.bflow[3] == vector<int>({0, 0}));
        REQUIRE(graph.bflow[4] == vector<int>({0, 0}));
        REQUIRE(graph.bflow[5] == vector<int>({0, 0}));
    }

    cleanup(graph,net);
}


TEST_CASE("2 PULSE: basic net test"){
    network net = parse("data/basic-net.max");
    omp_set_num_threads(4);
    resgraph graph = setup(&net);
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
    SECTION("aflows"){
        REQUIRE(graph.aflow[0] == vector<int>({5,15}));
        REQUIRE(graph.aflow[1] == vector<int>({ 5,0 })); 
        REQUIRE(graph.aflow[2] == vector<int>({ 5,5 })); 
        REQUIRE(graph.aflow[3] == vector<int>({ 0 })); 
        REQUIRE(graph.aflow[4] == vector<int>({ 0 })); 
        REQUIRE(graph.aflow[5] == vector<int>({  })); 
    }
    SECTION("bflows"){
        REQUIRE(graph.bflow[0] == vector<int>({   }));
        REQUIRE(graph.bflow[1] == vector<int>({ -5 })); 
        REQUIRE(graph.bflow[2] == vector<int>({ -15})); 
        REQUIRE(graph.bflow[3] == vector<int>({ -5,-5 })); 
        REQUIRE(graph.bflow[4] == vector<int>({ 0,-5 })); 
        REQUIRE(graph.bflow[5] == vector<int>({ 0,0 })); 
    }

}

TEST_CASE("3 PULSE: basic net test"){
    network net = parse("data/basic-net.max");
    omp_set_num_threads(4);
    resgraph graph = setup(&net);
    omp_set_num_threads(1);
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
