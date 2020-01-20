#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <unistd.h>
#include <limits>
#include "mpi-data.h"
#include "mpi-pr.h"
#include "mpi-comm.h"
using namespace std;



int main(int argc, char **argv){

int rank;
int size;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

string test(argv[1]);
//string init_str = "/home/annabro/final/data/";
string init_str = "data/";
//string init_str = "parallel-max-flow/data/";
string infile = init_str+test+".max";
string outfile = init_str+test+"my.soln";
if (rank==0){
    cout << infile << endl;
}
//wait_for_debugger();

network net;
resgraph graph;
comm_data cds;

double prev_time=0, curr_time=0, parse_t=0, setup_t=0, setup_cd_t=0, algo_t;
int num_runs=100;
curr_time=MPI_Wtime();
prev_time = curr_time;

for (int i=0; i<num_runs; i++){
    net = parse(infile, rank, size);
    curr_time = MPI_Wtime();
    parse_t += curr_time - prev_time;
    prev_time = curr_time;

    graph = setup(&net, rank, size);
    curr_time = MPI_Wtime();
    setup_t += curr_time - prev_time;
    prev_time = curr_time;

    cds = setup_cd(&graph, rank, size);
    curr_time = MPI_Wtime();
    setup_cd_t += curr_time - prev_time;
    prev_time = curr_time;

    async_pr(&graph, &cds, rank,size);
    curr_time = MPI_Wtime();
    algo_t += curr_time - prev_time;
    cleanup(&graph,&net,&cds);
    curr_time =MPI_Wtime();
    prev_time = curr_time;
}

//for (int i=0; i<100; i++){
//    graph = setup (&net, rank,size);
//    cds = setup_cd(&graph,rank,size);
//    async_pr(&graph, &cds,rank,size);
//}

double preprocess = parse_t+setup_t+setup_cd_t;
double total = preprocess+ algo_t;

//while (graph.n_act > 0){
//    pulse(&graph);
//}

//output(graph, outfile, total);
double timing[6];
timing[0]=parse_t;
timing[1]=setup_t;
timing[2]=setup_cd_t;
timing[3]=algo_t;
timing[4]=preprocess;
timing[5]=total;

MPI_Allreduce(MPI_IN_PLACE, timing, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
for (int i=0; i<6; i++){
    timing[i] = timing[i] / (size*num_runs);
}

if (rank==0){
    printf("Avg elapsed time over %d runs, with %d processors:\n  Parse: %f sec, Setup: %f sec, CD Setup: %f sec,\n  Algo: %f sec, Preprocess: %f, Total: %f\n", num_runs, size, timing[0], timing[1], timing[2], timing[3], timing[4], timing[5]);
}
//print_total(&graph,rank);
//print_flow(&graph, rank);



//check(init_str+test+".soln", outfile);


MPI_Finalize();

return 0;

}
