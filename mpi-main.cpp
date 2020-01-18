#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <unistd.h>
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
cout << infile << endl;

//wait_for_debugger();
double start = MPI_Wtime();
network net = parse(infile, rank, size);

double parse_t = MPI_Wtime();
resgraph graph = setup(&net, rank, size);
double setup_t = MPI_Wtime();
comm_data cds = setup_cd(&graph, rank, size);
double setup_cd_t = MPI_Wtime();
async_pr(&graph, &cds, rank,size);

double end = MPI_Wtime();
double total_parse = parse_t-start;
double total_setup = setup_t - parse_t;
double total_cd = setup_cd_t - setup_t;
double total_algo = end-setup_cd_t;
double total = end-start;
//while (graph.n_act > 0){
//    pulse(&graph);
//}

//output(graph, outfile, total);

printf("Elapsed time on Processor %d:\n  Parse: %f sec, Setup: %f sec, \n  CD Setup: %f sec, Algo: %f sec\n", rank, total_parse, total_setup, total_cd, total_algo);

print_total(&graph,rank);
//print_flow(&graph, rank);

cleanup(&graph,&net,&cds);


//check(init_str+test+".soln", outfile);


MPI_Finalize();

return 0;

}
