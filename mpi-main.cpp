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
network net = parse(infile, rank, size);

double start = MPI_Wtime();
resgraph graph = setup(&net, rank, size);
comm_data cds = setup_cd(&graph, rank, size);
async_pr(&graph, &cds, rank,size);

double end = MPI_Wtime();
double total = end-start;

//while (graph.n_act > 0){
//    pulse(&graph);
//}

//output(graph, outfile, total);

printf("Elapsed time: %f sec on processor %d\n", end, rank);

print_flow(&graph, rank);

cleanup(&graph,&net);


//check(init_str+test+".soln", outfile);


MPI_Finalize();

return 0;

}
