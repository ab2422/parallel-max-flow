#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include "mpi-pr.h"
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

async_pr(&graph, rank,size);

double end = MPI_Wtime();

//while (graph.n_act > 0){
//    pulse(&graph);
//}

//double end = omp_get_wtime();
//double total = end - start;
//output(graph, outfile, total);

//cout << "Elapsed time: " << total << endl;
printf("Elapsed time: %d sec on processor %d\n", end-start, rank);

cleanup(&graph,&net);


//check(init_str+test+".soln", outfile);


MPI_Finalize();


return 0;

}
