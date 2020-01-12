#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include "mpi-pr.h"
using namespace std;



int main(int argc, char **argv){

printf("starting main\n");

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
printf("about to go into parse, p%d\n", rank);
network net = parse(infile, rank, size);
printf("done parse, about to setup, p%d\n",rank);

//double start = omp_get_wtime();
resgraph graph = setup(&net, rank, size);
printf("done setup, p%d\n", rank);
//while (graph.n_act > 0){
//    pulse(&graph);
//}

//double end = omp_get_wtime();
//double total = end - start;
//output(graph, outfile, total);

//cout << "Elapsed time: " << total << endl;


cleanup(&graph,&net);


//check(init_str+test+".soln", outfile);

printf("about to finalize, p%d\n",rank);

MPI_Finalize();

printf("Totally done!\n");

return 0;

}
