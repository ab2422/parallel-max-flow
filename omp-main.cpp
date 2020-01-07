#include <iostream>
#include <omp.h>
#include <vector>
#include "omp-pr-track-active.h"
using namespace std;



int main(int argc, char **argv){

cout << "Max threads: " << omp_get_max_threads()<<endl;

string test(argv[1]);
//string init_str = "/home/annabro/final/data/";
string init_str = "data/";
string infile = init_str+test+".max";
string outfile = init_str+test+"my.soln";

network net = parse(infile);


double start = omp_get_wtime();
resgraph graph;
setup(&net, &graph);
while (graph.n_act > 0){
    pulse(&graph);
}

double end = omp_get_wtime();
double total = end - start;
output(graph, outfile, total);

cout << "Elapsed time: " << total << endl;

//cleanup(graph,net);

check(init_str+test+".soln", outfile);

return 1;

}
