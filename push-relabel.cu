#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
using namespace std;

/*
* data[5*v : 5*v+4] = [d, d', e, delta e, deg]
* nbhr[3*v : 3*(v+1)-1] = 
*           [w_1, c(v,w_1), f(v,w_1), w_2, c(v,w_2), f(v,w_2), ... ] 
*           for w_i adj to v
*/
__global__
void pulse(int *data, int *nbhr, int n, int mDeg, int *active){
    int v = blockDim.x*blockIdx.x + threadIdx.x; 
    int e = data[5*v+2];
    int d = data[5*v];
    int d_prime = 4*n*n; // an upper bound to # pulses
    int deg = data[5*v+4];
    int i=0;
    int w_ind=0;
    int r_temp=0;
    int delta=0;
    // active means e>0
    if ((v<n) && (e>0) ){
        atomicAdd(active, 1);
        // stage 1 & d' comp for stage 2
        while ((e>0) && (i<deg) ) {
            w_ind = mDeg*3*v+3*i;
            r_temp=nbhr[w_ind+1] - nbhr[w_ind+2];
            if (r_temp>0){
                d_prime = min(d_prime, data[5*nbhr[w_ind]]+1); 
                if ((data[5*nbhr[w_ind]] == d-1) && (r_temp>0) ) {
                    delta = min(e,r_temp);
                    e = e+delta;
                    cudaAtomicAdd(&data[5*nbhr[w_ind]+3], -delta);
                }
            }
            i++;
        }
        // stage 3 prep
        if (e>0) {
            data[5*v+1]=d_prime;
        } else{
            data[5*v+1]=d;
        }
    }
}


__global__
void fin_pulse(int *data, int *nbhr, n, mDeg){
    int v = blockDim.x*blockIdx.x + threadIdx.x; 
    // stage 3
    if (data[5*v+2]>0){
        data[5*v] = data[5*v+1];
    }
    data[5*v+2] += data[5*v+3];
}


/*
Parses a file filename provided in DIMACS netflow format. 
* Stores #verts in *n, #edges in *m.
* nbhr[v] = {w, c(v,w), w', c(v,w'), ... } for nbhrs w
* deg[v] = degree of v
* Stores src ind in *src, sink ind in *sink
* Returns 1 if parse was successful, 0 else
*/
bool parse(const char *filename, int *n, int *m, int **deg, int **nbhr, int *src, int *sink){
    char line[20];
    ifstream file;
    file.open(filename, ios::in);
    int i=0;
    if (file.is_open()) {
        while (file.getline(line)){
            if (line[0]=='p'){
                if (!((line[1]==' ')&&(line[2]=='m')&&(line[3]=='a')&&(line[4]=='x')&&(line[5]==' '))) {
                    return 0;
                } else {
                    i=6;
                    *n = atoi( &(line[6]));
                    while ((alnum(line[i]))){
                        i++;
                    }
                    *m = atoi(&(line[i]));
                }
            } else if (line[0]=='n'){
                i=2;
                while ((alnum(line[i]))){
                    i++;
                }
                if (line[i+1]=='s'){
                    *src = atoi(&(line[2]));
                } else if (line[i+1]=='t'){
                    *sink = atoi(&(line[2]));
                } else {
                    return 0;
                }
            } else if (line[0]=='a'){
                // deal w/ arcs
            } else if (line[0]!='c'){
                return 0;
            }
        }
        file.close();
    }
    return 1;
}

void main(int argc, char **argv){

int n = 100;
int mDeg = n;

int host_data[5*n] = {0};
for (int v=0; v<n; v++){
    host_data[5*v + 4] = 1; //TODO make this deg(v)
}
host_data[5*0+4]=n; //init d(s)=n. assumes s=0, t=n
// TODO init the edges??
int host_nbhr[3*n*mDeg] = {0};


cudaError_t cudaStatus = cudaSetDevice(0);
if (cudaStatus != cudaSuccess){
    cout << "Initialization of device failed" << endl;
}

int *dev_data;
cudaStatus = cudaMalloc(&dev_data, 5*n*sizeof(int));
if (cudaStatus != cudaSuccess){
    cout << "Data malloc failed" << endl;
}
int *dev_nbhr;
cuda Status = cudaMalloc(&dev_nbhr, 3*n*mDeg*sizeof(int));
if (cudaStatus != cudaSuccess){
    cout << "Nbhrs malloc failed" << endl;
}

cudaStatus = cudaMemcpy(dev_data, host_data, 5*n*sizeof(int), cudaMemcpyHostToDevice);
if (cudaStatus != cudaSuccess){
    cout << "Data memcpy failed" << endl;
}
cudaStatus = cudaMemcpy(dev_nbhr, host_nbhr, 3*n*mDeg*sizeof(int), cudaMemcpyHostToDevice);
if (cudaStatus != cudaSuccess){
    cout << "Nbhr memcpy failed" << endl;
}

int num_threads = 16;
int num_blocks = 16;



}
