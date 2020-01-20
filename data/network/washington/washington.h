#ifndef WASHINGTON_NETFLOWS_H
#define WASHINGTON_NETFLOWS_H

#include <stdio.h>

#define FAILURE    0
#define SUCCESS    1
#define FALSE      0
#define TRUE       1


#define MAX_N     20000  
/* #define MAX_N 60 */

#define MAX_CAP   100000000

/* Dimacs problem types */
#define UNDEFINED        0
#define MINCOSTFLOW      1
#define MAXFLOW          2
#define ASSIGNMENT       3


#define MAX_RANGE 1000000
#define MAX_DEGREE 20
#define VERY_BIG 1000000

typedef struct enode {
  struct enode *next;
  struct enode *mate;
  int c;
  int f;
  int h;
  int t;
  int flag;
} Edge;


typedef struct {
  Edge *A[MAX_N];
  int V[MAX_N];
  int size;
  int max_v;
} Graph;

typedef struct {
  int head, tail, size;
  int *data;
} Queue;


Graph* Mesh(int d1, int d2, int r);

Graph* RLevel(int d1, int d2, int r);

Graph *R2Level(int d1, int d2, int r);

Graph *Match(int n, int d);

Graph *SquareMesh(int d, int deg, int r);

Graph *BasicLine(int n,int m,int deg);

Graph *ExponentialLine(int n, int m, int deg);

Graph *DExponentialLine(int n, int m, int deg);

Graph *DinicBadCase(int n);

Graph *GoldBadCase(int n);

Graph *Cheryian(int n, int m, int c);

void Gadget(int a, int b, int n, int m, int c, Graph *G);

void Bridge(int a, int b, int n, Graph *G);

void Sink(int k, Graph *G);

int NewVertex(Graph *G);

void InitGraph(Graph *G);

Graph *CopyGraph(Graph *G1);

void AddVertex(int v, Graph *G);

void AddEdge(int v1, int v2, int a, Graph *G);

Edge *EdgeLookup(int v1, int v2, Graph* G);

void UEdgeArray(Edge *E[], int m, Graph *G);

int EdgeCount(Graph* G);

char *Alloc(int n);

void Barf(char* s);

int EOF_Test(FILE *f);

int SkipLine(FILE *f);

void Skip(FILE *f);

int GetString(FILE *f, char *buff);

void StrAppend(char *s1, char *s2, char *s3);

int GetInt(FILE *f);

void PutInt(int i, FILE *f);

char ReadChar(FILE *f);

int Min(int x, int y);

int Max(int x, int y);

int Abs(int x);

FILE *OpenFile(char *c);

Queue *MakeQueue(int n);

int Dequeue(Queue *Q);

void Enqueue(Queue *Q, int k);

int QSize(Queue *Q);

int QueueEmpty(Queue *Q);

void RandomPermutation(int perm[], int n);

void RandPerm(int perm[], int n);

int RandomInteger(int high, int low);

void InitRandom(int seed);

void RandomSubset(int low, int high, int n, int *x);

Graph *InputFlowGraph(FILE *f, int *s, int *t);

Graph *InputFlow(FILE *f, Graph *G, int *s);

int GraphOutput(Graph *G, FILE *f, int s, int t);

int OutputFlow(Graph *G, FILE *f, int s);

int PrintFlow(Graph *G, int s);

void PrintGraph(Graph *G);

int WriteVertex(int v, Graph *G, FILE *f);

int WriteVertex2(int v, Graph *G, FILE *f);

int WriteVertex3(int v, Graph *G, FILE *f);

#endif
