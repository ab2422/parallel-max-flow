#include<stdio.h>

/*  This program reads in points in d-dimensional space and */
/*  produces a sequence-graph, where point i is connected to */
/*  its k sequence neighbors. This is not a nearest-neighbor */
/*  graph because a neigbor is just the next in sequence.     */
/*  Writen by C. McGeoch at DIMACS, July 1991                 */ 


#define MAXDIM 20            /* max for d-space */ 
#define MAXK   10000         /* max number of neighbors */ 

main(argc, argv) 
int argc;
char *argv[];
{
   int k; 
   char str[5];
   int nodes, dimension; 
   int i, d, val; 
   int index, nodecount; 
   int max;
   int sum; 
   char c; 
   char control;
   char buf[50]; 
   int name[MAXK]; 
   int this[MAXDIM]; 
   int point[MAXK][MAXDIM];  
     
  if (argc == 1) k = 2;      /* default is two neighbors each */ 
  else  k = atoi(argv[1]);   /* get vertex degree from command line */ 

  if (k > MAXK) {
    printf("Too many neighbors--recompile with new MAXK\n");
    exit();
  }
  nodecount = 0;
  index = 0; 

  while (scanf("%c  ", &control )  != EOF) {
    if (control == 'v') { /* process a vertex */ 
               nodecount++;

               for (d = 1; d<= dimension; d++) { /*  read in the point */ 
                   scanf( " %d ", &val );
                   this[d] = val; 
		 }

	       if (nodecount > k) /* handle first k separately */ 
                   max = k;   
               else 
		   max = nodecount-1;
	       
	       for (i = 1; i <= max; i++) {
                 sum = 0;
		 for (d = 1; d <= dimension; d++) {
                      /* find manhattan distance */ 
                      sum += abs(this[d] - point[i][d]); 
		    } 
                 printf("e %d  %d  %d \n", nodecount, name[i], sum); 
	       }/* for each neighbor*/

	       /* save the point */ 
               index++; 
               if (index > k) index = 1;
	       for (d=1; d <= dimension; d++) {
		 point[index][d] = this[d];
                 name[index] = nodecount; 
	       }
	     
	     } /* if vertex */

    else if (control == 'p') { /* process problem line */
      scanf("%s  %d  %d", str, &nodes, &dimension);
    }/* else */
    else if (control == 'c') {  /* otherwise ignore a comment line */
      do {c = getchar(); } while ((c != '\n') && (c != EOF));
    }

  }/*while scanf */
 }/*main*/

