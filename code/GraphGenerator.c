#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "myRandom.h"



//**********************************************************************************//
//---------------------------------- VARIABLES -------------------------------------//
//**********************************************************************************//

//Structure for a RRG vertex
struct var {
    int deg, *edges;
    double *J;
} *graph;

char *D, *J;
int L, N, C;
int *GraphEdges;
FILE *fout;


//**********************************************************************************//
//---------------------------------- FUNCTIONS -------------------------------------//
//**********************************************************************************//

//Print errors
void error(char *string) {
    fprintf(stderr, "ERROR! %s\n", string);
    exit(EXIT_FAILURE);
}


//Load the makeRRG() function
#include "makeRRG.h"
#include "makeLRG.h"


// Function to print each pair i, j and their coupling value in the original format
void printGraphEdges(struct var *graph, FILE *output) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < graph[i].deg; j++) {
            fprintf(output,"%5d %5d %14.8f\n", i, graph[i].edges[j],  graph[i].J[j]);
        }
    }
    fprintf(output, "\n");
}

void printGraphNodes(struct var *graph, FILE *output) {
    for (int i = 0; i < N; i++) {
        fprintf(output, "ID:%i\tEDG:(", i);
        
        if (graph[i].deg > 0) {
            for (int j = 0; j < graph[i].deg-1; j++) fprintf(output, "%i, ", graph[i].edges[j]);
            fprintf(output, "%i)\n",graph[i].edges[graph[i].deg-1]);
        }else{
            fprintf(output, "NONE)\n");
        }
    }
    fprintf(output, "\n");
}

void graphINFO(struct var *graph) {
    int supC = 0, eqC = 0, infC = 0, TOTedges = 0;
    
    for (int i = 0; i < N; i++) {
        TOTedges += graph[i].deg;
        if(graph[i].deg == C) eqC++;
        else if(graph[i].deg > C) supC++;
        else if(graph[i].deg < C) infC++;
    }
//    fprintf(stderr, "\n\tGenerated graph with N = %i, and ", N);
//    if ( eqC == N ) fprintf(stderr, "Connectivity EXACTLY C = %i\n\n", C);
//    else {
//        fprintf(stderr, "AVERAGE Connectivity <C> = %.1f with:\n", (double)TOTedges/N);
//        fprintf(stderr, "\t - %i nodes with connectivity = C\n", eqC);
//        fprintf(stderr, "\t - %i nodes with connectivity > C\n", supC);
//        fprintf(stderr, "\t - %i nodes with connectivity < C\n\n", infC);
//    }
}

//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {
    
    char filename[100];
    int M, s1, s2;
    double coupling;
    
    if (argc != 5) {
        fprintf(stderr, "Usage: %s {D: 2/3/RRG} {J: Gaussian/PM1} <N> <C>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
        
    D = argv[1];
    J = argv[2];
    L = (int)atoi(argv[3]); //Number of nodes
//    L = 1 << ((int)atoi(argv[3])); //Number of nodes
    if( strcmp(D, "RRG") ) N = pow(L, (int)atoi(D));
    else N = L;
    C = (int)atoi(argv[4]);
    
    
    M = N * C / 2;
    
    //Initalize variables
    initRandom();
    
    //Allocate all
    GraphEdges = (int *)calloc(C*N, sizeof(int));
    graph = (struct var *)calloc(N, sizeof(struct var));

    
    //Make RRG
    if( !strcmp(D, "RRG") )
        makeRRG(N,C,GraphEdges);
    else
        makeLRG(L,N,(int)atoi(D),C,GraphEdges);

    
    //Initialize graph
    for(int i=0; i<N; i++){
        graph[i].edges = (int *)calloc(C, sizeof(int));
        graph[i].J = (double *)calloc(C, sizeof(double));
        graph[i].deg = 0;
    }
    
    //Assign vertices to each node
    for(int i=0; i<M; i++){
        s1 = GraphEdges[2*i];
        s2 = GraphEdges[2*i+1];
        
        graph[s1].edges[graph[s1].deg++] = s2;
        graph[s2].edges[graph[s2].deg++] = s1;
        
        if( strcmp(J,"Gaussian") == 0 ) coupling = gaussRan();
        else if( strcmp(J,"PM1") == 0 ) coupling = pm1;
        
        graph[s1].J[graph[s1].deg-1] = coupling;
        graph[s2].J[graph[s2].deg-1] = coupling;
    }
    

    //Print graph
    if( !strcmp(D, "RRG") ){
        sprintf(filename, "RRG_N%i_C%i.conf", N, C);
        fout = fopen(filename, "w");
        fprintf(fout, "#Random regular graph with connectivity C = %i\n", C);
    }else{
        sprintf(filename, "LRG%sd_N%i_C%i.conf", D, N, C);
        fout = fopen(filename, "w");
        fprintf(fout, "#Random lattice graph with connectivity C = %i\n", C);
    }
    
    fprintf(fout, "#Vertices:\tBonds:\n");
    fprintf(fout, "#\t%i\t\t%i\n", N, M);
    
    printGraphEdges(graph, fout);
//    printGraphNodes(graph, stderr);
    graphINFO(graph);

    fclose(fout);
    return 0;
}
