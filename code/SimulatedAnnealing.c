//Copyright (c) 2026 Paolo Rissone and Federico Ricci-Tersenghi


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "myRandom.h"


//**********************************************************************************//
//---------------------------------- VARIABLES -------------------------------------//
//**********************************************************************************//

#define MAX_LINE_LENGTH 1024
#define NSTEPS 1


//Structure for a RRG vertex
struct var {
    int ID, deg, spin, *edges, flip;
    double *J;
} *v, *bestv;

struct results {
    int iter;
    double beta, E;
} *optimal;


int N, M, MAXMEAS, MAXITERS, MAXRUNS;
double C=0., betaf, dCut=0.;
int *graph, *list;
double *coup;
FILE *fout, *fout2;


//**********************************************************************************//
//---------------------------------- FUNCTIONS -------------------------------------//
//**********************************************************************************//

//Print errors
void error(char *string) {
    fprintf(stderr, "\nERROR! %s\n\n", string);
    exit(EXIT_FAILURE);
}


//------------------------ Import Graph ------------------------

int is_blank_line(const char *line) {
    while (*line) {
        if (*line != ' ' && *line != '\t' && *line != '\n' && *line != '\r') {
            return 0; // Found non-whitespace character
        }
        line++;
    }
    return 1; // Line is blank or contains only whitespace
}

// Function to check if an edge exists
int edge_exists(int i1, int i2, int *graph, int edge_count) {
    for (int i = 0; i < edge_count; i++)
        if ((graph[2 * i] == i1 && graph[2 * i + 1] == i2) ||
            (graph[2 * i] == i2 && graph[2 * i + 1] == i1)) return 1;
    return 0;
}

// Function to read the v data from a file with bi-directional edge addition
void importGraph(const char *filename){
    
    char line[MAX_LINE_LENGTH];
    int current_line = 0, i, i1, i2, d1, d2, setIDX0=0;
    double coupling;
    FILE *fgraph;
    
    fgraph = fopen(filename, "r");
    if(!fgraph) error("Cannot open file!\n");
    
    // Extract N and M from the third line
    while( current_line<4 ){
        fgets(line, sizeof(line), fgraph);
        current_line++;
        if( current_line == 3 ) sscanf(line + 1, "%d %d", &N, &M);
        else if( current_line == 4 ){ //Make sure sites are numbered [0,N-1]
            sscanf(line, "%d %d %lf", &i1, &i2, &coupling);
            if( i1 == 1 ) setIDX0 = 1;
        }
    }
    fclose(fgraph);
        
    // Read lines up until the last one for edge data
    graph = (int *)calloc(2*M, sizeof(int));
    coup = (double *)calloc(M, sizeof(double));
    v = (struct var *)calloc(N, sizeof(struct var));
    bestv = (struct var *)calloc(N, sizeof(struct var));
    list = (int *)calloc(N, sizeof(int));
    
    optimal = (struct results *)calloc(MAXRUNS * MAXMEAS, sizeof(struct results));

    fgraph = fopen(filename, "r");
    if(!fgraph) error("Cannot open file!\n");
    
    int edge_count = 0;
    while ( fgets(line, sizeof(line), fgraph) ) {
        if (line[0] == '#') continue;
        sscanf(line, "%i %i %lf", &i1, &i2, &coupling);
        
        i1 -= setIDX0;
        i2 -= setIDX0;
               
        // Add edge if it does not already exist
        if ( !edge_exists(i1, i2, graph, edge_count) || strstr(filename, "BPSP") != NULL ) {
            graph[2 * edge_count] = i1;
            graph[2 * edge_count + 1] = i2;
            coup[edge_count] = -coupling;
            v[i1].deg++;
            v[i2].deg++;
            edge_count++;
        }
    }
    fclose(fgraph);

    //Definition and population of the v array
    for (i = 0; i < N; i++) {
        v[i].edges = (int *)calloc(v[i].deg, sizeof(int));
        v[i].J = (double *)calloc(v[i].deg, sizeof(double));
        v[i].deg = 0;
        list[i] = i;
    }
    for (i = 0; i < M; i++) {
        i1 = graph[2*i];
        i2 = graph[2*i+1];
        d1 = v[i1].deg;
        d2 = v[i2].deg;
        v[i1].edges[d1] = i2;
        v[i2].edges[d2] = i1;
        v[i1].J[d1] = coup[i];
        v[i2].J[d2] = coup[i];
        v[i1].deg++;
        v[i2].deg++;
    }
}

//Print v info
void graphINFO(void){
    for(int i=0; i<N; i++){
        C += v[i].deg;
        for(int j=0; j<v[i].deg; j++) dCut += v[i].J[j];
    }
//    fprintf(stderr, "\nImported graph with N = %d, M = %d, and <C> = %.1f\n", N, M, C/N);
//    fprintf(stderr, "#MAXCUT to ISING GS conversion: dE = %.3f\n", 0.5 * dCut);
}




//------------------------ Graph Operations ------------------------

//Function to copy all values from src to dest
void InitGraph(struct var *dest, const struct var *src) {
    for (int i = 0; i < N; i++) {
        dest[i].deg = src[i].deg;
        dest[i].edges = (int *)calloc(src[i].deg, sizeof(int));
        dest[i].J = (double *)calloc(src[i].deg, sizeof(double));
                
        memcpy(dest[i].edges, src[i].edges, src[i].deg * sizeof(int));
        memcpy(dest[i].J, src[i].J, src[i].deg * sizeof(double));
    }
}

void CopySpin(struct var *dest, const struct var *src) {
    for (int i = 0; i < N; i++) dest[i].spin = src[i].spin;
}

// Function to print each pair i, j and their coupling value in the original format
void printGraphEdges(struct var *v, FILE *output) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < v[i].deg; j++) {
            // Print each pair (i, neighbor) and coupling only if i < neighbor to avoid duplicates
//            if (i < neighbor) {
            fprintf(output,"%5d %5d %14.8f\n", i, v[i].edges[j],  v[i].J[j]);
//            fprintf(output,"%5d %5d %5d %5d %14.8f\n", i, v[i].spin, v[i].edges[j], v[v[i].edges[j]].spin,  v[i].J[j]);

//            }
        }
    }
    fprintf(output, "\n");
}

void printGraphNodes(struct var *v, FILE *output) {
    for (int i = 0; i < N; i++) {
        fprintf(output, "ID:%i\tS:%i\tEDG:(", v[i].ID, v[i].spin);
        
        if (v[i].deg > 0) {
            for (int j = 0; j < v[i].deg-1; j++) fprintf(output, "%i, ", v[i].edges[j]);
            fprintf(output, "%i)\n",v[i].edges[v[i].deg-1]);
        }else{
            fprintf(output, "NONE)\n");
        }
    }
    fprintf(output, "\n");
}


//------------------------ Spin ------------------------

//Random initialization of the spins
void initSpins(void) {
    for (int i = 0; i < N; i++) {
        v[i].spin = pm1;
        v[i].flip = 0;
    }
}

void printSpin(struct var *v, FILE *output) {
    fprintf(output, "# ");
    for(int i=0; i<N; i++){
        fprintf(output, "%i ", v[i].spin);
    }
    fprintf(output, "\n");
}


//------------------- Energy --------------------

//Compute total energy
double Energy(void){
    double E=0.;
    for(int i=0; i<N; i++){
        for(int j=0; j<v[i].deg; j++){
            int s = v[i].edges[j];
            E += v[i].spin * v[s].spin * v[i].J[j];
        }
    }
    E *= -0.5;
    return E;
}

double dEnergy(int i) {
    int s;
    double dE = 0.;
    for(int j=0; j<v[i].deg; j++){
        s = v[i].edges[j];
        dE += v[i].spin * v[s].spin * v[i].J[j]; //Note the inversion of spin i to update the energy
    }
    return dE;
}


//------------------- Monte Carlo --------------------

//MC step
void MCstep(double beta) {
    double dE;
    
    for (int i = 0; i < N; i++) {
        dE = dEnergy(i);
        if( dE <= 0.0 || beta * dE <= - log(FRANDOM) ) {
            v[i].spin = -v[i].spin;
            v[i].flip++;
        }
    }
}


//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {
    
    char *input, filename[100];
    int iter, run, m, r, meas, bestIter;
    double E, bestE, minE = 0., beta0, dbeta, beta;
    clock_t start, end;
    
    if (argc != 5) {
        fprintf(stderr, "\n\tUsage: %s beta <log_2(iters)> RUNS 'graph.dat'\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    betaf = (double)atof(argv[1]);
    MAXMEAS = (int)atof(argv[2]);
    MAXITERS = 1 << MAXMEAS; //Max number of iterations
    MAXRUNS = (int)atof(argv[3]);
    input = argv[4];
    
    initRandom();
    importGraph(input);
    graphINFO();
    
    
    //----------------------------- SA ----------------------------------
    
    InitGraph(bestv, v);
    beta0 = 0.;
    dbeta = (betaf - beta0) / (MAXITERS / NSTEPS);
    run = 0;
    
    start = clock(); //For execution time
    
    while( run < MAXRUNS ){
        
        //Initialize variables
        initSpins();
        iter = 1;
        m = 1;
        meas = 2;
        bestE = 1.e20;
        beta = beta0;
        
        while( iter <= MAXITERS ) {
            
            MCstep(beta);
            E = Energy();
            
            //Update best values
            if ( E < bestE ){
                bestE = E;
                bestIter = iter;
                
                if (bestE < minE){
                    minE = bestE;
                    CopySpin(bestv, v);
                }
            }
            
            if ( iter % NSTEPS == 0 ) beta += dbeta;
            
            if ( iter % meas == 0 ){
                optimal[m + run * MAXMEAS].E = bestE;
                optimal[m + run * MAXMEAS].iter = bestIter;
                optimal[m + run * MAXMEAS].beta = beta;
                
                meas *= 2;
                m++;
            }
            iter++;
        }
        
        run++;
    }
    end = clock();

    
    //---------------------------- Output ----------------------------
    
    sprintf(filename, "N%i_M%i_betaf%.2f_it2^%i.sa", N, M, betaf, (int)atoi(argv[2]) );
    fout = fopen(filename, "w");
    fprintf(fout, "#maxIters\tbeta\tbestIter\tbestE\n");
    
    for(m = 6; m <= MAXMEAS; m++){
        for(r = 0; r < run; r++){
            if ( optimal[m + r * MAXMEAS].E != 0. )
            fprintf(fout, "%8i\t%8.1f\t%5i\t%10.6f\n", (int)pow(2,m), optimal[m + r * MAXMEAS].beta, optimal[m + r * MAXMEAS].iter, optimal[m + r * MAXMEAS].E);
        }
        fprintf(fout, "\n");
    }
    
    fprintf(fout, "\n#Execution time: %.6f seconds\n", ((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(fout, "#No. runs: %i\n", MAXRUNS);
    fprintf(fout, "#Optimal configuration has E = %.6f\n", minE);
    fprintf(fout, "#MAXCUT to ISING GS conversion: dE = %.3f\n", 0.5 * dCut);
    fprintf(fout, "#Final spin configuration:\n");
    printSpin(bestv, fout);
    fclose(fout);

    return 0;
}
