#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "myRandom.h"


//**********************************************************************************//
//---------------------------------- VARIABLES -------------------------------------//
//**********************************************************************************//

#define MAX_LINE_LENGTH 1024

#define LARGE 1e4
#define SMALL 1.e-4
#define EPS 1e-8 //BP acceptance threshold

//#define MAXSIZE 1e7 //Max subgraph size
#define MAXBPITERS 1e2 //Max bp iters per subgraph
#define ITERS4FIELD 1e2 //Max BP iters for external field activation

//#define MAXMEAS 16 //Max iters (subgraph)
//#define MAXITERS 1 << MAXMEAS //Max iters (subgraph)
#define WAITINGTIME 1e20 //Max wating time (dt) between consecutive energy drops

//#define MAXRUNS 1e2 //Max number of runs
#define VARIABILITY 1e3 //Max number of consecutive runs with no new enery min

#define NSTEPS 1 //Resampling frequency
//#define mu0 0.05 //Mutation rate


//Structure for a v node
struct var {
    int ID, deg, *edges, dummySpin, flip, isBoundary, isActive, isVisited;
    double *J, **pm, *m, sumU, mag, field;
} *v;

struct results {
    int iter, MAXiter, accepted, failed;
    double E, BPitersAVRG, sizeS, sizeB;
} *optimal;

char *MAXFRUS;
int N, M, MAXMEAS, MAXITERS, MAXRUNS, sizeS, sizeB, BPiter, BPitersAVRG, MAXiters, q, edgesB, edgesI;
int *graph, *S, *B;
double C = 0., beta, dEising = 0., R;
double *coup;
FILE *fout1, *fout2;




//**********************************************************************************//
//---------------------------------- FUNCTIONS -------------------------------------//
//**********************************************************************************//

//-------------------- Generic Functions ----------------------

//Print errors
void error(char *string) {
    fprintf(stderr, "ERROR! %s\n", string);
    exit(EXIT_FAILURE);
}

//Compute the sign of a number
int sign(double x){
    return (x > 0) - (x < 0);
}


//------------------------ Graph Operations ------------------------

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
    optimal = (struct results *)calloc(MAXRUNS * (MAXMEAS + 1), sizeof(struct results));
    
    fgraph = fopen(filename, "r");
    if(!fgraph) error("Cannot open file!\n");
    
    int edge_count = 0;
    while ( fgets(line, sizeof(line), fgraph) ) {
        if (line[0] == '#' || is_blank_line(line) ) continue;
        sscanf(line, "%i %i %lf", &i1, &i2, &coupling);
                
        i1 -= setIDX0;
        i2 -= setIDX0;
                       
        // Add edge if it does not already exist. An exeption is made for BPSP graph, where repetitions are allowed in some cases.
        if ( !edge_exists(i1, i2, graph, edge_count) || strstr(filename, "BPSP") != NULL ) {
            graph[2 * edge_count] = i1;
            graph[2 * edge_count + 1] = i2;
            coup[edge_count] = -coupling; //WARNING: SET TO '-COPULING' FOR GSET
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
        v[i].pm = (double **)calloc(v[i].deg, sizeof(double*));
        v[i].m = (double *)calloc(v[i].deg, sizeof(double));
        v[i].sumU = 0.;
        v[i].deg = 0;
        v[i].isActive = 0;
        v[i].isBoundary = 0;
        v[i].isVisited = 0;
    }
    for (i = 0; i < M; i++) {
        i1 = graph[2*i];
        i2 = graph[2*i+1];
        d1 = v[i1].deg;
        d2 = v[i2].deg;
        v[i1].ID = i1;
        v[i2].ID = i2;
        v[i1].edges[d1] = i2;
        v[i2].edges[d2] = i1;
        v[i1].J[d1] = coup[i];
        v[i2].J[d2] = coup[i];
        v[i1].pm[d1] = &(v[i2].m[d2]);
        v[i2].pm[d2] = &(v[i1].m[d1]);
        v[i1].deg++;
        v[i2].deg++;
    }
}

//Print v info
void ComputeConnectivity(void){
    for(int i=0; i<N; i++) {
        C += v[i].deg;
        for(int j=0; j<v[i].deg; j++) 
            dEising += v[i].J[j];
    }
    C /= N;
}

// Function to print each pair i, j and their coupling value in the original format
void printGraphEdges(struct var *graph, FILE *output) {
    for(int i = 0; i < N; i++){
        for(int j = 0; j < graph[i].deg; j++){
            fprintf(output,"%5d %5d %14.8f\n", i, graph[i].edges[j],  graph[i].J[j]);
        }
    }
    fprintf(output, "\n");
}

//Print graph nodes
void printGraphNodes(int *spins, struct var *graph, FILE *output) {
    for (int i = 0; i < N; i++) {
        fprintf(output, "ID:%i\tS:%i\tEDG:(", i, spins[i]);
        
        if (graph[i].deg > 0) {
            for (int j = 0; j < graph[i].deg-1; j++) fprintf(output, "%i, ", graph[i].edges[j]);
            fprintf(output, "%i)\n",graph[i].edges[graph[i].deg-1]);
        }else{
            fprintf(output, "NONE)\n");
        }
    }
    fprintf(output, "\n");
}



//------------------------ Build Cluster ------------------------

#include "ClusterMake.h"


//----------------------------------- Spins --------------------------------------//

//Random initialization of the spins
void initSpins(int *spins) {
    for(int i = 0; i < N; i++){
        spins[i] = pm1;
        v[i].mag = 0.;
        v[i].flip = 0;
    }
}

//Update cluster spins
void UpdateSpins(int *spins){
    q = 0;
    for(int i=0; i<sizeS; i++){
        int p = S[i];
        if( v[p].mag == 0 ) v[p].mag = pm1;
        if( spins[p] * v[p].mag == - 1 ){
            spins[p] = v[p].mag;
            q++; //Compute overlap (note that if mag = 0, the spin is not assigned).
        }
    }
}

//Random mutatation of the spin configuration
void MutateSpins( int *spins, double mutation_rate ) {
    for (int i = 0; i < sizeS; i++){
        int n = S[i];
        if( FRANDOM < mutation_rate ) spins[n] = pm1;
    }
}

//Print spin values
void printSpin(int *spins, FILE *output) {
    fprintf(output, "# ");
    for(int i=0; i<N; i++)
        fprintf(output, "%i ", spins[i]);
    fprintf(output, "\n");
}


//-------------------------------- Monte Carlo ---------------------------------//

double MCdEnergy(int i, int *spins) {
    int s;
    double dE = 0.;
    for(int j=0; j<v[i].deg; j++){
        s = v[i].edges[j];
        dE += v[i].J[j] * spins[i] * spins[s]; //Note the inversion of spin i to update the energy
    }
    return dE;
}

void MCstepT0(int N, int *spins) {
    double dE;
    
    for (int i = 0; i < N; i++) {
        dE = MCdEnergy(i, spins);
        if( dE <= 0.0 ) spins[i] = -spins[i];
    }
}


//----------------------------------- Energy -------------------------------------//

//Compute total energy
double Energy(int *spins){
    double E = 0.;
    for(int i=0; i<N; i++){
        for(int j=0; j<v[i].deg; j++){
            int s = v[i].edges[j];
            E += v[i].J[j] * spins[i] * spins[s];
        }
    }
    E *= -0.5;
    return E;
}


//------------------------------ Belief Propagation -------------------------------//

#include "ClusterRunBP.h"


//------------------------------ Frustration Control ------------------------------//

int InitFrustration(int *spins, struct var *graph){
    int F0, F = 0;
    if( strcmp(MAXFRUS, "AUTO") != 0 ) F = (int)atoi(MAXFRUS);
    else {
        for(int i=0; i<10; i++){
            F0 = 5;
            do{
                BuildClusterRND(v, F0);
                oneClusterUpdate(spins);
                EraseCluster(v);
                F0 *= 2;
            }
            while( BPiter < MAXBPITERS && R > 0.25 );
            if( F0 > F && R > 0.25 ) F = (int) (0.5 * F0);
        }
    }
    return F;
}

void UpdateFrustration(int *F, int iter, int ACCEPTED){
    if( strcmp(MAXFRUS, "AUTO") == 0 /*&& iter % NSTEPS == 0*/ ){
        if( (double)ACCEPTED/iter > 0.55 && R > 0.25 ) *F += 5;
        else if( *F > 5 && ( (double)ACCEPTED/iter < 0.45 || R <= 0.25 ) ) *F -= 10;
    }
}




//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {
    
    char *infile, filename[100];
    int F, run, dr, m, meas, iter, dt, bestIter, TOTiters = 0, CONVERGED, ACCEPTED, FAILED, totACCEPTED = 0, totFAILED = 0;
    double E, Enext, dE, bestE, minE, AVRGsizeS, AVRGsizeB, dEi[2], dEb[2];
    clock_t start, end;
        
    if (argc != 6) {
        fprintf(stderr, "\n\tUsage: %s MAXFRUS beta <log_2(iters)> RUNS 'graph.conf'\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    MAXFRUS = argv[1];
    beta  = (double)atof(argv[2]);
    MAXMEAS = (int)atoi(argv[3]);
    MAXITERS = 1 << MAXMEAS;
    MAXRUNS = (int)atoi(argv[4]);
    infile = argv[5];

    //Initalize variables
    initRandom();
    importGraph(infile);
    ComputeConnectivity();


    if ( beta >= 1.e3 ){
        sprintf(filename, "N%i_M%i_f%s_betaInf.bp", N, M, MAXFRUS);
        fout1 = fopen(filename, "w");
        
        sprintf(filename, "N%i_M%i_f%s_betaInfEVO.bp", N, M, MAXFRUS);
        fout2 = fopen(filename, "w");
    }else{
        sprintf(filename, "N%i_M%i_f%s_beta%.2f.bp", N, M, MAXFRUS, beta);
        fout1 = fopen(filename, "w");
        
        sprintf(filename, "N%i_M%i_f%s_beta%.2fEVO.bp", N, M, MAXFRUS, beta);
        fout2 = fopen(filename, "w");
    }
    
    
    //----------------------------- BP ----------------------------------
    
    fprintf(fout2, "#%s %6s %6s %6s %8s %10s %10s %6s %16s\n","iter", "F", "itBP", "|S|", "R", "dEi", "dEb", "q", "E");

    run = 0;
    dr = 0;
    minE = 1e20;
    int spins[N], bestSpins[N];
    start = clock(); //For execution time
    
    while( dr < VARIABILITY && run < MAXRUNS ){
        
        initSpins(spins);
        E = Energy(spins);
        
        //Monte Carlo initialization (stops when gets to a min)
        Enext = E;
        do{
            E = Enext;
            MCstepT0(N, spins);
            Enext = Energy(spins);
        }while( Enext < E );
        
        //Running BP
        memcpy(bestSpins, spins, sizeof(bestSpins));

        iter = 1;
        dt = 0;
        m = 1;
        meas = 2; //Minimum number of steps
        bestE = E;
        bestIter = iter;
        BPitersAVRG = 0.;
        CONVERGED = 0;
        AVRGsizeS = 0.;
        AVRGsizeB = 0.;
        ACCEPTED = 0;
        FAILED = 0;
        F = InitFrustration(spins, v);

        fprintf(fout2, "#run %i\n", run + 1);
        
        while( dt < WAITINGTIME && iter <= MAXITERS ){
            
            //Build cluster
            UpdateFrustration(&F,iter,ACCEPTED);
            BuildClusterRND(v, F);
//            printCluster(v);

            AVRGsizeS += sizeS;
            AVRGsizeB += sizeB;
                        
            //Run BP on the cluster
            oneClusterUpdate(spins);
            
            if ( BPiter < MAXBPITERS && R > 0 && sizeS < N ){
                BPitersAVRG += BPiter;
                CONVERGED++;
                EnergyCluster(dEi, dEb, spins);
                dE = dEi[0] + dEb[0];
                
                if ( dE <= 0.0 || beta * dE <= - log(FRANDOM) ) {
                    
                    UpdateSpins(spins);
                    ACCEPTED++;
                    
                    if ( q > 0 ){
                        E = Energy(spins);
                        
                        //Update best param
                        if (E < bestE){ //Update minumum of the run
                            bestE = E;
                            bestIter = iter;
                            dt = 0.;
                            
                            if (bestE < minE){ //Update absolute minimum
                                minE = bestE;
                                memcpy(bestSpins, spins, sizeof(bestSpins));
                                dr = 0;
                            }
                        }
                        fprintf(fout2, "%6i %6i %6i %6i %8.3f %10.3f %10.3f %6i %16.6f\n", iter, F, BPiter, sizeS, R, dEi[0], dEb[0], q, E);
                        fflush(fout2);
                    }
                }
                else FAILED++;
            }
            else FAILED++;
            
//            if ( iter % NSTEPS == 0 ) MutateSpins( spins, mu0 );

            if ( iter % meas == 0 ){
                optimal[m + run * MAXMEAS].E = bestE;
                optimal[m + run * MAXMEAS].iter = bestIter;
                optimal[m + run * MAXMEAS].sizeS = AVRGsizeS / iter;
                optimal[m + run * MAXMEAS].sizeB = AVRGsizeB / iter;
                optimal[m + run * MAXMEAS].BPitersAVRG = (double)BPitersAVRG / CONVERGED;
                optimal[m + run * MAXMEAS].accepted = ACCEPTED;
                optimal[m + run * MAXMEAS].failed = FAILED;
                
                meas *= 2;
                m++;
            }
            
            EraseCluster(v);
            iter++;
            dt++;
        }
        
        fprintf(fout2, "\n\n");
        fflush(fout2);
        
        optimal[run * MAXMEAS].MAXiter = iter - 1;
        TOTiters += iter;
        totACCEPTED += ACCEPTED;
        totFAILED += FAILED;
        
        run++;
        dr++;

    }
    end = clock();

    //Output
    fprintf(fout1, "#Graph with N = %d, M = %d, and <C> = %.1f\n", N, M, C);
    fprintf(fout1, "#%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%10s\n","Iter", "totIters", "accepted", "failed", "<BPiters>", "<|S|>", "<|B|>", "bestIter", "bestE");
    
    for(m = 6; m <= MAXMEAS; m++){
        for(int r = 0; r < run; r++){
            if ( optimal[m + r * MAXMEAS].E != 0. )
                fprintf(fout1, "%8i\t%8i\t%8i\t%8i\t%8.1f\t%8.1f\t%8.1f\t%8i\t%10.4f\n", (int)pow(2,m), optimal[r * MAXMEAS].MAXiter, optimal[m + r * MAXMEAS].accepted, optimal[m + r * MAXMEAS].failed, optimal[m + r * MAXMEAS].BPitersAVRG, optimal[m + r * MAXMEAS].sizeS, optimal[m + r * MAXMEAS].sizeB, optimal[m + r * MAXMEAS].iter, optimal[m + r * MAXMEAS].E);
        }
        fprintf(fout1, "\n");
    }
    
    fprintf(fout1, "\n#Execution time: %.6f seconds\n", ((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(fout1, "#No. runs: %i\n", run);
    fprintf(fout1, "#Average iterations per run: %.1f\n", (double)(TOTiters - 1) / run);
    fprintf(fout1, "#Average accepted moves: %.1f\n", (double)totACCEPTED / run);
    fprintf(fout1, "#Average rejected moves: %.1f\n\n", (double)totFAILED / run);

    fprintf(fout1, "#MAXCUT to ISING GS conversion: dE = %.3f\n", 0.5 * dEising);
    fprintf(fout1, "#Optimal configuration has E = %.6f\n", minE);
    fprintf(fout1, "#Final spin configuration:\n");
    printSpin(bestSpins, fout1);
    
    fclose(fout1);
    fclose(fout2);
    
    return 0;
}
