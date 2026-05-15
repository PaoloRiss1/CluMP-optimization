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

//Structure for a RRG vertex
struct var {
    int ID, deg, spin, *edges, flips;
    double fitness, *J;
} *v, **tree, *bestv;

struct results {
    int iter;
    double E;
} *optimal;


int N, M, MAXITERS, MAXMEAS, MAXRUNS;
double C=0., B, Gamma = 0.01;
double *cumulant;
int *graph;
double *coup;
FILE *fout, *fout2;


//**********************************************************************************//
//---------------------------------- FUNCTIONS -------------------------------------//
//**********************************************************************************//

//Print errors
void error(char *string) {
    fprintf(stderr, "ERROR! %s\n", string);
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
    cumulant = (double *)calloc(N, sizeof(double));
    tree = (struct var **)calloc(N, sizeof(struct var *));
    
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
        v[i].deg = 0;
        v[i].ID = i;
        v[i].fitness=0.;
        v[i].flips=0;
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
    for(int i=0; i<N; i++) C += v[i].deg;
//    fprintf(stderr, "\nImported graph with N = %d, M = %d, and <C> = %.1f\n", N, M, C/N);
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

void CopyVariables(struct var *dest, const struct var *src) {
    for (int i = 0; i < N; i++){
        dest[i].spin = src[i].spin;
        dest[i].fitness = src[i].fitness;
    }
}

// Function to print each pair i, j and their coupling value in the original format
void printGraphEdges(struct var *v, FILE *output) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < v[i].deg; j++) {
            // Print each pair (i, neighbor) and coupling only if i < neighbor to avoid duplicates
//            if (i < neighbor) {
            fprintf(output,"%5d %5d %14.8f\n", i, v[i].edges[j],  v[i].J[j]);
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
            fprintf(output, "%i)\tF:%f\n",v[i].edges[v[i].deg-1], v[i].fitness);
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
        v[i].flips = 0;
    }
}

void flipSpin(int i) {
    v[i].spin = -v[i].spin;
    v[i].flips++;
}

void printSpin(struct var *v, FILE *output) {
    fprintf(output, "# ");
    for(int i=0; i<N; i++){
        fprintf(output, "%i ", v[i].spin);
    }
    fprintf(output, "\n");
}


//------------------- Fitness & Energy --------------------

//Initialize fitness
void initFitness(void) {
    int s;
    for(int i=0; i<N; i++){
        v[i].fitness = 0.;
        for(int j=0; j<v[i].deg; j++){
            s = v[i].edges[j];
            v[i].fitness += v[i].J[j] * v[s].spin;
        }
        v[i].fitness *= v[i].spin;
    }
}

void updateFitness(int i) {
    int s;
    v[i].fitness -= Gamma * ( v[i].flips - 1 ); //fitness + 'aging' term
    v[i].fitness = -v[i].fitness;
    v[i].fitness += Gamma * v[i].flips;

    for(int j=0; j<v[i].deg; j++){
        s = v[i].edges[j];
        
        if(v[s].spin == v[i].spin)
            v[s].fitness += 2. * v[i].J[j];
        else
            v[s].fitness -= 2. * v[i].J[j];
        

    }
}


//Compute total energy
double InitEnergy(void){
    double E=0.;
    for(int i=0; i<N; i++){
        E -= 0.5 * v[i].fitness;
    }
    return E;
}

double UpdateEnergy(int i, double E) {
    int s;
    for(int j=0; j<v[i].deg; j++){
        s = v[i].edges[j];
        E += 2. * (-v[i].spin) * v[s].spin * v[i].J[j]; //Note the inversion of spin i to update the energy
    }
    return E;
}


//--------------------- Select rank ------------------------

double tau(void){
    double tau = 1. + 1./log(N);
    return tau;
}

void initRank(void){
    
    int k;
    double norm = 0., sum = 0., pdf[N], t;
    
    t = tau();
    
    for(k=1; k<=N; k++){
        pdf[k-1] = 1./pow(k, t);
        norm += pdf[k-1];
    }
    
    for(k=1; k<=N; k++){
        sum += pdf[k-1];
        cumulant[k-1] = sum / norm;
    }
}

int Rank(void) {
    
    double x = FRANDOM;
    for (int k = 0; k < N; k++) {
        if (x < cumulant[k]) {
            return k;
        }
    }
    return N-1;
}


//------------------------- Sorting Tree ------------------------

void initTree(void){
    int i;
    for(i=0; i<N; i++){
        tree[i] = &v[i];
    }
}

void updateTree(void){
    short int n,leaf=0,node=0;
    struct var *help;
    
    for(;;){
        for(n=0; n<2; n++){
            if(++leaf < N){
                if(tree[node]->fitness > tree[leaf]->fitness){
//                    printf("PRE: %d\t%f\t%d\t%f\n", node, tree[node]->fitness, leaf, tree[leaf]->fitness);
                    help = tree[node];
                    tree[node] = tree[leaf];
                    tree[leaf] = help;
//                    printf("SWAP: %d\t%f\t%d\t%f\n", node, tree[node]->fitness, leaf, tree[leaf]->fitness);
                }//end if p
            }else
                return;
        }//end for n
        ++node;
    }//end for()
}

void printTree(FILE *output){
    for(int i=0; i<N; i++){
        fprintf(output, "ID:%i\tS:%i\tEDG:(", tree[i]->ID, tree[i]->spin);
        if (tree[i]->deg > 0) {
            for(int j=0; j<tree[i]->deg-1; j++){
                fprintf(output, "%i, ",tree[i]->edges[j]);
            }
            fprintf(output, "%i)\tF:%f\n",tree[i]->edges[tree[i]->deg-1], tree[i]->fitness);
        }else{
            fprintf(output, "NONE)\n");
        }
    }
    fprintf(output, "\n");
}







//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {
    
    char *input, filename[100];
    int meas, m, iter, run, r, bestIter, k, ID;
    double E, bestE, minE = 0.;
    clock_t start, end;
    
    if (argc != 4) {
        fprintf(stderr, "\n\tUsage: %s <log_2(iters)> RUNS 'graph.dat'\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    MAXMEAS = (int)atoi(argv[1]);
    MAXITERS = 1 << MAXMEAS; //Max number of iterations
    MAXRUNS = (int)atoi(argv[2]);
    input = argv[3];
    
    //Initalize variables
    initRandom();
    importGraph(input);
    
    //Compute average connectivity
    for(int i=0; i<N; i++) C += v[i].deg;

    //----------------------------- EO ----------------------------------
        
    InitGraph(bestv, v);
    
    int OptimalIter[MAXRUNS * (MAXMEAS-6+1)];
    memset(OptimalIter, 0, sizeof(OptimalIter));
    double OptimalE[MAXRUNS * (MAXMEAS-6+1)];
    memset(OptimalE, 0, sizeof(OptimalE));

    run = 0;
    start = clock(); //For execution time
    
    while( run < MAXRUNS ){
        
        //Initialize variables
        initSpins();
        initFitness();
        E = InitEnergy();
        
        //Init ranking
        initTree();
        updateTree();
        initRank();
        
        bestE = 1.e10;
        bestIter = iter;
        iter = 1;
        m = 1;
        meas = 2;
        
        while( iter <= MAXITERS ){
            
            k = Rank();
            
            ID  = tree[k]->ID;
            flipSpin(ID);
            updateFitness(ID);
            E = UpdateEnergy(ID, E);
            
            if (E < bestE){
                bestE = E;
                bestIter = iter;
                
                if (bestE < minE){
                    minE = bestE;
                    CopyVariables(bestv, v);
                }
            }
    
            if ( iter % meas == 0 ){
                optimal[m + run * MAXMEAS].E = bestE;
                optimal[m + run * MAXMEAS].iter = bestIter;
                
                meas *= 2;
                m++;
            }
            updateTree();
            iter++;
        }
        run++;
    }
    end = clock();
           
    
    //---------------------------- Output ----------------------------

    sprintf(filename, "N%i_M%i_it2^%i.eo", N, M, (int)atoi(argv[1]) );
    fout = fopen(filename, "w");
    fprintf(fout, "#maxIters\tbestIter\tbestE\n");
    
    for(m = 6; m <= MAXMEAS; m++){
        for(r = 0; r < run; r++){
            if ( optimal[m + r * MAXMEAS].E != 0. )
            fprintf(fout, "%8i\t%5i\t%10.6f\n", (int)pow(2,m), optimal[m + r * MAXMEAS].iter, optimal[m + r * MAXMEAS].E);
        }
        fprintf(fout, "\n");
    }
    
    fprintf(fout, "\n#Execution time: %.6f seconds\n", ((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(fout, "#No. runs: %i\n", MAXRUNS);
    fprintf(fout, "#Optimal configuration has E = %.6f\n", minE);
    fprintf(fout, "#Final spin configuration:\n");
    printSpin(bestv, fout);
    fclose(fout);
    
    return 0;
}
