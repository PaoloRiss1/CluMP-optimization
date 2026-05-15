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
#define NSTEPS 100

//Structure for a RRG vertex
struct var {
    int ID, deg, *edges;
    double *J;
} *v, *bestv;

typedef struct {
    int *spins, iter;
    double E;
} replica;

typedef struct {
    double average;
    double variance;
    double median;
    double percentile_25;
    double percentile_75;
    double percentile_90;
    double percentile_95;
    double min_energy;
    double max_energy;
    double fraction_at_global_min;
    int count_at_global_min;
} ReplicaStats;

struct results {
    int iter, MAXiter;
    double beta, E;
    ReplicaStats stats;
} *optimal;

int N, M, MAXMEAS, MAXITERS, MAXRUNS, MAXREPLICAS;
double betaf, C=0, dEising=0.;
int *graph, *list;
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
    list = (int *)calloc(N, sizeof(int));
    optimal = (struct results *)calloc(MAXRUNS * (MAXMEAS + 1), sizeof(struct results));

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
        for(int j=0; j<v[i].deg; j++) dEising += v[i].J[j];
    }
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

//------------------------ Spin ------------------------

//Random initialization of the spins
void initSpins(int *spins) {
    for(int i = 0; i < N; i++) spins[i] = pm1;
}

//Print spin values
void printSpin(replica *v, FILE *output) {
    fprintf(output, "# ");
    for(int i=0; i<N; i++)
        fprintf(output, "%i ", v->spins[i]);
    fprintf(output, "\n");
}

//------------------- Energy --------------------

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

double dEnergy(int i, int *spins) {
    int s;
    double dE = 0.;
    for(int j=0; j<v[i].deg; j++){
        s = v[i].edges[j];
        dE += v[i].J[j] * spins[i] * spins[s];
    }
    return dE;
}

//------------------------ Replica Operations ------------------------//

void CopyReplica(replica *dest, const replica *src) {
    if (dest->spins != NULL) free(dest->spins);  // Free existing memory
    dest->spins = (int *)calloc(N, sizeof(int));
    memcpy(dest->spins, src->spins, N * sizeof(int));
    dest->E = src->E;
    dest->iter = src->iter;
}

void FreeReplica(replica *r) {
    if (r->spins != NULL) {
        free(r->spins);
        r->spins = NULL;
    }
}

// Resampling function
void ResampleReplicas(double dbeta, replica *population, replica *new_population) {
    int i, j;
    double *cum_weights = (double*)malloc(MAXREPLICAS * sizeof(double));
    double total_weight = 0.0;
    
    for (i = 0; i < MAXREPLICAS; i++) {
        double weight = exp( dbeta * population[i].E );
        total_weight += weight;
        cum_weights[i] = total_weight;
    }
    
    for (i = 0; i < MAXREPLICAS; i++) {
        j = 0;
        double target = FRANDOM * total_weight;
        while (j < MAXREPLICAS - 1 && target > cum_weights[j]) j++;
        CopyReplica(&new_population[i], &population[j]);
    }
    
    for (i = 0; i < MAXREPLICAS; i++)
        CopyReplica(&population[i], &new_population[i]);
    
    for (i = 0; i < MAXREPLICAS; i++) {
        free(new_population[i].spins);
        new_population[i].spins = NULL;
    }
    
    free(cum_weights);
}

//------------------------ Replica Stats ------------------------//

int compare_doubles(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// Function to compute percentile from sorted array
double compute_percentile(double *sorted_energies, int n, double percentile) {
    if (n == 0) return 0.0;
    if (n == 1) return sorted_energies[0];
    
    double index = (percentile / 100.0) * (n - 1);
    int lower_index = (int)floor(index);
    int upper_index = (int)ceil(index);
    
    if (lower_index == upper_index) {
        return sorted_energies[lower_index];
    } else {
        double weight = index - lower_index;
        return sorted_energies[lower_index] * (1.0 - weight) +
               sorted_energies[upper_index] * weight;
    }
}

// Main function to compute replica statistics
void ComputeReplicaStats(replica *population, double global_min_energy, ReplicaStats *stats) {
    
    double *energies = (double*)malloc(MAXREPLICAS * sizeof(double));
    double sum = 0.0;
    double sum_squared = 0.0;
    int count_at_min = 0;
    stats->min_energy = population[0].E;
    stats->max_energy = population[0].E;
    
    for (int i = 0; i < MAXREPLICAS; i++) {
        energies[i] = population[i].E;
        sum += energies[i];
        sum_squared += energies[i] * energies[i];
        
        // Update min/max
        if (energies[i] < stats->min_energy) {
            stats->min_energy = energies[i];
        }
        if (energies[i] > stats->max_energy) {
            stats->max_energy = energies[i];
        }
        
        // Count replicas at global minimum (with small tolerance for floating point)
        if (fabs(energies[i] - global_min_energy) < 1e-10) {
            count_at_min++;
        }
    }
    
    // Compute average
    stats->average = sum / MAXREPLICAS;
    
    // Compute variance
    stats->variance = (sum_squared / MAXREPLICAS) - (stats->average * stats->average);
    
    // Compute fraction at global minimum
    stats->count_at_global_min = count_at_min;
    stats->fraction_at_global_min = (double)count_at_min / MAXREPLICAS;
    
    // Sort energies for percentile calculations
    qsort(energies, MAXREPLICAS, sizeof(double), compare_doubles);
    
    // Compute percentiles
    stats->median = compute_percentile(energies, MAXREPLICAS, 50.0);
    stats->percentile_25 = compute_percentile(energies, MAXREPLICAS, 25.0);
    stats->percentile_75 = compute_percentile(energies, MAXREPLICAS, 75.0);
    stats->percentile_90 = compute_percentile(energies, MAXREPLICAS, 90.0);
    stats->percentile_95 = compute_percentile(energies, MAXREPLICAS, 95.0);
    
    // Clean up
    free(energies);
}

//------------------- Monte Carlo --------------------

//MC step
void MCstep(double beta, int *spins) {
    double dE;
    
    for (int i = 0; i < N; i++) {
        dE = dEnergy(i, spins);
        if (dE <= 0.0 || beta * dE <= - log(FRANDOM)) spins[i] = -spins[i];
    }
}

//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {
    
    char *input, filename[100];
    int run, iter, meas, r, m, currIter, bestIter;
    double E, bestE, minE = 0., beta, beta0, dbeta;
    clock_t start, end;
    
    if (argc != 6) {
        fprintf(stderr, "\nUsage: %s betaf R <log_2(iters)> RUNS 'graph.dat'\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    betaf = (double)atof(argv[1]);
    MAXREPLICAS = (int)atoi(argv[2]);
    MAXMEAS = (int)atoi(argv[3]);
    MAXITERS = 1 << MAXMEAS;
    MAXRUNS = (int)atoi(argv[4]);
    input = argv[5];

    initRandom();
    importGraph(input);
    graphINFO();
    
    //----------------------------- PA ----------------------------------
    
    replica *population = (replica *)calloc(MAXREPLICAS, sizeof(replica));
    replica *new_population = (replica *)calloc(MAXREPLICAS, sizeof(replica));
    replica bestRun, bestReplica;
    replica *bestReplicas = (replica *)calloc(MAXREPLICAS, sizeof(replica));
    
    // Initialize replica structures
    for (r = 0; r < MAXREPLICAS; r++) {
        population[r].spins = NULL;
        new_population[r].spins = NULL;
        bestReplicas[r].spins = NULL;
    }
    bestRun.spins = NULL;
    bestReplica.spins = NULL;
    
    InitGraph(bestv, v);
    run = 0;
    bestRun.E = 1e10;  // Initialize to large value
    
    beta0 = 0.;
    dbeta = (betaf - beta0) / MAXITERS;

    start = clock();
    
    while( run < MAXRUNS ){
        
        // ======= Replicas Init =======
        for (r = 0; r < MAXREPLICAS; r++) {
            
            if (population[r].spins != NULL) free(population[r].spins);
            population[r].spins = (int *)calloc(N, sizeof(int));
            
            initSpins(population[r].spins);
            population[r].E = Energy(population[r].spins);  // Initialize energy
            population[r].iter = 0;
        }
        
        iter = 1;
        m = 1;
        meas = 2;
        bestReplica.E = 1.e10;  // Initialize to large value
        beta = beta0;
        
        while( iter <= MAXITERS ) {
            
            // ======= Annealing Loop =======
            for (r = 0; r < MAXREPLICAS; r++) {
                MCstep(beta, population[r].spins);
                population[r].E = Energy(population[r].spins);
                population[r].iter = iter;
                                
                //Update best values at iter t over all replicas
                if (population[r].E < bestReplica.E){
                    CopyReplica(&bestReplica, &population[r]);
                    
                    if (bestReplica.E < bestRun.E){
                        CopyReplica(&bestRun, &bestReplica);
                     
                        for (int n = 0; n < MAXREPLICAS; n++)
                            CopyReplica(&bestReplicas[n], &population[n]);
                    }
                    ComputeReplicaStats(population, bestReplica.E, &optimal[m + run * MAXMEAS].stats);
                }
            }
            
            //Resampling
            if ( iter % NSTEPS == 0 ) ResampleReplicas( dbeta, population, new_population );
            
            //Update temperature
            beta += dbeta;

            //Update results
            if ( iter % meas == 0 ){
                optimal[m + run * MAXMEAS].E = bestReplica.E;
                optimal[m + run * MAXMEAS].iter = bestReplica.iter;
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

    sprintf(filename, "N%i_M%i_betaf%.2f_R%i_it2^%i.pa", N, M, betaf, MAXREPLICAS, MAXMEAS );
    fout = fopen(filename, "w");
    fprintf(fout, "#maxIters\tbeta\tEminR\tEmaxR\tAvrgR\tVarR\tMedR\tp25R\tp75R\tp90R\tp95R\tbestEfrac\tbestIter\tbestE\n");
    
    for(m = 6; m <= MAXMEAS; m++){
        for(r = 0; r < run; r++){
            ReplicaStats *current_stats = &optimal[m + r * MAXMEAS].stats;

            fprintf(fout, "%8i\t%6.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%10.3f\t%6.3f\t%8i\t%10.6f\n", (int)pow(2,m), optimal[m + r * MAXMEAS].beta, current_stats->min_energy, current_stats->max_energy, current_stats->average, current_stats->variance, current_stats->median, current_stats->percentile_25, current_stats->percentile_75, current_stats->percentile_90, current_stats->percentile_95, current_stats->fraction_at_global_min, optimal[m + r * MAXMEAS].iter, optimal[m + r * MAXMEAS].E);
        }
        fprintf(fout, "\n");
    }
    
    fprintf(fout, "\n#Execution time: %.6f seconds\n", ((double) (end - start)) / CLOCKS_PER_SEC);
    fprintf(fout, "#No. runs: %i\n", run);
    fprintf(fout, "#No. replicas per run: %i\n", MAXREPLICAS);
    fprintf(fout, "#MAXCUT to ISING GS conversion: dE = %.3f\n", 0.5 * dEising);
    fprintf(fout, "#Optimal configuration has E = %.6f\n", bestRun.E);
    fprintf(fout, "#Final spin configuration:\n");
    if (bestRun.spins != NULL) printSpin(&bestRun, fout);
    
    fclose(fout);
    
    //Last iter replicas
    sprintf(filename, "N%i_M%i_betaf%.2f_R%i_it2^%iReplicas.pa", N, M, betaf, MAXREPLICAS, MAXMEAS );
    fout = fopen(filename, "w");
    fprintf(fout, "#Replicas at iteration %i (max iter)\n", MAXITERS);
    fprintf(fout, "#Replica\tE\n\n");
    
    for(r = 0; r < run; r++){
        fprintf(fout, "#run %i\n", r + 1);
        for(int n = 0; n < MAXREPLICAS; n++)
            fprintf(fout, "%8i\t%10.6f\n", n + 1, population[n].E);
        fprintf(fout, "\n");
    }
    fclose(fout);
    
    // Clean up memory
    for (r = 0; r < MAXREPLICAS; r++) {
        FreeReplica(&population[r]);
        FreeReplica(&new_population[r]);
        FreeReplica(&bestReplicas[r]);
    }
    FreeReplica(&bestRun);
    FreeReplica(&bestReplica);
    free(population);
    free(new_population);
    free(bestReplicas);
    
    return 0;
}
