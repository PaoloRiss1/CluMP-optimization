#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>


//**********************************************************************************//
//---------------------------------- VARIABLES -------------------------------------//
//**********************************************************************************//

#define MAX_LINE_LENGTH 1024
#define FNORM (2.3283064365e-10)
#define RANDOM ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)

//Structure for a v node
struct var {
    int ID, deg, *edges, dummySpin, flip, isBoundary, isActive, isVisited;
    double *J, **pm, *m, sumU, mag, field;
} *v;

int use_fixed_seed = 0;
unsigned myrand, ira[256], fixed_seed = 0;
unsigned char ip, ip1, ip2, ip3;

int N, M, sizeS, sizeB, sizeC, edgesB, edgesI;
int MAXFRUS, *graph, *S, *B;
double *coup, R;
FILE *fout;




//**********************************************************************************//
//---------------------------------- FUNCTIONS -------------------------------------//
//**********************************************************************************//

unsigned rand4init(void) {
    unsigned long long y;
    y = (myrand * 16807LL);
    myrand = (y & 0x7fffffff) + (y >> 31);
    if (myrand & 0x80000000)
        myrand = (myrand & 0x7fffffff) + 1;
    return myrand;
}

void initRandom(void) {
    unsigned i;
    ip = 128;
    ip1 = ip - 24;
    ip2 = ip - 55;
    ip3 = ip - 61;
    for (i = ip3; i < ip; i++)
        ira[i] = rand4init();
}




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


//------------------------ Build Cluster ------------------------

#include "ClusterMake.h"




//**********************************************************************************//
//-------------------------------------- MAIN --------------------------------------//
//**********************************************************************************//

int main(int argc, char *argv[]) {

    char *infile, filename[100];

    if (argc < 3) {
        fprintf(stderr, "\n\tUsage: %s MAXFRUS seed[optional] 'graph.conf'\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    MAXFRUS = (int)atoi(argv[1]);
    infile = argv[2];

    // ---- SEED HANDLING ----
    if (argc == 4) {
        use_fixed_seed = 1;
        fixed_seed = (unsigned)atoi(argv[3]);
        myrand = fixed_seed;
    } else {
        use_fixed_seed = 0;
        FILE *devran = fopen("/dev/random", "r");
        fread(&myrand, 4, 1, devran);
    }

    initRandom();
    importGraph(infile);

    sprintf(filename, "N%i_M%i_f%i_cluster.conf", N, M, MAXFRUS);
    fout = fopen(filename, "w");
    fprintf(fout,"#i\tj\tJ_ij\tCluster\tBoundary\n");
    
    //Build & print cluster
    BuildClusterRND(v, MAXFRUS);
    printCluster(v, fout);
    EraseCluster(v);
    
    return 0;
}
