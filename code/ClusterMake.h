//Find external boundary nodes
void BoundaryNodes(struct var *graph, int *sizeB, int sizeS){
    double edgesC = 0.;
    edgesI = 0, edgesB = 0, *sizeB = 0;
    for(int i=0; i<sizeS; i++){
        int node = S[i];
        for(int j = 0; j < graph[node].deg; j++){
            int neigh = graph[node].edges[j];
            if ( !graph[neigh].isActive ){
                edgesC += 2;
                edgesB++;
                if ( !graph[neigh].isBoundary ){
                    B[(*sizeB)++] = neigh;
                    graph[neigh].isBoundary = 1;
                }
            }
            else edgesC++;
        }
    }
    edgesC *= 0.5;
    edgesI = edgesC - edgesB;
    R = edgesB / edgesC;
}


//
void BuildClusterRND(struct var *graph, int F){
    int sizeScoda = 0, numFrus = 0;
    int Scoda[M];
    
    S = (int *)calloc(N, sizeof(int));    // Cluster membership array
    B = (int *)calloc(N, sizeof(int));    // Array of the boundary nodes

    sizeS = 0;
    sizeB = 0;
    
    // Step 1: Start with a random spin
    int startNode = (int)(FRANDOM * N);
    S[sizeS++] = startNode;                 // Add start node to the cluster
    graph[startNode].dummySpin = pm1;       // Initialize its spin
    graph[startNode].isActive = 1;
    
    // Add neighbors of startNode to Scoda
    for(int i = 0; i < graph[startNode].deg; i++)
        Scoda[sizeScoda++] = graph[startNode].edges[i];

    // Step 2: While there are nodes in Scoda
    while( sizeScoda > 0 && numFrus < F ){
        
        int randIndex = (int)(FRANDOM * sizeScoda);
        int currentNode = Scoda[randIndex];
        Scoda[randIndex] = Scoda[--sizeScoda]; //Remove currentNode from Scoda by replacing it with the last element
        
        if ( graph[currentNode].isActive ) continue; //If candidate node already belongs to S -> skip
        
        int sum = 0, sumSJ = 0;
        double sumJ = 0.0;
        
        for(int i = 0; i < graph[currentNode].deg; i++){
            
            int neighbor = graph[currentNode].edges[i];
            
            if( graph[neighbor].isActive ){
                sum++;
                sumSJ += graph[neighbor].dummySpin * sign(graph[currentNode].J[i]);
                sumJ += graph[neighbor].dummySpin * graph[currentNode].J[i];
            }
        }
        
        if ( numFrus + (sum - abs(sumSJ)) / 2. <= F ){
            
            // Add currentNode to the cluster
            S[sizeS++] = currentNode;
            graph[currentNode].isActive = 1;
            graph[currentNode].dummySpin = sign(sumJ);
            
            // Add neighbors of currentNode to Scoda if they are not already in the cluster
            for(int i = 0; i < graph[currentNode].deg; i++){
                int neighbor = graph[currentNode].edges[i];
                if( !graph[neighbor].isActive )
                    Scoda[sizeScoda++] = neighbor;
            }
            
            numFrus += (sum - abs(sumSJ)) / 2;
        }
    }
    BoundaryNodes(graph, &sizeB, sizeS);
}



//Build the most compact graph possible (pseudo-deterministic with FIFO neighbors extraction)
void BuildClusterFIFO(struct var *graph, int F) {
    int i, iRead = 0, iWrite = 0, offset, idx, numFrus = 0;
    struct var *Scoda[M], *ps, *pv, *neighbor;

    // Allocate memory
    S = (int *)calloc(N, sizeof(int));
    B = (int *)calloc(N, sizeof(int));

    sizeS = 0;
    sizeB = 0;

    // Step 1: Start with a random spin
    int startIndex = (int)(FRANDOM * N);
    struct var *startNode = &graph[startIndex];

    S[sizeS++] = startIndex;
    startNode->dummySpin = pm1;
    startNode->isActive = 1;

    // Add neighbors of startNode to Scoda (FIFO queue)
    offset = (int)(FRANDOM * startNode->deg);
    for (i = 0; i < startNode->deg; i++) {
        idx = (offset + i) % startNode->deg;
        neighbor = &graph[startNode->edges[i]];
        Scoda[iWrite++] = neighbor;
    }

    // Step 2: Cluster growth
    while (iRead < iWrite /*&& sizeS < MAXSIZE*/) {
        ps = Scoda[iRead++];  // Pop first node

        if (ps->isActive) continue;  // Already in cluster → skip

        int sum = 0, sumSJ = 0;
        double sumJ = 0.0;

        for (i = 0; i < ps->deg; i++) {
            pv = &graph[ps->edges[i]];
            if (pv->isActive) {
                sum++;
                sumSJ += pv->dummySpin * sign(ps->J[i]);
                sumJ += pv->dummySpin * ps->J[i];
//                corr -= graph[ps->ID].dummySpin * graph[i].dummySpin * graph[ps->ID].J[i];
            }
        }
//        corr /= sum;

        // Case 1: No frustration
        if ( numFrus + (sum - abs(sumSJ)) / 2. <= F /*&& FRANDOM < fmin(1., fmax(0.,corr))*/ ) {
            ps->isActive = 1;
            ps->dummySpin = sign(sumJ);
            S[sizeS++] = ps - graph;  // Convert pointer back to index

            offset = (int)(FRANDOM * ps->deg);
            for (i = 0; i < ps->deg; i++) {
                idx = (offset + i) % ps->deg;
                pv = &graph[ps->edges[idx]];
                if (!pv->isActive)
                    Scoda[iWrite++] = pv;
            }
            
            numFrus += (sum - abs(sumSJ)) / 2;
        }
    }

    // Final step: build boundary
    BoundaryNodes(graph, &sizeB, sizeS);
}


//Compute cluster energy
void EnergyCluster(double *dEi, double *dEb, int *spins){
    
    int p, s;
    double Ei0 = 0., Eb0 = 0., Eis = 0., Ebs = 0., Eit = 0., Ebt = 0.;

    for(int i = 0; i < sizeS; i++){
        p = S[i];
        for (int j = 0; j < v[p].deg; j++){
            s = v[p].edges[j];
            if( v[s].isBoundary ) {
                Eb0 += 2. * v[p].J[j] * spins[p] * spins[s];
                Ebs += 2. * v[p].J[j] * v[p].mag * spins[s];
                Ebt += 2. * v[p].J[j] * v[p].dummySpin * v[s].dummySpin;
            }
            else if( v[s].isActive ) {
                Ei0 += v[p].J[j] * spins[p] * spins[s];
                Eis += v[p].J[j] * v[p].mag * v[s].mag;
                Eit += v[p].J[j] * v[p].dummySpin * v[s].dummySpin;
            }
        }
    }
    dEi[0] = - 0.5 * (Eis - Ei0);
    dEb[0] = - 0.5 * (Ebs - Eb0);
    dEi[1] = - 0.5 * (Eit - Ei0);
    dEb[1] = - 0.5 * (Ebt - Eb0);
}


//Erase cluster
void EraseCluster(struct var *graph){
    for(int i=0; i<N; i++){
        graph[i].isBoundary = 0;
        graph[i].isActive = 0;
        graph[i].isVisited = 0;
    }
    free(S);
    S = NULL;
    free(B);
    B = NULL;
}


//Print cluster configuration
void printCluster(struct var *graph, FILE *output){
    struct var *pv;
    for(int i=0; i<N; i++){
        pv = graph + i;
        for(int j=0; j<pv->deg; j++){
            fprintf(output, "%5d\t%5d\t%8.5f\t%2d\t%2d\n", pv->ID, pv->edges[j], pv->J[j], pv->isActive, pv->isBoundary);
        }
    }
}
