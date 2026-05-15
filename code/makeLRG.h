void makeLRG(int L, int N, int D, int C, int *list) {
    
    int i, index, neighbor, edge_idx=0, node, count;
    int neighbors[2 * D];
    int *degree, *node_list;
    
    degree = (int*)calloc(N, sizeof(int));
    int *adj_matrix = (int*)calloc(N * N, sizeof(int));
    int *adj_matrixDILU = (int*)calloc(N * N, sizeof(int));


    if (C > 2 * D) error("C must be smaller than 2*D");
    if (N != pow(L, D)) error("N must be L^D for a valid lattice");

    // Step 1: Build full 2D-regular lattice
    for (node = 0; node < N; node++) {
        int coords[D];
        int neighbor_coords[D];
        int count = 0;
        
        // Convert node index to D-dimensional coordinates
        int temp_node = node;
        for (int d = 0; d < D; d++) {
            coords[d] = temp_node % L;
            temp_node /= L;
        }
        
        // Generate 2D nearest neighbors with periodic boundaries
        for (int d = 0; d < D; d++) {
            for (int delta = -1; delta <= 1; delta += 2) {
                for (int k = 0; k < D; k++) neighbor_coords[k] = coords[k];
                neighbor_coords[d] = (coords[d] + delta + L) % L;
                
                int neighbor_index = 0;
                for (int k = D - 1; k >= 0; k--) {
                    neighbor_index = neighbor_index * L + neighbor_coords[k];
                }
                
                neighbor = neighbor_index;
                
                if (!adj_matrix[node * N + neighbor]) {
                    adj_matrix[node * N + neighbor] = 1;
                    adj_matrix[neighbor * N + node] = 1;
                    adj_matrixDILU[node * N + neighbor] = 1;
                    adj_matrixDILU[neighbor * N + node] = 1;
                    degree[node]++;
                    degree[neighbor]++;
                }
            }
        }
    }
    
    // Graph dilution (if C < 2 * D)
    if (C < 2 * D) {

        //Random remove of edges until each node has degree C
        for (node = 0; node < N; node++) {
            if (degree[node] == C) continue;
            
            // Step 1: Collect current neighbors
            int neighbor_list[C], neighbor_count = 0;
            for (int j = 0; j < N; j++) {
                if (adj_matrix[node * N + j]) {
                    neighbor_list[neighbor_count++] = j;
                }
            }

            while (degree[node] > C && neighbor_count > 0) {
                int index = (int)(FRANDOM * neighbor_count);
                int neighbor = neighbor_list[index];
                neighbor_list[index] = neighbor_list[--neighbor_count];
                
                if (degree[neighbor] > C) {
                    adj_matrixDILU[node * N + neighbor] = 0;
                    adj_matrixDILU[neighbor * N + node] = 0;
                    degree[node]--;
                    degree[neighbor]--;
                }
            }
        }
    }

    //Copy edges from adj_matrix to list
    edge_idx = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if (adj_matrixDILU[i * N + j]) {
                list[edge_idx++] = i;
                list[edge_idx++] = j;
                
            }
        }
    }

    free(degree);
    free(adj_matrix);
    free(adj_matrixDILU);
}
