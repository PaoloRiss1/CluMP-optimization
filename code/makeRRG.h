//Copyright (c) 2026 Paolo Rissone and Federico Ricci-Tersenghi


void makeRRG(int N, int C, int * list) {
    int M, i, j, k, l, flag, tmp, i1, i2, e1, e2;
    int *deg, **edge;
    
    if (N * C % 2) error("RRG with both C and N odd do not exist!");
    if (C >= N) error("C must be smaller than N");
    M = (C * N) / 2;
    
    deg = (int*)calloc(N, sizeof(int));
    edge = (int**)calloc(N, sizeof(int*));
    edge[0] = (int*)calloc(C*N, sizeof(int));
    for (i = 1; i < N; i++)
        edge[i] = edge[i-1] + C;
    // fill list
    for (i = 0; i < N; i++)
        for (j = 0; j < C; j++)
            list[i*C+j] = i;
    // randomize list
    for (i = 0; i < C*N; i++) {
        j = (int)(FRANDOM * C*N);
        tmp = list[i];
        list[i] = list[j];
        list[j] = tmp;
    }
    // remove self-links
    do {
        flag = 0;
        for (i = 0; i < M; i++)
            if (list[2*i] == list[2*i+1]) {
                j = (int)(FRANDOM * C*N);
                tmp = list[2*i];
                list[2*i] = list[j];
                list[j] = tmp;
                flag = 1;
            }
    } while(flag);
    // fill edge
    for (i = 0; i < M; i++) {
        i1 = list[2*i];
        i2 = list[2*i+1];
        edge[i1][deg[i1]++] = i;
        edge[i2][deg[i2]++] = i;
    }
    // check deg
    for (i = 0; i < N; i++)
        if (deg[i] != C) error("in deg");
    // check self-links
    for (i = 0; i < M; i++)
        if (list[2*i] == list[2*i+1]) error("self-loop");
    // remove double-links
    do {
        flag = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < C-1; j++)
                for (k = j+1; k < C; k++) {
                    e1 = edge[i][j];
                    e2 = edge[i][k];
                    if ((list[2*e1] == list[2*e2] && list[2*e1+1] == list[2*e2+1]) ||
                        (list[2*e1] == list[2*e2+1] && list[2*e1+1] == list[2*e2])) {
                        do {
                            e2 = (int)(FRANDOM * M);
                        } while (list[2*e1] == list[2*e2] ||
                                 list[2*e1] == list[2*e2+1] ||
                                 list[2*e1+1] == list[2*e2] ||
                                 list[2*e1+1] == list[2*e2+1]);
                        i1 = list[2*e1];
                        i2 = list[2*e2];
                        list[2*e1] = i2;
                        list[2*e2] = i1;
                        l = 0;
                        while (edge[i1][l] != e1) l++;
                        if (l >= C) error("primo l");
                        edge[i1][l] = e2;
                        l = 0;
                        while (edge[i2][l] != e2) l++;
                        if (l >= C) error("secondo l");
                        edge[i2][l] = e1;
                        flag = 1;
                    }
                }
        }
    } while(flag);
    free(edge[0]);
    free(edge);
    free(deg);
}
