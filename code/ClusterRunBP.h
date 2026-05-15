//Initialize messages in cluster
void initMessages(int *spins){
    int i, j, s;
    struct var *pv;
    
    for(i = 0; i < sizeS; i++){
        pv = v + S[i];
        pv->field = 0.0;
//        pv->field = SMALL * gaussRan();
        for (j = 0; j < pv->deg; j++){
            pv->m[j] = 0.0;
//            pv->m[j] = - 1. + 2. * gaussRan();
            s = pv->edges[j];
            if ( v[s].isBoundary ) {
                if ( fabs(pv->J[j]) == 1. )
                    *(pv->pm[j]) = spins[s] * (1. - FRANDOM * SMALL);
                else
                    *(pv->pm[j]) = spins[s] * LARGE;
//                *(pv->pm[j]) = spins[s] * LARGE;
            }
        }
    }
}


double u_J(double x, double J) {
    if (x * x < J * J) {
        if (J > 0.)
            return x;
        else
            return -x;
    } else {
        if (x > 0.)
            return J;
        else
            return -J;
    }
}

double oneBPstep(int num) {
    int index, i, temp;
    double maxDiff=0.0, sum;
    struct var *pv;
    
    while (num) {
        index = (int)(FRANDOM * num);
        pv = v + S[index];
        
        // Reshuffle S to exclude last extracted
        temp = S[index];
        S[index] = S[--num];
        S[num] = temp;
        
        sum = pv->field;
        for (i = 0; i < pv->deg; i++)
            sum += u_J(*(pv->pm[i]), pv->J[i]);
        
        if( fabs(pv->sumU - sum) > maxDiff )
            maxDiff = fabs(pv->sumU - sum);
        pv->sumU = sum;
        pv->mag = sign(sum);

        for (i = 0; i < pv->deg; i++)
            pv->m[i] = sum - u_J(*pv->pm[i], pv->J[i]);
    }
    return maxDiff;
}

void oneClusterUpdate(int *spins) {
    
    BPiter = 0;
    double maxDiff;

    initMessages(spins);
    do {
        maxDiff = oneBPstep(sizeS);
        BPiter++;
    } while (maxDiff > EPS && BPiter < MAXBPITERS);
}




