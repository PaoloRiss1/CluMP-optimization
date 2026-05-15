//Copyright (c) 2026 Paolo Rissone and Federico Ricci-Tersenghi


#include<stdint.h>

#define LCG_A 1664525UL
#define LCG_C 1013904223UL
#define LCG_M 0xFFFFFFFFUL

#define FNORM (2.3283064365e-10)
#define RANDOM getRandomValue()
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

// Improved LCG step function
unsigned rand4init(void) {
    myrand = (LCG_A * myrand + LCG_C) & LCG_M;
    return myrand;
}

// Function to get a high-quality seed
unsigned getHighQualitySeed(void) {
    unsigned seed = 0;
    FILE *entropy_file;
    
    // Try to read from system entropy
    entropy_file = fopen("/dev/urandom", "rb");
    if (entropy_file) {
        if (fread(&seed, sizeof(seed), 1, entropy_file) == 1) {
            fclose(entropy_file);
            return seed;
        }
        fclose(entropy_file);
    }
    
    seed = (unsigned)time(NULL);
    seed ^= (unsigned)clock();
    seed ^= (unsigned)(uintptr_t)&seed;
    seed ^= (unsigned)(uintptr_t)&entropy_file;
    
    // Hash-like mixing
    seed ^= seed >> 16;
    seed *= 2654435761U;
    seed ^= seed >> 16;
    
    return seed ? seed : 123456789U;
}

// Initialization function
void initRandom(void) {
    unsigned i, warmup_rounds;
    
    // Get high-quality seed
    myrand = getHighQualitySeed();
//    myrand = 123346;
    
    // Initialize all array positions
    for (i = 0; i < 256; i++) {
        ira[i] = rand4init();
    }
    
    // Set up indices with proper spacing
    ip = 128;
    ip1 = (ip - 24) & 255;  // Use bitwise AND for wraparound
    ip2 = (ip - 55) & 255;
    ip3 = (ip - 61) & 255;
    
    // Warm up the generator by running it several times
    // This helps eliminate any initial correlations
    warmup_rounds = 1000 + (myrand & 1023); // 1000-2023 warmup iterations
    
    for (i = 0; i < warmup_rounds; i++) {
        // Use the RANDOM macro to warm up
        ira[ip] = ira[ip1] + ira[ip2];
        ira[ip] ^= ira[ip3];
        
        ip = (ip + 1) & 255;
        ip1 = (ip1 + 1) & 255;
        ip2 = (ip2 + 1) & 255;
        ip3 = (ip3 + 1) & 255;
    }
}

// Alternative simpler macro
unsigned getRandomValue(void) {
    unsigned result;
    
    result = (ira[ip1] + ira[ip2]) ^ ira[ip3];
    ira[ip] = result;
    
    ip = (ip + 1) & 255;
    ip1 = (ip1 + 1) & 255;
    ip2 = (ip2 + 1) & 255;
    ip3 = (ip3 + 1) & 255;
    
    return result;
}

// Generate Gaussian distributed random number
double gaussRan(void) {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    
    if (iset == 0) {
        do {
            v1 = 2.0 * FRANDOM - 1.0;
            v2 = 2.0 * FRANDOM - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}
