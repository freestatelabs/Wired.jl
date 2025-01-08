/*  Computational kernel for Wired.jl
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports Float32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Convenience function to create a 2D array of zeros
float** zeros(int Nrows, int Ncols) {

    float** array = calloc(Nrows, sizeof(float*));
    for (int i=0; i<Nrows; i++) {
        array[i] = calloc(Ncols, sizeof(float));
    }

    return array;
}

// De-allocate a 2D array
void free2darray(float** array, int Ncols) {

    for (int i=0; i<Ncols; i++) {
        free(array[i]);
    }

    free(array);
}

// Norm of a 3xN matrix along direction N
void normcols(float* b, float** A, int N) {
    for (int i=0; i<N; i++) {
        b[i] = sqrt(pow(A[0][i],2) +pow(A[1][i],2) +pow(A[2][i],2));
    }
}

// Cross the columns of A with vector b
void crosscols(float** C, float** A, float* b, int N) {

    for (int i=0; i<N; i++) {
        C[0][i] = A[1][i] * b[2] - A[2][i] * b[1];
        C[1][i] = A[2][i] * b[0] - A[0][i] * b[2];
        C[2][i] = A[0][i] * b[1] - A[1][i] * b[0];
    }
}


//
void bfield_wires(float* Bx, float* By, float* Bz, float* x, float* y, float* z, 
            float* a0x, float* a0y, float* a0z, float* a1x, float* a1y, float* a1z, 
            float* I, float* R, int Nn, int Nw, float mu_r)
{

    float d; 
    float a[3];
    float* b = aligned_alloc(32, 32*Nn);
    float* c = aligned_alloc(32, 32*Nn);

    // Outer loop over sources 
    for (int i=0; i<Nw; i++) {

        d = (1e-7) * I[i];
        a[0] = a1x[i] - a0x[i];
        a[1] = a1y[i] - a0y[i];
        a[2] = a1z[i] - a0z[i];

        // Inner loop over nodes 
    }

    free(b);
    free(c);
}

// Solve with optimization level 2: auto-vectorization using the compiler
void solve2(float** bfield, float** nodes, float** sources, int Nn, int Ns) {

    float* a = calloc(Nn, sizeof(float)); 
    float** b = zeros(3,Nn);
    float** c = zeros(3,Nn);
    float** cxa = zeros(3,Nn);
    float* norm_cxa = calloc(Nn, sizeof(float)); 
    float* dot_ac = calloc(Nn, sizeof(float)); 
    float* dot_ab = calloc(Nn, sizeof(float)); 
    float* norm_c = calloc(Nn, sizeof(float)); 
    float* norm_b = calloc(Nn, sizeof(float)); 
    float* s = calloc(Nn, sizeof(float)); 
    float d;

    // Outer loop over sources
    for (int j=0; j<Ns; j++) {

        d = (1e-7) * sources[j][6];
        a[0] = sources[j][3] - sources[j][0];
        a[1] = sources[j][4] - sources[j][1];
        a[2] = sources[j][5] - sources[j][2];

        // Inner loop over nodes
        for (int i=0; i<Nn; i++) {
            b[0][i] = sources[j][0] - nodes[0][i];
            b[1][i] = sources[j][1] - nodes[1][i];
            b[2][i] = sources[j][2] - nodes[2][i];

            c[0][i] = sources[j][3] - nodes[0][i];
            c[1][i] = sources[j][4] - nodes[1][i];
            c[2][i] = sources[j][5] - nodes[2][i];
        }

        crosscols(cxa, c, a, Nn);
        normcols(norm_cxa, cxa, Nn);
        normcols(norm_b, b, Nn);
        normcols(norm_c, c, Nn);

        for (int i=0; i<Nn; i++) {
            dot_ac[i] = a[0]*c[0][i] + a[1]*c[1][i] + a[2]*c[2][i];
        }
            
        for (int i=0; i<Nn; i++) {
            dot_ab[i] = a[0]*b[0][i] + a[1]*b[1][i] + a[2]*b[2][i];
        }

        for (int i=0; i<Nn; i++) {
            s[i] = d * (dot_ac[i]/norm_c[i] - dot_ab[i]/norm_b[i]) * pow(norm_cxa[i],-2);
        }

        for (int i=0; i<Nn; i++) {
            bfield[0][i] += cxa[0][i] * s[i]; 
            bfield[1][i] += cxa[1][i] * s[i]; 
            bfield[2][i] += cxa[2][i] * s[i];
        }
    }    

    free(a); 
    free2darray(b, 3);
    free2darray(c, 3);
    free2darray(cxa, 3);
    free(norm_cxa);
    free(dot_ac);
    free(dot_ab);
    free(norm_b);
    free(norm_c);
    free(s);
}
