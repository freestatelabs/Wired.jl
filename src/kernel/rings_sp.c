/*  Computational kernel for Wired.jl - Ring Sources
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports float32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#define ITMAX 100 
#define ERRMAX 1e-12
const float pi = M_PI;        // for readability

/*
    Ring 
A circular current-carrying ring with height above XY plane H, major radius R, 
minor radius r, and total current I.
*/
typedef struct {
    float H;
    float R; 
    float r;
    float I;
} Ring;


/*
    float ellipK(float k2) 

Calculate the complete elliptic integral of the first kind (numerically)
*/
static inline double ellipK(double k2) {

    if (fabs(k2 - 1.0) <= ERRMAX) {
        return -1;
    }

    int it = 1; 
    double err = 2*ERRMAX; 
    double a0 = 1.0; 
    double g0 = sqrt(1 - k2);

    double a1;
    double g1;

    while ((it < ITMAX) && (err > ERRMAX)) {
        a1 = 0.5*(a0+g0);
        g1 = sqrt(a0*g0);
        a0 = a1;
        g0 = g1;

        err = fabs(a0 - g0);
        it += 1;
    }

    return pi/(2*a0);
}

/*
    double ellipE(double k2)

Calculate the complete elliptic integral of the second kind (numerically)
*/
static inline double ellipE(double k2) {

    double E = 0; 
    if (fabs(k2 - 1.0) <= ERRMAX) {
        E = 1.0;
    }
    else if (fabs(k2 - 0.0) <= ERRMAX) {
        E = pi/2;
    }
    else {
        double err = 2*ERRMAX;
        int n = 0; 
        double an = 1; 
        double gn = sqrt(1 - k2);
        double cn = fabs(an*an - gn*gn);
        double esumn = cn * pow(2.0, (double)(n-1));
        double an1, gn1, cn1, esumn1;

        while ((n < ITMAX) && (err > ERRMAX)) {
            n += 1;

            an1 = (an + gn) / 2.0;
            gn1 = sqrt(an * gn);
            cn1 = fabs(an1*an1 - gn1*gn1);
            esumn1 = esumn + cn1*(pow(2, (double)(n-1)));

            err = fabs(esumn1 - esumn);

            an = an1; 
            gn = gn1; 
            esumn = esumn1;
        }

        E = (1 - esumn)*pi/(2*an);
    }

    return E;
}

int bfield_rings(float* restrict Bx, float* restrict By, float* restrict Bz, float* restrict x, float* restrict y, float* restrict z, 
                Ring* restrict rings, int Nn, int Nr, float mu_r, int check_inside)
{
    float* rho = aligned_alloc(32, 32*Nn);
    float* rho2 = aligned_alloc(32, 32*Nn);
    float* r2 = aligned_alloc(32, 32*Nn);
    float* alpha2 = aligned_alloc(32, 32*Nn);
    float* beta = aligned_alloc(32, 32*Nn);
    float* beta2 = aligned_alloc(32, 32*Nn);
    float* k2 = aligned_alloc(32, 32*Nn);
    float* K = aligned_alloc(32, 32*Nn);
    float* E = aligned_alloc(32, 32*Nn);
    float* _Bx = aligned_alloc(32, 32*Nn);
    float* _By = aligned_alloc(32, 32*Nn); 
    float* _Bz = aligned_alloc(32, 32*Nn);
    float* jc; 
    if (check_inside > 0) jc = aligned_alloc(32, 32*Nn);
    float C, R, R2, a2;

    // Calculate the node variables first
    for (int j=0; j<Nn; j++) {
        rho2[j] = x[j]*x[j] + y[j]*y[j];
        rho[j] = sqrt(rho2[j]);
    }
    for (int j=0; j<Nn; j++) {
        r2[j] = x[j]*x[j] + y[j]*y[j] + z[j]*z[j];
    }

    for (int i=0; i<Nr; i++) {

        R = rings[i].R;
        R2 = R*R;
        C = mu_r * (4e-7) * rings[i].I;

        // Calculuate alpha, beta, k2, and elliptic integrals now 
        for (int j=0; j<Nn; j++) {
            alpha2[j] = R2 + r2[j] - 2*R*rho[j];
        }
        for (int j=0; j<Nn; j++) {
            beta2[j] = R2 + r2[j] + 2*R*rho[j];     // todo opt based on alpha2?
            beta[j] = sqrt(beta2[j]);
        }
        for (int j=0; j<Nn; j++) {
            k2[j] = 1 - alpha2[j]/beta2[j];
        }
        for (int j=0; j<Nn; j++) {
            K[j] = ellipK(k2[j]);
        }
        for (int j=0; j<Nn; j++) {
            E[j] = ellipE(k2[j]);
        }

        // Now we have everything we need to calculate B
        // Bx first 
        float var1, var2, var3, var4;
        for (int j=0; j<Nn; j++) {
            _Bx[j] = ((C * x[j] * z[j]) / (2*alpha2[j]*beta[j]*rho2[j])) * ((R2 + r2[j]) * E[j] - alpha2[j]*K[j]); 
        }

        for (int j=0; j<Nn; j++)
        {
            _By[j] = (y[j]/x[j]) * _Bx[j];
        }

        for (int j=0; j<Nn; j++) {
            _Bz[j] = (C / (2*alpha2[j]*beta[j])) * ((R2 - r2[j]) * E[j] + alpha2[j]*K[j]); 
        }

        // Current density correction
        if (check_inside > 0) {
            // Calculate current density correction 
            for (int j=0; j<Nn; j++) {
                jc[j] = (alpha2[j] < a2) ? ((alpha2[j] > 0) ? alpha2[j]/a2 : 0.0) : 1.0;
            }

            for (int j=0; j<Nn; j++) {
                _Bx[j] *= jc[j]; 
                _By[j] *= jc[j]; 
                _Bz[j] *= jc[j]; 
            }
        }


        // Copy to output array
        for (int j=0; j<Nn; j++) {
            Bx[j] += _Bx[j];
            By[j] += _By[j];
            Bz[j] += _Bz[j];
        }
    }

    free(rho); free(rho2); free(r2); free(alpha2); free(beta2); free(k2); 
    free(K); free(E); free(_Bx); free(_By); free(_Bz);

    return 0;
}

#define NUMRINGS 1000
#define NUMNODES 1000
#define NUMIT 100

int main() {
    const int Nn = NUMNODES;
    const int Nr = NUMRINGS;
    const float mu_r = 1.0;
    const int check_inside = 1;
    const float H = 0; 
    const float R = 2; 
    const float r = 0.1;

    float Bx[NUMNODES] = {0};
    float By[NUMNODES] = {0};
    float Bz[NUMNODES] = {0};
    float x[NUMNODES] = {0};
    float y[NUMNODES] = {0};
    float z[NUMNODES] = {0};
    float I = 10000 / NUMRINGS;
    float totaltime = 0;

    for (int i=0; i<NUMNODES; i++) {
        x[i] = 0.1;
        y[i] = 0.2;
        z[i] = 0.3;
    }

    Ring rings[NUMRINGS];
    
    for (int i=0; i<NUMRINGS; i++) 
    {
        rings[i] = (Ring) {.H=H, .R=R, .r=r, .I=I};
    }

    clock_t start, stop; 
    clock();
    for (int i=0; i<NUMIT; i++) {

        for (int j=0; j<NUMNODES; j++) {
            Bx[j] = 0;
            By[j] = 0; 
            Bz[j] = 0;
        }

        start = clock();
        int val = bfield_rings(Bx, By, Bz, x, y, z, rings, Nn, Nr, mu_r, check_inside);
        stop = clock();

        totaltime += (float)(stop - start)/CLOCKS_PER_SEC;
    }

    totaltime /= NUMIT;
    
    printf("result is: \n");
    for (int i=0; i<10; i++) {
        printf("%f %f %f\n", Bx[i], By[i], Bz[i]);
    }

    printf("Elapsed time: %f\n", totaltime);

    return 0;

}