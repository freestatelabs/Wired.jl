/*  Computational kernel for Wired.jl - Wire Sources
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports Float32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>


// Testing @ccall from Julia
void test(float* a, float* b) {
    printf("Hello from C. Your number is %f", a[0]);
    b[0] = a[0];
}

typedef struct {
    float a0[3];
    float a1[3];
    float I; 
    float R;
} Wire;


static inline float mag3(float x, float y, float z){
    return pow(x*x + y*y + z*z, 0.5);
}

static inline float dot3(float a1, float a2, float a3, float b1, float b2, float b3) {
    float res1 = a1*b1;
    float res2 = a2*b2;
    float res3 = a3*b3;
    float result = res1 + res2;
    result += res3;
    return result;
}

int bfield_wires(float* _Bx, float* _By, float* _Bz, const float* x, const float* y, const float* z, 
           const Wire* wires, int Nn, int Nw, float mu_r, int check_inside)
{

    float d; 
    float mag;
    float ac;
    float ab;
    float ac_ab;
    float a[3] = {0};
    float* cx = aligned_alloc(32, 32*Nn);
    float* cy = aligned_alloc(32, 32*Nn);
    float* cz = aligned_alloc(32, 32*Nn);
    float* bx = aligned_alloc(32, 32*Nn);
    float* by = aligned_alloc(32, 32*Nn);
    float* bz = aligned_alloc(32, 32*Nn);
    float* Bx = aligned_alloc(32, 32*Nn);
    float* By = aligned_alloc(32, 32*Nn);
    float* Bz = aligned_alloc(32, 32*Nn);
    float* g = aligned_alloc(32, 32*Nn);
    float* jc = NULL;
    float* r = NULL;
    if (check_inside > 0) {
        // Only allocate if we have to
        jc = aligned_alloc(32, 32*Nn);
        r = aligned_alloc(32, 32*Nn);
    }

    // exit if any of the inputs don't exist
    if (!(x && y && z && wires)) {
        printf("error!\n");
        return 1;
    }

    // Outer loop over sources (wires)
    for (int i=0; i<Nw; i++) {

        d = mu_r * (1e-7) * wires[i].I;
        a[0] = wires[i].a1[0] - wires[i].a0[0];
        a[1] = wires[i].a1[1] - wires[i].a0[1];
        a[2] = wires[i].a1[2] - wires[i].a0[2]; 

        // Calculate the b and c vectors, which point from the node 
        //  to the start and end of the Wire
        for (int j=0; j<Nn; j++) {
            bx[j] = wires[i].a0[0] - x[j];
            cx[j] = wires[i].a1[0] - x[j];
        }
        for (int j=0; j<Nn; j++) {
            by[j] = wires[i].a0[1] - y[j];
            cy[j] = wires[i].a1[1] - y[j];
        }
        for (int j=0; j<Nn; j++) {
            bz[j] = wires[i].a0[2] - z[j];
            cz[j] = wires[i].a1[2] - z[j];
        } 

        // cross rows of c with a, store in B
        // B = cxa
        for (int j=0; j<Nn; j++) {
            Bx[j] = cy[j]*a[2] - cz[j]*a[1];    //   cy*az  -  cz*ay
            By[j] = cz[j]*a[0] - cx[j]*a[2];    // -(cx*az) +  ax*cz
            Bz[j] = cx[j]*a[1] - cy[j]*a[0];    //   cx*ay  -  cy*ax
        }

        // At the same time, calculate the distance from the center of the 
        //  conductor to the node and determine correction factor, but only if
        //  requested.
        // r = |c x a| / |a|
        if (check_inside > 0) {
            mag = 1/mag3(a[0], a[1], a[2]);               // |a|
            for (int j=0; j<Nn; j++) {
                r[j] = mag*mag3(Bx[j], By[j], Bz[j]);     // |c x a| 
            }

            mag = 1/(wires[i].R * wires[i].R);
            for (int j=0; j<Nn; j++) {
                jc[j] = mag*(r[j]*r[j]);
            } 
        }
 
        // Divide B by the magnitude squared and include constant d
        // B = cxa --> B = d*cxa/mag(cxa)^2
        // mag = sqrt(Bx^2 + By^2 + Bz^2)
        // mag^2 = Bx^2 + By^2 + Bz^2
        for (int j=0; j<Nn; j++) {
            // TODO: check if divide by zero
            // assert((pow(Bx[j], 2) + pow(By[j], 2) + pow(Bz[j], 2)) > 1e-8);
            float denom = ((Bx[j]*Bx[j]) + (By[j]*By[j]) +(Bz[j]*Bz[j]) );
            g[j] = d/denom;
        }
        for (int j=0; j<Nn; j++) {     
            Bx[j] *= g[j];
            By[j] *= g[j]; 
            Bz[j] *= g[j];    
        }

        // Apply correction factor if requested
        if (check_inside > 0) {
            for (int j=0; j<Nn; j++) {
                jc[j] = (r[j] < wires[i].R) ? ((r[j] > 0) ? jc[j] : 0.0) : 1.0;
            }
            for (int j=0; j<Nn; j++) {
                Bx[j] *= jc[j];
                By[j] *= jc[j];  
                Bz[j] *= jc[j];    
            }
        } 
        
        // Calculate the dot products and store in B 
        // B *= (a*c/mag(c) - a*b/mag(b))
        for (int j=0; j<Nn; j++) {
            // TODO: check if divide by zero
            ac_ab = dot3(a[0], a[1], a[2], cx[j], cy[j], cz[j]) / mag3(cx[j], cy[j], cz[j]);
            ac_ab -= dot3(a[0], a[1], a[2], bx[j], by[j], bz[j]) / mag3(bx[j], by[j], bz[j]);
            Bx[j] *= ac_ab;
            By[j] *= ac_ab;
            Bz[j] *= ac_ab;
        }

        // copy to output array 
        for (int j=0; j<Nn; j++) {
            _Bx[j] += Bx[j];
            _By[j] += By[j];
            _Bz[j] += Bz[j];
        }  
    }

    free(bx); free(by); free(bz); free(cx); free(cy); free(cz); free(g);

    if (jc) free(jc);
    if (r) free(r);

    return 0;
} 

#define NUMWIRES 1000
#define NUMNODES 1000
#define NUMIT 1000

int main() {
    const int Nn = NUMNODES;
    const int Nw = NUMWIRES;
    const float mu_r = 1.0;
    const int check_inside = 1;
    float Bx[NUMNODES] = {0};
    float By[NUMNODES] = {0};
    float Bz[NUMNODES] = {0};
    float x[NUMNODES] = {0};
    float y[NUMNODES] = {0};
    float z[NUMNODES] = {0};
    float I = 1000 / NUMWIRES;
    float totaltime = 0;

    for (int i=0; i<NUMNODES; i++) {
        x[i] = 1.0;
    }

    Wire wires[NUMWIRES];
    
    for (int i=0; i<NUMWIRES; i++) 
    {
        wires[i] = (Wire) {.a0 = {0, 0, -1e9}, .a1 = {0, 0, 1e9}, .I=I, .R=0.1};
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
        int val = bfield_wires(Bx, By, Bz, x, y, z, wires, Nn, Nw, mu_r, check_inside);
        stop = clock();

        totaltime += (float)(stop - start)/CLOCKS_PER_SEC;


    }

    totaltime /= NUMIT;
    
    printf("Wires:\n");
    for (int i=0; i<10; i++) {
        for (int j=0; j<3; j++) {
            printf("%f ", wires[i].a0[j]);
        }
        for (int j=0; j<3; j++) {
            printf("%f ", wires[i].a1[j]);
        }
        printf("%f %f", wires[i].I, wires[i].R);
        printf("\n");
    }
    printf("result is: \n");
    for (int i=0; i<10; i++) {
        printf("%f %f %f\n", Bx[i], By[i], Bz[i]);
    }

    printf("Elapsed time: %f\n", totaltime);

    return 0;
}