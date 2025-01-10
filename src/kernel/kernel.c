/*  Computational kernel for Wired.jl
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports Float32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define SMALLEST_FLOAT 1e-8

typedef struct {
    float a0[3];
    float a1[3];
    float I; 
    float R;
} Wire;

// Testing @ccall from Julia
void test(float* a, float* b) {
    printf("Hello from C. Your number is %f", a[0]);
    b[0] = a[0];
}

static float mag3(float x, float y, float z){
    return pow(x*x + y*y + z*z, 0.5);
}

static float dot3(float a1, float a2, float a3, float b1, float b2, float b3) {
    return a1*b1 + a2*b2 + a3*b3;
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
            by[j] = wires[i].a0[1] - y[j];
            bz[j] = wires[i].a0[2] - z[j];
        }
        for (int j=0; j<Nn; j++) {
            cx[j] = wires[i].a1[0] - x[j];
            cy[j] = wires[i].a1[1] - y[j];
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
        //  requested
        if (check_inside > 0) {
            for (int j=0; j<Nn; j++) {
                r[j] = mag3(Bx[j], By[j], Bz[j]);     // |c x a| 
            }
            mag = 1/mag3(a[0], a[1], a[2]);            // 1 / |a|
            for (int j=0; j<Nn; j++) {
                r[j] *= mag;
            }
            mag = pow(wires[i].R, -2);
            for (int j=0; j<Nn; j++) {
                jc[j] = pow(r[j],2);;
            }
            for (int j=0; j<Nn; j++) {
                jc[j] *= mag;
            } 
        }
 
        // Divide B by the magnitude squared and include constant d
        // B = cxa --> B = d*cxa/mag(cxa)^2
        // mag = sqrt(Bx^2 + By^2 + Bz^2)
        // mag^2 = Bx^2 + By^2 + Bz^2
        for (int j=0; j<Nn; j++) {
            // TODO: check if divide by zero
            // assert((pow(Bx[j], 2) + pow(By[j], 2) + pow(Bz[j], 2)) > 1e-8);
            g[j] = d / (pow(Bx[j], 2) + pow(By[j], 2) + pow(Bz[j], 2));
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


int main() {
    const int Nn = 5;
    const int Nw = 5;
    const float mu_r = 1.0;
    const int check_inside = 1;
    float Bx[5] = {0};
    float By[5] = {0};
    float Bz[5] = {0};
    float x[5] = {1e-9}; //{1, 1, 1, 1, 1};
    float y[5] = {0};
    float z[5] = {0};

    Wire wires[5] = {
        {{0, 0, -1e6}, {0, 0, 1e6}, 200, 0.1},
        {{0, 0, -1e6}, {0, 0, 1e6}, 200, 0.1},
        {{0, 0, -1e6}, {0, 0, 1e6}, 200, 0.1},
        {{0, 0, -1e6}, {0, 0, 1e6}, 200, 0.1},
        {{0, 0, -1e6}, {0, 0, 1e6}, 200, 0.1}
    };

    int val = bfield_wires(Bx, By, Bz, x, y, z, wires, Nn, Nw, mu_r, check_inside);


    printf("Wires:\n");
    for (int i=0; i<Nw; i++) {
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
    for (int i=0; i<Nn; i++) {
        printf("%f %f %f\n", Bx[i], By[i], Bz[i]);
    }

    return 0;
}