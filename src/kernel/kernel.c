/*  Computational kernel for Wired.jl
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports Float32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Testing @ccall from Julia
void test(float* a, float* b) {
    printf("Hello from C. Your number is %f", a[0]);
    b[0] = a[0];

}

static inline float mag3(float x, float y, float z){
    return pow(pow(x,2) + pow(y,2) + pow(z,2), 0.5);
}

static inline float dot3(float a1, float a2, float a3, float b1, float b2, float b3) {
    return a1*b1 + a2*b2 + a3*b3;
}

typedef struct wire {
    float a0[3];
    float a1[3];
    float I; 
    float R;
} Wire;


void bfield_wires(float* Bx, float* By, float* Bz, float* x, float* y, float* z, 
            Wire* wires[], int Nn, int Nw, float mu_r)
{

    float d; 
    float a[3];
    float* cx = aligned_alloc(32, 32*Nn);
    float* cy = aligned_alloc(32, 32*Nn);
    float* cz = aligned_alloc(32, 32*Nn);
    float* bx = aligned_alloc(32, 32*Nn);
    float* by = aligned_alloc(32, 32*Nn);
    float* bz = aligned_alloc(32, 32*Nn);
    float mag;
    float ac;
    float ab;
    float ac_ab;

    // Outer loop over sources 
    for (int i=0; i<Nw; i++) {

        d = mu_r * (1e-7) * wires[i]->I;
        a[0] = wires[i]->a1[0] - wires[i]->a0[0];
        a[1] = wires[i]->a1[1] - wires[i]->a0[1];
        a[2] = wires[i]->a1[2] - wires[i]->a0[2];

        // Calculate the c-vector
        for (int j=0; j<Nn; j++) {
            bx[j] = wires[i]->a0[0] - x[j];
            by[j] = wires[i]->a0[1] - y[j];
            bz[j] = wires[i]->a0[2] - z[j];
            cx[j] = wires[i]->a1[0] - x[j];
            cy[j] = wires[i]->a1[1] - y[j];
            cz[j] = wires[i]->a1[2] - z[j];
        }

        // cross rows of c with a, store in B
        // B = cxa
        for (int j=0; j<Nn; j++) {
            Bx[j] = cy[j]*a[2] - cz[j]*a[1];    // cy*az - cz*ay
            By[j] = cz[j]*a[0] - cx[j]*a[2];    // -(cx*az) + ax*cz
            Bz[j] = cx[j]*a[1] - cy[j]*a[0];    // cx*ay - cy*ax
        }

        // Divide B by the magnitude squared and include constant d
        // B = d*cxa/mag(cxa)^2
        // mag = sqrt(Bx^2 + By^2 + Bz^2)
        // mag^2 = Bx^2 + By^2 + Bz^2
        for (int j=0; j<Nn; j++) {
            mag = d/(pow(Bx[j],2) + pow(By[j],2) + pow(Bz[j],2));
            Bx[j] *= mag;
            By[j] *= mag;    
            Bz[j] *= mag;    
        }

        // Calculate the dot products and store in B 
        // B *= (a*c/mag(c) - a*b/mag(b))
        for (int j=0; j<Nn; j++) {
            ac_ab = dot3(a[0], a[1], a[2], cx[j], cy[j], cz[j]) / mag3(cx[j], cy[j], cz[j]) - dot3(a[0], a[1], a[2], bx[j], by[j], bz[j]) / mag3(bx[j], by[j], bz[j]);
            Bx[j] *= ac_ab;
            By[j] *= ac_ab;
            Bz[j] *= ac_ab;
        }
    }

    free(bx); free(by); free(bz); free(cx); free(cy); free(cz);
}

int main() {
    int Nn = 1000;
    int Nw = 2; 

    float* Bx = aligned_alloc(32, 32*Nn); 
    float* By = aligned_alloc(32, 32*Nn); 
    float* Bz = aligned_alloc(32, 32*Nn); 
    float* x = aligned_alloc(32, 32*Nn); 
    float* y = aligned_alloc(32, 32*Nn); 
    float* z = aligned_alloc(32, 32*Nn); 
    Wire** wires; 
    float mu_r = 1.0; 

    bfield_wires(Bx, By, Bz, x, y, z, wires, Nn, Nw, mu_r);

}
