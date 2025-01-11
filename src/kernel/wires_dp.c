/*  Computational kernel for Wired.jl - Wire Sources
    Runs 2-3x faster than equivalent Julia code when compiled natively

    Notes
    - Supports double32's only
    - Threads are managed by Julia
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

// Testing @ccall from Julia
void test(double* a, double* b) {
    printf("Hello from C. Your number is %f", a[0]);
    b[0] = a[0];
}

// Match the Wire definition in Julia
typedef struct {
    double a0[3];
    double a1[3];
    double I; 
    double R;
} Wire;

// Magnitude of a 3-length vector
// mag(x,y,z) = sqrt(x^2 + y^2 + z^2)
static inline double mag3(double x, double y, double z){
    return pow(x*x + y*y + z*z, 0.5);
}

// Dot product of two 3-length vectors 
// dot3(a,b) = a1*b1 + a2*b2 + a3*b3
static inline double dot3(double a1, double a2, double a3, double b1, double b2, double b3) {
    return a1*b1 + a2*b2 + a3*b3;
}

// Calculate the Bfield generated a sequence of node points (x,y,z) by a series
//   of Wire objects
int bfield_wires(double* Bx, double* By, double* Bz, 
                const double* x, const double* y, const double* z, 
           const Wire* wires, int Nn, int Nw, double mu_r, int check_inside)
{

    double d; 
    double mag;
    double ac;
    double ab;
    double ac_ab;
    double a[3] = {0};
    double* cx = aligned_alloc(32, 32*Nn);
    double* cy = aligned_alloc(32, 32*Nn);
    double* cz = aligned_alloc(32, 32*Nn);
    double* bx = aligned_alloc(32, 32*Nn);
    double* by = aligned_alloc(32, 32*Nn);
    double* bz = aligned_alloc(32, 32*Nn);
    double* _Bx = aligned_alloc(32, 32*Nn);
    double* _By = aligned_alloc(32, 32*Nn);
    double* _Bz = aligned_alloc(32, 32*Nn);
    double* g = aligned_alloc(32, 32*Nn);
    double* jc = NULL;
    double* r = NULL;
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

        // Calculate characteristics of the Wire
        // d = mu_r * mu0 * I / (4pi)
        // a is a vector which points from start to end of Wire
        d = mu_r * (1e-7) * wires[i].I;
        a[0] = wires[i].a1[0] - wires[i].a0[0];
        a[1] = wires[i].a1[1] - wires[i].a0[1];
        a[2] = wires[i].a1[2] - wires[i].a0[2]; 

        // Calculate the b and c vectors, which point from the Node 
        //    to the start and end of the Wire.
        // Group the 3 loops by xj, yj, zj to keep them in cache (?)
        // Appears that it might be slightly effective
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
        // Splitting this into three loops caused it to run slower
        for (int j=0; j<Nn; j++) {
            _Bx[j] = cy[j]*a[2] - cz[j]*a[1];    //   cy*az  -  cz*ay
            _By[j] = cz[j]*a[0] - cx[j]*a[2];    // -(cx*az) +  ax*cz
            _Bz[j] = cx[j]*a[1] - cy[j]*a[0];    //   cx*ay  -  cy*ax
        }

        // At the same time, calculate the distance from the center of the 
        //  conductor to the node and determine correction factor, but only if
        //  requested
        if (check_inside > 0) {
            mag = 1/mag3(a[0], a[1], a[2]); 
            
            for (int j=0; j<Nn; j++) {
                r[j] = mag*mag3(_Bx[j], _By[j], _Bz[j]);     // |c x a| 
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
            // Splitting this up may improve performance
            double denom = ((_Bx[j]*_Bx[j]) + (_By[j]*_By[j]) +(_Bz[j]*_Bz[j]));
            g[j] = d/denom;
        }
        for (int j=0; j<Nn; j++) {     
            _Bx[j] *= g[j];
            _By[j] *= g[j]; 
            _Bz[j] *= g[j];    
        }

        // Apply correction factor if requested
        if (check_inside > 0) {
            for (int j=0; j<Nn; j++) {
                jc[j] = (r[j] < wires[i].R) ? ((r[j] > 0) ? jc[j] : 0.0) : 1.0;
            }
            for (int j=0; j<Nn; j++) {
                _Bx[j] *= jc[j];
                _By[j] *= jc[j];  
                _Bz[j] *= jc[j];    
            }
        } 
        
        // Calculate the dot products and store in B 
        // B *= (a*c/mag(c) - a*b/mag(b))
        for (int j=0; j<Nn; j++) {
            ac_ab = dot3(a[0], a[1], a[2], cx[j], cy[j], cz[j]) / mag3(cx[j], cy[j], cz[j]);
            ac_ab -= dot3(a[0], a[1], a[2], bx[j], by[j], bz[j]) / mag3(bx[j], by[j], bz[j]);
            _Bx[j] *= ac_ab;
            _By[j] *= ac_ab;
            _Bz[j] *= ac_ab;
        }

        // copy to output array 
        for (int j=0; j<Nn; j++) {
            Bx[j] += _Bx[j];
            By[j] += _By[j];
            Bz[j] += _Bz[j];
        }  
    }

    free(bx); free(by); free(bz); free(cx); free(cy); free(cz); free(g);

    if (jc) free(jc);
    if (r) free(r);

    return 0;
} 


// Define a test case for checking the code and for profiling speed
// Expected result: By = 0.0002 T
#define NUMWIRES 1000
#define NUMNODES 1000
#define NUMIT 1000

int main() {
    const int Nn = NUMNODES;
    const int Nw = NUMWIRES;
    const double mu_r = 1.0;
    const int check_inside = 1;
    double Bx[NUMNODES] = {0};
    double By[NUMNODES] = {0};
    double Bz[NUMNODES] = {0};
    double x[NUMNODES] = {0};
    double y[NUMNODES] = {0};
    double z[NUMNODES] = {0};
    double I = 1000 / NUMWIRES;
    double totaltime = 0;

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
        totaltime += (double)(stop - start)/CLOCKS_PER_SEC;

    }

    totaltime /= NUMIT;
    
    printf("result is: \n");
    for (int i=0; i<10; i++) {
        printf("%f %f %f\n", Bx[i], By[i], Bz[i]);
    }

    printf("Elapsed time: %f\n", totaltime);

    return 0;
}