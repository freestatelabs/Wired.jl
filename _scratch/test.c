
#include <stdlib.h>
#include <stdio.h>

typedef struct wire {
    float a0[3];
    float a1[3];
    float I; 
    float R;
} Wire;

void testwires(Wire wires[], int Nwires, float a[]){
    // continually overwrite the values of a
    for (int i=0; i<Nwires; i++)
    {
        a[0] = wires[i].a0[0];
        a[1] = wires[i].a0[1];
        a[2] = wires[i].a0[2];
        a[3] = wires[i].a1[0];
        a[4] = wires[i].a1[1];
        a[5] = wires[i].a1[2];
        a[6] = wires[i].I;
        a[7] = wires[i].R;
    }
} 

int main() {
    int Nwires = 10; 
    int size_float = sizeof(float);
    int size_wire = sizeof(Wire);
    float* a = aligned_alloc(size_float, 8*size_float);
    Wire* wires = (Wire *)aligned_alloc(size_wire, Nwires*size_wire);
    //Wire* wires = (Wire *)malloc(Nwires * sizeof(Wire));

    for (int i=0; i<Nwires; i++) {
        wires[i].a0[0] = 1.0;
        wires[i].a0[1] = 2.0;
        wires[i].a0[2] = 3.0;
        wires[i].a1[0] = 4.0;
        wires[i].a1[1] = 5.0;
        wires[i].a1[2] = 6.0;
        wires[i].I = 7.0;
        wires[i].R = 8.0;
    }
    
    testwires(wires, Nwires, a);
    for (int i=0; i<8; i++) {
        printf("a[%i] = %f\n", i, a[i]);
    }
    free(a);
    free(wires);
    return 0;
}


void testwire(Wire* wire, float* a){
    a[0] = wire->a0[0];
    a[1] = wire->a0[1];
    a[2] = wire->a0[2];
    a[3] = wire->a1[0];
    a[4] = wire->a1[1];
    a[5] = wire->a1[2];
    a[6] = wire->I;
    a[7] = wire->R;
} 
void mmult(float* A, float* B, float* C, int N) {

    for (int i=0; i<N; i++) {
        C[i] = A[i]*B[i];
    }
    
}

void vdot(float* a, float* b, float* c, int N) {
    for (int i=0; i<N; i++) {
        c[i] = a[i]*b[i];
    }
}