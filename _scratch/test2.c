#include <stdlib.h>
#include <stdio.h>

typedef struct wire {
    float a0[3];
    float a1[3];
    float I; 
    float R;
} Wire;

void testwires(Wire** wires, int Nwires, float* a){
    // continually overwrite the values of a
    for (int i = 0; i < Nwires; i++) {
        a[0] = wires[i]->a0[0];
        a[1] = wires[i]->a0[1];
        a[2] = wires[i]->a0[2];
        a[3] = wires[i]->a1[0];
        a[4] = wires[i]->a1[1];
        a[5] = wires[i]->a1[2];
        a[6] = wires[i]->I;
        a[7] = wires[i]->R;
    }
} 

int main() {
    int Nwires = 10; 
    int size_float = sizeof(float);
    int size_wire = sizeof(Wire); // Correct size for a Wire struct
    float* a = aligned_alloc(size_float, 8 * size_float); // Allocate memory for 8 floats

    // Allocate memory for an array of Wire pointers
    Wire** wires = aligned_alloc(sizeof(Wire*), Nwires * sizeof(Wire*)); 

    float I = 7; 
    float R = 8;
    for (int i = 0; i < Nwires; i++) {
        // Dynamically allocate memory for each Wire struct
        wires[i] = malloc(size_wire); // Allocate memory for a Wire struct
        if (wires[i] == NULL) {
            // Handle allocation failure
            perror("Failed to allocate memory for Wire");
            return 1;
        }
        // Initialize the Wire struct
        for (int j = 0; j < 3; j++) {
            wires[i]->a0[j] = 1.0f + j; // {1, 2, 3}
            wires[i]->a1[j] = 4.0f + j; // {4, 5, 6}
        }
        wires[i]->I = I;
        wires[i]->R = R;
    }
    
    // Call the function to populate the array a
    testwires(wires, Nwires, a);

    // Free allocated memory
    for (int i = 0; i < Nwires; i++) {
        free(wires[i]);
    }
    free(wires);
    free(a);

    return 0;
}
