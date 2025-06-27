#include <stdio.h>
#include "ar_c.h"

int fortran_index(int n, int i, int j) {
    // This function converts a 2D index to a 1D index
    // Fortran uses column-major order, so we need to adjust the index accordingly
    return i  +  n * j;

}

void ArC(
    int *id, int *nc, double *n, double *V, double *T,
    double *Ar,
    double *ArV, double *ArT, double *ArTV, double *ArV2, double *ArT2,
    double *Arn, double *ArVn, double *ArTn, double *Arn2
) {
    // This function is just an example and does not perform any operations.
    // In a real implementation, you would include the logic to compute the values
    // based on the input parameters.

    int idx;

    for (int i = 0; i < *nc; i++) {
        *Ar = *Ar + n[i] + *V + *T; // Just an example operation
    }

    // 2D arrays must be interfaced like 1D arrays, we provide the 
    // fortran_index function to ensure correct indexing.
    for (int i = 0; i < *nc; i++) {
        for (int j = 0; j < *nc; j++) {
            idx = fortran_index(*nc, i, j);
            Arn2[idx] = idx;
        }
    }
};

int main() {
    // The id specifier is just used in case the user wants to have a 
    // particular id on the C-side, but it is not used in the Fortran code.
    int const id = 1;

    // The number of components is fixed to 2 for this example.
    int const nc = 2;

    // Critical constants are required for ArModels
    double Tc[nc], Pc[nc], w[nc];

    Tc[0] = 304.2; Tc[1] = 369.8;
    Pc[0] = 73.8; Pc[1] = 45.0;
    w[0] = 0.225; w[1] = 0.185;

    // The ArModel is a function pointer to the model function that will be used.
    set_yaeos_model(id, nc, Tc, Pc, w, *ArC);
    return 0;
}