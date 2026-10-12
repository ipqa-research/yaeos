#include <stdbool.h>

// How the C implementation should be used
typedef void (*ar_func_t)(
    int nc, 
    const double *n, double v, double t,
    bool calc_Ar,   double *Ar,
    bool calc_ArV,  double *ArV,
    bool calc_ArT,  double *ArT,
    bool calc_ArTV, double *ArTV,
    bool calc_ArV2, double *ArV2,
    bool calc_ArT2, double *ArT2,
    bool calc_Arn,  double *Arn,   // vector nc
    bool calc_ArVn, double *ArVn,  // vector nc
    bool calc_ArTn, double *ArTn,  // vector nc
    bool calc_Arn2, double *Arn2   // matrix nc x nc (column-major!)
);

// Factory exposed to Fortran
extern void make_c_model(ar_func_t f, int nc, void **model_ptr);

// Example
void my_ar(
    int nc, const double *n, double v, double t,
    bool calc_Ar,   double *Ar,
    bool calc_ArV,  double *ArV,
    bool calc_ArT,  double *ArT,
    bool calc_ArTV, double *ArTV,
    bool calc_ArV2, double *ArV2,
    bool calc_ArT2, double *ArT2,
    bool calc_Arn,  double *Arn,
    bool calc_ArVn, double *ArVn,
    bool calc_ArTn, double *ArTn,
    bool calc_Arn2, double *Arn2
) {
    if (calc_Ar)  *Ar  = -1.0;
    if (calc_ArV) *ArV = -2.0;
    // ...
    if (calc_Arn) {
        for (int i = 0; i < nc; i++) Arn[i] = n[i] * (-1.0);
    }
    if (calc_Arn2) {
        // COLUMN-MAJOR: Arn2[i + j*nc]
        for (int i = 0; i < nc; i++)
            for (int j = 0; j < nc; j++)
                Arn2[i + j*nc] = 0.0;
    }
}