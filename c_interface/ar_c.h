extern void set_yaeos_model(
    int id,
    int nc,
    double Tc[nc], double Pc[nc], double w[nc],
    void (*)(int *, int *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *, double *)
);

void ArC(
    int *id, int *nc, double *n, double *V, double *T, 
    double *Ar, 
    double *ArV, double *ArT, double *ArTV, double *ArV2, double *ArT2, 
    double *Arn, double *ArVn, double *ArTn, double *Arn2
);