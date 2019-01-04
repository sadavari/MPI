#include "z.update.cpp"

void FLUX ( double **U, double **Fr, double **Fz, double *z, double D, int Mrp1, int Mzp1, double dr, double dz, double time, int myID, int nWRs);
void PDE ( double **U, double **Fr, double **Fz, double *r, int Mrp1, int Mzp1, double dr, double dz, double dt);
