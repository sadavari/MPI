#include <math.h>

void FLUX ( double **U, double **Fr, double **Fz, double *z, double D, int Mrp1, int Mzp1, double dr, double dz, double time, int myID, int nWRs)
{
  int i, j;
  for ( i=0; i<=Mrp1; i++)
  {
    if ( myID == 1)
    {
      U[0][i] = 0.0;
      Fz[1][i] = -(2.0*D/dz)*(U[1][i]-U[0][i]);	
    }
    else
    {
      Fz[1][i] = -(D/dz)*(U[1][i]-U[0][i]);			
    }
    if (myID == nWRs)
    {
      U[Mzp1][i] = 0.0;
      Fz[Mzp1][i] = -(2.0*D/dz)*(U[Mzp1][i]-U[Mzp1-1][i]);
    }
    else
    {
      Fz[Mzp1][i] = -(D/dz)*(U[Mzp1][i]-U[Mzp1-1][i]);
    }
  }
  for (j=0;j<=Mzp1;j++)
  {
    U[j][0] = 0.0;
    U[j][Mrp1] = exp(-time)*log(2.0)*sin(z[j]);
    Fr[j][1] = -(2.0*D/dr)*(U[j][1]-U[j][0]);
    Fr[j][Mrp1] = -(2.0*D/dr)*(U[j][Mrp1]-U[j][Mrp1-1]);
  }
  for(j=0;j<=Mzp1;j++)
  {
    for (i=2;i<Mrp1;i++)
    {
      Fr[j][i] = -(D/dr)*(U[j][i]-U[j][i-1]);
    }
  }
  for(i=0;i<=Mrp1;i++)
  {
    for (j=2;j<Mzp1;j++)
    {
      Fz[j][i] = -(D/dz)*(U[j][i]-U[j-1][i]);
    }
  }
}
void PDE ( double **U, double **Fr, double **Fz, double *r, int Mrp1, int Mzp1, double dr, double dz, double dt)
{
  double dr2 = 0.5*dr, dz2 = 0.5*dz,V;
  int i,j;
  double pi = 4*atan(1.0);
  for (i=1;i<Mrp1;i++)
  {
    for (j=1;j<Mzp1;j++)
    {
      V = 2.0*pi*r[i]*dr*dz;
      U[j][i] = U[j][i]-(2.0*pi*dt/V)*((Fr[j][i+1]*(r[i]+dr2)*dz)-(Fr[j][i]*(r[i]-dr2)*dz)+(r[i]*dr)*(Fz[j+1][i]-Fz[j][i]) );
    }
  }
}

