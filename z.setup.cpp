#include <math.h>

void MESH ( double *r, double *z, double Rin, double Rout, double Z0, double Z, double dr, double dz, int Mrp1, int Mzp1, int myID, int nWRs)
{
  int i,j;
  r[0] = Rin;
  r[1] = Rin + 0.5*dr;
  if (myID ==1 || myID ==0)
  {
    z[0] = Z0;
  }
  else
  {
    z[0] = Z0 - 0.5*dz;	
  }
  z[1] = Z0 + 0.5*dz;
  for (i=2;i<Mrp1;i++)
  {
    r[i] = r[i-1] + dr;
  }
  r[Mrp1] = Rout;
  for (i=2;i<Mzp1;i++)
  {
    z[i] = z[i-1] + dz;
  }
  if (myID ==nWRs || myID ==0)
  {
    z[Mzp1] = Z;
  }
  else
  {
    z[Mzp1] = Z + 0.5*dz;
  }
}
void INIT ( double *r, double *z, double **U, int Mrp1, int Mzp1)
{
  int i,j;
  for (i=0;i<=Mrp1;i++)
  {
    for ( j=0;j<=Mzp1;j++)
    {
      U[j][i] = log(r[i])*sin(z[j]);
    }
  }
}

