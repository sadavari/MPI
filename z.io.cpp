#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;
void output ( ofstream *Uprofile, double *r, double *z, double **U, double time, int Mrp1, int Mzp1)
{
  int i,j;
  *Uprofile << "# time: " << time << endl;
  for ( j=0; j<=Mzp1; j++)
  {
    for (i=0; i<=Mrp1; i++)
    {
      *Uprofile << U[j][i] << " ";
    }
    *Uprofile << endl;
  }
  *Uprofile << endl;
}
