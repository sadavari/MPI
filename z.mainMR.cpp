// mainMR.cpp 
//

#include <math.h>
#include "io.h"
#include "messaging.h"
#include <mpi.h>

using namespace std;

#define MAX(a,b) (((a)>(b))?(a):(b))

double comparison (ofstream *ferr, double *r, double *z, double **U, double time,int Mrp1, int Mzp1);

void master(int nWRs, int myID)
{
  double Rin,Rout,Z,Z0,D,dtout,factor,time,tout,tend;
  double dr,dz,dtexpl,dt;
  int Mr,Mz,Mz_tot,nsteps,Nmax,Niparams,Nparams,ierr;
  int Mrp1, Mzp1;
  double MMr,MMz;
  FILE *fptr;
  ofstream   Uprofile("Uprofile.dat", ios::out);//File for U-profile
  ofstream   ferr("error.dat", ios::out);//File for errors
  fptr=fopen("INPUT.dat","r+");
  fscanf(fptr,"Rin=%lf Rout=%lf Z=%lf MMr=%lf Mz_tot=%d D=%lf factor=%lf dtout=%lf tend=%lf",&Rin,&Rout,&Z,&MMr,&Mz_tot,&D,&factor,&dtout,&tend); 
  fclose(fptr);     
  Nparams = 11;// # of double params
  double *params = new double[Nparams];
  Z0 = 0;
	//Rin=1.0; Rout=2.0; Z=3.1415926535; MMr=32; MMz=10.185917; D=1.0; /*k=1.0; rho=1.0; Cp=1.0;*/ factor=0.95; tend=1.0; dtout=0.1; 
  Mr = (Rout - Rin) * MMr;
  dr = 1.0/(double)MMr;
  MMz=Mz_tot/Z;//total # of nodes
  Mz = Mz_tot / nWRs;// # of mesh nodes at each process
  dz = 1.0/(double)MMz;
  dtexpl = 0.5*(dr*dz)*(dr*dz)/(D*((dr*dr)+(dz*dz)));
  dt = factor * dtexpl;
  Mrp1 = Mr+1;
  Mzp1 = Mz_tot+1;
  double *r = new double [Mrp1+1];
  double *z = new double [Mzp1+1];
	// Allocating contiguous memory
  double **U = new (nothrow)double* [Mzp1+1];
  double **Fr = new (nothrow)double* [Mzp1+1];
  double **Fz = new (nothrow)double* [Mzp1+1];
  if(U==0 || Fr==0 || Fz==0)
    printf("Allocation Error\n");
  else
  {
    U[0] = new double[(Mzp1+1)*(Mrp1+1)];
    Fr[0] = new double[(Mzp1+1)*(Mrp1+1)];
    Fz[0] = new double[(Mzp1+1)*(Mrp1+1)];
    for (int i=1; i<=Mzp1; i++)
    {
      U[i] = U[i-1]+Mrp1+1;
      Fr[i] = Fr[i-1]+Mrp1+1;
      Fz[i] = Fz[i-1]+Mrp1+1;
    }
  }
  MESH( r, z, Rin, Rout, Z0, Z, dr, dz, Mrp1, Mzp1, myID, nWRs);// initializes x vector for the main (big) problem.
  params[0] = Rin; params[1] = Rout; params[2] = Z; params[3] = D; params[4] = factor; params[5] = tend; params[6] = dtout; params[7] = MMr; params[8] = Mz_tot;
  /*for (int i=1;i<=8;i++)
    cout<<params[i]<<endl;
  for (int i=0;i<=Mzp1;i++)
    cout<<z[i]<<endl;*/
  ierr = MPI_Bcast(params, Nparams, MPI_DOUBLE, 0, MPI_COMM_WORLD);// sends double params to all workers
  time = 0;
  tout = dtout;
  RECV_output_MPI( nWRs, Mr, Mz, U);//receives U at time = 0 from all workers
  output ( &Uprofile, r, z, U, time, Mrp1, Mzp1);// writes the profile at time =0 to the file
  do
  {
    ierr = MPI_Barrier(MPI_COMM_WORLD);// synchronizes all workers			
    time += dt;
    if ( time >= tout)
    {
      RECV_output_MPI( nWRs, Mr, Mz, U);//receives U at time = tout from all workers
      output ( &Uprofile, r, z, U, time, Mrp1, Mzp1);// writes the profile at time = tout to the file
      tout += dtout;
    }	
  }while(time<=tend);
  ierr = MPI_Barrier(MPI_COMM_WORLD);// synchronizes all workers
  RECV_output_MPI( nWRs, Mr, Mz, U);//receives final U from all workers
  cout << "Maximum error is: " << comparison (&ferr, r, z, U, time, Mrp1, Mzp1) << endl;
  Uprofile.close();
  ferr.close();
}
double comparison (ofstream *ferr, double *r, double *z, double **U, double time,int Mrp1, int Mzp1)
{
  int i,j;
  double U_exact, error, maxerr=0;
  for (i =0;i<=Mrp1;i++)
  {
    for (j=0;j<=Mzp1;j++)
    {
      U_exact = exp(-time) * log(r[i]) * sin(z[j]);
      error = U_exact - U[j][i];
      maxerr = MAX (maxerr , fabs(error));
      *ferr << error << "     ";
    }
    *ferr << endl;
  }
  return maxerr;
}
