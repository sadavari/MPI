// mainWR.cpp 
//

#include <math.h>
#include "update.h"

using namespace std;

void worker(int nWRs, int myID)
{
  double Rin,Rout,Z,D,dtout,factor,time,tout,tend;
  double dr,dz,dtexpl,dt;
  int Mr,Mz,Mz_tot,nsteps,Nmax,Niparams,Nparams,ierr;
  double MMr,MMz;
  double my_Z0,my_Z;
  int NodeUP,NodeDN;
  int Mrp1, Mzp1;
  Nparams = 11;// # of double params
  double *params = new double[Nparams];
  ierr = MPI_Bcast(params, Nparams, MPI_DOUBLE, 0, MPI_COMM_WORLD);// receives double params from master
  Rin = params[0]; Rout = params[1]; Z = params[2]; D = params[3]; factor = params[4]; tend = params[5]; dtout = params[6]; MMr = params[7]; Mz_tot = params[8];
  NodeUP = myID + 1;//calculat the up node's ID
  NodeDN = myID -1;//calculat the down node's ID
  Mr = (Rout - Rin) * MMr;
  dr = 1.0/(double)MMr;
  MMz=Mz_tot/Z;//total # of nodes
  Mz = Mz_tot / nWRs;// # of mesh nodes at this process
  dz = 1.0/(double)MMz;
  dtexpl = 0.5*(dr*dz)*(dr*dz)/(D*((dr*dr)+(dz*dz)));
  dt = factor * dtexpl;	
  Mrp1 = Mr+1;
  Mzp1 = Mz+1;
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
  my_Z0 = Mz * (myID-1) * dz;
  my_Z = Mz * myID * dz;
  if (myID == nWRs)
    my_Z = Z;
  MESH( r, z, Rin, Rout, my_Z0, my_Z, dr, dz, Mrp1, Mzp1, myID, nWRs);// initialize local z vector for this process.
	// Initialization
  nsteps = 0;
  time = 0;
  tout = dtout;
  Nmax = tend/dt+1;	
  INIT( r, z, U, Mrp1, Mzp1);//set initial values
  SEND_output_MPI ( myID, Mr, Mz, U);//send U at time = 0 to master
  do
  {
    ierr = MPI_Barrier(MPI_COMM_WORLD);//synchronize with other workers		
    EXCHANGE_bry_MPI( nWRs, myID, NodeUP, NodeDN, Mr, Mz, U);// exchange boundary values
    FLUX( U, Fr, Fz, z, D, Mrp1, Mzp1, dr, dz, time, myID, nWRs);// calculate fluxes at this time step
    PDE( U, Fr, Fz, r, Mrp1, Mzp1, dr, dz, dt);// calculate U at this time step
    time += dt;
    if ( time >= tout)
    {
      SEND_output_MPI ( myID, Mr, Mz, U);//send U at time = tout to master
      tout += dtout;
    }
  }while (time<=tend);
  EXCHANGE_bry_MPI( nWRs, myID, NodeUP, NodeDN, Mr, Mz, U);// exchange final boundary values
  ierr = MPI_Barrier(MPI_COMM_WORLD);//synchronize with other workers
  SEND_output_MPI ( myID, Mr, Mz, U);//send final U to master 
}
