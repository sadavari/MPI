// lab9.cpp : Defines the entry point for the console application.
//

#include <math.h>
#include "setup.h"
#include "mainMR.h"
#include "mainWR.h"
#include <mpi.h>

using namespace std;
int main(int argc, char* argv[])
{
  int nPROC, nWRs, mster, myID, ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &nPROC);//returns nPROC
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &myID);//assignes myID
  mster = 0;//master has rank = 0
  nWRs = nPROC -1;//=number of workers
  if( myID == mster)
  {
    double tt0 = MPI_Wtime();// start CPU timer on master
    master(nWRs, myID);
    double tt1 = MPI_Wtime(); //end timer
    cout << ">>main>> MR timing = " << tt1-tt0 << " sec on " << nWRs << " WRs" << endl;
  }
  else
  {
    worker( nWRs, myID);  
    if ( ierr !=0 )
    {
      cout <<">>>>worker: " << myID << " ended with ierr = " << ierr << endl;
    }
  }
  ierr = MPI_Finalize();
  return 0;
}
