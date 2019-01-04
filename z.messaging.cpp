#include <mpi.h>

void RECV_output_MPI( int nWRs, int Mr, int Mz, double **U)
{
	int I2, J2,ierr, Ime, msgtag;
	MPI_Status status;
	MPI_Datatype mytype;
	I2 = Mr +2;//# of elements to be received
	ierr = MPI_Type_contiguous(I2, MPI_DOUBLE, &mytype);
	ierr = MPI_Type_commit ( &mytype);

	J2 = Mz + 2;
	
	for ( int i=1; i<= nWRs; i++)//receive local U from every worker
	{
		Ime = (i-1) * Mz;
		msgtag = 1000 + i;
		MPI_Recv (&U[Ime][0], J2, mytype,  MPI_ANY_SOURCE, msgtag, MPI_COMM_WORLD, &status);
	}

	ierr = MPI_Type_free(& mytype);
}

void SEND_output_MPI( int myID, int Mr, int Mz, double **U)
{
	
        int I2, J2,ierr, msgtag;
		MPI_Datatype mytype;
        I2 = Mr +2;//# of elements to be received from each worker
        ierr = MPI_Type_contiguous(I2, MPI_DOUBLE, &mytype);
        ierr = MPI_Type_commit ( &mytype);

        int mster = 0;
		J2 = Mz + 2;
        msgtag = 1000 + myID;

        MPI_Send (&U[0][0], J2, mytype, mster, msgtag, MPI_COMM_WORLD);//send local U from the worker to master

        ierr = MPI_Type_free(& mytype);
}
void EXCHANGE_bry_MPI( int nWRs, int myID, int NodeUP, int NodeDN, int Mr, int Mz, double **U)
{
	int I2,ierr, Ime, msgtag, msgUP, msgDN;
	MPI_Status status;
	MPI_Datatype mytype;
	int Iup = Mz;
	int Iup1 = Iup + 1;
	msgUP = 10;
	msgDN = 20;
	I2 = Mr + 2 ;//# of boundary values to be sent to neighbers
	ierr = MPI_Type_contiguous(I2, MPI_DOUBLE, &mytype);
	ierr = MPI_Type_commit ( &mytype);


	if ( myID != 1 )
	{
		msgtag = msgDN;
		ierr = MPI_Send (&U[1][0], 1, mytype, NodeDN, msgtag, MPI_COMM_WORLD);//sned boundary value to the lower node
	}
	if ( myID != nWRs )
	{
		msgtag = msgDN;
		ierr = MPI_Recv (&U[Iup1][0], 1, mytype, NodeUP, msgtag, MPI_COMM_WORLD, &status);//receive boundary value from the right node
	}

	if ( myID != nWRs )
	{
		msgtag = msgUP;
		ierr = MPI_Send (&U[Iup][0], 1, mytype, NodeUP, msgtag, MPI_COMM_WORLD);//sned boundary value to the right node
	}
	if ( myID != 1 )
	{
		msgtag = msgUP;
		ierr = MPI_Recv (&U[0][0], 1, mytype, NodeDN, msgtag, MPI_COMM_WORLD, &status);//receive boundary value from the left node
	}

	ierr = MPI_Type_free(& mytype);

}