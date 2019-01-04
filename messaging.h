#include "z.messaging.cpp"

void RECV_output_MPI( int nWRs, int Mr, int Mz, double **U);
void SEND_output_MPI( int myID, int Mr, int Mz, double **U);
void EXCHANGE_bry_MPI( int nWRs, int myID, int NodeUP, int NodeDN, int Mr, int Mz, double **U);