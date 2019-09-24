#include <stdio.h>
#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
  const int ndata = 5;

  MPI::Init(argc, argv);
  int numproc = MPI::COMM_WORLD.Get_size();
  int myid    = MPI::COMM_WORLD.Get_rank();

  // Fill up the array with data to send to the destination node. Note
  // that the contents of the array will
  int * sendbuf = new int[ndata * numproc];
  int * recvbuf = new int[ndata];
  
  int root = 0;

  if (myid == root) {
    cout << "SEND " << myid << " : ";
    for (int i = 0; i < ndata * numproc; ++i) {
      sendbuf[i] = i;
      cout << sendbuf[i] << " ";
    }
    cout << endl;
  }

  MPI::COMM_WORLD.Scatter(sendbuf, ndata, MPI_INT, recvbuf, ndata, MPI_INT, root);

  cout << "RECV " << myid << " : ";
  for (int i = 0; i < ndata; ++i) {
    cout << recvbuf[i] << " ";
  }
  cout << endl;

  MPI::Finalize();

  delete[] sendbuf;
  delete[] recvbuf;
}
