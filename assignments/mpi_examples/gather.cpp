#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{

  MPI::Init(argc, argv);
  int numproc = MPI::COMM_WORLD.Get_size();
  int myid    = MPI::COMM_WORLD.Get_rank();

  // Fill up the array with data to send to the destination node. Note
  // that the contents of the array will depend on the process ID
  int ndata = 3 * numproc;
  int * sendarray = new int[ndata];
  cout << "SEND " << myid << " : ";
  for (int i = 0; i < ndata; ++i) {
    sendarray[i] = 100 * myid + i;
    cout << sendarray[i] << " ";
  }
  cout << endl;

  // Need to allocate memory to receive the data - note that it needs
  // to be big enough to receive from all processors
  int * recvarray = new int[ndata * numproc];

  int dest = 1;
  MPI::COMM_WORLD.Gather(sendarray, ndata, MPI_INT, recvarray, ndata, MPI_INT, dest);

  // Look at what has been gathered - note the difference depending on
  // the node
  cout << "RECV " << myid << " : ";
  for (int i = 0; i < ndata * numproc; ++i) {
    cout << recvarray[i] << " ";
  }
  cout << endl;

  MPI::Finalize();

  delete[] sendarray;
}
