#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[])
{
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI::Init(argc,argv);
  int numproc = MPI::COMM_WORLD.Get_size();
  int myid    = MPI::COMM_WORLD.Get_rank();

  // If master process, get N from user input
  int N = 0;
  if (myid == 0) N = atoi(argv[1]);
  cout << "id=" << myid << " N=" << N << " before broadcast" << endl;

  // Broadcast from master to all slaves. If this is the master
  // prcess, then this routine is sending, if this is a slave process,
  // then the routine is receiving
  MPI::COMM_WORLD.Bcast(&N, 1, MPI_INT, 0);
  cout << "id=" << myid << " N=" << N << " after broadcast" << endl;

  int sum = 0;
  int first = N / numproc * myid + 1;
  int last  = N / numproc * (myid + 1);
  cout << "id=" << myid << " first=" << first << " last=" << last << endl;
  for (int i=first; i <= last; ++i) sum +=i;

  // Combine the sums computed on each node, into result. If this is
  // one of the slaves, this routine will send sum to the master. If
  // this is the master process, , then it will receive the sums and
  // combine them.
  int result = 0;
  cout << "id=" << myid << " before reduce: sum=" << sum << " result=" << result << endl;
  MPI::COMM_WORLD.Reduce(&sum, &result, 1, MPI_INT, MPI_SUM, 0);
  cout << "id=" << myid << " after reduce: sum=" << sum << " result=" << result << endl;

  if (myid==0) cout << "Sum= " << result << endl;

  MPI::Finalize();
}
