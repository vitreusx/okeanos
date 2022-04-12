#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int numProcesses, myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  srand(time(NULL));

  if (myRank == 0) {
    int rank, data;
    for (int node = 0; node < numProcesses - 1; ++node) {
      MPI_Recv((void *)&rank, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, NULL);
      MPI_Recv((void *)&data, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD,
               NULL);
      printf("[node %d] received %d from node #%d\n", myRank, data, rank);
    }
  } else {
    int data = rand();
    printf("[node %d] generated %d\n", myRank, data);
    MPI_Send((const void *)&myRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send((const void *)&data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  }

  return 0;
}