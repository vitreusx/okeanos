#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int numProcesses, myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  for (size_t s = (size_t)1 << 10; s < ((size_t)1 << 28); s <<= 1) {
    void *buf = malloc(s);
    int from = (myRank - 1 + numProcesses) % numProcesses;
    int next = (myRank + 1) % numProcesses;

    if (myRank == 0) {
      double startTime = MPI_Wtime();
      MPI_Send((const void *)buf, s, MPI_BYTE, next, 0, MPI_COMM_WORLD);
      MPI_Recv(buf, s, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      double endTime = MPI_WTime();

      double execTime = endTime - startTime;
      printf("n = %d, s = %d, t = %f\n", numProcesses, s, execTime);
    } else {
      MPI_Recv(buf, s, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Send((const void *)buf, s, MPI_BYTE, next, 0, MPI_COMM_WORLD);
    }
  }

  return 0;
}