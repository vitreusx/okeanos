#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int numProcesses, myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  printf("numProcesses = %d, myRank = %d\n", numProcesses, myRank);

  return 0;
}