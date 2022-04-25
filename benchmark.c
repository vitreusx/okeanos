#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_REPEATS 36

int compare_double(const void *x, const void *y) {
  double xval = *(const double *)x, yval = *(const double *)y;
  if (xval < yval)
    return -1;
  else if (xval > yval)
    return 1;
  else
    return 0;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int numProcesses, myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  for (int s = 1 << 10; s < (1 << 28); s <<= 1) {
    void *buf = malloc((size_t)s);
    int from = (myRank - 1 + numProcesses) % numProcesses;
    int next = (myRank + 1) % numProcesses;

    double execTimes[NUM_REPEATS];
    for (int rep = 0; rep < NUM_REPEATS; ++rep) {
      if (myRank == 0) {
        double startTime = MPI_Wtime();
        MPI_Send((const void *)buf, s, MPI_BYTE, next, 0, MPI_COMM_WORLD);
        MPI_Recv(buf, s, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        double endTime = MPI_Wtime();

        execTimes[rep] = endTime - startTime;
      } else {
        MPI_Recv(buf, s, MPI_BYTE, from, MPI_ANY_TAG, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        MPI_Send((const void *)buf, s, MPI_BYTE, next, 0, MPI_COMM_WORLD);
      }
    }

    if (myRank == 0) {
      qsort(execTimes, NUM_REPEATS, sizeof(double), compare_double);
      int start = 2, end = NUM_REPEATS - 2, n = end - start;

      double mean = 0.0;
      for (int rep = start; rep < end; ++rep)
        mean += execTimes[rep];
      mean /= (double)n;

      double var = 0.0;
      for (int rep = start; rep < end; ++rep)
        var += (execTimes[rep] - mean) * (execTimes[rep] - mean);
      var /= (double)(n - 1);

      double std = sqrt(var);

      printf("%d %d %f %f\n", numProcesses, s, mean, std);
    }
  }

  return 0;
}