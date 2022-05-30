/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include "laplace-common.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <sys/time.h>
#include <vector>

#define OPTION_VERBOSE "--verbose"
#define ASYNC_COMM
#define ALLREDUCE

static void printUsage(char const *progName) {
  std::cerr
      << "Usage:" << std::endl
      << "    " << progName << " [--verbose] <N>" << std::endl
      << "Where:" << std::endl
      << "   <N>         The number of points in each dimension (at least 4)."
      << std::endl
      << "   " << OPTION_VERBOSE << "   Prints the input and output systems."
      << std::endl;
}

static InputOptions parseInput(int argc, char *argv[], int numProcesses) {
  int numPointsPerDimension = 0;
  bool verbose = false;
  int errorCode = 0;

  if (argc < 2) {
    std::cerr << "ERROR: Too few arguments!" << std::endl;
    printUsage(argv[0]);
    errorCode = 1;
    MPI_Finalize();
  } else if (argc > 3) {
    std::cerr << "ERROR: Too many arguments!" << std::endl;
    printUsage(argv[0]);
    errorCode = 2;
    MPI_Finalize();
  } else {
    int argIdx = 1;

    if (argc == 3) {
      if (strncmp(argv[argIdx], OPTION_VERBOSE, strlen(OPTION_VERBOSE)) != 0) {
        std::cerr << "ERROR: Unexpected option '" << argv[argIdx] << "'!"
                  << std::endl;
        printUsage(argv[0]);
        errorCode = 3;
        MPI_Finalize();
      }
      verbose = true;
      ++argIdx;
    }

    numPointsPerDimension = std::strtol(argv[argIdx], nullptr, 10);

    if ((numPointsPerDimension < 4) ||
        (numProcesses > numPointsPerDimension / 2)) {
      /* If we had a smaller grid, we could use the sequential version. */
      std::cerr << "ERROR: The number of points, '" << argv[argIdx]
                << "', should be an iteger greater than or equal to 4; and at "
                   "least 2 points per process!"
                << std::endl;
      printUsage(argv[0]);
      MPI_Finalize();
      errorCode = 4;
    }
  }

  return {numPointsPerDimension, verbose, errorCode};
}

static std::tuple<int, double> performAlgorithm(int myRank, int numProcesses,
                                                GridFragment *frag,
                                                double omega, double epsilon) {
  int startRowIncl = frag->firstRowIdxIncl + (myRank == 0 ? 1 : 0);
  int endRowExcl = frag->lastRowIdxExcl - (myRank == numProcesses - 1 ? 1 : 0);
  int numTotalRows = frag->lastRowIdxExcl - frag->firstRowIdxIncl + 2;

  double maxDiff = 0;
  int numIterations = 0;

  /* TODO: change the following code fragment */
  /* Implement asynchronous communication of neighboring elements */
  /* and computation of the grid */
  /* the following code just recomputes the appropriate grid fragment */
  /* but does not communicate the partial results */
  while (true) {
    maxDiff = 0.0;

    for (int color = 0; color < 2; ++color) {
      int comm_color = 1 - color;

      auto *low_extra = frag->data[comm_color][0],
           *low = frag->data[comm_color][1],
           *high = frag->data[comm_color][numTotalRows - 2],
           *high_extra = frag->data[comm_color][numTotalRows - 1];
      auto low_extra_n =
               frag->getNumColorPointsInRow(startRowIncl - 1, comm_color),
           low_n = frag->getNumColorPointsInRow(startRowIncl, comm_color),
           high_n = frag->getNumColorPointsInRow(endRowExcl - 1, comm_color),
           high_extra_n = frag->getNumColorPointsInRow(endRowExcl, comm_color);

      if (myRank % 2 == 0) {
        if (myRank - 1 >= 0) {
          MPI_Send(low, low_n, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD);
          MPI_Recv(low_extra, low_extra_n, MPI_DOUBLE, myRank - 1, MPI_ANY_TAG,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (myRank + 1 < numProcesses) {
          MPI_Send(high, high_n, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD);
          MPI_Recv(high_extra, high_extra_n, MPI_DOUBLE, myRank + 1,
                   MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      } else {
        if (myRank + 1 < numProcesses) {
          MPI_Recv(high_extra, high_extra_n, MPI_DOUBLE, myRank + 1,
                   MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(high, high_n, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD);
        }
        if (myRank - 1 >= 0) {
          MPI_Recv(low_extra, low_extra_n, MPI_DOUBLE, myRank - 1, MPI_ANY_TAG,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Send(low, low_n, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD);
        }
      }

      for (int rowIdx = startRowIncl; rowIdx < endRowExcl; ++rowIdx) {
        for (int colIdx = 1 + (rowIdx % 2 == color ? 1 : 0);
             colIdx < frag->gridDimension - 1; colIdx += 2) {
          double tmp =
              (GP(frag, rowIdx - 1, colIdx) + GP(frag, rowIdx + 1, colIdx) +
               GP(frag, rowIdx, colIdx - 1) + GP(frag, rowIdx, colIdx + 1)) /
              4.0;
          double diff = GP(frag, rowIdx, colIdx);
          GP(frag, rowIdx, colIdx) = (1.0 - omega) * diff + omega * tmp;
          diff = fabs(diff - GP(frag, rowIdx, colIdx));

          if (diff > maxDiff) {
            maxDiff = diff;
          }
        }
      }
    }

    ++numIterations;
#ifndef ALLREDUCE
    double globalMaxDiff = 0.0;
    for (int rank = 0; rank < numProcesses; ++rank) {
      double rankMaxDiff = maxDiff;
      MPI_Bcast(&rankMaxDiff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
      globalMaxDiff = std::max(globalMaxDiff, rankMaxDiff);
    }
#else
    double globalMaxDiff = 0.0;
    MPI_Allreduce(&maxDiff, &globalMaxDiff, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
#endif

    if (globalMaxDiff <= epsilon)
      break;
  }

  /* no code changes beyond this point should be needed */

  return std::make_tuple(numIterations, maxDiff);
}

static std::tuple<int, double>
performAlgorithmAsync(int myRank, int numProcesses, GridFragment *frag,
                      double omega, double epsilon) {
  int startRowIncl = frag->firstRowIdxIncl + (myRank == 0 ? 1 : 0);
  int endRowExcl = frag->lastRowIdxExcl - (myRank == numProcesses - 1 ? 1 : 0);
  int numTotalRows = frag->lastRowIdxExcl - frag->firstRowIdxIncl + 2;

  double maxDiff = 0;
  int numIterations = 0;

  MPI_Request requests[4];

  /* TODO: change the following code fragment */
  /* Implement asynchronous communication of neighboring elements */
  /* and computation of the grid */
  /* the following code just recomputes the appropriate grid fragment */
  /* but does not communicate the partial results */
  while (true) {
    maxDiff = 0.0;

    for (int color = 0; color < 2; ++color) {
      int comm_color = 1 - color;

      auto *low_extra = frag->data[comm_color][0],
           *low = frag->data[comm_color][1],
           *high = frag->data[comm_color][numTotalRows - 2],
           *high_extra = frag->data[comm_color][numTotalRows - 1];
      auto low_extra_n =
               frag->getNumColorPointsInRow(startRowIncl - 1, comm_color),
           low_n = frag->getNumColorPointsInRow(startRowIncl, comm_color),
           high_n = frag->getNumColorPointsInRow(endRowExcl - 1, comm_color),
           high_extra_n = frag->getNumColorPointsInRow(endRowExcl, comm_color);

      if (myRank % 2 == 0) {
        if (myRank - 1 >= 0) {
          MPI_Isend(low, low_n, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD,
                    &requests[0]);
          MPI_Irecv(low_extra, low_extra_n, MPI_DOUBLE, myRank - 1, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &requests[1]);
        }
        if (myRank + 1 < numProcesses) {
          MPI_Isend(high, high_n, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD,
                    &requests[2]);
          MPI_Irecv(high_extra, high_extra_n, MPI_DOUBLE, myRank + 1,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &requests[3]);
        }
      } else {
        if (myRank + 1 < numProcesses) {
          MPI_Irecv(high_extra, high_extra_n, MPI_DOUBLE, myRank + 1,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &requests[0]);
          MPI_Isend(high, high_n, MPI_DOUBLE, myRank + 1, 0, MPI_COMM_WORLD,
                    &requests[1]);
        }
        if (myRank - 1 >= 0) {
          MPI_Irecv(low_extra, low_extra_n, MPI_DOUBLE, myRank - 1, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &requests[2]);
          MPI_Isend(low, low_n, MPI_DOUBLE, myRank - 1, 0, MPI_COMM_WORLD,
                    &requests[3]);
        }
      }

      for (int rowIdx = startRowIncl + 1; rowIdx < endRowExcl - 1; ++rowIdx) {
        for (int colIdx = 1 + (rowIdx % 2 == color ? 1 : 0);
             colIdx < frag->gridDimension - 1; colIdx += 2) {
          double tmp =
              (GP(frag, rowIdx - 1, colIdx) + GP(frag, rowIdx + 1, colIdx) +
               GP(frag, rowIdx, colIdx - 1) + GP(frag, rowIdx, colIdx + 1)) /
              4.0;
          double diff = GP(frag, rowIdx, colIdx);
          GP(frag, rowIdx, colIdx) = (1.0 - omega) * diff + omega * tmp;
          diff = fabs(diff - GP(frag, rowIdx, colIdx));

          if (diff > maxDiff) {
            maxDiff = diff;
          }
        }
      }

      if (myRank % 2 == 0) {
        if (myRank - 1 >= 0) {
          MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
          MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
        }
        if (myRank + 1 < numProcesses) {
          MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
          MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        }
      } else {
        if (myRank + 1 < numProcesses) {
          MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
          MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
        }
        if (myRank - 1 >= 0) {
          MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
          MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
        }
      }

      for (int rowIdx : {startRowIncl, endRowExcl - 1}) {
        for (int colIdx = 1 + (rowIdx % 2 == color ? 1 : 0);
             colIdx < frag->gridDimension - 1; colIdx += 2) {
          double tmp =
              (GP(frag, rowIdx - 1, colIdx) + GP(frag, rowIdx + 1, colIdx) +
               GP(frag, rowIdx, colIdx - 1) + GP(frag, rowIdx, colIdx + 1)) /
              4.0;
          double diff = GP(frag, rowIdx, colIdx);
          GP(frag, rowIdx, colIdx) = (1.0 - omega) * diff + omega * tmp;
          diff = fabs(diff - GP(frag, rowIdx, colIdx));

          if (diff > maxDiff) {
            maxDiff = diff;
          }
        }
      }
    }

    ++numIterations;
#ifndef ALLREDUCE
    double globalMaxDiff = 0.0;
    for (int rank = 0; rank < numProcesses; ++rank) {
      double rankMaxDiff = maxDiff;
      MPI_Bcast(&rankMaxDiff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);
      globalMaxDiff = std::max(globalMaxDiff, rankMaxDiff);
    }
#else
    double globalMaxDiff = 0.0;
    MPI_Allreduce(&maxDiff, &globalMaxDiff, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
#endif

#ifdef ASYNC_COMM

#endif

    if (globalMaxDiff <= epsilon)
      break;
  }

  /* no code changes beyond this point should be needed */

  return std::make_tuple(numIterations, maxDiff);
}

int main(int argc, char *argv[]) {
  int numProcesses;
  int myRank;
  struct timeval startTime {};
  struct timeval endTime {};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  auto inputOptions = parseInput(argc, argv, numProcesses);

  if (inputOptions.getErrorCode() != 0) {
    return inputOptions.getErrorCode();
  }

  auto numPointsPerDimension = inputOptions.getNumPointsPerDimension();
  auto isVerbose = inputOptions.isVerbose();

  double omega = Utils::getRelaxationFactor(numPointsPerDimension);
  double epsilon = Utils::getToleranceValue(numPointsPerDimension);

  auto gridFragment =
      new GridFragment(numPointsPerDimension, numProcesses, myRank);
  gridFragment->initialize();

  if (gettimeofday(&startTime, nullptr)) {
    gridFragment->free();
    std::cerr << "ERROR: Gettimeofday failed!" << std::endl;
    MPI_Finalize();
    return 6;
  }

  /* Start of computations. */

#ifdef ASYNC_COMM
  auto result =
      performAlgorithmAsync(myRank, numProcesses, gridFragment, omega, epsilon);
#else
  auto result =
      performAlgorithm(myRank, numProcesses, gridFragment, omega, epsilon);
#endif

  /* End of computations. */

  if (gettimeofday(&endTime, nullptr)) {
    gridFragment->free();
    std::cerr << "ERROR: Gettimeofday failed!" << std::endl;
    MPI_Finalize();
    return 7;
  }

  double duration =
      ((double)endTime.tv_sec + ((double)endTime.tv_usec / 1000000.0)) -
      ((double)startTime.tv_sec + ((double)startTime.tv_usec / 1000000.0));

  std::cerr << "Statistics: duration(s)=" << std::fixed << std::setprecision(10)
            << duration << " #iters=" << std::get<0>(result)
            << " diff=" << std::get<1>(result) << " epsilon=" << epsilon
            << std::endl;

  if (isVerbose) {
    gridFragment->printEntireGrid(myRank, numProcesses);
  }

  gridFragment->free();
  MPI_Finalize();
  return 0;
}
