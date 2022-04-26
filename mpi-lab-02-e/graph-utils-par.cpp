/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include "graph-base.h"
#include "graph-utils.h"
#include <cassert>
#include <mpi.h>

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {
  return (numVertices * myRank) / numProcesses;
}

Graph *createAndDistributeGraph(int numVertices, int numProcesses, int myRank) {
  assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);

  auto begin = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank);
  auto end = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1);

  auto graph = allocateGraphPart(numVertices, begin, end);

  if (graph == nullptr) {
    return nullptr;
  }

  assert(graph->numVertices > 0 && graph->numVertices == numVertices);
  assert(graph->firstRowIdxIncl >= 0 &&
         graph->lastRowIdxExcl <= graph->numVertices);

  for (int row = begin; row < end; ++row)
    initializeGraphRow(graph->data[row - begin], row, numVertices);

  return graph;
}

void collectAndPrintGraph(Graph *graph, int numProcesses, int myRank) {
  assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
  assert(graph->numVertices > 0);
  assert(graph->firstRowIdxIncl >= 0 &&
         graph->lastRowIdxExcl <= graph->numVertices);

  auto *buf = new int[graph->numVertices];

  for (int rank = 0; rank < numProcesses; ++rank) {
    auto begin =
        getFirstGraphRowOfProcess(graph->numVertices, numProcesses, rank);
    auto end =
        getFirstGraphRowOfProcess(graph->numVertices, numProcesses, rank + 1);

    for (int row = begin; row < end; ++row) {
      void *src = (rank == myRank) ? graph->data[row - begin] : buf;
      MPI_Bcast(src, graph->numVertices, MPI_INT, rank, MPI_COMM_WORLD);
      printGraphRow((int *)src, row, graph->numVertices);
    }
  }

  delete[] buf;
}

void destroyGraph(Graph *graph, int numProcesses, int myRank) {
  freeGraphPart(graph);
}
