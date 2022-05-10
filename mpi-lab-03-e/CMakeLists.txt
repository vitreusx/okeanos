cmake_minimum_required(VERSION 3.10)
SET(CMAKE_CXX_COMPILER CC)   # use CC for okeanos

project(mpi_lab_03 CXX)

set(CMAKE_CXX_STANDARD 14)


include_directories(.)

add_executable(laplace-seq.exe
        laplace-common.h
        laplace-common.cpp
        laplace-seq.cpp)

add_executable(laplace-par.exe
        laplace-common.h
        laplace-common.cpp
        laplace-par.cpp)

add_executable(ring-nonblocking.exe
        ring-nonblocking.cpp)

add_executable(scatter-gather-max.exe
        scatter-gather-max.cpp)

find_package(MPI REQUIRED)

target_link_libraries(laplace-par.exe PRIVATE MPI::MPI_CXX)
target_link_libraries(ring-nonblocking.exe PRIVATE MPI::MPI_CXX)
target_link_libraries(scatter-gather-max.exe PRIVATE MPI::MPI_CXX)
