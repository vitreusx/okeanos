benchmark: benchmark.cpp
	mpicc -std=c++17 -O2 -o benchmark benchmark.cpp