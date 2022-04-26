benchmark: benchmark.c
	cc -O2 -lm -o $@ $<

floyd: floyd.cpp
	cc -O2 -o $@ $<