benchmark: benchmark.c
	cc -O2 -lm -o $@ $<

floyd: floyd.c
	cc -O2 -o $@ $<