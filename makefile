convolve: main.c convolve.c
	gcc -std=c99 -Wall -O3 -msse3 -mtune=core2 -o convolve main.c convolve.c -I.
