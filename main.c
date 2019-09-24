/**
 * File containing the main function
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "convolve.h"



int main(int argc, char **argv){
    if(argc != 3){
        printf("Wrong number of arguments supplied (2 required).\n");
        return -1;
    }
    if(argc == 3){
        printf("Rows: %s Cols: %s\n", argv[1], argv[2]);
    }
    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);
    short Dx[rows * cols];
    short Dy[rows * cols];
    convolve(rows, cols, &Dx[0], &Dy[0]);
    return 0;
}
