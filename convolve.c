
/*File containing convolution functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include "convolve.h"
#include "xmmintrin.h"
#include "pmmintrin.h"

#define DBG 0

#define DBG_PRINTF(...) do{ if(DBG) printf(__VA_ARGS__);} while(0)


#define KERNEL_SIZE (3)

#define INT_MAX (0x7FFFFFFF)
#define INT_MIN (0x80000000)

union { __m128i ret4; 
            int ret[4]; 
        } helper_un;

// static inline __m128i char_to_int(unsigned char* char4){
//     helper_un.ret[0] = (int)(unsigned int)char4[0];
//     helper_un.ret[1] = (int)(unsigned int)char4[1];
//     helper_un.ret[2] = (int)(unsigned int)char4[2];
//     helper_un.ret[3] = (int)(unsigned int)char4[3];
//     return helper_un.ret4;
// }

//linear search for finding smallest
inline int find_smallest(int rows, int cols, short* D){
	int smallest = INT_MAX;
	int i, total;
	total = rows * cols; 
	for(i = 0; i < total; i++){
		if(D[i] < smallest){
			smallest = D[i];
		}
	}
	return smallest;
}

inline int find_largest(int rows, int cols, short* D){
	int largest = INT_MIN;
	int i, total;
	total = rows * cols; 
	for(i = 0; i < total; i++){
		if(D[i] > largest){
			largest = D[i];
		}
	}
	return largest;
}

//convolution for horizontal, keeping track of reused terms to avoid memory accesses
//only one memory access per element
inline void convolve_horizontal(int rows, int cols, short* Dx, unsigned char* M){
	//unsigned char prev, curr, next;
    int index, first_ind, i, j;
    for(i = 0; i < rows; i++){
        first_ind = i * cols;
        // prev = M[first_ind];
        // curr = M[first_ind + 1];     
        Dx[first_ind] = -M[first_ind + 1];
        for(j = 1; j < cols - 1; j++){
            index = first_ind + j;
            //next = ;
            Dx[index] = M[index - 1] - M[index + 1]; 
        }
        Dx[first_ind + cols - 1] = M[first_ind + cols - 2];
    }
}

inline void convolve_vertical(int rows, int cols, short* Dy, unsigned char* M){
	int i, j, first_ind;
    for(j = 0; j < cols; j++){
        Dy[j] = (short)(unsigned short)M[cols + j];
    }
	for(i = 1; i < rows - 1; i++){
        first_ind = i * cols;
        for(j = 0; j < cols; j++){
            Dy[first_ind + j] = (short)(unsigned short)M[(i + 1) * cols + j] - (short)(unsigned short)M[(i - 1) * cols + j];  
        }
    }
    int last_ind = (rows - 1) * cols;
    int last_last_ind = (rows - 2) * cols;
    if(rows - 2 < 0){
        last_last_ind = 0;
    }
    for(j = 0; j < cols; j++){
        Dy[last_ind + j] = -(short)(unsigned short)M[last_last_ind + j];
    }
}

//assumed the array is zero padded with thickness of one
void convolve(int rows, int cols, short* Dx, short* Dy){
    FILE *file1;
    FILE *file2;
    file1 = fopen("./horizontal_matrix.txt", "w");
    file2 = fopen("./vertical_matrix.txt", "w");
    double time_taken;
    clock_t start, end;

    //populate and print the array
    srand(time(0));
    unsigned char M[rows * cols];
    int i, j;
    int first_ind = 0;
    start = clock();
    fprintf(file1, "Original Matrix:\n");
    fprintf(file2, "Original Matrix:\n");
    for(i = 0; i < rows; i++){
        first_ind = i * cols;
        for(j = 0; j < cols; j++){
            M[first_ind + j] = rand();
            DBG_PRINTF("%03d ", (int)(unsigned int)M[first_ind + j]);
            fprintf(file1, "%03d ", (int)M[first_ind + j]);
            fprintf(file2, "%03d ", (int)M[first_ind + j]);
        }
        fprintf(file1, "\n");
        fprintf(file2, "\n");
        DBG_PRINTF("\n");
    }
    fprintf(file1, "\n");
    fprintf(file2, "\n");
    end = clock();
    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
    printf("Matrix Population Time: %fms\n", time_taken);

    //convolve horizontally
	start = clock(); 
    convolve_horizontal(rows, cols, Dx, &M[0]);
    end = clock();
    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
    printf("Horizontal Convolution Time: %fms\n", time_taken);
 
    //vertical convolution
    start = clock();
    convolve_vertical(rows, cols, Dy, &M[0]);
    end = clock();
    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
    printf("Vertical Convolution Time: %fms\n", time_taken);

    //store matrix to file
    fprintf(file1, "Horizontal Matrix:\n");
    fprintf(file2, "Vertical Matrix\n");
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            fprintf(file1, "%03d ", (int)Dx[i * cols + j]);
            fprintf(file2, "%03d ", (int)Dy[i * cols + j]);
        }
        fprintf(file1, "\n");
        fprintf(file2, "\n");
    }

    //Dx smallest
    start = clock();
    int Dx_small = find_smallest(rows, cols, Dx);
    int Dx_large = find_largest(rows, cols, Dx);
    end = clock();
    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
    printf("Dx Min: %d Max: %d\n", Dx_small, Dx_large);
    printf("Dx Min/Max Time: %fms\n", time_taken);

    //Dy smallest
    start = clock();
    int Dy_small = find_smallest(rows, cols, Dy);
    int Dy_large = find_largest(rows, cols, Dy);
    end = clock();
    time_taken = ((double)(end - start)) / CLOCKS_PER_SEC * 1000;
    printf("Dy Min: %d Max: %d\n", Dy_small, Dy_large);
    printf("Dy Min/Max Time: %fms\n", time_taken);
    

}
