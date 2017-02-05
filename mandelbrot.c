/* sample hello world program  *
 *  C Michael Hanke 2006-10-19 */

#include <stdio.h>
#include <string.h>
#include <mpi.h>

#define M 2048
#define N 2048

typedef struct complex {
    double real;
    double imag;
} Complex;

int mandelbrot(Complex c) {
    // implement the mandelbrot algorithm.
}

main(int argc, char **argv) {
    int rank, size, tag, rc, i, j;
    MPI_Status status;
    unsigned char color[M*N];
    FILE *fp;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // argv = nx, ny
    // check if nx * ny = size or throw exceptions
    // check if divided or throw exceptions
    // allocate memories for each calculations
    // DATA_TYPE??

    tag = 100;
    if (rank == 0) {
        // write to file.
        // firstly write rank 0's result.
    	fp = fopen('color.txt', 'w');
        for (j = 0; j < N; j++) {
            for (i = 0; i < M; i++) {
                fprintf(fp, '%hhu', color[i+j*M]);
            }
            fprintf(fp, '\n');
        }
        fclose(fp);
    } else {
        rc = MPI_Send(message, /* size */, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }
	
    
    rc = MPI_Finalize();
}
