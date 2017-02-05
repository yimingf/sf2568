/* sample hello world program  *
 *  C Michael Hanke 2006-10-19 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#define M 2048
#define N 2048

typedef struct complex {
    double real;
    double imag;
} Complex;

int mb(Complex c) {
    int count, num_iter;
    Complex z;
    double foo;

    count = 0;
    num_iter = 256;
    z.real = 0;
    z.imag = 0;

    while ((z.real*z.real+z.imag*z.imag<4.0) && (count<num_iter)) { // we know if |z|<2 then z is considered converge
        foo = z.real*z.real-z.imag*z.imag+c.real;
        z.imag = 2*z.real*z.imag+c.imag;
        z.real = foo;
        count++;
    }
    return count;
}

main(int argc, char **argv) {
    int rank, size, tag, rc, i, j, k, num_rows, start, end;
    double r, rx, ry, foo;
    int *data, *data_foo;
    Complex c;
    MPI_Status status;
    FILE *fp;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // test if divided
    assert(N%size == 0);
    num_rows = N/size;
    // allocate memories for each calculations
    data = (int*)malloc(num_rows*M*sizeof(int));
    data_foo = data; // a pointer for further use.

    start = rank*num_rows;
    end = start+num_rows-1;

    r = 0.01; // display range and radius.
    rx = 0.25;
    ry = -0.01;

    for (i = start; i < end; i++) {
        c.imag = i/((double)M)*(2*r)+ry; // map array position to 2d-position
        for (j = 0; j < M; j++) {
            c.real = j/((double)N)*(2*r)+rx;
            foo = mb(c);
            *data++ = foo;
        }
    }
    data = data_foo;

    tag = 100;
    if (rank == 0) {
        // write to file.
        fp = fopen("color.txt", "w");
        printf("num_rows %d\n", num_rows);
        for (i = 0; i < num_rows; i++) {
            for (j = 0; j < M; j++) {
                fprintf(fp, "%hhu ", (unsigned char)data[j+i*M]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
        // then others' results.
        for (k = 1; k < size; k++) {
            MPI_Recv(data, num_rows*N, MPI_INTEGER, k, tag, MPI_COMM_WORLD, &status);
            printf("received message from process %d\n", k);
            fp = fopen("color.txt", "a");
            for (i = 0; i < num_rows; i++) {
                for (j = 0; j < M; j++) {
                    fprintf(fp, "%hhu ", (unsigned char)data[j+i*M]);
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }
    } else {
        rc = MPI_Send(data, num_rows*N, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
    }
    
    rc = MPI_Finalize();
}