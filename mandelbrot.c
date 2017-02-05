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

int mb(Complex c) {
    int count, num_iter;
    Complex z;
    double foo;

    count = 0;
    num_iter = 256;
    z.real = 0;
    z.imag = 0;

    while ((z.real*z.real+z.imag*z.imag<4.0) && (count<num_iter)) {
        foo = z.real*z.real-z.imag*z.imag+c.real;
        z.imag = 2*z.real*z.imag+c.imag;
        z.real = foo;
        count++;
    }
    return count;
}

main(int argc, char **argv) {
    int rank, size, tag, rc, i, j, num_rows, start, end;
    double r = 2.0;
    double foo;
    int *data, *data_foo;
    Complex c;
    MPI_Status status;
    FILE *fp;

    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // assume divided.
    num_rows = M/size;
    // allocate memories for each calculations
    data = (int *) malloc(num_rows*N*sizeof(int));
    data_foo = data; // a pointer for further use.

    start = rank*num_rows;
    end = start+num_rows-1;

    for (i = start; i < end; i++) {
        c.real = i/((double)M)*(2*r)-r;
        for (j = 0; j < N; j++) {
            c.imag = j/((double)N)*(2*r)-r;
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
        fwrite(data, num_rows*N, sizeof(int), fp);
        fclose(fp);
        // then others' results.
        for (i = 1; i < size; i++) {
            MPI_Recv(data, num_rows*N, MPI_INTEGER, i, tag, MPI_COMM_WORLD, &status);
            printf("received message from process %d\n", i);
            fp = fopen("color.txt", "a");
            fwrite(data, num_rows*N, sizeof(int), fp);
            fclose(fp);
        }
    } else {
        rc = MPI_Send(data, num_rows*N, MPI_INTEGER, 0, tag, MPI_COMM_WORLD);
    }
    
    rc = MPI_Finalize();
}