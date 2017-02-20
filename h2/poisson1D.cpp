#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "mpi.h"

using namespace std;

// coefficient functions, analytical answer is u = sin(5pi*x)
double r(const double x) {
  return -3*x;
}
double f(const double x) {
  return exp(-3*x);
}

int main(int argc, char** argv) {
  const int n = 1000, iter_max = 1000000;
  const double h = 1.0/(n+1);
  int rank, size, tag = 100;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int L = n/size, R = n%size;
  int I = (n+size-rank-1)/size;
  int first = rank*L+min(rank, R)+1;

  typedef vector<double> Vec;
  Vec rr(I+2), ff(I+2);
  double x_i;

  for(int i = 0; i<ff.size(); ++i) {
    x_i = (first-1+i)*h;
    rr[i] = r(x_i);
    ff[i] = f(x_i);
  }
  
  Vec u(I+2, 0), u_new(I+2, 0);
  int right = rank+1, left = rank-1;
  if(rank == size-1) {
    right = MPI_PROC_NULL;
  } else if(rank == 0) {
    left = MPI_PROC_NULL;
  }

  for(int step = 0; step < iter_max; ++step) {
     if(rank%2 == 0) { // red
      MPI_Send(&u[I], 1, MPI_DOUBLE, right, tag, MPI_COMM_WORLD);
      MPI_Recv(&u.back(), 1, MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&u[1], 1, MPI_DOUBLE, left, tag, MPI_COMM_WORLD);
      MPI_Recv(&u.front(), 1, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &status);
    } else { // black
      MPI_Recv(&u.front(), 1, MPI_DOUBLE, left, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&u[1], 1, MPI_DOUBLE, left, tag, MPI_COMM_WORLD);
      MPI_Recv(&u.back(), 1,MPI_DOUBLE, right, tag, MPI_COMM_WORLD, &status);
      MPI_Send(&u[I], 1,MPI_DOUBLE, right, tag, MPI_COMM_WORLD);
    }
    
    // Jacobi iteration over inner points
    for(int i = 1; i < u.size() - 1; ++i) {
      u_new[i] = (u[i-1]+u[i+1]-h*h*ff[i])/(2.0-h*h*rr[i]);
    }

    u.swap(u_new);  // runs in O(1)
  }

  const char file[] = "out";
  ofstream fout;

  MPI_Recv(0, 0, MPI_CHAR, left, tag, MPI_COMM_WORLD, &status);

  if(rank == 0) {
    fout.open(file); fout << 0 << ' ';
  } // left boundary
  else fout.open(file, ios::app);
  
  for(int i = 1; i < u.size()-1; ++i)
    fout << u[i] << ' ';

  if(rank == size-1) fout << 0;
  fout.close();

  MPI_Send(0, 0, MPI_CHAR, right, tag, MPI_COMM_WORLD);
  
  MPI_Finalize();
  return 0;
}