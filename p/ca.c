 /* Cellular automata - Sequential version */
 
 #include <stdio.h>
 #include <stdlib.h>
 #define N 62
 #define MAX_ITERATION N //Number of iteration
 
 void init_mat (int *A, int length);
 void writemat (int *A);
 int f (int x1,int x2,int x3,int x4);
 
 int main (int argc, char **argv) {
   int* A=(int*)malloc(N*N*sizeof(int));
   int* B=(int*)malloc((N-2)*(N-2)*sizeof(int));
   int i,j,k=0;
   init_mat (A,N);
   while(k++<MAX_ITERATION) {
   //Calculation of the new matrix in a temporary matrix
     for (i=0;i<N-2;i++)
       for (j=0;j<N-2;j++)
         *(B+i*N+j)=f(*(A+(i+1)*N+j),*(A+(i+1)*N+j+2),*(A+i*N+j+1),*(A+(i+2)*N+j+1));
     //Copy of the new matrix
     for (i=1;i<N-1;i++)
       for (j=1;j<N-1;j++)
         *(A+i*N+j)=*(B+(i-1)*N+j-1);
     //Update of the border
     for (i=1;i<N-1;i++) {
       *(A+i)=*(A+(N-2)*N+i);
       *(A+(N-1)*N+i)=*(A+N+i);
       *(A+i*N)=*(A+i*N+N-2);
       *(A+i*N+N-1)=*(A+i*N+1);
     }
   }
   writemat(A);
   return 0;
 }
 
 void init_mat (int *A, int length) {
   int i,j;
   for (i=0;i<length;i++)
     for (j=0;j<length;j++)
       *(A+i*length+j)=i+j>length/2&&i+j<length*3/2;
   //Each border of the matrix will have the same value than the border of his opposite side
   for (i=1;i<length-1;i++) {
     *(A+i)=*(A+(length-2)*length+i);
     *(A+(length-1)*length+i)=*(A+length+i);
     *(A+i*length)=*(A+i*length+length-2);
     *(A+i*length+length-1)=*(A+i*length+1);
   }
 }
 
 void writemat (int *A) {
   int i,j;
   for (i = 0; i < N; i++) {
     for (j = 0; j < N; j++)
       printf ("%d ", *(A + i*N +j));
     printf("\n");
   }
 }
 
 int f (int x1,int x2,int x3,int x4) {
   return x1+x2+x3+x4>2;
 }