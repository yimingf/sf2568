#include <stdio.h>    

struct pt {
  int numFox; // foxes
  int numRab; // rabbits
  int numVeg; // vegetation.
};

int main() {
  struct pt p;      
  printf("%d %d %d\n", sizeof(p.numFox), sizeof(p), sizeof(char));
  return(0);
}