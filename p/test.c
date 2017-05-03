#include <stdio.h>    
#include <stdlib.h>

int main() {
  int i = 0;
  while(1) {
    int j = rand()%1460;
    if (j == 0) {
      printf("%d\n", i);
      return(0);
    }
    i++;
  }
  return(0);
}