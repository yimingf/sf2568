#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdbool.h>

typedef enum {
  N_UPDATE,
  W_UPDATE,
  E_UPDATE,
  S_UPDATE,
  BOUNDS_MESSAGE,
  GENERATION_MESSAGE,
  NEIGHBORLIST_MESSAGE,
  BOARD_MESSAGE,
  PARTITION_MESSAGE
} tagType;

struct partition {
  int x0;
  int y0;
  int x1;
  int y1;
};
struct pt {
  int numFox; // foxes
  int numRab; // rabbits
  int numVeg; // vegetation.
};
struct partition* generateBoard (int w, int l, int* processes);
int* neighborList (int);

#endif