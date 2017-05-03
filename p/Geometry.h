#ifndef GEOMETRY_H
#define GEOMETRY_H

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
struct partition* generateBoard (int width, int length, int* processes);
int* neighborList (int);

#endif