#ifndef GEOMETRY_H
#define GEOMETRY_H

struct partition {
  int startX;
  int startY;
  int lengthX;
  int lengthY;
};
struct pt {
  int numFox; // foxes
  int numRab; // rabbits
  int numVeg; // vegetation.
};
struct partition* generateBoard (int width, int length, int* processes);
int* neighborList (int);

#endif