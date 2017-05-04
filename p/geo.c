/*
  kth sf2568 parpro-17 (parallel computations for large-scale problems) project.
  this skeleton from: https://github.com/wallstop/MPI-Game-of-Life
*/

#include "geo.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int wx, wy, sideSmall, sideLarge, numDivX, numDivY, extraPartitions;
bool isSmallestSideVertical;
int numParts;
struct partition* partitions;

void compareSides (void);

void setDivisions (int foo, int bar) {
  bool inPlace;
  int tempFactor;

  if (((sideLarge%foo == 0) || (sideSmall%bar == 0)) && (foo<=sideLarge && bar<=sideSmall)) {
    inPlace = true;
  } else if (((sideSmall%foo == 0) || (sideLarge%bar == 0)) && (foo<=sideSmall && bar<=sideLarge)) {
    inPlace = false;
  } else {
    inPlace = true;
  }

  if ((isSmallestSideVertical && !inPlace) || (!isSmallestSideVertical && inPlace)) {
    tempFactor = foo;
    foo = bar;
    bar = tempFactor;
  }

  numDivX = foo;
  numDivY = bar;
}

void determineFactors (void) {
  int i;
  int root;
  int total;
  int highestFactor;
  int otherFactor;

  total = numParts+extraPartitions;
  root = sqrt(total);

  for (i = 1; i <= root; i++) {
    if (total % i == 0) {
      highestFactor = i;
    }
  }
  otherFactor = total/highestFactor;

  if ((highestFactor > sideLarge || otherFactor > sideSmall) && (highestFactor > sideSmall || otherFactor > sideLarge)) {
    if (sideSmall > total/sideSmall)
      setDivisions(sideSmall, total/sideSmall);
    else
      setDivisions(total/sideSmall, sideSmall);
  } else if (otherFactor > numParts) {
    setDivisions(numParts, 1);
  } else {
    setDivisions(highestFactor, otherFactor);
  }
}

void allocatepartitions () {
  partitions = malloc(sizeof(struct partition)*numParts);
}

void deallocatepartitions () {
  free(partitions);
}

int findMissingYSpace (int row) {
  int i;
  int result;
  result = 0;

  for(i = 0; (row - numDivX * i) >= 0; i++) {
    result += partitions[row - numDivX * i].y1;
  }
  result = wy - result;
  return result;
}

void split (void) {
  int i, xLength, yLength, extraX, extraY, tempExtraX, totalX;
  bool foo = false;

  allocatepartitions();
  compareSides();
  extraPartitions = 0;

  if (numParts > sideLarge) {
    while((numParts % sideSmall != 0) && (numParts % sideLarge != 0))
      numParts--;
  }

  determineFactors();

  xLength = wx/numDivX;
  yLength = wy/numDivY;
  extraX = wx%numDivX;
  extraY = wy%numDivY;

  tempExtraX = extraX;
  totalX = 0;

  for (i = 0; i<numParts; i++) {
    if (i == 0) {
      partitions[i].x0 = 0;
      partitions[i].y0 = 0;
      partitions[i].isRed = false; // initial
    } else if (totalX % wx == 0) {
      tempExtraX = extraX;
      partitions[i].x0 = 0;
      partitions[i].y0 = partitions[i-1].y1 + partitions[i-1].y0;
      foo = !foo; // switch between red/black
      partitions[i].isRed = foo;
      if (extraY > 0) {
        extraY--;
      }
    } else {
      partitions[i].x0 = partitions[i-1].x0 + partitions[i-1].x1;
      partitions[i].y0 = partitions[i-1].y0;
      partitions[i].isRed = !partitions[i-1].isRed;
    }

    totalX += xLength;
    partitions[i].x1 = xLength;
    partitions[i].y1 = yLength;

    if (extraY > 0)
      partitions[i].y1 += 1;

    if (tempExtraX > 0) {
      partitions[i].x1 += 1;
      totalX++;
      tempExtraX--;
    }

    printf("part %d isRed? %d\n", i, partitions[i].isRed);
  }
}

void compareSides (void) {
  if (wx < wy) {
    sideSmall = wx;
    sideLarge = wy;
    isSmallestSideVertical = false;
  } else {
    sideSmall = wy;
    sideLarge = wx;
    isSmallestSideVertical = true;
  }
}

int* neighborList (int partitionNumber) {
  int* list;
  int positionX;
  int positionY;
  int currentPosition;

  currentPosition = 0;
  list = malloc(sizeof(int)*4); // 1. 0W 1E 2N 3S
  positionX = partitionNumber%numDivX;
  positionY = partitionNumber/numDivX;

  for (int i = -1; i <= 1; i++) { // 1X
    if (i != 0) {
      if ((positionX+i<0) || (positionX+i>=numDivX)) {
        list[currentPosition] = -1;
      } else {
        list[currentPosition] = positionX+i+positionY*numDivX;
      }
      currentPosition++;
    }
  }
  for (int i = -1; i <= 1; i++) { // 1Y
    if (i != 0) {
      if ((positionY+i<0) || (positionY+i>=numDivY)) {
        list[currentPosition] = -1;
      } else {
        list[currentPosition] = positionX+(positionY+i)*numDivX;
      }
      currentPosition++;
    }
  }
  return list;
}

struct partition *generateBoard (int w, int l, int *processes) {
  wx = w;
  wy = l;
  if (*processes>(w*l)) {
    *processes=w*l;
  }
  numParts = *processes;

  split();
  *processes = numParts;
  return partitions;
}