#include "Geometry.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int widthX;
int widthY;
int smallestSide;
int largestSide;
int numXDivisions;
int numYDivisions;
int extraPartitions;
bool isSmallestSideVertical;
int numParts;
struct partition* partitions;

void compareSides (void);

void setDivisions (int factor1, int factor2) {
  bool inPlace;
  int tempFactor;

  if(((largestSide % factor1 == 0) || (smallestSide % factor2 == 0)) && (factor1 <= largestSide && factor2 <= smallestSide)) {
    inPlace = true;
  } else if(((smallestSide % factor1 == 0) || (largestSide % factor2 == 0)) && (factor1 <= smallestSide && factor2 <= largestSide)) {
    inPlace = false;
  } else {
    inPlace = true;
  }

  if((isSmallestSideVertical && !inPlace) || (!isSmallestSideVertical && inPlace)) {
    tempFactor = factor1;
    factor1 = factor2;
    factor2 = tempFactor;
  }

  numXDivisions = factor1;
  numYDivisions = factor2;
}

void determineFactors (void) {
  int i;
  int root;
  int total;
  int highestFactor;
  int otherFactor;

  total = numParts + extraPartitions;
  root = sqrt(total);

  for(i = 1; i <= root; i++)
    if(total % i == 0)
      highestFactor = i;

  otherFactor = total / highestFactor;

  if ((highestFactor > largestSide || otherFactor > smallestSide) && (highestFactor > smallestSide || otherFactor > largestSide)) {
    if(smallestSide > total/smallestSide)
      setDivisions(smallestSide, total/smallestSide);
    else
      setDivisions(total/smallestSide, smallestSide);
  } else if (otherFactor > numParts) {
    setDivisions(numParts, 1);
  } else {
    setDivisions(highestFactor, otherFactor);
  }
}

void allocatepartitions() {
  partitions = malloc(sizeof(struct partition)*numParts);
}

void deallocatepartitions() {
  free(partitions);
}

int findMissingYSpace(int row) {
  int i;
  int result;
  result = 0;

  for(i = 0; (row - numXDivisions * i) >= 0; i++) {
    result += partitions[row - numXDivisions * i].y1;
  }
  result = widthY - result;
  return result;
}

int findNumberOfCellsInRow(int row) {
  int i;
  int result;
  result = 0;

  for(i = 0; (row - numXDivisions * i) >= 0; i++) {
    result++;
  }
  return result;
}

void splitGeometry(void) {
  int i;
  int xLength;
  int yLength;
  int extraX;
  int extraY;
  int tempExtraX;
  int totalX;

  allocatepartitions();
  compareSides();
  extraPartitions = 0;

  if (numParts > largestSide) {
    while((numParts % smallestSide != 0) && (numParts % largestSide != 0))
      numParts--;
  }

  determineFactors();

  xLength = widthX / numXDivisions;
  yLength = widthY / numYDivisions;
  extraX = (widthX % numXDivisions);
  extraY = (widthY % numYDivisions);

  tempExtraX = extraX;
  totalX = 0;

  for (i = 0; i < numParts ; i++) {
    if (i == 0) {
      partitions[i].x0 = 0;
      partitions[i].y0 = 0;

    } else if (totalX % widthX == 0) {
      tempExtraX = extraX;
      partitions[i].x0 = 0;
      partitions[i].y0 = partitions[i-1].y1 + partitions[i-1].y0;

      if(extraY > 0) {
        extraY--;
      }
    } else {
      partitions[i].x0 = partitions[i-1].x0 + partitions[i-1].x1;
      partitions[i].y0 = partitions[i-1].y0;
    }

    totalX += xLength;
    partitions[i].x1 = xLength;
    partitions[i].y1 = yLength;

    if(extraY > 0)
      partitions[i].y1 += 1;

    if(tempExtraX > 0) {
      partitions[i].x1 += 1;
      totalX++;
      tempExtraX--;
    }
  }
}

void compareSides(void) {
  int smallResult;
  int largeResult;

  if(widthX < widthY) {
    smallResult = widthX;
    largeResult = widthY;
    isSmallestSideVertical = false;
  } else {
    smallResult = widthY;
    largeResult = widthX;
    isSmallestSideVertical = true;
  }

  smallestSide = smallResult;
  largestSide = largeResult;
}

int* neighborList(int partitionNumber) {
  int* list;
  int positionX;
  int positionY;
  int currentPosition;

  currentPosition = 0;
  list = malloc(sizeof(int)*4); // 1. 0W 1E 2N 3S
  positionX = partitionNumber%numXDivisions;
  positionY = partitionNumber/numXDivisions;

  for (int i = -1; i <= 1; i++) { // 1X
    if (i != 0) {
      if ((positionX+i<0) || (positionX+i>=numXDivisions)) {
        list[currentPosition] = -1;
      } else {
        list[currentPosition] = positionX+i+positionY*numXDivisions;
      }
      currentPosition++;
    }
  }
  for (int i = -1; i <= 1; i++) { // 1Y
    if (i != 0) {
      if ((positionY+i<0) || (positionY+i>=numYDivisions)) {
        list[currentPosition] = -1;
      } else {
        list[currentPosition] = positionX+(positionY+i)*numXDivisions;
      }
      currentPosition++;
    }
  }
  return list;
}

struct partition *generateBoard(int width, int length, int *processes) {
  widthX = width;
  widthY = length;
  if(*processes > (width * length))
    *processes = width * length;
  numParts = *processes;

  compareSides();
  splitGeometry();

  *processes = numParts;
  return partitions;
}