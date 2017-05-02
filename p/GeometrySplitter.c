#include "GeometrySplitter.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int widthX; //Width of board
int widthY; //Length of board

int smallestSide;   //Smallest side of the board
int largestSide;    //Largest side of the board (not used)

int numXDivisions;  //Number of partitions along the x axis (typical case)
int numYDivisions;  //Number of partitions along the y axis (typical case)

int extraPartitions; // Number of extra partitions need to be dealt with after factoring
bool isSmallestSideVertical;    //Which axis is the smallest side
int numberOfPartitions; //Number of actual partitions
struct partition* partitions;

void compareSides (void);    //Function prototyping

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

/* Given the number of partitions, determines the two closest factors after partition size is (potentially) massaged upwards to create a nicely divisible partition number */
void determineFactors (void) {
  int i;
  int root;
  int total;
  int highestFactor;
  int otherFactor;

  total = numberOfPartitions + extraPartitions;

  root = sqrt(total);

  for(i = 1; i <= root; i++)
    if(total % i == 0)
      highestFactor = i;

  otherFactor = total / highestFactor;

  //If these factors are unable to map to the sides
  if((highestFactor > largestSide || otherFactor > smallestSide) && (highestFactor > smallestSide || otherFactor > largestSide)) {
    if(smallestSide > total/smallestSide)
      setDivisions(smallestSide, total/smallestSide);
    else
      setDivisions(total/smallestSide, smallestSide);
  } else if (otherFactor > numberOfPartitions) {
    setDivisions(numberOfPartitions, 1);
  } else {
    setDivisions(highestFactor, otherFactor);
  }
}

void allocatepartitions() {
  partitions = malloc(sizeof(struct partition)*numberOfPartitions);
}

void deallocatepartitions() {
  free(partitions);
}

int findMissingYSpace(int row) {
  int i;
  int result;
  result = 0;

  for(i = 0; (row - numXDivisions * i) >= 0; i++) {
    result += partitions[row - numXDivisions * i].lengthY;
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
  extraPartitions = 0; //Range of 0 - (smallestSide - 1)

  if (numberOfPartitions > largestSide) {
    while((numberOfPartitions % smallestSide != 0) && (numberOfPartitions % largestSide != 0))
      numberOfPartitions--;
  }

  determineFactors();

  xLength = widthX / numXDivisions;
  yLength = widthY / numYDivisions;
  extraX = (widthX % numXDivisions);
  extraY = (widthY % numYDivisions);

  tempExtraX = extraX;
  totalX = 0;

  for (i = 0; i < numberOfPartitions ; i++) {
    if (i == 0) {
      partitions[i].startX = 0;
      partitions[i].startY = 0;

    } else if (totalX % widthX == 0) {
      tempExtraX = extraX;
      partitions[i].startX = 0;
      partitions[i].startY = partitions[i-1].lengthY + partitions[i-1].startY;

      if(extraY > 0) {
        extraY--;
      }
    } else {
      partitions[i].startX = partitions[i-1].startX + partitions[i-1].lengthX;
      partitions[i].startY = partitions[i-1].startY;
    }

    totalX += xLength;
    partitions[i].lengthX = xLength;
    partitions[i].lengthY = yLength;

    if(extraY > 0)
      partitions[i].lengthY += 1;

    if(tempExtraX > 0) {
      partitions[i].lengthX += 1;
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

/* Partitions the board and returns a struct array of the partitions. Modifies processes as it sees fit */
struct partition *generateBoard(int width, int length, int *processes) {
  widthX = width;
  widthY = length;
  if(*processes > (width * length))
    *processes = width * length;
  numberOfPartitions = *processes;

  compareSides();
  splitGeometry();

  *processes = numberOfPartitions;
  return partitions;
}