#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "GeometrySplitter.h"

int identity;
int actualPartitions;
struct partition *partitionArray;
char* localBoard;
char* nextGenBoard;
char* masterBoard;

int masterBoard_columns;
int masterBoard_rows;

struct partition myCoords;

int numberOfGenerations;
int myNeighborIDs[4]; // 1
MPI_Status lastStatus;
MPI_Request lastRequest;

int localBoard_Size;

char** allocatedMemory;
int numberOfMemoryAllocations;

typedef enum {
  NW_UPDATE,
  N_UPDATE,
  NE_UPDATE,
  W_UPDATE,
  E_UPDATE,
  SW_UPDATE,
  S_UPDATE,
  SE_UPDATE,
  BOUNDS_MESSAGE,
  GENERATION_MESSAGE,
  NEIGHBORLIST_MESSAGE,
  BOARD_MESSAGE,
  PARTITION_MESSAGE
} tagType;

bool isAlive(int x, int y);
void swapBoards();

char getArray(int x, int y) {
  return localBoard[x+(y*(myCoords.lengthX+4))]; // 1
}

void setNextArray(int x, int y, int value) {
  nextGenBoard[x+(y*(myCoords.lengthX+4))] = value; // 1
}

void freeMemory() {
  for(int i = 0; i < numberOfMemoryAllocations; i++)
    free(allocatedMemory[i]);
}

/* Parses the input file character by character. Assumes proper file format of ITERATIONS\bCOLUMNS\bROWS\bARRAY_STUFF. Sets the resulting data to *currentGenBoard */
void parseFile(int argc, char ** argv) {
  int numberOfProcesses;

  //Keeps track of the total number of *s and .s in the file
  int counter;

  //Keeps track of the current character being read from the file
  char currentChar;

  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
  counter = 0;

  FILE * filePtr = fopen(argv[1], "r");

  //Breaks if the user pointed to a file that doesn't exist
  if(filePtr == NULL) {
    printf("Could not find file! Please restart and retry.");
    exit(1);
  }

  fscanf(filePtr, "%d", &numberOfGenerations);
  fscanf(filePtr, "%d", &masterBoard_columns);
  fscanf(filePtr, "%d", &masterBoard_rows);

  actualPartitions = numberOfProcesses;

  // TODO: change to row*col*(struct)
  masterBoard = malloc(masterBoard_rows*masterBoard_columns*sizeof(char));
  // TODO: change to fox-rabbit configuration
  while(counter < (masterBoard_rows * masterBoard_columns)) {
    if(fscanf(filePtr, "%c", &currentChar) == EOF) {
      printf("Your file specification's jacked up, might want to check it out.");
      exit(1);
    }
    if(currentChar == 42) {
      masterBoard[counter] = 1; // TODO: new box (num_fox, num_rabbit);
      counter++;
    }
    else if(currentChar == 46) {
      masterBoard[counter] = 0;
      counter++;
    }
  }

  fclose(filePtr);
}

/* Called when the final board configurations have been calculated, all processes submit their sections for gather */
void finalizeBoard() {
  int size;

  if(!identity) {
    char* incomingBoard;

    for(int i = 0; i < actualPartitions; i++) {
      size = (partitionArray[i].lengthX+4)*(partitionArray[i].lengthY+4); // 1
      incomingBoard = malloc(sizeof(char) * size);

      MPI_Recv(incomingBoard, size, MPI_CHAR, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);

      for(int k = 0; k < partitionArray[i].lengthY; k++) {
        for(int j = 0; j < partitionArray[i].lengthX; j++) {
          int clientEquivLocation = ((k+2)*(partitionArray[i].lengthX+4)+(j+2)); // 1

          masterBoard[j+partitionArray[i].startX+((k+partitionArray[i].startY)*masterBoard_columns)] = incomingBoard[clientEquivLocation];
        }
      }

      free(incomingBoard);
    }

    printf("\nFinal board configuration: \n");

    for(int i = 0; i < masterBoard_columns*masterBoard_rows; i++) {
      if(masterBoard[i] == 1)
        printf("*");
      else
        printf(".");
      if((i + 1) % (masterBoard_columns) == 0)
        printf("\n");
    }

    freeMemory();
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Finalize();
}

void initializeBoard(int argc, char** argv) {
  int* curLoopNeighbors;
  int curLoopBoard_Size;
  char* curLoopBoard;
  int l;
  int numberOfProcessors;

  parseFile(argc, argv);

  partitionArray = generateBoard(masterBoard_columns, masterBoard_rows, &actualPartitions);

  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcessors);

  printf("Forcing %d partitions\n", actualPartitions);

  numberOfMemoryAllocations = actualPartitions;
  // TODO: change to struct pointers.
  allocatedMemory = malloc(sizeof(char*)*(numberOfMemoryAllocations+8));

  for(int i = 0; i < actualPartitions; i++) {
    // TODO: change neighbor sturcture.
    curLoopNeighbors = neighborList(i);
    MPI_Isend(&partitionArray[i], sizeof(struct partition), MPI_BYTE, i, BOUNDS_MESSAGE, MPI_COMM_WORLD, &lastRequest);
    MPI_Isend(curLoopNeighbors, 4, MPI_INT, i, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 1

    curLoopBoard_Size = (partitionArray[i].lengthX+4)*(partitionArray[i].lengthY+4);
    curLoopBoard = malloc(sizeof(char) *curLoopBoard_Size);

    l = 0;

    for(int k = partitionArray[i].startY-2; k <= (partitionArray[i].startY+partitionArray[i].lengthY+1); k++) { // 1
      for(int j = partitionArray[i].startX-2; j <= (partitionArray[i].startX+partitionArray[i].lengthX+1); j++) {
        if((k<0)||(j<0)||(k>=masterBoard_rows)||(j>=masterBoard_columns))
          curLoopBoard[l++] = 0; // boundary
        else
          curLoopBoard[l++] = masterBoard[j+(k*masterBoard_columns)];
      }
    }

    allocatedMemory[i] = curLoopBoard;
    MPI_Isend(curLoopBoard, curLoopBoard_Size, MPI_CHAR, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest);
  }

  printf("\nInitial board: \n");
  for(int i = 0; i < masterBoard_columns * masterBoard_rows; i++) {
    if(masterBoard[i] == 1)
      printf("*");
    else
      printf(".");
    if((i + 1) % (masterBoard_columns) == 0)
      printf("\n");
  }

  for (int i = 0; i < numberOfProcessors; i++)
    MPI_Isend(&actualPartitions, 1, MPI_INT, i, PARTITION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
  for (int i = 0; i < numberOfProcessors; i++)
    MPI_Isend(&numberOfGenerations, 1, MPI_INT, i, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
}

void calculateBoard() {
  char* memoryArray[4]; // 1
  int memoryUsed;

  char* sendEdge;
  int sendSize;
  int recvSize;
  int tag;

  char* recvEdge;
  int currentNeighbor;

  while (numberOfGenerations-- > 0) { // each generation.
    memoryUsed = 0;

    if (identity < actualPartitions) {
      for (int j = 0; j < 4; j++) { // 1
        if (myNeighborIDs[j] > -1) {
          switch (j) {
            case 0: // W
              sendEdge = malloc(2*sizeof(char)*myCoords.lengthY);
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k] = localBoard[(2*myCoords.lengthX+10)+k*(myCoords.lengthX+4)];
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+myCoords.lengthY] = localBoard[(2*myCoords.lengthX+11)+k*(myCoords.lengthX+4)];
              }
              sendSize = 2*myCoords.lengthY;
              tag = E_UPDATE;
              break;
            case 1: // E
              sendEdge = malloc(2*sizeof(char)*myCoords.lengthY);
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k] = localBoard[(3*myCoords.lengthX+8)+k*(myCoords.lengthX+4)];
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+myCoords.lengthY] = localBoard[(3*myCoords.lengthX+9)+k*(myCoords.lengthX+4)];
              }
              sendSize = 2*myCoords.lengthY;
              tag = W_UPDATE;
              break;
            case 2:
              sendEdge = malloc(2*sizeof(char)*myCoords.lengthX);
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k] = localBoard[(myCoords.lengthX+4)*2+2+k];
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*3+2+k];
              }
              sendSize = 2*myCoords.lengthX;
              tag = S_UPDATE;
              break;
            case 3:
              sendEdge = malloc(2*sizeof(char)*myCoords.lengthX);
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k] = localBoard[(myCoords.lengthX+4)*myCoords.lengthY+2+k];
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+1)+2+k];
              }
              sendSize = 2*myCoords.lengthX;
              tag = N_UPDATE;
              break;
          }

          memoryArray[memoryUsed++] = sendEdge;
          MPI_Isend(sendEdge, sendSize, MPI_CHAR, myNeighborIDs[j], tag, MPI_COMM_WORLD, &lastRequest);
        }
      }

      for (int j = 0; j < 4; j++) { // 1
        currentNeighbor = myNeighborIDs[j];

        if (currentNeighbor > -1) {
          switch(j) {
            case 0:
              recvEdge = malloc(sizeof(char)*2*myCoords.lengthY);
              recvSize = 2*myCoords.lengthY;
              tag = W_UPDATE;
              break;
            case 1:
              recvEdge = malloc(sizeof(char)*2*myCoords.lengthY);
              recvSize = 2*myCoords.lengthY;
              tag = E_UPDATE;
              break;
            case 2:
              recvEdge = malloc(sizeof(char)*2*myCoords.lengthX);
              recvSize = 2*myCoords.lengthX;
              tag = N_UPDATE;
              break;
            case 3:
              recvEdge = malloc(sizeof(char)*2*myCoords.lengthX);
              recvSize = 2*myCoords.lengthX;
              tag = S_UPDATE;
              break;
          }

          MPI_Recv(recvEdge, recvSize, MPI_CHAR, currentNeighbor, tag, MPI_COMM_WORLD, &lastStatus);

          switch(tag) {
            case N_UPDATE:
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[2+i] = recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[myCoords.lengthX+6+i] = recvEdge[myCoords.lengthX+i];
              }
              break;
            case W_UPDATE:
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(2+i)*(myCoords.lengthX+4)] = recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(2+i)*(myCoords.lengthX+4)+1] = recvEdge[myCoords.lengthY+i];
              }
              break;
            case E_UPDATE:
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(3+i)*(myCoords.lengthX+4)-2] = recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(3+i)*(myCoords.lengthX+4)-1] = recvEdge[myCoords.lengthY+i];
              }
              break;
            case S_UPDATE:
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+2)+2+i] = recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+3)+2+i] = recvEdge[myCoords.lengthX+i];
              }
              break;
          }

          free(recvEdge);
        }
      }

      for(int j = 0; j < (myCoords.lengthY + 4); j++) { // 1
        for(int i = 0; i < (myCoords.lengthX + 4); i++) {
          if(isAlive(i, j))
            setNextArray(i, j, 1);
          else
            setNextArray(i, j, 0);
        }
      }

      swapBoards();

      for(int i = 0; i < memoryUsed; i++)
        free(memoryArray[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(identity < actualPartitions)
    MPI_Isend(localBoard, (myCoords.lengthX+4)*(myCoords.lengthY+4), MPI_CHAR, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 1

  MPI_Barrier(MPI_COMM_WORLD);
}

bool isAlive(int x, int y) {
  bool alive;
  int numNeighbors;
  alive = false;
  numNeighbors = 0;

  for(int j = -1; j <= 1; j++)
    for(int i = -1; i <= 1; i++) {
      if(!(((x+i)<0) || ((x+i)>=(myCoords.lengthX+4)) || ((y+j)<0) || ((y+j)>=(myCoords.lengthY+4)) || (i==0 && j==0))) // 1
        numNeighbors += getArray(x + i, y + j);
    }

  if((getArray(x,y) == 1) && (numNeighbors == 2 || numNeighbors == 3))
    alive = true;
  else if((getArray(x,y) == 0) && (numNeighbors == 3))
    alive = true;

  return alive;
}

void swapBoards(void) {
  char* tempBoard;
  tempBoard = localBoard;
  localBoard = nextGenBoard;
  nextGenBoard = tempBoard;
}

void initMPI(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &identity);

  numberOfMemoryAllocations = 0; // garbage colllection

  if (!identity) {
    MPI_Comm_size(MPI_COMM_WORLD, &actualPartitions);
    initializeBoard(argc, argv);
  }
  else
    allocatedMemory = malloc(sizeof(char*) * 8);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Recv(&actualPartitions, 1, MPI_INT, 0, PARTITION_MESSAGE, MPI_COMM_WORLD, &lastStatus);
  MPI_Recv(&numberOfGenerations, 1, MPI_INT, 0, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastStatus);

  if(identity < actualPartitions) {
    MPI_Recv(&myCoords, sizeof(struct partition), MPI_BYTE, 0, BOUNDS_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    MPI_Recv(myNeighborIDs, 4, MPI_INT, 0, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastStatus); // 1
    localBoard_Size = (myCoords.lengthX+4)*(myCoords.lengthY+4); // 1
    localBoard = malloc(sizeof(char) * localBoard_Size);
    MPI_Recv(localBoard, localBoard_Size, MPI_CHAR, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    nextGenBoard = malloc(sizeof(char) * localBoard_Size);
  }
}

int main(int argc, char ** argv) {
  initMPI(argc, argv);

  calculateBoard();

  finalizeBoard();
  return 0;
}

