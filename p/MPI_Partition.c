#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "GeometrySplitter.h"

int id;
int numPartitions;
struct partition* partitions;
struct pt* localBoard; // 2
struct pt* nextGenBoard;
struct pt* masterBoard;

int masterBoardCol;
int masterBoardRow;
struct partition myCoords;

int numberOfGenerations;
int myNeighborIDs[4]; // 1
MPI_Status lastStatus;
MPI_Request lastRequest;

int localBoardSize;

struct pt** allocatedMemory; // 2
int numMalloc;

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

bool isAlive(int x, int y);
void swapBoards();

struct pt getArray(int x, int y) { // 2
  return localBoard[x+(y*(myCoords.lengthX+4))]; // 1
}

void setReproduce(int x, int y, int f, int r, int v) { // 2
  nextGenBoard[x+(y*(myCoords.lengthX+4))].numFox = f; // 1
  nextGenBoard[x+(y*(myCoords.lengthX+4))].numRab = r;
  nextGenBoard[x+(y*(myCoords.lengthX+4))].numVeg = v;
}

void setMigration(int x, int y, int f, int r) { // 2
  nextGenBoard[x+(y*(myCoords.lengthX+4))].numFox += f;
  nextGenBoard[x+(y*(myCoords.lengthX+4))].numRab += r;
}

void freeMemory() {
  for(int i = 0; i < numMalloc; i++)
    free(allocatedMemory[i]);
}

void parseFile(int argc, char** argv) {
  int numProcesses;

  int counter;
  int f, r, v;

  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  counter = 0;

  FILE * filePtr = fopen(argv[1], "r");
  if(filePtr == NULL) {
    printf("Could not find file! Please restart and retry.");
    exit(1);
  }

  fscanf(filePtr, "%d", &numberOfGenerations);
  fscanf(filePtr, "%d", &masterBoardCol);
  fscanf(filePtr, "%d", &masterBoardRow);

  numPartitions = numProcesses;

  masterBoard = malloc(masterBoardRow*masterBoardCol*sizeof(struct pt)); // 2
  while(counter < (masterBoardRow*masterBoardCol)) { // 2
    fscanf(filePtr, "%d", &f);
    fscanf(filePtr, "%d", &r);
    fscanf(filePtr, "%d", &v);
    masterBoard[counter].numFox = f;
    masterBoard[counter].numRab = r;
    masterBoard[counter].numVeg = v;

    counter++;
  }
  fclose(filePtr);
}

void finalizeBoard() {
  int size;

  if(!id) {
    struct pt* incomingBoard; // 2

    for(int i = 0; i < numPartitions; i++) {
      size = (partitions[i].lengthX+4)*(partitions[i].lengthY+4); // 1
      incomingBoard = malloc(sizeof(struct pt)*size); // 2
      MPI_Recv(incomingBoard, sizeof(struct pt)*size, MPI_BYTE, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);

      for(int k = 0; k < partitions[i].lengthY; k++) {
        for(int j = 0; j < partitions[i].lengthX; j++) {
          int equivLocation = ((k+2)*(partitions[i].lengthX+4)+(j+2)); // 1
          masterBoard[j+partitions[i].startX+((k+partitions[i].startY)*masterBoardCol)] = incomingBoard[equivLocation];
        }
      }
      free(incomingBoard);
    }

    printf("\nFinal board configuration: \n"); // 2
    for(int i = 0; i < masterBoardCol*masterBoardRow; i++) {
      printf("%d ", masterBoard[i].numFox);
      if((i + 1) % (masterBoardCol) == 0)
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
  int* neighbors;
  int currentBoardSize;
  struct pt* currentBoard; // 2
  int l;
  int numProcessors;

  parseFile(argc, argv);
  partitions = generateBoard(masterBoardCol, masterBoardRow, &numPartitions);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcessors);
  printf("Forcing %d partitions\n", numPartitions);
  numMalloc = numPartitions;
  allocatedMemory = malloc(sizeof(struct pt*)*(numMalloc+8)); // 2

  for(int i = 0; i < numPartitions; i++) {
    neighbors = neighborList(i);
    MPI_Isend(&partitions[i], sizeof(struct partition), MPI_BYTE, i, BOUNDS_MESSAGE, MPI_COMM_WORLD, &lastRequest);
    MPI_Isend(neighbors, 4, MPI_INT, i, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 1
    currentBoardSize = (partitions[i].lengthX+4)*(partitions[i].lengthY+4);
    currentBoard = malloc(sizeof(struct pt) *currentBoardSize);

    l = 0;
    for (int k = partitions[i].startY-2; k<=(partitions[i].startY+partitions[i].lengthY+1); k++) { // 1
      for(int j = partitions[i].startX-2; j<=(partitions[i].startX+partitions[i].lengthX+1); j++) {
        if((k<0)||(j<0)||(k>=masterBoardRow)||(j>=masterBoardCol)) { // 2
          currentBoard[l].numFox = 0; // boundary
          currentBoard[l].numRab = 0;
          currentBoard[l].numVeg = 0;
          l++;
        } else { // 2
          currentBoard[l].numFox = masterBoard[j+(k*masterBoardCol)].numFox;
          currentBoard[l].numRab = masterBoard[j+(k*masterBoardCol)].numRab;
          currentBoard[l].numVeg = masterBoard[j+(k*masterBoardCol)].numVeg;
          l++;
        }
      }
    }

    allocatedMemory[i] = currentBoard;
    MPI_Isend(currentBoard, sizeof(struct pt)*currentBoardSize, MPI_BYTE, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 2
  }

  printf("\nInitial board: \n"); // 2
  for(int i = 0; i < masterBoardCol * masterBoardRow; i++) {
    printf("%d ", masterBoard[i].numFox);
    if((i + 1) % (masterBoardCol) == 0)
      printf("\n");
  }

  for (int i = 0; i < numProcessors; i++)
    MPI_Isend(&numPartitions, 1, MPI_INT, i, PARTITION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
  for (int i = 0; i < numProcessors; i++)
    MPI_Isend(&numberOfGenerations, 1, MPI_INT, i, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
}

void calculateBoard() {
  int* memoryArray[4]; // 1
  int memoryUsed;

  int* sendEdge; // 2
  int sendSize;
  int recvSize;
  int tag;

  int* recvEdge; // 2
  int currentNeighbor;

  while (numberOfGenerations-- > 0) { // each generation.
    memoryUsed = 0;

    // TODO: implement calculations and generate a migration array.

    // communication.
    if (id < numPartitions) {
      for (int j = 0; j < 4; j++) { // 1
        if (myNeighborIDs[j] > -1) {
          switch (j) {
            case 0: // W
              sendEdge = malloc(3*sizeof(int)*myCoords.lengthY); // 2
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k] = localBoard[(2*myCoords.lengthX+10)+k*(myCoords.lengthX+4)].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+myCoords.lengthY] = localBoard[(2*myCoords.lengthX+11)+k*(myCoords.lengthX+4)].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+2*myCoords.lengthY] = localBoard[(2*myCoords.lengthX+11)+k*(myCoords.lengthX+4)].numRab; // r
              }
              sendSize = 3*myCoords.lengthY;
              tag = E_UPDATE;
              break;
            case 1: // E
              sendEdge = malloc(3*sizeof(int)*myCoords.lengthY);
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k] = localBoard[(3*myCoords.lengthX+8)+k*(myCoords.lengthX+4)].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+myCoords.lengthY] = localBoard[(3*myCoords.lengthX+8)+k*(myCoords.lengthX+4)].numRab; // r
              }
              for (int k = 0; k < myCoords.lengthY; k++) {
                sendEdge[k+2*myCoords.lengthY] = localBoard[(3*myCoords.lengthX+9)+k*(myCoords.lengthX+4)].numFox; // r
              }
              sendSize = 3*myCoords.lengthY;
              tag = W_UPDATE;
              break;
            case 2: // N
              sendEdge = malloc(3*sizeof(int)*myCoords.lengthX);
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k] = localBoard[(myCoords.lengthX+4)*2+2+k].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*3+2+k].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+2*myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*3+2+k].numRab; // r
              }
              sendSize = 3*myCoords.lengthX;
              tag = S_UPDATE;
              break;
            case 3: // S
              sendEdge = malloc(3*sizeof(int)*myCoords.lengthX);
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k] = localBoard[(myCoords.lengthX+4)*myCoords.lengthY+2+k].numFox; // f
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*myCoords.lengthY+2+k].numRab; // r
              }
              for (int k = 0; k < myCoords.lengthX; k++) {
                sendEdge[k+2*myCoords.lengthX] = localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+1)+2+k].numFox; // f
              }
              sendSize = 3*myCoords.lengthX;
              tag = N_UPDATE;
              break;
          }

          memoryArray[memoryUsed++] = sendEdge;
          MPI_Isend(sendEdge, sendSize, MPI_INT, myNeighborIDs[j], tag, MPI_COMM_WORLD, &lastRequest);
        }
      }

      for (int j = 0; j < 4; j++) { // 1
        currentNeighbor = myNeighborIDs[j];

        if (currentNeighbor > -1) {
          switch(j) {
            case 0:
              recvEdge = malloc(sizeof(int)*3*myCoords.lengthY);
              recvSize = 3*myCoords.lengthY;
              tag = W_UPDATE;
              break;
            case 1:
              recvEdge = malloc(sizeof(int)*3*myCoords.lengthY);
              recvSize = 3*myCoords.lengthY;
              tag = E_UPDATE;
              break;
            case 2:
              recvEdge = malloc(sizeof(int)*3*myCoords.lengthX);
              recvSize = 3*myCoords.lengthX;
              tag = N_UPDATE;
              break;
            case 3:
              recvEdge = malloc(sizeof(int)*3*myCoords.lengthX);
              recvSize = 3*myCoords.lengthX;
              tag = S_UPDATE;
              break;
          }

          MPI_Recv(recvEdge, recvSize, MPI_INT, currentNeighbor, tag, MPI_COMM_WORLD, &lastStatus);
          switch(tag) {
            case W_UPDATE:
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(2+i)*(myCoords.lengthX+4)].numFox += recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(2+i)*(myCoords.lengthX+4)].numRab += recvEdge[myCoords.lengthY+i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(2+i)*(myCoords.lengthX+4)+1].numFox += recvEdge[2*myCoords.lengthY+i];
              }
              break;
            case E_UPDATE:
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(3+i)*(myCoords.lengthX+4)-2].numFox += recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(3+i)*(myCoords.lengthX+4)-1].numFox += recvEdge[myCoords.lengthY+i];
              }
              for (int i = 0; i < myCoords.lengthY; i++) {
                localBoard[(3+i)*(myCoords.lengthX+4)-1].numRab += recvEdge[2*myCoords.lengthY+i];
              }
              break;
            case N_UPDATE:
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[2+i].numFox += recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[2+i].numRab += recvEdge[myCoords.lengthX+i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[myCoords.lengthX+6+i].numFox += recvEdge[2*myCoords.lengthX+i];
              }
              break;
            case S_UPDATE:
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+2)+2+i].numFox += recvEdge[i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+3)+2+i].numFox += recvEdge[myCoords.lengthX+i];
              }
              for (int i = 0; i < myCoords.lengthX; i++) {
                localBoard[(myCoords.lengthX+4)*(myCoords.lengthY+3)+2+i].numRab += recvEdge[2*myCoords.lengthX+i];
              }
              break;
          }

          free(recvEdge);
        }
      }

      for(int j = 2; j < (myCoords.lengthY + 2); j++) { // 1
        for(int i = 2; i < (myCoords.lengthX + 2); i++) { // 2
          setReproduce(i, j, localBoard[i+j*(myCoords.lengthY+4)].numFox, localBoard[i+j*(myCoords.lengthY+4)].numRab, localBoard[i+j*(myCoords.lengthY+4)].numVeg);
        }
      }
      swapBoards();

      for(int i = 0; i < memoryUsed; i++)
        free(memoryArray[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(id < numPartitions)
    MPI_Isend(localBoard, sizeof(struct pt)*(myCoords.lengthX+4)*(myCoords.lengthY+4), MPI_BYTE, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 1 // 2

  MPI_Barrier(MPI_COMM_WORLD);
}

// bool isAlive(int x, int y) {
//   bool alive;
//   int numNeighbors;
//   alive = false;
//   numNeighbors = 0;

//   for(int j = -1; j <= 1; j++)
//     for(int i = -1; i <= 1; i++) {
//       if(!(((x+i)<0) || ((x+i)>=(myCoords.lengthX+4)) || ((y+j)<0) || ((y+j)>=(myCoords.lengthY+4)) || (i==0 && j==0))) // 1
//         numNeighbors += getArray(x + i, y + j);
//     }

//   if((getArray(x,y) == 1) && (numNeighbors == 2 || numNeighbors == 3))
//     alive = true;
//   else if((getArray(x,y) == 0) && (numNeighbors == 3))
//     alive = true;

//   return alive;
// }

void swapBoards(void) {
  struct pt* tempBoard; // 2
  tempBoard = localBoard;
  localBoard = nextGenBoard;
  nextGenBoard = tempBoard;
}

void initMPI(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  numMalloc = 0; // garbage colllection

  if (!id) {
    MPI_Comm_size(MPI_COMM_WORLD, &numPartitions);
    initializeBoard(argc, argv);
  }
  else
    allocatedMemory = malloc(sizeof(struct pt*) * 8);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Recv(&numPartitions, 1, MPI_INT, 0, PARTITION_MESSAGE, MPI_COMM_WORLD, &lastStatus);
  MPI_Recv(&numberOfGenerations, 1, MPI_INT, 0, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastStatus);

  if(id < numPartitions) {
    MPI_Recv(&myCoords, sizeof(struct partition), MPI_BYTE, 0, BOUNDS_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    MPI_Recv(myNeighborIDs, 4, MPI_INT, 0, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastStatus); // 1
    localBoardSize = (myCoords.lengthX+4)*(myCoords.lengthY+4); // 1
    localBoard = malloc(sizeof(struct pt)*localBoardSize);
    MPI_Recv(localBoard, sizeof(struct pt)*localBoardSize, MPI_BYTE, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    nextGenBoard = malloc(sizeof(struct pt)*localBoardSize);
  }
}

int main(int argc, char ** argv) {
  initMPI(argc, argv);
  calculateBoard();
  finalizeBoard();
  return 0;
}