#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include "Geometry.h"

int id;
int numPartitions;
struct partition* partitions;
struct pt* myBoard; // 2
struct pt* nextBoard;
struct pt* masterBoard;

int masterBoardCol;
int masterBoardRow;
struct partition conf;

int numIter;
int myNeighborIDs[4]; // 1
MPI_Status lastStatus;
MPI_Request lastRequest;

int myBoardSize;

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
  return myBoard[x+y*conf.lengthX];
}

void setReproduce(int x, int y, int f, int r, int v) { // 2
  nextBoard[x+y*conf.lengthX].numFox = f;
  nextBoard[x+y*conf.lengthX].numRab = r;
  nextBoard[x+y*conf.lengthX].numVeg = v;
}

void setMigration(int x, int y, int f, int r) { // 2
  nextBoard[x+y*conf.lengthX].numFox += f;
  nextBoard[x+y*conf.lengthX].numRab += r;
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

  fscanf(filePtr, "%d", &numIter);
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
      size = partitions[i].lengthX*partitions[i].lengthY;
      incomingBoard = malloc(sizeof(struct pt)*size); // 2
      MPI_Recv(incomingBoard, sizeof(struct pt)*size, MPI_BYTE, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);

      for(int k = 0; k < partitions[i].lengthY; k++) {
        for(int j = 0; j < partitions[i].lengthX; j++) {
          int equivLocation = j+k*partitions[i].lengthX; // 1
          masterBoard[j+partitions[i].startX+((k+partitions[i].startY)*masterBoardCol)] = incomingBoard[equivLocation];
        }
      }
      free(incomingBoard);
    }

    printf("\nFinal board confuration: \n"); // 2
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
    MPI_Isend(neighbors, 4, MPI_INT, i, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastRequest);
    currentBoardSize = partitions[i].lengthX*partitions[i].lengthY;
    currentBoard = malloc(sizeof(struct pt) *currentBoardSize);

    l = 0;
    for (int k = partitions[i].startY; k < partitions[i].startY+partitions[i].lengthY; k++) {
      for(int j = partitions[i].startX; j < partitions[i].startX+partitions[i].lengthX; j++) {
        currentBoard[l].numFox = masterBoard[j+k*masterBoardCol].numFox;
        currentBoard[l].numRab = masterBoard[j+k*masterBoardCol].numRab;
        currentBoard[l].numVeg = masterBoard[j+k*masterBoardCol].numVeg;
        l++;
      }
    }

    allocatedMemory[i] = currentBoard;
    MPI_Isend(currentBoard, sizeof(struct pt)*currentBoardSize, MPI_BYTE, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 2
  }

  printf("\nInitial board: \n"); // 2
  for(int i = 0; i < masterBoardCol*masterBoardRow; i++) {
    printf("%d ", masterBoard[i].numFox);
    if((i+1)%(masterBoardCol) == 0)
      printf("\n");
  }

  for (int i = 0; i < numProcessors; i++)
    MPI_Isend(&numPartitions, 1, MPI_INT, i, PARTITION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
  for (int i = 0; i < numProcessors; i++)
    MPI_Isend(&numIter, 1, MPI_INT, i, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastRequest);
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

  while (numIter-- > 0) { // each generation.
    memoryUsed = 0;

    // TODO: implement calculations and generate a migration array.

    // communication.
    if (id < numPartitions) {
      for (int j = 0; j < 4; j++) { // 1
        if (myNeighborIDs[j] > -1) {
          switch (j) {
            case 0: // W
              sendEdge = malloc(3*sizeof(int)*conf.lengthY); // 2
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k] = myBoard[k*conf.lengthX].numFox; // f
              }
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k+conf.lengthY] = myBoard[1+k*conf.lengthX].numFox; // f
              }
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k+2*conf.lengthY] = myBoard[1+k*conf.lengthX].numRab; // r
              }
              sendSize = 3*conf.lengthY;
              tag = E_UPDATE;
              break;
            case 1: // E
              sendEdge = malloc(3*sizeof(int)*conf.lengthY);
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k] = myBoard[(k+1)*conf.lengthX-2].numFox;
              }
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k+conf.lengthY] = myBoard[(k+1)*conf.lengthX-2].numRab; // r
              }
              for (int k = 0; k < conf.lengthY; k++) {
                sendEdge[k+2*conf.lengthY] = myBoard[(k+1)*conf.lengthX-1].numFox; // r
              }
              sendSize = 3*conf.lengthY;
              tag = W_UPDATE;
              break;
            case 2: // N
              sendEdge = malloc(3*sizeof(int)*conf.lengthX);
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k] = myBoard[k].numFox; // f
              }
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k+conf.lengthX] = myBoard[conf.lengthX+k].numFox; // f
              }
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k+2*conf.lengthX] = myBoard[conf.lengthX+k].numRab; // r
              }
              sendSize = 3*conf.lengthX;
              tag = S_UPDATE;
              break;
            case 3: // S
              sendEdge = malloc(3*sizeof(int)*conf.lengthX);
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k] = myBoard[(conf.lengthY-2)*conf.lengthX+k].numFox; // f
              }
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k+conf.lengthX] = myBoard[(conf.lengthY-2)*conf.lengthX+k].numRab; // r
              }
              for (int k = 0; k < conf.lengthX; k++) {
                sendEdge[k+2*conf.lengthX] = myBoard[(conf.lengthY-1)*conf.lengthX+k].numFox; // f
              }
              sendSize = 3*conf.lengthX;
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
              recvEdge = malloc(sizeof(int)*3*conf.lengthY);
              recvSize = 3*conf.lengthY;
              tag = W_UPDATE;
              break;
            case 1:
              recvEdge = malloc(sizeof(int)*3*conf.lengthY);
              recvSize = 3*conf.lengthY;
              tag = E_UPDATE;
              break;
            case 2:
              recvEdge = malloc(sizeof(int)*3*conf.lengthX);
              recvSize = 3*conf.lengthX;
              tag = N_UPDATE;
              break;
            case 3:
              recvEdge = malloc(sizeof(int)*3*conf.lengthX);
              recvSize = 3*conf.lengthX;
              tag = S_UPDATE;
              break;
          }

          MPI_Recv(recvEdge, recvSize, MPI_INT, currentNeighbor, tag, MPI_COMM_WORLD, &lastStatus);
          switch(tag) {
            case W_UPDATE:
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[i*conf.lengthX].numFox += recvEdge[i];
              }
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[i*conf.lengthX].numRab += recvEdge[conf.lengthY+i];
              }
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[1+i*conf.lengthX].numFox += recvEdge[2*conf.lengthY+i];
              }
              break;
            case E_UPDATE:
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[(i+1)*conf.lengthX-2].numFox += recvEdge[i];
              }
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[(i+1)*conf.lengthX-1].numFox += recvEdge[conf.lengthY+i];
              }
              for (int i = 0; i < conf.lengthY; i++) {
                myBoard[(i+1)*conf.lengthX-1].numRab += recvEdge[2*conf.lengthY+i];
              }
              break;
            case N_UPDATE:
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[i].numFox += recvEdge[i];
              }
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[i].numRab += recvEdge[conf.lengthX+i];
              }
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[i+conf.lengthX].numFox += recvEdge[2*conf.lengthX+i];
              }
              break;
            case S_UPDATE:
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[(conf.lengthY-2)*conf.lengthX+i].numFox += recvEdge[i];
              }
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[(conf.lengthY-1)*conf.lengthX+i].numFox += recvEdge[conf.lengthX+i];
              }
              for (int i = 0; i < conf.lengthX; i++) {
                myBoard[(conf.lengthY-1)*conf.lengthX+i].numRab += recvEdge[2*conf.lengthX+i];
              }
              break;
          }

          free(recvEdge);
        }
      }

      for(int j = 0; j < conf.lengthY; j++) { // 1
        for(int i = 0; i < conf.lengthX; i++) { // 2
          setReproduce(i, j, myBoard[i+j*conf.lengthX].numFox, myBoard[i+j*conf.lengthX].numRab, myBoard[i+j*conf.lengthX].numVeg);
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
    MPI_Isend(myBoard, sizeof(struct pt)*conf.lengthX*conf.lengthY, MPI_BYTE, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest); // 1 // 2

  MPI_Barrier(MPI_COMM_WORLD);
}

// bool isAlive(int x, int y) {
//   bool alive;
//   int numNeighbors;
//   alive = false;
//   numNeighbors = 0;

//   for(int j = -1; j <= 1; j++)
//     for(int i = -1; i <= 1; i++) {
//       if(!(((x+i)<0) || ((x+i)>=(conf.lengthX+4)) || ((y+j)<0) || ((y+j)>=(conf.lengthY+4)) || (i==0 && j==0))) // 1
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
  tempBoard = myBoard;
  myBoard = nextBoard;
  nextBoard = tempBoard;
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
  MPI_Recv(&numIter, 1, MPI_INT, 0, GENERATION_MESSAGE, MPI_COMM_WORLD, &lastStatus);

  if(id < numPartitions) {
    MPI_Recv(&conf, sizeof(struct partition), MPI_BYTE, 0, BOUNDS_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    MPI_Recv(myNeighborIDs, 4, MPI_INT, 0, NEIGHBORLIST_MESSAGE, MPI_COMM_WORLD, &lastStatus); // 1
    myBoardSize = conf.lengthX*conf.lengthY; // 1
    myBoard = malloc(sizeof(struct pt)*myBoardSize);
    MPI_Recv(myBoard, sizeof(struct pt)*myBoardSize, MPI_BYTE, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);
    nextBoard = malloc(sizeof(struct pt)*myBoardSize);
  }
}

int main(int argc, char ** argv) {
  initMPI(argc, argv);
  calculateBoard();
  finalizeBoard();
  return 0;
}