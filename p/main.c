/*
  kth sf2568 parpro-17 (parallel computations for large-scale problems) project.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <math.h>

#include "geo.h"
#include "alg.h"

int id, numIter, numPartitions, myBoardSize, myNeighborIDs[4];
int f, r, v, f_p, r_p, lifeRabbit;
int masterBoardCol, masterBoardRow, numMalloc;
struct partition *partitions, conf;
struct pt *myBoard, *nextBoard, *masterBoard, **allocatedMemory;
MPI_Status lastStatus;
MPI_Request lastRequest;

int reproduceRateRabbit = 63;
int reproduceRateFox = 180;
double growthRateVeg = 1.1;
int lifeFox = 460; // 1460 (4 yrs) is too long for long simulation.

void swapBoards();

void setBeforeMig(int x, int y, int f, int r, int v) { // 2
  if (f<0) { // boundary check.
    nextBoard[x+y*conf.x1].numFox = 0;
  } else if (f>65535) {
    nextBoard[x+y*conf.x1].numFox = 65535;
  } else {
    nextBoard[x+y*conf.x1].numFox = f;
  }

  if (r<0) {
    nextBoard[x+y*conf.x1].numRab = 0;
  } else if (r>65535) {
    nextBoard[x+y*conf.x1].numRab = 65535;
  } else {
    nextBoard[x+y*conf.x1].numRab = r;
  }
  
  nextBoard[x+y*conf.x1].numVeg = v;
}

void setMig(int x, int y, int f, int r) { // 2
  nextBoard[x+y*conf.x1].numFox += f;
  nextBoard[x+y*conf.x1].numRab += r;
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
    printf("file not found");
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
      size = partitions[i].x1*partitions[i].y1;
      incomingBoard = malloc(sizeof(struct pt)*size);
      
      MPI_Recv(incomingBoard, sizeof(struct pt)*size, MPI_BYTE, i, BOARD_MESSAGE, MPI_COMM_WORLD, &lastStatus);
      for(int k = 0; k < partitions[i].y1; k++) {
        for(int j = 0; j < partitions[i].x1; j++) {
          int equivLocation = j+k*partitions[i].x1; // 1
          masterBoard[j+partitions[i].x0+((k+partitions[i].y0)*masterBoardCol)] = incomingBoard[equivLocation];
        }
      }
      free(incomingBoard);
    }
    
    printf("\nFinal board config:\n"); // 2
    for(int i = 0; i < masterBoardCol*masterBoardRow; i++) {
      printf("%d ", masterBoard[i].numFox);
      if((i + 1) % (masterBoardCol) == 0)
        printf("\n");
    }
    for(int i = 0; i < masterBoardCol*masterBoardRow; i++) {
      printf("%d ", masterBoard[i].numRab);
      if((i + 1) % (masterBoardCol) == 0)
        printf("\n");
    }
    for(int i = 0; i < masterBoardCol*masterBoardRow; i++) {
      printf("%d ", masterBoard[i].numVeg);
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
    currentBoardSize = partitions[i].x1*partitions[i].y1;
    currentBoard = malloc(sizeof(struct pt) *currentBoardSize);

    l = 0;
    for (int k = partitions[i].y0; k < partitions[i].y0+partitions[i].y1; k++) {
      for(int j = partitions[i].x0; j < partitions[i].x0+partitions[i].x1; j++) {
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

void calculateBoard (void) {
  int currentNeighbor, day, mem;
  int *migrationFox, *migrationRabbit;
  int *sendEdge, *recvEdge, sendSize, recvSize, tag;
  int *memArray[4];
  bool birthFlagRabbit, birthFlagFox;
  
  migrationFox = malloc(sizeof(int)*(conf.x1+4)*(conf.y1+4));
  migrationRabbit = malloc(sizeof(int)*(conf.x1+2)*(conf.y1+2));
  for (day = 0; day<numIter; day++) {
    for (int j = 0; j<(conf.y1+4); j++) {
      for (int i = 0; i<(conf.x1+4); i++) {
        migrationFox[i+j*(conf.x1+4)] = 0;
      }
    }
    for (int j = 0; j<(conf.y1+2); j++) {
      for (int i = 0; i<(conf.x1+2); i++) {
        migrationRabbit[i+j*(conf.x1+2)] = 0;
      }
    }
    mem = 0;
    birthFlagRabbit = false;
    birthFlagFox = false;

    // 1. birthing day of rabbits
    if (day%reproduceRateRabbit == 0) {
      birthFlagRabbit = true;
    }
    // 2. birthing day of foxes.
    if (day%reproduceRateFox == 0) {
      birthFlagFox = true;
    }

    for(int j = 0; j < conf.y1; j++) { // before migration.
      for(int i = 0; i < conf.x1; i++) {
        f = myBoard[i+j*conf.x1].numFox;
        r = myBoard[i+j*conf.x1].numRab;
        v = myBoard[i+j*conf.x1].numVeg;

        if (birthFlagRabbit) { // 1.
          r += birthRabbit(r, v);
        }
        if (birthFlagFox) { // 2.
          f += birthFox(f, r);
        }

        // 3. vegetation.
        v = (int)(floor(growthRateVeg*v))-r;
        if (v<100) {
          v = 100;
        }
        if (v>1000) {
          v = 1000;
        }

        // 4. predator (foxes)
        if (v<600) {
          r_p = (int)(floor(f*4/7));
          if (r_p>r) {
            r_p = (int)(floor(f*2/7));
          }
        } else {
          r_p = (int)(floor(f*2/7));
        }
        if (r_p>r || r==0) {
          r = 0;
          f = (int)(floor(0.9*f)); // starve
        } else {
          r -= r_p;
        }
        // 5. foxes' natural death
        if (f>0) {
          if (rand()%lifeFox == 0) {
            f--;
          }
        }

        // 6. rabbits' starve and natural death
        if (r>0) {
          if (v>350) {
            lifeRabbit = 216;
          } else if (v>250) {
            lifeRabbit = 144;
          } else if (v>150) {
            lifeRabbit = 72;
          } else {
            lifeRabbit = 36;
          }

          if (rand()%lifeRabbit == 0) {
            r--;
          }
        }
        setBeforeMig(i, j, f, r, v);
      }
    }

    for (int j = 0; j<conf.y1; j++) {
      for (int i = 0; i<conf.x1; i++) {
        f = myBoard[i+j*conf.x1].numFox;
        r = myBoard[i+j*conf.x1].numRab;

        f_p = r_p = 0;
        r_p = ((int)floor(r*0.1));
        if (r_p != 0) {
          r_p = rand()%r_p;
        }
        f_p = ((int)floor(f*0.12));
        if (f_p != 0) {
          f_p = rand()%f_p;
        }

        int foo = 0;
        for (int k = -2; k<3; k++) { // fox
          if (k!=0) {
            if ((conf.x0+i+k>=0) && (conf.x0+i+k<masterBoardCol)) { // x
              migrationFox[(i+2+k)+(j+2)*(conf.x1+4)] += f_p;
              foo++;
            }
            if ((conf.y0+i+k>=0) && (conf.y0+i+k<masterBoardRow)) { // y
              migrationFox[(i+2)+(j+2+k)*(conf.x1+4)] += f_p;
              foo++;
            }
          }
        }
        f_p *= (-foo);

        foo = 0;
        for (int k = -1; k<2; k++) { // rabbit
          if (k!=0) {
            if ((conf.x0+i+k>=0) && (conf.x0+i+k<masterBoardCol)) { // x
              migrationRabbit[(i+1+k)+(j+1)*(conf.x1+2)] += r_p;
              foo++;
            }
            if ((conf.y0+i+k>=0) && (conf.y0+i+k<masterBoardRow)) { // y
              migrationRabbit[(i+1)+(j+1+k)*(conf.x1+2)] += r_p;
              foo++;
            }
          }
        }
        r_p *= (-foo);
        setMig(i, j, f_p, r_p);
      }
    }

    if (id < numPartitions) {
      for (int j = 0; j < 4; j++) { // 1
        if (myNeighborIDs[j] > -1) {
          switch (j) {
            case 0: // W
              sendEdge = malloc(3*sizeof(int)*conf.y1); // 2
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k] = migrationFox[(k+2)*(conf.x1+4)]; // f
              }
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k+conf.y1] = migrationFox[1+(k+2)*(conf.x1+4)]; // f
              }
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k+2*conf.y1] = migrationRabbit[(k+1)*(conf.x1+2)]; // r
              }
              sendSize = 3*conf.y1;
              tag = E_UPDATE;
              break;
            case 1: // E
              sendEdge = malloc(3*sizeof(int)*conf.y1);
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k] = migrationFox[(k+3)*(conf.x1+4)-2];
              }
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k+conf.y1] = migrationRabbit[(k+2)*(conf.x1+2)-1]; // r
              }
              for (int k = 0; k < conf.y1; k++) {
                sendEdge[k+2*conf.y1] = migrationFox[(k+3)*(conf.x1+4)-2]; // f
              }
              sendSize = 3*conf.y1;
              tag = W_UPDATE;
              break;
            case 2: // N
              sendEdge = malloc(3*sizeof(int)*conf.x1);
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k] = migrationFox[k+2]; // f
              }
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k+conf.x1] = migrationFox[k+6+conf.x1]; // f
              }
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k+2*conf.x1] = migrationRabbit[k+1]; // r
              }
              sendSize = 3*conf.x1;
              tag = S_UPDATE;
              break;
            case 3: // S
              sendEdge = malloc(3*sizeof(int)*conf.x1);
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k] = migrationFox[(k+2)+(conf.y1+2)*(conf.x1+4)]; // f
              }
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k+conf.x1] = migrationRabbit[(k+1)+(conf.y1+1)*(conf.x1+2)]; // r
              }
              for (int k = 0; k < conf.x1; k++) {
                sendEdge[k+2*conf.x1] = migrationFox[(k+2)+(conf.y1+3)*(conf.x1+4)]; // f
              }
              sendSize = 3*conf.x1;
              tag = N_UPDATE;
              break;
          }

          memArray[mem++] = sendEdge;
          MPI_Isend(sendEdge, sendSize, MPI_INT, myNeighborIDs[j], tag, MPI_COMM_WORLD, &lastRequest);
        }
      }

      for (int j = 0; j < 4; j++) { // 1
        currentNeighbor = myNeighborIDs[j];

        if (currentNeighbor > -1) {
          switch(j) {
            case 0:
              recvEdge = malloc(sizeof(int)*3*conf.y1);
              recvSize = 3*conf.y1;
              tag = W_UPDATE;
              break;
            case 1:
              recvEdge = malloc(sizeof(int)*3*conf.y1);
              recvSize = 3*conf.y1;
              tag = E_UPDATE;
              break;
            case 2:
              recvEdge = malloc(sizeof(int)*3*conf.x1);
              recvSize = 3*conf.x1;
              tag = N_UPDATE;
              break;
            case 3:
              recvEdge = malloc(sizeof(int)*3*conf.x1);
              recvSize = 3*conf.x1;
              tag = S_UPDATE;
              break;
          }

          MPI_Recv(recvEdge, recvSize, MPI_INT, currentNeighbor, tag, MPI_COMM_WORLD, &lastStatus);
          switch(tag) {
            case W_UPDATE:
              for (int i = 0; i < conf.y1; i++) {
                migrationFox[2+i*(conf.x1+4)] += recvEdge[i];
              }
              for (int i = 0; i < conf.y1; i++) {
                migrationRabbit[1+i*(conf.x1+2)] += recvEdge[conf.y1+i];
              }
              for (int i = 0; i < conf.y1; i++) {
                migrationFox[3+i*(conf.x1+4)] += recvEdge[2*conf.y1+i];
              }
              break;
            case E_UPDATE:
              for (int i = 0; i < conf.y1; i++) {
                migrationFox[(i+3)*conf.x1-4] += recvEdge[i];
              }
              for (int i = 0; i < conf.y1; i++) {
                migrationFox[(i+3)*conf.x1-3] += recvEdge[conf.y1+i];
              }
              for (int i = 0; i < conf.y1; i++) {
                migrationRabbit[(i+2)*conf.x1-2] += recvEdge[2*conf.y1+i];
              }
              break;
            case N_UPDATE:
              for (int i = 0; i < conf.x1; i++) {
                migrationFox[i+2+2*(conf.x1+4)] += recvEdge[i];
              }
              for (int i = 0; i < conf.x1; i++) {
                migrationRabbit[i+1+(conf.x1+2)] += recvEdge[conf.x1+i];
              }
              for (int i = 0; i < conf.x1; i++) {
                migrationFox[i+2+3*(conf.x1+4)] += recvEdge[2*conf.x1+i];
              }
              break;
            case S_UPDATE:
              for (int i = 0; i < conf.x1; i++) {
                migrationFox[i+2+conf.y1*(conf.x1+4)] += recvEdge[i];
              }
              for (int i = 0; i < conf.x1; i++) {
                migrationFox[i+2+(conf.y1+1)*(conf.x1+4)] += recvEdge[conf.x1+i];
              }
              for (int i = 0; i < conf.x1; i++) {
                migrationRabbit[i+1+conf.y1*(conf.x1+2)] += recvEdge[2*conf.x1+i];
              }
              break;
          }

        }
      }

      for (int j = 0; j<conf.y1; j++) {
        for (int i = 0; i<conf.x1; i++) {
          setMig(i, j, migrationFox[(i+2)+(j+2)*(conf.x1+4)], migrationRabbit[(i+1)+(j+1)*(conf.x1+2)]);
        }
      }
      swapBoards();

      if (id== 0 && day%1000 == 0) {
        printf("%d days has gone.\n", day);
      }

      for(int i = 0; i < mem; i++) {
        free(memArray[i]);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  free(migrationFox);
  free(migrationRabbit);
  if(id < numPartitions) {
    MPI_Isend(myBoard, sizeof(struct pt)*conf.x1*conf.y1, MPI_BYTE, 0, BOARD_MESSAGE, MPI_COMM_WORLD, &lastRequest);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void swapBoards (void) {
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
    myBoardSize = conf.x1*conf.y1; // 1
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