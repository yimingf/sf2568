/*
  kth sf2568 parpro-17 (parallel computations for large-scale problems) project.
*/

#include "alg.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

int birthRabbit (int r, int v) {
  if (r<2) {
    return 0;
  } else if (v<200) {
    if (r>700) {
      return 2;
    } else {
      return 3;
    }
  } else if (v<500) {
    if (r>700) {
      return 3;
    } else {
      return 4;
    }
  } else if (v<800) {
    if (r>700) {
      return 4;
    } else if (r>200) {
      return 5;
    } else {
      return 6;
    }
  } else {
    if (r>5000) {
      return 5;
    } else if (r>700) {
      return 7;
    } else if (r>200) {
      return 8;
    } else {
      return 9;
    }
  }
}

int birthFox (int f, int r) {
  if (f<2) {
    return 0;
  } else if (r<3) {
    if (f>100) {
      return 0;
    } else if (f>50) {
      return 1;
    } else {
      return 2;
    }
  } else if (r<10) {
    if (f>100) {
      return 1;
    } else if (f>50) {
      return 2;
    } else {
      return 3;
    }
  } else if (r<40) {
    if (f>100) {
      return 2;
    } else if (f>10) {
      return 3;
    } else {
      return 4;
    }
  } else {
    if (f>50) {
      return 3;
    } else if (f>10) {
      return 4;
    } else {
      return 5;
    }
  }
}

int lifeSpanRabbit (int v) {
  if (v>350) {
    return 144;
  } else if (v>250) {
    return 72;
  } else if (v>150) {
    return 36;
  } else {
    return 6;
  }
}