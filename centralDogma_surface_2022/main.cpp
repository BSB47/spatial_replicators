#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <iostream>
#include <iterator>
#include <random>
#include <vector>

void singleRun() {
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);
  randomizedGrid.visualize();
}

int main() {
  singleRun();
  return 0;
}
