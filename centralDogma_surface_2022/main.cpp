#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <utility>
#include <vector>

void singleRun() {
  long init_time{0};
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);

  for (long t{init_time}; t <= 2000; t++) {
    /* if (t % Para::display_interval == 0) { */
    /* randomizedGrid.visualize(); */
    /* } */
    randomizedGrid.writeFile(t, 0, &newCA::testComplex);
    randomizedGrid.update_squares();
  }
}

int main() {
  singleRun();
  return 0;
}
