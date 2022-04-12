#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <utility>
#include <vector>

void singleRun() {
  long init_time{0};
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);

  for (long t{init_time}; t <= Para::max_time; t++) {
    if (t % Para::display_interval == 0) {
      randomizedGrid.visualize();
    }
    randomizedGrid.writeFile(t, 0, &newCA::testComplex);
    randomizedGrid.update_squares();
    /* if (t == 50) { */
    /*   randomizedGrid.testDensity(t); */
    /* } */
  }
  char kill{};
  std::cout << "type q+enter to quit: "; // Needs to go!
  std::cin >> kill;
  if (kill == 'q') {
    std::cout << "display_p is dead! quitting now...";
  }
}

int main() {
  singleRun();
  /* newCA testguy{10, 10}; */
  /* Molecule tester{Molecule::p, Molecule::free}; */
  /* std::cout << tester.getTypeReplicator(); */
  /* testguy.decayRoll(tester); */
  /* std::cout << tester.getTypeReplicator(); */
  return 0;
}
