#include "cash-display.hpp"
#include "molecule.h"
#include "para.h"

#include <iostream>
#include <iterator>
#include <random>
#include <vector>

namespace DiceRoller {
std::random_device rd;
std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
std::mt19937 twister{ss};
}; // namespace DiceRoller

void singleRun() {

  std::vector<Molecule> plane(
      Para::sys_nrow *
      Para::sys_ncol); // Using a 1D array to represent the 2D array in order to
                       // avoid double-dereferencing when passing array to
                       // function. The size of the 1D array is (nrow * ncol).
  std::uniform_int_distribution typeInitializer{1, 3};
  for (int i = 0; i < std::size(plane); i++) {
    switch (typeInitializer(DiceRoller::twister)) {
    case 1:
      plane[i] = Molecule(i, Molecule::p);
      break;
    case 2:
      plane[i] = Molecule(i, Molecule::q);
      break;
    default:
      plane[i] = Molecule(i, Molecule::s);
      break;
    }
    std::cout << plane[i].getTypeReplicator() << 'z';
  }

  delete display_p;
}

int main() {
  singleRun();
  return 0;
}
