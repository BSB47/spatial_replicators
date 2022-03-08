#include "molecule.h"
#include "para.h"
#include <iostream>
#include <iterator>
#include <random>

namespace Random {
std::random_device rd;
std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
std::mt19937 twister{ss};
}; // namespace Random

void singleRun() {
  Molecule
      plane[(Para::sys_nrow * Para::sys_ncol) -
            1]{}; // Using a 1D array to represent the 2D array in order to
                  // avoid double-dereferencing when passing array to function.
                  // The size of the 1D array is (nrow * ncol) - 1.
  std::uniform_int_distribution typeInitializer{1, 3};
  for (int i = 0; i < std::size(plane); i++) {
    switch (typeInitializer(Random::twister)) {
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
    std::cout << plane[i].getTypeReplicator() << ' ';
  }
}

int main() { void singleRun(); }
