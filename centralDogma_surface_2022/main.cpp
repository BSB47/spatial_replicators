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
  randomizedGrid.visualize(); // double free error from this call
  randomizedGrid.update_self_replication();
  char kill{};
  std::cout << "type q+enter to quit: ";
  std::cin >> kill;
  if (kill == 'q') {
    std::cout << "display_p is dead! quitting now...";
  }
}

int main() {
  singleRun();
  return 0;
}
