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
  std::vector<CashPanelInfo> panel_info{
      1}; // 1 panel for visualization of the CA.
  int win_row{Para::sys_nrow + 2};
  int win_col{Para::sys_ncol + 2};
  panel_info[0].n_row = Para::sys_nrow;
  panel_info[0].n_col = Para::sys_ncol;
  panel_info[0].o_row = Para::cash_margin;
  panel_info[0].o_col = Para::cash_margin;

  CashDisplay *display_p;
  display_p = new CashDisplay(win_row, win_col, panel_info);

  display_p->color_rgb(0, 0, 0, 0);
  display_p->color_rgb(1, 255, 255, 255);
  display_p->color_rgb(2, 20, 20, 20);
  display_p->color_rgb(3, 255, 255, 255);
  display_p->color_rgb(4, 0, 0, 255);     // blue
  display_p->color_rgb(5, 0, 200, 0);     // dark green
  display_p->color_rgb(6, 255, 0, 0);     // red
  display_p->color_rgb(7, 188, 143, 143); // brown
  display_p->color_rgb(8, 255, 255, 255);
  display_p->color_rgb(9, 0, 0, 0);
  display_p->color_rgb(10, 150, 150, 150);

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
}

int main() {
  singleRun();
  return 0;
}
