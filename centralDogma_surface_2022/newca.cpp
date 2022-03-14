#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"
#include "para.h"
#include "random.h"

#include <iostream>
#include <random>

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  for (unsigned row = 0; row < nrow; row++) {
    for (unsigned col = 0; col < ncol; col++) {
      switch (DiceRoller::typeInitializer(DiceRoller::twister)) {
      case 1:
        plane.cell(row, col).setTypeRep(Molecule::p);
        break;
      case 2:
        plane.cell(row, col).setTypeRep(Molecule::q);
        break;
      default:
        plane.cell(row, col).setTypeRep(Molecule::s);
        break;
      }
      /* std::cout << plane.cell(row, col).getTypeReplicator() << ' '; tested,
       * creation and population of grid/cells seems to work. However, seems to
       * require inclusion of newCA.cpp in testCA.cpp (Why? Would we have to do
       * this for main.cpp?) -- RESOLVED By using inline variables in random.h*/
    }
  }
  for (unsigned row = 0; row < 10; row++) {
    for (unsigned col = 0; col < 10; col++) {
      std::cout << plane.cell(row, col).getTypeReplicator() << ' ';
    }
  }

  std::vector<CashPanelInfo> panel_info(1);
  panel_info[0].n_row = nrow;
  panel_info[0].n_col = ncol;
  panel_info[0].o_row = 0;
  panel_info[0].o_col = 0;

  display_p = new CashDisplay( // display_p is a pointer to the CASH window
      Para::sys_nrow, Para::sys_ncol,
      panel_info); // window_row/col are sys_nrow/col for testing.
  display_p->color_rgb(blue, 0, 0, 255);
  display_p->color_rgb(red, 255, 0, 0);
  display_p->color_rgb(white, 255, 255, 255);

  display_p->open_window();
  display_p->open_png();
}

void newCA::visualize() {
  plane_to_display();
  display_p->draw_window();
  display_p->draw_png();
}

void newCA::plane_to_display() {
  for (unsigned row{0}; row < nrow; row++) {
    for (unsigned col{0}; col < ncol; col++) {
      switch (plane.cell(row, col).getTypeReplicator()) {
      case Molecule::p:
        display_p->put_pixel(1, row, col, blue);
        break;
      case Molecule::q:
        display_p->put_pixel(1, row, col, red);
        break;
      case Molecule::s:
        display_p->put_pixel(1, row, col, white);
        break;
      }
    }
  }
}

void newCA::update_self_replication() {
  unsigned row{};
  unsigned col{};
  for (int i{0}; i < Para::grid_size; i++) {
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);

    /* If some replicator is non-s, choose a neighbour at random; then change
     * this neighbour's type to be the same as the chosen
     * replicator's. */
    if (plane.cell(row, col).getTypeReplicator() != Molecule::s) {
      Molecule &someNei{plane.neigh_wrap(
          row, col, DiceRoller::randomNei(DiceRoller::twister))};
      if (someNei.getTypeReplicator() == Molecule::s) {
        someNei.setTypeRep(plane.cell(row, col).getTypeReplicator());
      }
      /* else do nothing */
    } else
      continue;
  }
}

newCA::~newCA() {
  if (display_p)
    delete display_p;
  display_p = nullptr;
}
