#include "cash-display.hpp"
#include "cellular-automata.hpp"
#include "molecule.h"

#ifndef NEWCA_H
#define NEWCA_H

class newCA {
private:
  const unsigned nrow{};
  const unsigned ncol{};

  CA2D<Molecule> plane;

  CashDisplay *display_p{};

  enum cashColors {
    blue,  // p
    red,   // q
    white, // s
  };

public:
  newCA(const unsigned a_nrow, const unsigned a_ncol);
  void visualize();
  void plane_to_display();

  ~newCA();
};

#endif
