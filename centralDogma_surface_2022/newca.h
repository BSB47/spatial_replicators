#include "cellular-automata.hpp"
#include "molecule.h"

#ifndef NEWCA_H
#define NEWCA_H

class newca {
private:
  CA2D<Molecule> plane;

  enum cashColors {
    blue,  // p
    red,   // q
    white, // s
  };

public:
  newca(const unsigned nrow, const unsigned ncol);
  ~newca();
};

#endif
