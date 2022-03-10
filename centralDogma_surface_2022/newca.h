#include "cellular-automata.hpp"
#include "molecule.h"

#ifndef NEWCA_H
#define NEWCA_H

class newCA {
private:
  const unsigned nrow;
  const unsigned ncol;

  CA2D<Molecule> plane;

  enum cashColors {
    blue,  // p
    red,   // q
    white, // s
  };

public:
  newCA(const unsigned a_nrow, const unsigned a_ncol);

  ~newCA();
};

#endif
