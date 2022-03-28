#include "cash-display.hpp"
#include "cellular-automata.hpp"
#include "molecule.h"

#include <fstream>

#ifndef NEWCA_H
#define NEWCA_H

class newCA {
private:
  const unsigned nrow{};
  const unsigned ncol{};

  std::uint32_t pqDensity{};
  CA2D<Molecule> plane;

  CashDisplay *display_p{nullptr};

  std::ofstream output{"output/output.txt"};

  enum cashColors {
    blue,  // p
    red,   // q
    white, // s
  };

  enum Panels {
    CA,
    proportion_p,
  };

public:
  /* display and outupt */
  newCA(const unsigned a_nrow, const unsigned a_ncol);
  void visualize(const long t); // initializing and opening window/png
  void writeFile(const long t); // writing density of replicators to a file
  void testDensity(const long t);
  void printDensity();
  /* function to test if the pqDensity
   * incre-/decrements are working.*/
  void plane_to_display(); // letting display_p put pixels into plane

  /* actual simulation */
  void decay(Molecule &mole);
  void diffuse(Molecule &mole, unsigned row, unsigned col);
  void complexFormation(Molecule &mole);
  void update_squares();

  ~newCA();
};

#endif
