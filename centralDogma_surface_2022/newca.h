#include "cash-display.hpp"
#include "cellular-automata.hpp"

#include <fstream>
#include <memory>

#ifndef NEWCA_H
#define NEWCA_H

class Molecule;

class newCA {
private:
  const unsigned nrow{};
  const unsigned ncol{};

  CA2D<std::unique_ptr<Molecule>> plane;

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
  void visualize(const long t);  // initializing and opening window/png
  void writeFile(const long t);  // writing density of replicators to a file
  int testDensity(const long t); // counts number of p's and q's
  int testComplex(const long t, int c); // counters number of a certain complex
  void plane_to_display(); // letting display_p put pixels into plane

  /* actual simulation */
  void decay(Molecule *mole);
  void diffuse(std::unique_ptr<Molecule> &mole, unsigned row, unsigned col);
  int determineComplex(Molecule *mole, Molecule *someNei, int mole_type);
  void formingComplex(int complex, Molecule *mole, int neiNum,
                      Molecule *someNei);
  void update_squares();

  ~newCA();
};

#endif
