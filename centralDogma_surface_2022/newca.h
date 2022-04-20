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

  CashDisplay *display_p{
      nullptr}; // need to disable copy and assignment for this class!!!
                // Otherwise will violate rule of three

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
  void visualize(); // initializing and opening window/png
  void writeFile(const long t, int c,
                 int (newCA::*fcn)(int)); // writing density of a certain
                                          // replicator/complex to a file

  int testSimple(int c);   // counts number of simple p's and q's
  int testDensity(int c);  // counts number of p's and q's indiscriminately
  int testComplex(int c);  // counts number of a certain complex
  void plane_to_display(); // letting display_p put pixels into plane

  /* actual simulation */
  void decay(Molecule *mole);
  void diffuse(unsigned row, unsigned col);
  int determineComplex(const double myFate, double &cumuProb, Molecule *mole,
                       Molecule *someNei, const int mole_type);
  void formingComplex(int complex, Molecule *mole, int neiNum,
                      Molecule *someNei);
  void formingComplex(int complex, Molecule *mole,
                      Molecule *someNeiWM); // overload for well-mixed for now
  void update_squares();

  ~newCA();
};

#endif
