#include "cash-display.hpp"
#include "cellular-automata.hpp"

#include <fstream>
#include <memory>
#include <string_view>

#ifndef NEWCA_H
#define NEWCA_H

void singleRun();
class Molecule;

class newCA {
private:
  const unsigned nrow{};
  const unsigned ncol{};

  CA2D<std::unique_ptr<Molecule>> plane;

  CashDisplay *display_p{
      nullptr}; // need to disable copy and assignment for this class!!!
                // Otherwise will violate rule of three

  std::ofstream density{"output/density.txt"};
  std::ofstream field{"output/field.txt"};
  std::ofstream kParam{"output/catalytic_param.txt"};
  std::ofstream kDistr{"output/distribution.txt"};

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
  void writeDensity(const long t,
                    int (newCA::*fcn)(int, int,
                                      int)); // writes density of all possible
                                             // replicators/complexes to a file

  void writeField(const long t); //
  void readField(
      std::string_view t, unsigned lowerB,
      unsigned upperB); // reads field into CA2D object, pass a lower and upper
                        // bound to only read part of a field (the partial
                        // field is bound by two empty columns either side)
  int testSimple(char type);  // counts number of simple p's and q's
  int testDensity(char type); // counts number of p's and q's indiscriminately
  int testComplex(int type1, int type2,
                  int cType); // counts number of a certain complex
  /* void testDummy(int x = 0); */
  void writeAverageK(const int k, const long t);
  void writeDistribution(
      const long t); // writes frequency distribution of all 8 k parameters
  void plane_to_display(); // letting display_p put pixels into plane
  void reallycba();

  /* actual simulation */
  void decay(Molecule *mole, unsigned dis = 0); // dis=1 for dissociation
  void diffuse(unsigned row, unsigned col);
  int determineComplex(const double myFate, double &cumuProb, Molecule *mole,
                       Molecule *someNei, const int mole_type);
  void formingComplex(int complex, Molecule *mole, int neiNum,
                      Molecule *someNei);
  void formingComplex(int complex, Molecule *mole,
                      Molecule *someNeiWM); // overload for well-mixed for now
  void replication(Molecule *mole, Molecule *someNei);
  void linearMutation(Molecule *someNei);
  void exponentialMutation(Molecule *someNei);
  void update_squares();

  friend void singleRun();

  ~newCA();
};

#endif
