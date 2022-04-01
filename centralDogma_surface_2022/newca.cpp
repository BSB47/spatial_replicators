#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"
#include "para.h"
#include "random.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <random>

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  for (unsigned row = 1; row <= nrow; row++) {
    for (unsigned col = 1; col <= ncol; col++) {
      switch (DiceRoller::typeInitializer(DiceRoller::twister)) {
      case 1:
        plane.cell(row, col).setTypeRep(Molecule::p);
        ++pqDensity;
        break;
      case 2:
        plane.cell(row, col).setTypeRep(Molecule::q);
        ++pqDensity;
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
  /* std::cout << "pq denstiy is currently: " << pqDensity << std::endl; */

  std::vector<CashPanelInfo> panel_info(1);
  panel_info[CA].n_row = nrow;
  panel_info[CA].n_col = ncol;
  panel_info[CA].o_row = 0;
  panel_info[CA].o_col = 0;

  display_p = new CashDisplay( // display_p is a pointer to the CASH window
      Para::sys_nrow, Para::sys_ncol, panel_info,
      Para::scale); // window_row/col are sys_nrow/col for testing.
  display_p->color_rgb(blue, 0, 0, 255);
  display_p->color_rgb(red, 255, 0, 0);
  display_p->color_rgb(white, 255, 255, 255);
}

void newCA::visualize(const long t) {
  display_p->open_window();
  display_p->open_png();

  if (t % Para::display_interval == 0) {
    plane_to_display();
    display_p->draw_window();
    display_p->draw_png();
  } else
    return;
}

void newCA::writeFile(const long t) {
  if (!output.is_open()) {
    output.open("output/output.txt", std::ios::app);
    /* std::cout << t << "opening output.txt" << std::endl; */
  }

  if (!output) {
    std::cerr << "Could not access output.txt\n";
  }

  output << t * Para::alpha << ' '
         << static_cast<double>(pqDensity) / Para::grid_size << '\n';
  output.close();
}

void newCA::testDensity(const long t) {
  std::uint32_t testDensity{};
  for (unsigned row{1}; row <= Para::sys_nrow; row++) {
    for (unsigned col{1}; col <= Para::sys_ncol; col++) {
      if (plane.cell(row, col).getTypeReplicator() != Molecule::s)
        ++testDensity;
    }
  }
  /* if (!output.is_open()) { */
  /*   output.open("output.txt", std::ios::app); */
  /*   output << "testing now!!! " << t << '\t' << testDensity; */
  /* } */
  std::cout << testDensity << std::endl;
}

void newCA::printDensity() { std::cout << pqDensity << std::endl; }

void newCA::plane_to_display() { // Does not paint display! Just transports
                                 // plane data to it.
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      switch (plane.cell(row, col).getTypeReplicator()) {
      case Molecule::p:
        display_p->put_pixel(CA, row, col, blue);
        break;
      case Molecule::q:
        display_p->put_pixel(CA, row, col, red);
        break;
      case Molecule::s:
        display_p->put_pixel(CA, row, col, white);
        break;
      }
    }
  }
}

void newCA::decay(Molecule &mole) {
  mole.setTypeRep(Molecule::s);
  --pqDensity;
}

void newCA::diffuse(Molecule &mole, unsigned row, unsigned col) {
  Molecule &someNei{
      plane.neigh_wrap(row, col, DiceRoller::randomNei(DiceRoller::twister))};

  Molecule tmp{std::move_if_noexcept(someNei)}; // make a copy of chosen nei
  someNei = std::move_if_noexcept(mole);
  mole = std::move_if_noexcept(tmp);
}

/* The following function returns a number which can be used to decipher if and
what complex forms between the chosen molecule (passed by reference to this
function) and a molecule in its neighbourhood (randomly chosen through this
function).
There are 3 kinds of numbers this function can return:
0-7: complex has formed and mole is the catalyst, someNei is the template
100-107: complex has formed and someNei is the catalyst, mole is the template
>=1000: complex did not form or someNei is S */
int newCA::determineComplex(unsigned row, unsigned col, Molecule &mole,
                            int mole_type) {
  const double complexProb{DiceRoller::probabilityGen(DiceRoller::twister)};
  double cumuProb{};

  Molecule &someNei{
      plane.neigh_wrap(row, col, DiceRoller::randomNei(DiceRoller::twister))};

  switch (mole_type) {
  case 0:
    switch (someNei.getTypeReplicator()) { // if mole & someNei are both P
    case 0:
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += mole.m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += someNei.m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1000; // complex did not form between two P molecules
    case 1:        // if mole is P and someNei is Q
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += mole.m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += someNei.m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1001; // complex did not form between P and Q
    default:
      return 1005; // someNei is a S
    }
  case 1:
    switch (someNei.getTypeReplicator()) { // if mole is Q and someNei is P
    case 0:
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += mole.m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += someNei.m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1002; // complex did not form between Q and P
    case 1:        // if both are Q
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += mole.m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += someNei.m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1003; // complex did not form between 2 Qs
    default:
      return 1005; // someNei is a S
    }
  default:
    return 1010; // something is broken

    /* switch (mole_type) { // 0=p, 1=q, 2=s */
    /* case 0: */
    /*   switch (someNei.getTypeReplicator()) { // 0=p, 1=q, 2=s */
    /*   case 0: // p1p2, where p1 is mole and the catalyst */
    /*     cumuProb += Para::beta * mole.m_kppp; */
    /*     if (complexProb <= cumuProb && mole.m_typeComp == 0 && */
    /*         someNei.m_typeComp == 0) { */
    /*       mole.m_typeComp = Molecule::cata; */
    /*       someNei.m_typeComp = Molecule::temp; */
    /*       someNei.m_repOutcome = Molecule::p; */
    /*       return; */
    /*     } else if (complexProb <= (cumuProb + mole.m_kppq) && */
    /*                mole.m_typeComp == 0 && someNei.m_typeComp == 0) { */
    /*       mole.m_typeComp = Molecule::cata; */
    /*       someNei.m_typeComp = Molecule::temp; */
    /*       someNei.m_repOutcome = Molecule::q; */
    /*     } */
    /*   } */
    /* } */
  }
}

void newCA::update_squares() {
  unsigned row{};
  unsigned col{};

  for (int i{1}; i <= Para::grid_size; i++) { // # of times is counted from 1
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);

    /* Molecule *mole{&(plane.cell(row, col))}; */
    Molecule &mole{(plane.cell(row, col))};
    auto mole_type{mole.getTypeReplicator()};
    const double myFate{DiceRoller::probabilityGen(DiceRoller::twister)};

    if (mole_type != Molecule::s) {
      if (myFate <= (Para::alpha * Para::decay_probability)) {
        decay(mole);
      } else if (myFate <= (Para::alpha * (Para::decay_probability +
                                           Para::diffusion_probability))) {
        diffuse(mole, row, col);
      } else {
        determineComplex(row, col, mole, mole_type);
      }
    }

    /* If chosen replicator is non-s, choose a neighbour at random; then
     * change this neighbour's type to be the same as the chosen
     * replicator's. */
    /* Molecule &someNei{plane.neigh_wrap( */
    /*     row, col, DiceRoller::randomNei(DiceRoller::twister))}; */
    /* if (someNei.getTypeReplicator() == Molecule::s) { */
    /*   someNei.setTypeRep(mole_type); */
    /* else do nothing */
    /* } */
  }
}

newCA::~newCA() {
  if (display_p)
    delete display_p;
  display_p = nullptr;
}
