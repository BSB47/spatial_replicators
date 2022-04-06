#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"
#include "para.h"
#include "random.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  for (unsigned row = 1; row <= nrow; row++) {
    for (unsigned col = 1; col <= ncol; col++) {
      auto tmp{std::make_unique<Molecule>()};
      plane.cell(row, col) = std::move(tmp);
      switch (DiceRoller::typeInitializer(DiceRoller::twister)) {
      case 1:
        plane.cell(row, col)->setTypeRep(Molecule::p);
        break;
      case 2:
        plane.cell(row, col)->setTypeRep(Molecule::q);
        break;
      default:
        plane.cell(row, col)->setTypeRep(Molecule::s);
        break;
      }
      /* std::cout << plane.cell(row, col).getTypeReplicator() << ' '; tested,
       * creation and population of grid/cells seems to work. However, seems
       * to require inclusion of newCA.cpp in testCA.cpp (Why? Would we have
       * to do this for main.cpp?) -- RESOLVED By using inline variables in
       * random.h*/
    }
  }

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
  display_p->open_window();
  display_p->open_png();
}

void newCA::visualize(const long t) {
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
         << static_cast<double>(testDensity(t)) / Para::grid_size << '\n';
  output.close();
}

int newCA::testDensity(const long t) {
  std::uint32_t testDensity{};
  for (unsigned row{1}; row <= Para::sys_nrow; row++) {
    for (unsigned col{1}; col <= Para::sys_ncol; col++) {
      if (plane.cell(row, col)->getTypeReplicator() != Molecule::s)
        ++testDensity;
    }
  }
  return testDensity;
}

void newCA::plane_to_display() { // Does not paint display! Just transports
                                 // plane data to it.
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      switch (plane.cell(row, col)->getTypeReplicator()) {
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

void newCA::decay(Molecule *mole) { mole->setTypeRep(Molecule::s); }

void newCA::diffuse(std::unique_ptr<Molecule> &mole_ptr, unsigned row,
                    unsigned col) {
  mole_ptr.swap(
      plane.neigh_wrap(row, col, DiceRoller::randomNei(DiceRoller::twister)));
}

/* The following function returns a number which can be used to decipher if and
what complex forms between the chosen molecule (passed by reference to this
function) and a molecule in its neighbourhood (randomly chosen through this
function).
There are 3 kinds of numbers this function can return:
0-7: complex has formed and mole is the catalyst, someNei is the template
100-107: complex has formed and someNei is the catalyst, mole is the template
>=1000: complex did not form or someNei is S */
int newCA::determineComplex(Molecule *mole, Molecule *someNei, int mole_type) {
  const double complexProb{DiceRoller::probabilityGen(DiceRoller::twister)};
  double cumuProb{};

  switch (mole_type) {
  case 0:
    switch (someNei->getTypeReplicator()) { // if mole & someNei are both P
    case 0:
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += mole->m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += someNei->m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1000; // complex did not form between two P molecules
    case 1:        // if mole is P and someNei is Q
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += mole->m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += someNei->m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1001; // complex did not form between P and Q
    default:
      return 1005; // someNei is a S
    }
  case 1:
    switch (someNei->getTypeReplicator()) { // if mole is Q and someNei is P
    case 0:
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += mole->m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += someNei->m_rateList[i];
        if (complexProb <= cumuProb)
          return (i + 100);
      }
      /* std::cout << complexProb << ' ' << cumuProb << std::endl; */
      return 1002; // complex did not form between Q and P
    case 1:        // if both are Q
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += mole->m_rateList[i];
        if (complexProb <= cumuProb)
          return i;
      }
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += someNei->m_rateList[i];
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

void newCA::formingComplex(int complex, Molecule *mole, int neiNum,
                           Molecule *someNei) {
  mole->bon_nei = neiNum;
  someNei->bon_nei = plane.neighbor_converter(neiNum);
  switch (complex) {
  case 0:
  case 3:
  case 4:
  case 7:
    mole->m_typeComp = Molecule::cata;
    someNei->m_typeComp = Molecule::tempP;
    std::cout << "complex formed" << std::endl;
    return;
  case 1:
  case 2:
  case 5:
  case 6:
    mole->m_typeComp = Molecule::cata;
    someNei->m_typeComp = Molecule::tempQ;
    std::cout << "complex formed" << std::endl;
    return;
  case 100:
  case 103:
  case 104:
  case 107:
    someNei->m_typeComp = Molecule::cata;
    mole->m_typeComp = Molecule::tempP;
    std::cout << "complex formed" << std::endl;
    return;
  case 101:
  case 102:
  case 105:
  case 106:
    someNei->m_typeComp = Molecule::cata;
    mole->m_typeComp = Molecule::tempQ;
    std::cout << "complex formed" << std::endl;
    return;
  case 1000 ... 1005:
    std::cout << "complex didn't form" << std::endl;
    return;
  case 1010:
    std::cout << "broken" << std::endl;
    throw -1;
  }
}

void newCA::update_squares() {
  unsigned row{};
  unsigned col{};

  for (int i{1}; i <= Para::grid_size; i++) { // # of times is counted from 1
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);

    std::unique_ptr<Molecule> &mole_ptr{(plane.cell(row, col))};
    Molecule *mole{plane.cell(row, col).get()};
    auto mole_type{mole->getTypeReplicator()};

    unsigned neiNum{
        mole->bon_nei
            ? plane.neigh_7_select(
                  DiceRoller::randomNei(DiceRoller::randomNei),
                  mole->bon_nei) // exclude nei if mole is already in a complex
            : static_cast<unsigned int>(
                  DiceRoller::randomNei(DiceRoller::twister))};
    Molecule *someNei{plane.neigh_wrap(row, col, neiNum).get()};

    const double myFate{DiceRoller::probabilityGen(DiceRoller::twister)};

    // sanity check: if not bonded to a neighbour, mole should be free, vice
    // versa
    assert((mole->bon_nei == 0 && mole->m_typeComp == Molecule::free) ||
           (mole->bon_nei != 0 && mole->m_typeComp != Molecule::free));

    if (mole_type != Molecule::s) {
      if (myFate <= (Para::alpha * Para::decay_probability)) {
        decay(plane.cell(row, col).get());
        // get() returns a raw pointer to the molecule owned
        // by the unique_ptr at the row, col
      } else if (myFate <= (Para::alpha * (Para::decay_probability +
                                           Para::diffusion_probability))) {
        diffuse(mole_ptr, row, col); // cannot use get() here because need to
                                     // swap the underlying unique_ptr I think
      } else if (myFate <= (Para::alpha * (Para::decay_probability +
                                           Para::diffusion_probability +
                                           Para::complex_probability)) &&
                 mole->m_typeComp == Molecule::free &&
                 someNei->m_typeComp == Molecule::free) {
        try {
          formingComplex(determineComplex(mole, someNei, mole_type), mole,
                         neiNum, someNei);
        } catch (signed int) {
          std::cerr << "Exception detected in formingComplex()";
        }
      }
    }
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

newCA::~newCA() {
  if (display_p)
    delete display_p;
  display_p = nullptr;
}
