#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"
#include "para.h"
#include "random.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <string_view>
#include <vector>

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  for (unsigned row = 1; row <= nrow; row++) {
    for (unsigned col = 1; col <= ncol; col++) {
      plane.cell(row, col) = std::make_unique<Molecule>();
      switch (DiceRoller::typeInitializer(DiceRoller::twister)) {
      case 1:
        plane.cell(row, col)->setTypeRep(Molecule::p);
        break;
      case 2:
        plane.cell(row, col)->setTypeRep(Molecule::q);
        break;
      case 3:
        plane.cell(row, col)->setTypeRep(Molecule::s);
        break;
      }
    }
  }

  std::vector<CashPanelInfo> panel_info(1);
  panel_info[CA].n_row = nrow;
  panel_info[CA].n_col = ncol;
  panel_info[CA].o_row = 0;
  panel_info[CA].o_col = 0;

  // display_p is a pointer to the CASH window
  display_p = new CashDisplay(
      Para::sys_nrow, Para::sys_ncol, panel_info,
      Para::scale); // window_row/col are sys_nrow/col for testing.
  display_p->color_rgb(blue, 0, 0, 255);
  display_p->color_rgb(white, 255, 255, 255);
  /* display_p->open_window(); */
  display_p->open_png();
  /* std::cout << testDensity(0) << '\n'; */
}

void newCA::visualize() {
  plane_to_display();
  /* display_p->draw_window(); */
  display_p->draw_png();
  return;
}

using newcaFcn = int (newCA::*)(int);
void newCA::writeFile(const long t, int c, newcaFcn fcn) {
  if (!output) {
    std::cerr << "Could not access output.txt\n";
  }

  output << t * Para::alpha << ' '
         << static_cast<double>((this->*fcn)(c)) / Para::grid_size << ' '
         << static_cast<double>(testSimple(0)) / Para::grid_size << ' '
         << static_cast<double>(testSimple(1)) / Para::grid_size << ' '
         << static_cast<double>(testDensity(1)) / Para::grid_size << '\n';
  if (t % 100 == 0)
    output.flush();
}

int newCA::testSimple(char type) {
  std::uint32_t simpleDensity{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      if (plane.cell(row, col)->getTypeReplicator() == type &&
          plane.cell(row, col)->m_typeComp == Molecule::free)
        ++simpleDensity;
    }
  }
  return simpleDensity;
}

int newCA::testDensity(char type) {
  std::uint32_t testDensity{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      if (plane.cell(row, col)->getTypeReplicator() == type)
        ++testDensity;
    }
  }
  return testDensity;
}

int newCA::testComplex(std::string_view compType) { // compType sanity checking
                                                    // is handled at entry point
  unsigned numb{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {

      /* if (compType == "PP") { */
      /*   if (plane.cell(row, col)->nei_ptr && */
      /*       plane.cell(row, col)->m_typeComp == Molecule::cata && */
      /*       plane.cell(row, col)->m_typeRep == Molecule::p && */
      /*       plane.cell(row, col)->nei_ptr->m_typeRep == Molecule::p) { */
      /*     assert(plane.cell(row, col)->nei_ptr->m_typeComp == Molecule::tempP
       * || */
      /*            plane.cell(row, col)->nei_ptr->m_typeComp ==
       * Molecule::tempQ); */
      /*     ++numb; */
      /*   } */
      /* } else if (compType == "PQ") { */
      /*   if (plane.cell(row, col)->nei_ptr && */
      /*       plane.cell(row, col)->m_typeComp == Molecule::cata && */
      /*       plane.cell(row, col)->m_typeRep == Molecule::p && */
      /*       plane.cell(row, col)->nei_ptr->m_typeRep == Molecule::p) { */
      /*     assert(plane.cell(row, col)->nei_ptr->m_typeComp == Molecule::tempP
       * || */
      /*            plane.cell(row, col)->nei_ptr->m_typeComp ==
       * Molecule::tempQ); */
      /*     ++numb; */
      /*   } */
      /* } */
    }
    std::cout << numb << std::endl;
    return numb;
  }
}

void newCA::plane_to_display() { // Does not paint display! Just
                                 // transports plane data to it.
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

void newCA::decay(Molecule *mole) {
  if (mole->nei_ptr != nullptr) { // dissociation by decay
    mole->bon_nei = 0;            // remove mole's index info on bonded nei
    mole->nei_ptr->bon_nei = 0;   // bonded nei loses its bonded nei (i.e. mole)
    mole->nei_ptr->nei_ptr = nullptr; // bonded nei loses nei_ptr (to mole)
    mole->nei_ptr->m_typeComp = Molecule::free; // bonded nei is now free
    mole->nei_ptr = nullptr; // mole loses pointer to bonded nei
  }
  mole->setTypeRep(Molecule::s); // mole decays but old bonded nei doesn't
  mole->m_typeComp = Molecule::free;
}

void newCA::diffuse(unsigned row, unsigned col) {
  unsigned neiNum{
      static_cast<unsigned int>(DiceRoller::randomNei(DiceRoller::twister))};
  auto &container_a{plane.cell(row, col)};
  auto &container_b{plane.neigh_wrap(row, col, neiNum)};

  unsigned nei_row{};
  unsigned nei_col{};
  plane.xy_neigh_wrap(row, col, neiNum, nei_row, nei_col);

  if (container_a->m_typeComp != Molecule::free &&
      container_b->m_typeComp != Molecule::free &&
      container_a->nei_ptr != container_b.get()) {
    plane.neigh_wrap(row, col, container_a->bon_nei)
        .swap(plane.neigh_wrap(nei_row, nei_col, container_b->bon_nei));
    unsigned tmp{container_a->nei_ptr->bon_nei};
    container_a->nei_ptr->bon_nei = container_b->nei_ptr->bon_nei;
    container_b->nei_ptr->bon_nei = tmp;
    tmp = container_a->bon_nei;
    container_a->bon_nei = container_b->bon_nei;
    container_b->bon_nei = tmp;
  } else if (container_a->m_typeComp != Molecule::free &&
             container_b->m_typeComp == Molecule::free) {
    container_b.swap(plane.neigh_wrap(row, col, container_a->bon_nei));
    container_a->nei_ptr->bon_nei = neiNum;
    container_a->bon_nei = plane.neighbor_converter(neiNum);
  } else if (container_a->m_typeComp == Molecule::free &&
             container_b->m_typeComp != Molecule::free) {
    plane.xy_neigh_wrap(row, col, neiNum, nei_row, nei_col);
    container_a.swap(plane.neigh_wrap(nei_row, nei_col, container_b->bon_nei));
    container_b->nei_ptr->bon_nei = plane.neighbor_converter(neiNum);
    container_b->bon_nei = neiNum;
  } else if (container_a->nei_ptr == container_b.get()) {
    assert(plane.neighbor_converter(neiNum) == container_b->bon_nei &&
           container_a->bon_nei == neiNum);
    container_a->bon_nei = plane.neighbor_converter(neiNum);
    container_b->bon_nei = neiNum;
  }
  container_a.swap(container_b);
}

/* The following function returns a number which can be used to decipher
if and what complex forms between the chosen molecule (passed by address to
this function) and a molecule in its neighbourhood (randomly chosen through
this function). There are 3 kinds of numbers this function can return: 0-7:
complex has formed and mole is the catalyst, someNei is the template
100-107: complex has formed and someNei is the catalyst, mole is the
template
>=1000: complex did not form or someNei is S */
int newCA::determineComplex(const double myFate, double &cumuProb,
                            Molecule *mole, Molecule *someNei,
                            const int mole_type) {

  Molecule::TypeReplicator nei_type{someNei->getTypeReplicator()};
  assert((someNei->m_typeComp == Molecule::free) &&
         (mole->getTypeReplicator() == mole_type) &&
         (mole_type != Molecule::s) &&
         (someNei->m_typeComp == Molecule::free) && (nei_type != Molecule::s));

  switch (mole_type) {
  case 0:
    switch (nei_type) {
    case 0: // if mole & someNei are both P
      /* cumuProb += Para::alpha * mole->m_rateList[0]; */
      /* if (myFate <= cumuProb) */
      /*   return 0; */
      /* cumuProb += Para::alpha * someNei->m_rateList[0]; */
      /* if (myFate <= cumuProb) */
      /*   return 100; */
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += Para::alpha * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += Para::alpha * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1000; // complex did not form between two P molecules
    case 1:        // if mole is P and someNei is Q
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += Para::alpha * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += Para::alpha * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1001; // complex did not form between P and Q
    default:
      return 1005; // someNei is a S
    }
  case 1:
    switch (someNei->getTypeReplicator()) { // if mole is Q and someNei is P
    case 0:
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += Para::alpha * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += Para::alpha * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1002; // complex did not form between Q and P
    case 1:        // if both are Q
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += Para::alpha * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += Para::alpha * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1003; // complex did not form between 2 Qs
    default:
      return 1005; // someNei is a S
    }
  default:
    return 1010; // something is broken
  }
}

void newCA::formingComplex(int complex, Molecule *mole, int neiNum,
                           Molecule *someNei) {
  mole->bon_nei = neiNum;
  mole->nei_ptr = someNei;
  someNei->bon_nei = plane.neighbor_converter(neiNum);
  someNei->nei_ptr = mole;
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
  case 1000 ... 1003:
    mole->bon_nei = 0;
    mole->nei_ptr = nullptr;
    someNei->bon_nei = 0;
    someNei->nei_ptr = nullptr;
    return;
  }
}

void newCA::formingComplex(int complex, Molecule *mole, Molecule *someNeiWM) {
  assert(mole->nei_ptr == nullptr && someNeiWM->nei_ptr == nullptr);
  mole->nei_ptr = someNeiWM;
  someNeiWM->nei_ptr = mole;
  switch (complex) {
  case 0:
  case 3:
  case 4:
  case 7:
    mole->m_typeComp = Molecule::cata;
    someNeiWM->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed mole-nei" << std::endl; */
    return;
  case 1:
  case 2:
  case 5:
  case 6:
    mole->m_typeComp = Molecule::cata;
    someNeiWM->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 100:
  case 103:
  case 104:
  case 107:
    someNeiWM->m_typeComp = Molecule::cata;
    mole->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed nei-mole" << std::endl; */
    return;
  case 101:
  case 102:
  case 105:
  case 106:
    someNeiWM->m_typeComp = Molecule::cata;
    mole->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 1000 ... 1003:
    mole->nei_ptr = nullptr;
    someNeiWM->nei_ptr = nullptr;
    return;
  }
}

void newCA::update_squares() {
  unsigned row{};
  unsigned col{};

  /* // DIFFUSION */
  /* for (int i{1}; i <= Para::grid_size; i++) { */
  /*   row = DiceRoller::randomRowOrCol(DiceRoller::twister); */
  /*   col = DiceRoller::randomRowOrCol(DiceRoller::twister); */
  /*   diffuse(row, col); */
  /* } */

  // COMPLEX FORMATION & REPLICATION
  for (int i{1}; i <= Para::grid_size; i++) { // # of times is counted from 1
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);

    // mole_ptr is the unique_ptr that owns a molecule; mole is the
    // pointer to the actual molecule, used for passing the molecule to
    // functions (obtained by unique_ptr.get())
    Molecule *mole{plane.cell(row, col).get()};
    auto mole_type{mole->getTypeReplicator()};

    /* unsigned neiNum{ */
    /*     mole->bon_nei */
    /*         ? plane.neigh_7_select( */
    /*               DiceRoller::randomNeiExcl(DiceRoller::twister), */
    /*               mole->bon_nei) // exclude nei if mole is already in a
     * complex */
    /*         : static_cast<unsigned int>( */
    /*               DiceRoller::randomNei(DiceRoller::twister))}; */
    /* Molecule *someNei{plane.neigh_wrap(row, col, neiNum).get()}; */
    Molecule *someNeiWM{
        plane
            .cell(DiceRoller::randomRowOrCol(DiceRoller::twister),
                  DiceRoller::randomRowOrCol(DiceRoller::twister))
            .get()};
    while (someNeiWM == mole) {
      someNeiWM = plane // Reassignment only needed in well-mixed system
                      .cell(DiceRoller::randomRowOrCol(DiceRoller::twister),
                            DiceRoller::randomRowOrCol(DiceRoller::twister))
                      .get();
    }

    const double myFate{DiceRoller::probabilityGen(DiceRoller::twister)};
    double cumuProb{};

    // sanity check: if not bonded to a neighbour, mole should be free,
    // vice versa
    /* assert((mole->bon_nei == 0 && mole->m_typeComp == Molecule::free)
     * || */
    /* (mole->bon_nei != 0 && mole->m_typeComp != Molecule::free)); */
    assert((mole->nei_ptr == nullptr && mole->m_typeComp == Molecule::free) ||
           (mole->nei_ptr != nullptr && mole->m_typeComp != Molecule::free));

    if (mole_type != Molecule::s) {
      cumuProb += Para::alpha * Para::decay_probability;
      if (myFate <= cumuProb) {
        decay(plane.cell(row, col).get());
      } else if ((mole->m_typeComp == Molecule::free) &&
                 (someNeiWM->m_typeRep != Molecule::s) &&
                 (someNeiWM->m_typeComp == Molecule::free)) {
        try {
          /* formingComplex(determineComplex(mole, someNeiWM, mole_type),
           * mole,
           */
          /*                neiNum, someNeiWM); */
          int myComplex{
              determineComplex(myFate, cumuProb, mole, someNeiWM, mole_type)};
          if (myComplex > 1003) // note that anything beyond 1003 means
                                // something in determineComplex is broken
                                // e.g. someNei is a s!
            throw -1;
          formingComplex(myComplex, mole, someNeiWM);
        } catch (int) {
          std::cerr << "Exception detected in determineComplex(); ";
        }
      }
      /* std::cout << cumuProb << std::endl; */
    }
  }
}

newCA::~newCA() {
  if (display_p)
    delete display_p;
  display_p = nullptr;
}
