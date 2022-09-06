#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"
#include "para.h"
#include "random.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  /* std::uniform_int_distribution test{1, 2}; */
  for (unsigned row = 1; row <= nrow; row++) {
    for (unsigned col = 1; col <= ncol; col++) {
      /* if (col % 2 == 0) */
      /*   plane.cell(row, col) = */
      /*       std::make_unique<Molecule>(Molecule::p, Molecule::tempP); */
      /* else */
      /*   plane.cell(row, col) = */
      /*       std::make_unique<Molecule>(Molecule::p, Molecule::cata); */
      /* plane.cell(row, col) = */
      /*     std::make_unique<Molecule>(Molecule::p, Molecule::free); */
      plane.cell(row, col) = std::make_unique<Molecule>();
      if (Para::read == false) {
        switch (DiceRoller::typeInitializer(DiceRoller::twister)) {
        case 1:
          plane.cell(row, col)->m_typeRep = Molecule::p;
          break;
        case 2:
          plane.cell(row, col)->m_typeRep = Molecule::q;
          break;
        case 3:
          plane.cell(row, col)->m_typeRep = Molecule::s;
          break;
        }
      }
    }
  }

  if (Para::read == true) {
    readField("2475", 0, 51);
  }

  if (Para::visualization == 1) {
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
    display_p->color_rgb(red, 255, 0, 0);
    display_p->color_rgb(white, 255, 255, 255);
    if (Para::movie == 1)
      display_p->open_window();
    display_p->open_png();
  }
}

void newCA::visualize() {
  plane_to_display();
  if (Para::movie == 1)
    display_p->draw_window();
  display_p->draw_png();
  return;
}

using newcaFcn = int (newCA::*)(int, int, int);
void newCA::writeDensity(const long t, newcaFcn fcn) {
  if (!density) {
    std::cerr << "Could not access density.txt\n";
  }

  density << t * Para::alpha << ' ';

  // loop through combinations of ppp (0, 0, 0) ... qqq (1, 1, 1); for numeric
  // representation of replicator types, see molecule.h enum
  for (int type1{0}; type1 <= 1; type1++) {
    for (int type2{0}; type2 <= 1; type2++) {
      for (int i{0}; i <= 1; i++) {
        density << static_cast<double>((this->*fcn)(type1, type2, i)) /
                       Para::grid_size
                << ' ';
      }
    }
  }

  for (int simpType{0}; simpType <= 2; simpType++) {
    density << static_cast<double>(testSimple(simpType)) / Para::grid_size
            << ' ';
    if (simpType == 2)
      density << '\n';
  }

  density.flush();
  /* std::cout << "flushed \n"; */
}

void newCA::writeField(const long t) {
  if (t == 0)
    field << DiceRoller::timeSeed << '\n';

  for (int x{1}; x <= nrow; x++)
    for (int y{1}; y <= ncol; y++)
      if (plane.cell(x, y)->m_typeRep != Molecule::s) {
        field << t * Para::alpha << ' ' << x << ' ' << y << ' '
              << plane.cell(x, y)->m_typeRep << ' '
              << plane.cell(x, y)->m_typeComp << ' ';
        if (plane.cell(x, y)->m_typeComp != Molecule::free) {
          unsigned neiX{};
          unsigned neiY{};
          plane.xy_neigh_wrap(x, y, plane.cell(x, y)->bon_nei, neiX, neiY);
          field << plane.cell(x, y)->bon_nei << ' ' << neiX << ' ' << neiY
                << ' ';
        }
        for (int i{0}; i < std::size(plane.cell(x, y)->m_rateList); i++)
          field << plane.cell(x, y)->m_rateList[i] << ' ';
        field << '\n';
      }
  field.flush();
  /* std::cout << "saved \n"; */
}

/* Because this function uses std::find, to find lines matching the exact t, one
 * should pass a time argument that includes a decimal point. This avoids the
 * situation of a later time point containing the string of the ealier time
 * point you've entered.
 */
void newCA::readField(std::string_view t, unsigned lowerB, unsigned upperB) {
  // this map is used for determining the m_myCataParam of the catalyst in
  // complexes we instantiate in a saved field, as that information is not
  // saved; for this purpose we need to know the m_typeComp of both complex
  // molecules; the TypeComplex enum is implicitly converted to numbers and used
  // as keys; the associated value is the index of the corresponding catalytic
  // parameter. For details see molecule.h.
  const std::map<std::string, unsigned> whichK{
      {"002", 0}, {"003", 1}, {"013", 2}, {"012", 3},
      {"102", 4}, {"103", 5}, {"113", 6}, {"112", 7}};
  std::ifstream fin{"field.txt"};
  if (!fin) {
    std::cerr << "Cannot find field.txt\n";
    exit(1);
  }
  std::string moleState{};
  auto isSpace{[](auto ch) {
    return (ch == ' ');
  }}; // lambda to test if element pointed by iterator is a space char
  std::string::iterator spaceAfterRow, spaceAfterCol, spaceAfterType,
      spaceAfterComp, spaceAfterBon, spaceAfterNeiX, spaceAfterNeiY;

  Molecule::TypeComplex compStatus;
  Molecule::TypeReplicator myType;

  while (std::getline(fin, moleState)) {
    unsigned myRow{};
    unsigned myCol{};
    unsigned neiX{};
    unsigned neiY{};
    int counter{};
    if (moleState.find(t) != std::string::npos) {
      for (auto ch{moleState.begin()}; ch != moleState.end(); ++ch) {
        if (*ch == ' ') {
          counter++;
          switch (counter) {
          case 1: {
            spaceAfterRow = std::find_if(ch + 1, moleState.end(), isspace);
            myRow = std::stoi(std::string(
                ch + 1, spaceAfterRow)); // std::string construction using
                                         // iterators makes a string from
                                         // [firstIt, lastIt);
            assert(myRow >= 1 || myRow <= Para::sys_nrow);
            continue;
          }
          case 2: {
            spaceAfterCol =
                std::find_if(spaceAfterRow + 1, moleState.end(), isSpace);
            myCol = std::stoi(std::string(spaceAfterRow + 1, spaceAfterCol));
            assert(myCol >= 1 || myCol <= Para::sys_ncol);
            continue;
          }
          case 3: {
            spaceAfterType =
                std::find_if(spaceAfterCol + 1, moleState.end(), isSpace);
            myType = static_cast<Molecule::TypeReplicator>(
                std::stoi(std::string(spaceAfterCol + 1, spaceAfterType)));
            assert(myType == 0 || myType == 1 || myType == 2);
            if (myCol >= lowerB && myCol <= upperB)
              plane.cell(myRow, myCol)->m_typeRep = myType;
            continue;
          }
          case 4: {
            spaceAfterComp =
                std::find_if(spaceAfterType + 1, moleState.end(), isSpace);
            compStatus = static_cast<Molecule::TypeComplex>(
                std::stoi(std::string(spaceAfterType + 1, spaceAfterComp)));
            assert(compStatus == 0 || compStatus == 1 || compStatus == 2 ||
                   compStatus == 3);
            if (myCol >= lowerB && myCol <= upperB)
              plane.cell(myRow, myCol)->m_typeComp = compStatus;
            continue;
          }
          case 5: {
            // the string after compStatus can be either bon_nei or k0
            // depending on compStatus; defining the following three variables
            // outside the following if statement is because I'd like to
            // conditionally initialize "start" (see below) and also for setting
            // m_myCataParam of catalyst
            spaceAfterBon =
                std::find_if(spaceAfterComp + 1, moleState.end(), isSpace);
            spaceAfterNeiX =
                std::find_if(spaceAfterBon + 1, moleState.end(), isSpace);
            spaceAfterNeiY =
                std::find_if(spaceAfterNeiX + 1, moleState.end(), isSpace);

            if (plane.cell(myRow, myCol)->m_typeComp > Molecule::free &&
                myCol >= lowerB && myCol <= upperB) {
              assert(plane.cell(myRow, myCol)->nei_ptr == nullptr &&
                     plane.cell(myRow, myCol)->bon_nei == 0);
              plane.cell(myRow, myCol)->bon_nei =
                  std::stoi(std::string(spaceAfterComp + 1, spaceAfterBon));

              // checking consistency between bon_nei and m_typeComp
              assert(plane.cell(myRow, myCol)->bon_nei != 0);

              neiX = std::stoi(std::string(spaceAfterBon + 1, spaceAfterNeiX));
              neiY = std::stoi(std::string(spaceAfterNeiX + 1, spaceAfterNeiY));

              plane.cell(myRow, myCol)->nei_ptr = plane.cell(neiX, neiY).get();
            }
            // assign catalytic parameters;
            if (myCol >= lowerB && myCol <= upperB) {
              auto &start(plane.cell(myRow, myCol)->m_typeComp == 0
                              ? spaceAfterComp
                              : spaceAfterNeiY);

              std::string buf;
              std::stringstream ss{std::string(start + 1, moleState.end())};
              std::vector<std::string> Ks;
              while (getline(ss, buf, ' '))
                Ks.push_back(buf);

              for (int i{0}; i <= 7; i++) {
                plane.cell(myRow, myCol)->m_rateList[i] = std::stod(Ks[i]);
              }
            }
            continue;
          }
          }
        }
      }
    }
  }

  /* for (int row{1}; row <= 512; row++) { */
  /*   if (plane.cell(row, lowerB)->m_typeComp != Molecule::free) { */
  /*     decay(plane.cell(row, lowerB)->nei_ptr); */
  /*     decay(plane.cell(row, lowerB).get()); */
  /*   } */
  /*   if (plane.cell(row, upperB)->m_typeComp != Molecule::free) { */
  /*     decay(plane.cell(row, upperB)->nei_ptr); */
  /*     decay(plane.cell(row, upperB).get()); */
  /*   } */
  /* } */

  for (int myRow{1}; myRow <= Para::sys_nrow; myRow++)
    for (int myCol{1}; myCol < upperB; myCol++) {
      if (plane.cell(myRow, myCol)->m_typeComp == 1) {
        std::string key{
            std::to_string(plane.cell(myRow, myCol)->m_typeComp) +
            std::to_string(plane.cell(myRow, myCol)->nei_ptr->m_typeRep) +
            std::to_string(plane.cell(myRow, myCol)->nei_ptr->m_typeComp)};
        plane.cell(myRow, myCol)->m_myCataParam = whichK.find(key)->second;
      }
    }

  for (int x{1}; x <= Para::sys_nrow; x++)
    for (int y{1}; y <= Para::sys_ncol; y++) {
      if (plane.cell(x, y)->m_typeComp != Molecule::free) {
        assert(plane.cell(x, y)->m_typeRep != Molecule::s &&
               plane.cell(x, y)->nei_ptr->m_typeRep != Molecule::s);
        assert(plane.cell(x, y)->nei_ptr->bon_nei ==
               plane.neighbor_converter(plane.cell(x, y)->bon_nei));
        assert(plane.cell(x, y)->nei_ptr->nei_ptr == plane.cell(x, y).get());
        assert((plane.cell(x, y)->m_myCataParam != 10) !=
               (plane.cell(x, y)->nei_ptr->m_myCataParam != 10));
      } else {
        assert(plane.cell(x, y)->bon_nei == 0 &&
               plane.cell(x, y)->nei_ptr == nullptr &&
               plane.cell(x, y)->m_myCataParam == 10);
      }
    }
}

/* void newCA::testDummy(int x) { */
/*   for (int i{0}; i <= 9; i++) { */
/*     int counts{}; */
/*     for (unsigned row{1}; row <= nrow; row++) { */
/*       for (unsigned col{1}; col <= ncol; col++) { */
/*         if (plane.cell(row, col)->dummy >= (i / 10.0) && */
/*             plane.cell(row, col)->dummy < ((i / 10.0) + 0.1)) */
/*           counts++; */
/*       } */
/*     } */
/*     std::cout << "0." << i << "--" << (i / 10.0) + 0.1 << '|'; */
/*     if (x == 1) */
/*       counts = log(counts); */
/*     for (int g{0}; g < counts; g++) */
/*       if (x == 0 && g % 10 == 0) */
/*         std::cout << '='; */
/*       else if (x == 1) */
/*         std::cout << '='; */
/*     std::cout << '\n'; */
/*   } */
/* } */

int newCA::testSimple(char type) {
  std::uint32_t simpleDensity{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      if (plane.cell(row, col)->m_typeRep == type &&
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
      if (plane.cell(row, col)->m_typeRep == type)
        ++testDensity;
    }
  }
  return testDensity;
}

int newCA::testComplex(int type1, int type2, int cType) {
  assert(cType == 0 || cType == 1);
  unsigned numb{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      if (plane.cell(row, col)->nei_ptr &&
          plane.cell(row, col)->getTypeReplicator() == type1 &&
          plane.cell(row, col)->nei_ptr->getTypeReplicator() == type2 &&
          plane.cell(row, col)->m_typeComp == Molecule::cata) {
        if (cType == 0) {
          if (plane.cell(row, col)->nei_ptr->m_typeComp == Molecule::tempP)
            ++numb;
        } else if (cType == 1) {
          if (plane.cell(row, col)->nei_ptr->m_typeComp == Molecule::tempQ)
            ++numb;
        }
      }
    }
  }
  return numb;
}

void newCA::writeAverageK(const int k, const long t) {
  if (k == 0)
    kParam << t * Para::alpha << ' ';

  long repCount{};
  double averageK{};
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      if (plane.cell(row, col)->m_typeRep != Molecule::s) {
        repCount += 1;
        averageK += plane.cell(row, col)->m_rateList[k];
      }
    }
  }

  if (repCount == 0)
    kParam << "0 ";
  else
    kParam << (averageK / repCount) << ' ';

  if (k == 7)
    kParam << '\n';

  kParam.flush();
}

void newCA::writeDistribution(const long t) {
  std::vector names{"kppp", "kppq", "kpqq", "kpqp",
                    "kqpp", "kqpq", "kqqq", "kqqp"};
  for (int k{0}; k < 8; k++) {
    std::map<int, int> histogram;
    for (unsigned row{1}; row <= nrow; row++)
      for (unsigned col{1}; col <= ncol; col++)
        if (plane.cell(row, col)->m_typeRep != Molecule::s) {
          ++histogram[static_cast<int>(plane.cell(row, col)->m_rateList[k] *
                                       1e04)];
        }
    kDistr << "Time " << t * Para::alpha << " Parameter " << names[k] << '\n';
    for (const auto &d : histogram)
      kDistr << d.first * 1e-4 << ' ' << d.second << '\n';
    kDistr.flush();
  }
}

void newCA::plane_to_display() { // Does not paint display! Just
                                 // transports plane data to it.
                                 // Also perhaps not a good idea to color
                                 // according to typeRep
  for (unsigned row{1}; row <= nrow; row++) {
    for (unsigned col{1}; col <= ncol; col++) {
      switch (plane.cell(row, col)->m_typeRep) {
      case Molecule::p:
        /* if (plane.cell(row, col)->m_typeComp == Molecule::cata) */
        display_p->put_pixel(CA, row, col, blue);
        /* else if (plane.cell(row, col)->m_typeComp == Molecule::tempP) */
        /*   display_p->put_pixel(CA, row, col, red); */
        /* else */
        /*   display_p->put_pixel(CA, row, col, white); */
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

void newCA::decay(Molecule *mole, unsigned dis) {
  if (mole->nei_ptr != nullptr) { // dissociation by decay
    mole->bon_nei = 0;            // remove mole's index info on bonded nei
    mole->nei_ptr->bon_nei = 0;   // bonded nei loses its bonded nei (i.e. mole)
    mole->nei_ptr->nei_ptr = nullptr; // bonded nei loses nei_ptr (to mole)
    mole->nei_ptr->m_typeComp = Molecule::free; // bonded nei is now free
    mole->nei_ptr->m_myCataParam = 10;
    mole->nei_ptr = nullptr; // mole loses pointer to bonded nei
  }

  assert(!dis || dis == 1);
  if (!dis) { // !dis for pure decay, dis for dissociation w/o decay
    mole->m_typeRep = Molecule::s; // mole decays but old bonded nei doesn't
  }
  mole->m_myCataParam = 10;
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
if and what complex forms between the chosen molecule (passed by address
to this function) and a molecule in its neighbourhood (randomly chosen
through this function). There are 3 kinds of numbers this function can
return: 0-7: complex has formed and mole is the catalyst, someNei is the
template 100-107: complex has formed and someNei is the catalyst, mole is
the template
>=1000: complex did not form or someNei is S */
int newCA::determineComplex(const double myFate, double &cumuProb,
                            Molecule *mole, Molecule *someNei,
                            const int mole_type) {
  Molecule::TypeReplicator nei_type{someNei->m_typeRep};
  assert((someNei->m_typeComp == Molecule::free) &&
         (mole->m_typeRep == mole_type) && (mole_type != Molecule::s) &&
         (someNei->m_typeComp == Molecule::free) && (nei_type != Molecule::s));

  switch (mole_type) {
  case 0:
    switch (nei_type) {
    case 0: // if mole & someNei are both P
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += Para::alpha * Para::beta * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{0}; i <= 1; i++) {
        cumuProb += Para::alpha * Para::beta * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1000; // complex did not form between two P molecules
    case 1:        // if mole is P and someNei is Q
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += Para::alpha * Para::beta * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += Para::alpha * Para::beta * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1001; // complex did not form between P and Q
    default:
      return 1005; // someNei is a S
    }
  case 1:
    switch (someNei->getTypeReplicator()) { // if mole is Q and
                                            /* someNei is */
                                            /* P *1/ */
    case 0:
      for (unsigned i{4}; i <= 5; i++) {
        cumuProb += Para::alpha * Para::beta * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{2}; i <= 3; i++) {
        cumuProb += Para::alpha * Para::beta * someNei->m_rateList[i];
        if (myFate <= cumuProb)
          return (i + 100);
      }
      /* std::cout << myFate << ' ' << cumuProb << std::endl; */
      return 1002; // complex did not form between Q and P
    case 1:        // if both are Q
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += Para::alpha * Para::beta * mole->m_rateList[i];
        if (myFate <= cumuProb)
          return i;
      }
      for (unsigned i{6}; i <= 7; i++) {
        cumuProb += Para::alpha * Para::beta * someNei->m_rateList[i];
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
    mole->m_myCataParam = complex;
    someNei->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 1:
  case 2:
  case 5:
  case 6:
    mole->m_typeComp = Molecule::cata;
    mole->m_myCataParam = complex;
    someNei->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 100:
  case 103:
  case 104:
  case 107:
    someNei->m_typeComp = Molecule::cata;
    someNei->m_myCataParam = (complex - 100);
    mole->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 101:
  case 102:
  case 105:
  case 106:
    someNei->m_typeComp = Molecule::cata;
    someNei->m_myCataParam = (complex - 100);
    mole->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 1000 ... 1003:
    mole->bon_nei = 0;
    mole->nei_ptr = nullptr;
    someNei->bon_nei = 0;
    someNei->nei_ptr = nullptr;
    return;
  }
}

// Overload for well-mixed; this function does not touch bon_nei
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
    mole->m_myCataParam = complex;
    someNeiWM->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed mole-nei" << std::endl; */
    return;
  case 1:
  case 2:
  case 5:
  case 6:
    mole->m_typeComp = Molecule::cata;
    mole->m_myCataParam = complex;
    someNeiWM->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 100:
  case 103:
  case 104:
  case 107:
    someNeiWM->m_typeComp = Molecule::cata;
    someNeiWM->m_myCataParam = (complex - 100);
    mole->m_typeComp = Molecule::tempP;
    /* std::cout << "complex formed nei-mole" << std::endl; */
    return;
  case 101:
  case 102:
  case 105:
  case 106:
    someNeiWM->m_typeComp = Molecule::cata;
    someNeiWM->m_myCataParam = (complex - 100);
    mole->m_typeComp = Molecule::tempQ;
    /* std::cout << "complex formed" << std::endl; */
    return;
  case 1000 ... 1003:
    mole->nei_ptr = nullptr;
    someNeiWM->nei_ptr = nullptr;
    return;
  }
}

void newCA::replication(Molecule *mole, Molecule *someNei) {
  assert(mole->nei_ptr && mole->nei_ptr->nei_ptr == mole);
  assert((mole->m_typeComp == Molecule::cata) !=
         (mole->nei_ptr->m_typeComp == Molecule::cata));
  Molecule *&templ{(mole->m_typeComp != Molecule::cata) ? mole : mole->nei_ptr};
  std::copy(std::begin(templ->m_rateList), std::end(templ->m_rateList),
            std::begin(someNei->m_rateList));
  for (int i{0}; i < std::size(templ->m_rateList); i++)
    assert(someNei->m_rateList[i] == templ->m_rateList[i]);
  assert((templ->m_typeComp == Molecule::tempP) !=
         (templ->m_typeComp == Molecule::tempQ));
  if (templ->m_typeComp == Molecule::tempP)
    someNei->m_typeRep = Molecule::p;
  else if (templ->m_typeComp == Molecule::tempQ)
    someNei->m_typeRep = Molecule::q;

  mole->bon_nei = 0;
  mole->nei_ptr->bon_nei = 0;
  mole->m_typeComp = Molecule::free;
  mole->nei_ptr->m_typeComp = Molecule::free;
  mole->m_myCataParam = 10;
  mole->nei_ptr->m_myCataParam = 10;
  mole->nei_ptr->nei_ptr = nullptr;
  mole->nei_ptr = nullptr;
  assert(!(mole->nei_ptr) && someNei->m_typeComp == Molecule::free);

  /* linearMutation(someNei); */
  exponentialMutation(someNei);
}

void newCA::linearMutation(Molecule *someNei) {
  if (DiceRoller::probabilityGen(DiceRoller::twister) <=
      Para::mutation_probability) {
    /* std::cout << "mutating\n"; */
    for (int i{0}; i < std::size(someNei->m_rateList); i++) {
      someNei->m_rateList[i] += DiceRoller::mutaGen(DiceRoller::twister);
      if (someNei->m_rateList[i] > 1)
        someNei->m_rateList[i] = 2 - someNei->m_rateList[i];
      else if (someNei->m_rateList[i] < 0)
        someNei->m_rateList[i] = 0;
      assert(someNei->m_rateList[i] <= 1 && someNei->m_rateList[i] >= 0);
    }
  }
}

void newCA::exponentialMutation(Molecule *someNei) {
  if (DiceRoller::probabilityGen(DiceRoller::twister) <=
      Para::mutation_probability) {
    for (int i{0}; i < std::size(someNei->m_rateList); i++) {
      someNei->m_rateList[i] *= exp(DiceRoller::mutaGen(DiceRoller::twister));
      if (someNei->m_rateList[i] > 1)
        someNei->m_rateList[i] = 2 - someNei->m_rateList[i];
      assert(someNei->m_rateList[i] <= 1 && someNei->m_rateList[i] >= 0);
    }
  }
}

void newCA::update_squares() {
  unsigned row{};
  unsigned col{};

  // DIFFUSION
  for (int i{1}; i <= Para::grid_size; i++) {
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);
    if (DiceRoller::probabilityGen(DiceRoller::twister) <=
        Para::alpha * Para::diffusion_probability)
      diffuse(row, col);
  }

  // COMPLEX FORMATION & REPLICATION
  for (int i{1}; i <= Para::grid_size; i++) { // # of times is counted from 1
    row = DiceRoller::randomRowOrCol(DiceRoller::twister);
    col = DiceRoller::randomRowOrCol(DiceRoller::twister);

    // mole_ptr is the unique_ptr that owns a molecule; mole is the
    // pointer to the actual molecule, used for passing the molecule to
    // functions (obtained by unique_ptr.get())
    Molecule *mole{plane.cell(row, col).get()};
    auto mole_type{mole->m_typeRep};

    unsigned neiNum{mole->bon_nei
                        ? plane.neigh_7_select(
                              DiceRoller::randomNeiExcl(DiceRoller::twister),
                              mole->bon_nei) // exclude nei if mole is already
                                             // in a complex
                        : static_cast<unsigned int>(
                              DiceRoller::randomNei(DiceRoller::twister))};
    Molecule *someNei{plane.neigh_wrap(row, col, neiNum).get()};
    /* Molecule *someNeiWM{ */
    /*     plane */
    /*         .cell(DiceRoller::randomRowOrCol(DiceRoller::twister), */
    /*               DiceRoller::randomRowOrCol(DiceRoller::twister)) */
    /*         .get()}; */
    /* while (someNeiWM == mole || someNeiWM == mole->nei_ptr) { */
    /*   someNeiWM = plane // Reassignment only needed in well-mixed system */
    /*                   .cell(DiceRoller::randomRowOrCol(DiceRoller::twister),
     */
    /*                         DiceRoller::randomRowOrCol(DiceRoller::twister))
     */
    /*                   .get(); */
    /* } */

    const double myFate{DiceRoller::probabilityGen(DiceRoller::twister)};
    double cumuProb{};

    // sanity check: if not bonded to a neighbour, mole should be free,
    // vice versa
    assert((mole->bon_nei == 0 && mole->m_typeComp == Molecule::free) ||
           (mole->bon_nei != 0 && mole->m_typeComp != Molecule::free));
    assert((mole->nei_ptr == nullptr && mole->m_typeComp == Molecule::free) ||
           (mole->nei_ptr != nullptr && mole->m_typeComp != Molecule::free));

    cumuProb += Para::alpha * Para::decay_probability;
    if (myFate <= cumuProb) {
      decay(mole);
    } else {
      if ((mole->m_typeComp == Molecule::free) &&
          (mole->m_typeRep != Molecule::s) &&
          (someNei->m_typeRep != Molecule::s) &&
          (someNei->m_typeComp == Molecule::free)) {
        assert(mole->m_myCataParam == 10 && someNei->m_myCataParam == 10);
        try {
          int myComplex{
              determineComplex(myFate, cumuProb, mole, someNei, mole_type)};
          if (myComplex > 1003) // note that anything beyond 1003 means
                                // something in determineComplex is broken
                                // e.g. someNei is a s!
            throw -1;
          formingComplex(myComplex, mole, neiNum, someNei);
          /* formingComplex(myComplex, mole, someNeiWM); */
          continue;
        } catch (int) {
          std::cerr << "Exception detected in determineComplex(); ";
          exit(2);
        }
      }
      /* cumuProb = 0; */
      assert(cumuProb == Para::alpha * Para::decay_probability);
      if ((mole->m_typeComp != Molecule::free)) {
        assert(mole->m_myCataParam <= 7 != mole->nei_ptr->m_myCataParam <= 7);
        double dissProb{
            mole->m_typeComp == Molecule::cata
                ? (1 - mole->m_rateList[mole->m_myCataParam])
                : (1 -
                   mole->nei_ptr->m_rateList[mole->nei_ptr->m_myCataParam])};
        /* assert(dissProb = 0.2); */
        cumuProb += Para::alpha * Para::beta * dissProb;
        if (myFate <= cumuProb)
        /* replication(mole, someNeiWM); */
        {
          /* std::cout << testSimple(0) << '\n'; */
          mole->m_typeComp = Molecule::free;
          mole->m_myCataParam = 10;
          mole->bon_nei = 0;
          mole->nei_ptr->m_typeComp = Molecule::free;
          mole->nei_ptr->m_myCataParam = 10;
          mole->nei_ptr->bon_nei = 0;
          mole->nei_ptr->nei_ptr = nullptr;
          mole->nei_ptr = nullptr;
          /* std::cout << testSimple(0) << '\n'; */
        } else if (myFate <= (cumuProb + Para::alpha * Para::gamma) &&
                   someNei->m_typeComp == Molecule::free &&
                   someNei->m_typeRep == Molecule::s)
          replication(mole, someNei);
      }
    }
  }
}

void newCA::reallycba() {
  int row{DiceRoller::randomRowOrCol(DiceRoller::twister)};
  int col{DiceRoller::randomRowOrCol(DiceRoller::twister)};
  double cumuProb{};
  double myFate{DiceRoller::probabilityGen(DiceRoller::twister)};

  Molecule *someNeiWM{plane
                          .cell(DiceRoller::randomRowOrCol(DiceRoller::twister),
                                DiceRoller::randomRowOrCol(DiceRoller::twister))
                          .get()};
  while (someNeiWM == plane.cell(row, col).get()) {
    someNeiWM = plane // Reassignment only needed in well-mixed system
                    .cell(DiceRoller::randomRowOrCol(DiceRoller::twister),
                          DiceRoller::randomRowOrCol(DiceRoller::twister))
                    .get();
  }

  if (plane.cell(row, col)->m_typeComp != Molecule::free) {
    if (myFate <= Para::alpha * Para::beta * 0.01) {
      /* assert(plane.cell(row, col)->m_typeComp != Molecule::free && */
      /*        plane.cell(row, col)->nei_ptr->m_typeComp !=
       * Molecule::free); */
      plane.cell(row, col)->m_typeComp = Molecule::free;
      plane.cell(row, col)->nei_ptr->m_typeComp = Molecule::free;
      plane.cell(row, col)->nei_ptr->nei_ptr = nullptr;
      plane.cell(row, col)->nei_ptr = nullptr;
    }
  } else if (plane.cell(row, col)->m_typeComp == Molecule::free &&
             someNeiWM->m_typeComp == Molecule::free) {
    if (myFate <= (cumuProb + Para::alpha * Para::beta * (1 - 0.01))) {
      plane.cell(row, col)->m_typeComp = Molecule::cata;
      someNeiWM->m_typeComp = Molecule::tempP;
      plane.cell(row, col)->nei_ptr = someNeiWM;
      someNeiWM->nei_ptr = plane.cell(row, col).get();
    } else if (myFate <=
               (cumuProb + Para::alpha * Para::beta * 2 * (1 - 0.01))) {
      plane.cell(row, col)->m_typeComp = Molecule::tempP;
      someNeiWM->m_typeComp = Molecule::cata;
      plane.cell(row, col)->nei_ptr = someNeiWM;
      someNeiWM->nei_ptr = plane.cell(row, col).get();
    }
  }
}

newCA::~newCA() {
  if (display_p)
    delete display_p;
  display_p = nullptr;
}
