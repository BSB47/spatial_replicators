#include "newca.h"
#include "para.h"

#include <iostream>

#ifndef MOLECULE_H
#define MOLECULE_H

class Molecule {
public:
  enum TypeReplicator {
    p,
    q,
    s,
  };

  enum TypeComplex {
    free,
    cata,
    tempP, // this means a template producing P
    tempQ,
  };

private:
  double m_rateList[8]{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  /* Below are the rates and their corresponding indices */
  /* 0: m_kppp    last two letters stand for template/product */
  /* 1: m_kppq */
  /* 2: m_kpqq */
  /* 3: m_kpqp */
  /* 4: m_kqpp */
  /* 5: m_kqpq */
  /* 6: m_kqqq */
  /* 7: m_kqqp */

  TypeReplicator m_typeRep{};
  TypeComplex m_typeComp{free};

  unsigned bon_nei{};

  double m_mutation_probability{};

public:
  Molecule() {}

  Molecule(TypeReplicator typeR, TypeComplex typeC)
      : m_typeRep{typeR}, m_typeComp{typeC} {
    std::cout << "I'm born\n";
  }

  ~Molecule() { std::cout << "I'm dead"; }

  /* Molecule(Molecule &&a) noexcept : m_ratePtr(a.m_ratePtr) { */
  /*   a.m_ratePtr = nullptr; */
  /* } */

  /* Molecule &operator=(Molecule &&a) noexcept { */
  /*   if (&a == this) // self-assignment check */
  /*     return *this; */

  /* m_ratePtr = a.m_ratePtr; */
  /* a.m_ratePtr = nullptr; */

  /* return *this; */
  /* } */

  const TypeReplicator &getTypeReplicator() const { return m_typeRep; }
  void setTypeRep(TypeReplicator myType) { m_typeRep = myType; }

  friend int newCA::determineComplex(Molecule *mole, Molecule *someNei,
                                     int mole_type);
  friend void newCA::formingComplex(int complex, Molecule *mole, int neiNum,
                                    Molecule *someNei);
  friend void newCA::update_squares();
};

#endif
