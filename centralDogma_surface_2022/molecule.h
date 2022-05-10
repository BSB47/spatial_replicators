#include "newca.h"
#include "para.h"

#include <iostream>
#include <iterator>

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
  double m_rateList[8]{0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  /* Below are the rates and their corresponding indices */
  /* 0: m_kppp    last two letters stand for template/product */
  /* 1: m_kppq */
  /* 2: m_kpqq */
  /* 3: m_kpqp */
  /* 4: m_kqpp */
  /* 5: m_kqpq */
  /* 6: m_kqqq */
  /* 7: m_kqqp */
  unsigned m_myCataParam{10};

  TypeReplicator m_typeRep{};
  TypeComplex m_typeComp{free};

  unsigned bon_nei{0};
  Molecule *nei_ptr{nullptr};

public:
  Molecule() {
    for (int i{0}; i < std::size(m_rateList); i++) {
      m_rateList[i] *= Para::beta;
    }
  }

  Molecule(TypeReplicator typeR, TypeComplex typeC)
      : m_typeRep{typeR}, m_typeComp{typeC} {
    std::cout << "I'm born\n";
    for (int i{0}; i < std::size(m_rateList); i++) {
      m_rateList[i] *= Para::beta;
    }
  }

  ~Molecule() {}

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
  void setTypeRep(TypeReplicator myType) {
    m_typeRep = myType;
  } // setters no longer useful since newCA is now a friend class!

  friend class newCA;
};

#endif
