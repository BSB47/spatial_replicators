#include "newca.h"
#include "para.h"

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
    /* occupied, */
    cata,
    temp,
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

  double *m_ratePtr{m_rateList};

  TypeReplicator m_typeRep{};
  TypeReplicator m_repOutcome{};
  TypeComplex m_typeComp{};

  Molecule *bon_nei;

  double m_mutation_probability{};

public:
  Molecule() = default;

  Molecule(TypeReplicator typeR, TypeComplex typeC)
      : m_typeRep{typeR}, m_typeComp{typeC} {}

  Molecule(Molecule &&a) noexcept : m_ratePtr(a.m_ratePtr) {
    a.m_ratePtr = nullptr;
  }

  Molecule &operator=(Molecule &&a) noexcept {
    if (&a == this) // self-assignment check
      return *this;

    m_ratePtr = a.m_ratePtr;
    a.m_ratePtr = nullptr;

    return *this;
  }

  const TypeReplicator &getTypeReplicator() const { return m_typeRep; }
  void setTypeRep(TypeReplicator myType) { m_typeRep = myType; }

  friend int newCA::determineComplex(unsigned row, unsigned col, Molecule &mole,
                                     int mole_type);
};

#endif
