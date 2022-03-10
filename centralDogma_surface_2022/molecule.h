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
    occupied,
    temp,
    cata,
  };

private:
  double m_kppp{}; // last two letters stand for template/product
  double m_kppq{};
  double m_kpqq{};
  double m_kpqp{};
  double m_kqpp{};
  double m_kqpq{};
  double m_kqqq{};
  double m_kqqp{};

  int m_cor{}; // This is coordinate in 1D for now.

  TypeReplicator m_typeRep{};
  TypeComplex m_typeComp{};

  double m_mutation_probability{};

public:
  Molecule() = default;

  Molecule(int cor, TypeReplicator typeR = s, TypeComplex typeC = free)
      : m_kppp{1}, m_kppq{1}, m_kpqq{1}, m_kpqp{1}, m_kqpp{1}, m_kqpq{1},
        m_kqqq{1}, m_kqqp{1}, m_cor{cor}, m_typeRep{typeR}, m_typeComp{typeC} {}

  const TypeReplicator &getTypeReplicator() { return m_typeRep; }
  int getCor() { return m_cor; }
  void setTypeRep(TypeReplicator &myType) { m_typeRep = myType; }
};

#endif
