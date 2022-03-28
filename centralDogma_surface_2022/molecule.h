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

  TypeReplicator m_typeRep{};
  TypeComplex m_typeComp{};

  double m_mutation_probability{};

public:
  Molecule() = default;

  Molecule(TypeReplicator typeR, TypeComplex typeC)
      : m_kppp{1}, m_kppq{1}, m_kpqq{1}, m_kpqp{1}, m_kqpp{1}, m_kqpq{1},
        m_kqqq{1}, m_kqqp{1}, m_typeRep{typeR}, m_typeComp{typeC} {}

  const TypeReplicator &getTypeReplicator() const { return m_typeRep; }
  void setTypeRep(TypeReplicator myType) { m_typeRep = myType; }
};

#endif
