#include "para.h"
#include <random>

#ifndef RANDOM_H
#define RANDOM_H

namespace DiceRoller {
inline std::random_device rd;
inline std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
inline std::mt19937 twister{ss};
inline std::uniform_int_distribution typeInitializer{1, 3};
inline std::uniform_int_distribution randomRowOrCol{1, Para::sys_nrow};

// Index should be from 1 to nrow/ncol!!! Boundaries are on
// index=0 and index=nrow+1/ncol+1

inline std::uniform_int_distribution randomNei{1, 8};
inline std::uniform_int_distribution randomNeiExcl{1, 7};
inline std::uniform_real_distribution<double> probabilityGen{0.0, 1.0};

} // namespace DiceRoller

#endif
