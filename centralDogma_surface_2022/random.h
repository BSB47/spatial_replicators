#include "para.h"
#include <chrono>
#include <random>

#ifndef RANDOM_H
#define RANDOM_H

namespace DiceRoller {
inline unsigned timeSeed{static_cast<unsigned int>(
    std::chrono::steady_clock::now().time_since_epoch().count())};
inline std::random_device rd;
inline std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
inline std::mt19937 twister{2998173820};
inline std::uniform_int_distribution typeInitializer{1, 3};
inline std::uniform_int_distribution randomRowOrCol{
    1, static_cast<int>(Para::sys_nrow)};

// Index should be from 1 to nrow/ncol!!! Boundaries are on
// index=0 and index=nrow+1/ncol+1

inline std::uniform_int_distribution randomNei{1, 8};
inline std::uniform_int_distribution randomNeiExcl{1, 7};
inline std::uniform_real_distribution<double> probabilityGen{0.0, 1.0};
inline std::uniform_real_distribution<double> mutaGen{
    -(Para::mutation_interval), Para::mutation_interval};

} // namespace DiceRoller

#endif
