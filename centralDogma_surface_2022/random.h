#include "para.h"
#include <random>

#ifndef RANDOM_H
#define RANDOM_H

namespace DiceRoller {
inline std::random_device rd;
inline std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
inline std::mt19937 twister{ss};
inline std::uniform_int_distribution typeInitializer{1, 3};
inline std::uniform_int_distribution randomRowOrCol{
    0, (Para::sys_nrow - 1)}; // sys_nrow = sys_ncol for now so this is ok
                              // namespace DiceRoller
inline std::uniform_int_distribution randomNei{1, 8};
inline std::uniform_real_distribution probabilityGen{0.0, 1.0};

} // namespace DiceRoller

#endif
