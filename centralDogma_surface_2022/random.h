#include <random>

#ifndef RANDOM_H
#define RANDOM_H

namespace DiceRoller {
inline std::random_device rd;
inline std::seed_seq ss{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
inline std::mt19937 twister{ss};
inline std::uniform_int_distribution typeInitializer{1, 3};
}; // namespace DiceRoller

#endif
