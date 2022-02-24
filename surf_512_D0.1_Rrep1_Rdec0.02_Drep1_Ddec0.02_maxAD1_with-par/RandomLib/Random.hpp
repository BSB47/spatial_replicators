/**
 * \file Random.hpp
 * \brief Header for Random.
 *
 * Use the random bits from RandomGenerator to produce random integers
 * of various sizes, random reals with various precisions, a random
 * probability, etc.
 *
 * Written by <a href="http://charles.karney.info/">Charles Karney</a>
 * <charles@karney.com> and licensed under the LGPL.  For more
 * information, see http://charles.karney.info/random/
 **********************************************************************/

#if !defined(RANDOM_H)
#define RANDOM_H

#define RCSID_RANDOM_H "$Id: Random.hpp 6384 2007-04-24 11:43:38Z ckarney $"

#if defined(_MSC_VER)
#define WINDOWS 1
// Disable throw, precedence, size_t->int warnings
#pragma warning (disable: 4290 4554 4267)
#else
#define WINDOWS 0
#endif

#if defined(__sparc)
#define SUN 1
#else
#define SUN 0
#endif

#if !defined(HAVE_BOOST_SERIALIZATION)
/**
 * Use boost serialization?
 **********************************************************************/
#define HAVE_BOOST_SERIALIZATION 0
#endif

/**
 * Use table, Power2::power2, for pow2?  This isn't necessary with g++
 * 4.x because calls to std::pow are optimized.
 **********************************************************************/
#if defined(__GNUC__) && __GNUC__ >= 4
#define RANDOM_POWERTABLE 0
#else
// otherwise use a lookup table
#define RANDOM_POWERTABLE 1
#endif

#if WINDOWS
#define RANDOM_LONGDOUBLEPREC 53
#elif SUN
#define RANDOM_LONGDOUBLEPREC 113
#else
/**
 * The precision of long doubles, used for sizing Power2::power2.
 **********************************************************************/
#define RANDOM_LONGDOUBLEPREC 64
#endif

#if !defined(HAVE_SSE2)
#define HAVE_SSE2 0
#endif

/**
 * Are denormalized reals of type RealType supported?
 **********************************************************************/
#define RANDOM_HASDENORM(RealType) 1

#if !defined(STATIC_ASSERT)
/**
 * A simple compile-time assert.
 **********************************************************************/
#define STATIC_ASSERT(cond) { enum{ STATIC_ASSERT_ENUM = 1/int(cond) }; }
#endif

// This includes RandomLib/RandomSeed.hpp and RandomLib/MT19937.hpp
#include "RandomLib/RandomCanonical.hpp"

namespace RandomLib {

#if defined(USE_SFMT) && USE_SFMT == 1
  typedef SRandomGenerator32 RandomGenerator32;
  typedef SRandomGenerator64 RandomGenerator64;
  typedef SRandom32          Random32;
  typedef SRandom64          Random64;
#else
  typedef MRandomGenerator32 RandomGenerator32;
  typedef MRandomGenerator64 RandomGenerator64;
  typedef MRandom32          Random32;
  typedef MRandom64          Random64;
#endif

#if defined(MTPREC64) && MTPREC64 == 1
  typedef RandomGenerator64 RandomGenerator;
#else
  /**
   * Point RandomGenerator to one of the specific MT19937 generators.
   **********************************************************************/
  typedef RandomGenerator32  RandomGenerator;
#endif

  /**
   * Hook Random to RandomGenerator
   **********************************************************************/
  typedef RandomCanonical<RandomGenerator> Random;
} // namespace RandomLib

#endif	// RANDOM_H
