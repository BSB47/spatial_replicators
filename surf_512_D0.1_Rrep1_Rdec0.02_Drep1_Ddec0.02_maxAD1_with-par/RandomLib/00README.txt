$Id: 00README.txt 6403 2007-05-21 14:58:59Z ckarney $

A random number library using the Mersenne Twister random number
generator.

Written by Charles Karney <charles@karney.com> and licensed under
the LGPL.  For more information, see

    http://charles.karney.info/random/

Files

    00README.txt -- this file
    Doxyfile -- Doxygen config file
    Random.doc -- main page of Doxygen documentation
    RandomSeed.hpp -- seed management
    MT19937.hpp -- implementation of Mersenne Twister
    SFMT19937.hpp -- implementation of SFMT generator
    RandomCanonical.hpp Power2.hpp -- Random integers, reals, booleans
    Random.hpp Random.cpp -- main include file plus implementation
    NormalDistribution.hpp -- sample from normal distribution
    ExponentialDistribution.hpp -- sample from exponential distribution
    RandomSelect.hpp -- sample from discrete distribution
    LeadingZeros.hpp -- count of leading zeros on random fraction
    ExponentialProb.hpp -- true with probability exp(-p)
    RandomNumber.hpp -- support for infinite precision randoms
    ExactExponential.hpp -- sample exactly from exponential distribution
    exphist.png exphist.pdf -- figures for documentation
    ExactPower.hpp -- sample exactly from power distribution
    powerhist.png powerhist.pdf -- figures for documentation
    Makefile -- Makefile (for Linux + gcc)
    RandomExample.cpp -- example code
    RandomPermutation.cpp -- prints a random permutation of integers
    Random.sln -- MS Visual C++ solution for library + examples
    RandomLib.vcproj -- MS Visual C++ project for library
    RandomExample.vcproj -- MS Visual C++ project for RandomExample
    RandomPermutation.vcproj -- MS Visual C++ project for RandomPermutation
    shuffle.sh -- shuffles the lines of a file

This is the 2007-05 version of the library.

Change between 2007-04 and 2007-05 versions:

 * Add SFMT19937 generators.

 * Introduce RandomGenerator::Name() to identify generator.

 * Change define used to make 64-bit generator the default.

 * Add RandomSelect::Weight.

 * Ensure portability to systems where RandomSeed::u32 is longer than 32
   bits.

Changes between 2006-12 and 2007-04 versions:

 * Add utilities RandomPermutation and shuffle.

 * Implement MSB ordering on binary I/O in a portable way.

Changes between 2006-11 and 2006-12 versions:

 * Add leapfrogging.  The output format needed to be changed to
   accommodate an extra word of data.  However, I/O routines can still
   read the 2006-11 version.

Changes between 2006-10 and 2006-11 versions:

 * Introduce RandomCanonical class which accepts the random generator
   as a template argument.

 * This allows the inclusion of 32-bit and 64-bit versions of mt19937.

 * Include checksum in I/O.

 * Include boost serialization.

Changes between 2006-09 and 2006-10 versions:

 * Make 64-bit ready so a 64-bit version of mt19937 can be dropped in.

 * Fix a bug in the seeding.  (This bug was trigged by seed length of
   624 or greater; so it was unlikely to have been encountered in
   practice.)

 * Stop the special case treatment for
   Random::IntegerC<T>(numeric_limits<T>::max()).  In some cases (e.g.,
   T = int) this now gives different (but equivalent) results.

Changes between 2006-08 and 2006-09 versions:

 * Add ExponentialProb, ExactExponential, ExactPower, and RandomNumber.

 * Fix weakness in the seeding algorithm.  A given seed now gives a
   random sequence different from previous version; so this is an
   incompatible change.

 * Restructure the documentation.

 * Allow constructors to accept vectors of any integral type and
   constructors with a pair of iterators.

Change between 2006-07 and 2006-08 versions:

 * Improve efficiency of Integer(n) where n is a power of two.
