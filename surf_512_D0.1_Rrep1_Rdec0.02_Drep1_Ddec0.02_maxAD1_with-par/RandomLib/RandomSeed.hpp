/**
 * \file RandomSeed.hpp
 * \brief Header for RandomSeed
 *
 * This provides a base class for random generators.
 *
 * Written by <a href="http://charles.karney.info/"> Charles Karney</a>
 * <charles@karney.com> and licensed under the LGPL.  For more
 * information, see http://charles.karney.info/random/
 **********************************************************************/

#if !defined(RANDOMSEED_H)
#define RANDOMSEED_H

#define RCSID_RANDOMSEED_H "$Id: RandomSeed.hpp 6407 2007-05-24 17:13:46Z ckarney $"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>		// For std::transform
#include <sstream>		// For VectorToString

namespace RandomLib {
  /**
   * \brief A base class for random generators
   *
   * This provides facilities for managing the seed and for converting the seed
   * into random generator state.
   *
   * The seed is taken to be a vector of unsigned longs of arbitrary length.
   * (Only the low 32 bit of each element of the vector are used.)  The class
   * provides several methods for setting the seed, static functions for
   * producing "random" and "unique" seeds, and facilities for converting the
   * seed to a string so that it can be printed easily.
   *
   * The seeding algorithms are those used by
   * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
   * MT19937</a> with some modifications to make all states accessible and to
   * minimize the likelihood of different seeds giving the same state.
   *
   * Finally some low-level routines are provided to facilitate the creation of
   * I/O methods for the random generator.
   *
   * A random generator class can be written based on this class.  The
   * generator class would use the base class methods for setting the seed and
   * for converting the seed into state.  It would provide the machinery for
   * advancing the state and for producing random data.  It is also responsible
   * for the routine to save and restore the generator state (including the
   * seed).
   *
   * Written by <a href="http://charles.karney.info/"> Charles Karney</a>
   * <charles@karney.com> and licensed under the LGPL.  The seeding algorithms
   * are adapted from those of
   * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
   * MT19937</a>.  For more information, see http://charles.karney.info/random/
   **********************************************************************/
  class RandomSeed {
  public:
    virtual ~RandomSeed() throw() = 0;
    /**
     * A type large enough to hold the seed words.  This is needs to hold 32
     * bits and is an unsigned long for portability.
     **********************************************************************/
    typedef unsigned long seed_type;

    /** \name Resetting the seed
     **********************************************************************/
    ///@{
    /**
     * Set the seed to a vector \a v.  Only the low \e 32 bits of each element
     * are used.
     **********************************************************************/
    template<typename IntType> void Reseed(const std::vector<IntType>& v)
      throw(std::bad_alloc) {
      // Reserve space for new seed so as not to change object on error.
      _seed.reserve(v.size());
      _seed.resize(v.size());
      std::transform(v.begin(), v.end(), _seed.begin(), SeedMask);
      Reset();
    }
    /**
     * Set the seed to [\a a, \a b) from a pair of iterators.  The iterator
     * must produce results which can be converted into seed_type.  Only the
     * low 32 bits of each element are used.
     **********************************************************************/
    template<typename InputIterator>
    void Reseed(InputIterator a, InputIterator b)
      throw(std::bad_alloc) {
      // Read new seed into temporary so as not to change object on error.
      std::vector<seed_type> t;
      std::transform(a, b, back_inserter(t), SeedMask);
      _seed.swap(t);
      Reset();
    }
    /**
     * Set the seed to [\a n].  Only the low 32 bits of \a n are used.
     **********************************************************************/
    void Reseed(seed_type n) throw(std::bad_alloc) {
      // Reserve space for new seed so as not to change object on error.
      _seed.reserve(1);
      _seed.resize(1);
      _seed[0] = SeedMask(n);
      Reset();
    }
    /**
     * Set the seed to [SeedWord()]
     **********************************************************************/
    void Reseed() throw(std::bad_alloc) { Reseed(SeedWord()); }
    /**
     * Set the seed from the string \e s using Random::StringToVector
     **********************************************************************/
    void Reseed(const std::string& s) throw(std::bad_alloc) {
      // Read new seed into temporary so as not to change object on error.
      std::vector<seed_type> t = StringToVector(s);
      _seed.swap(t);
      Reset();
    }
    ///@}

    /** \name Examining the seed
     **********************************************************************/
    ///@{
    /**
     * Return reference to the seed vector (read-only)
     **********************************************************************/
    const std::vector<seed_type>& Seed() const throw() { return _seed; }
    /**
     * Format the current seed suitable for printing.
     **********************************************************************/
    std::string SeedString() const throw(std::bad_alloc)
    { return VectorToString(_seed); }
    ///@}

    /** \name Resetting the random seed
     **********************************************************************/
    ///@{
    /**
     * Resets the sequence to its just-seeded state.  This needs to be declared
     * vitual here so that the Reseed functions can call it after saving the
     * seed.
     **********************************************************************/
    virtual void Reset() throw() = 0;
    ///@}

    /** \name Writing to and reading from a stream
     **********************************************************************/
    ///@{
    /**
     * Save the state of the Random object to an output stream.  Format is a
     * sequence of unsigned 32-bit integers written either in decimal (\a bin
     * false, text format) or in network order with most significant byte first
     * (\a bin true, binary format).  Data consists of:
     *  - version identifier (1 word)
     *  - count of following words (1 word)
     *  - data (count words content depends on random generator)
     *
     * Declared virtual so that operator<< can be defined.
     **********************************************************************/
    virtual void Save(std::ostream& os, bool bin = true) const
      throw(std::ios::failure) = 0;
    /**
     * Restore the state of the Random object from an input stream.  If \a bin,
     * read in binary, else use text format.  Declared virtual so that
     * operator>> can be defined.
     **********************************************************************/
    virtual void Load(std::istream& is, bool bin = true)
      throw(std::ios::failure, std::out_of_range, std::bad_alloc) = 0;
    ///@}

    /** \name Static functions for seed management
     **********************************************************************/
    ///@{
    /**
     * Return a 32 bits of data suitable for seeding the random generator.  The
     * result is obtained by combining data from /dev/urandom, gettimeofday,
     * time, and getpid to provide a reasonably "random" word of data.
     **********************************************************************/
    static seed_type SeedWord();
    /**
     * Return a vector of unsigned longs suitable for seeding the random
     * generator.  The vector is almost certainly unique; however, the results
     * of successive calls to Random::SeedVector() will be correlated.  If
     * several Random objects are required within a single program execution,
     * call Random::SeedVector once, print it out (!), push_back additional
     * data to identify the instance (e.g., loop index, thread ID, etc.), and
     * use the result to seed the Random object.
     **********************************************************************/
    static std::vector<seed_type> SeedVector();
    /**
     * Convert a vector into a string suitable for printing or as an argument
     * for Random::Reseed(const std::string& s).
     **********************************************************************/
    template<typename IntType>
    static std::string VectorToString(const std::vector<IntType>& v)
      throw(std::ios::failure) {
      std::ostringstream os;
      os << "[";
      for (typename std::vector<IntType>::const_iterator n = v.begin();
	   n != v.end(); ++n) {
	if (n != v.begin())
	  os << ",";
	// Normalize in case this is called by user.
	os << SeedMask(*n);
      }
      os << "]";
      return os.str();
    }
    /**
     * Convert a string into a vector of seed_type suitable for printing or as
     * an argument for Random::Reseed(const std::vector<seed_type>& v).  Reads
     * consecutive digits in string.  Thus "[1,2,3]" => [1,2,3]; "-0.123e-4" =>
     * [0,123,4], etc.  strtoul understands C's notation for octal and
     * hexadecimal, for example "012 10 0xa" => [10,10,10].  Reading of a
     * number stops at the first illegal character for the base.  Thus
     * "2006-04-08" => [2006,4,0,8] (i.e., 08 becomes two numbers).  Note that
     * input numbers greater than ULONG_MAX overflow to ULONG_MAX, which
     * probably will result in the number being interpreted as LONG_MASK.
     **********************************************************************/
    static std::vector<seed_type> StringToVector(const std::string& s)
      throw(std::bad_alloc);
    ///@}

    /** \name Basic I/O
     **********************************************************************/
    ///@{
    /**
     * Write the state of a generator to stream \a os as text
     **********************************************************************/
    friend std::ostream& operator<<(std::ostream& os, const RandomSeed& r) {
      r.Save(os, false);
      return os;
    }

    /**
     * Read the state of a generator from stream \a is as text
     **********************************************************************/
    friend std::istream& operator>>(std::istream& is, RandomSeed& r) {
      r.Load(is, false);
      return is;
    }
    ///@}

    /**
     * An unsigned type which holds at least 32 bits for Ran32().
     **********************************************************************/
    typedef unsigned u32;
    /**
     * An unsigned type which holds at least 64 bits for Ran64().
     **********************************************************************/
    typedef unsigned long long u64;

  protected:
    /**
     * The seed vector
     **********************************************************************/
    std::vector<seed_type> _seed;

    /**
     * Convert the seed vector into random generator state.  The state is an
     * array of \a n elements of integer type \a UIntType.  The lowest \a
     * BitWidth bits of each element are set.  An important property of this
     * method is that "close" seeds should produce "widely separated" states.
     * This allows the seeds to be set is some systematic fashion to produce a
     * set of uncorrelation random number sequences.
     *
     * Three instantiations are provided:
     *  - \a UIntType = u32, \a BitWidth = 32.  The method is that given by
     *    Brown, et al.,
     *    <a href="http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2079.pdf">
     *    Random Number Generation in C++0X: A Comprehensive Proposal, version
     *    3</a>.  This is an adaption of the array initialization procedure for
     *    <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
     *    MT19937</a>.  The modifications are
     *    - All \a n elements are initialized, allowing all states to be
     *      accessed.
     *    - A defect which allowed the same states to be given by close seeds
     *      is corrected.
     *    - The seed length can be zero.
     *  - \a IntType = u64, \a BitWidth = 64.  This uses the array seeding
     *    method for the
     *    <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html">
     *    MT19937-64</a> making the same adaptions as for the 32-bit version.
     *    The seed is a vector of 32-bit elements and these are packed into a
     *    64-bit word for use with the seeding procedure.
     *  - \a UIntType = u32, \a BitWidth = 128.  This is the
     *    <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html">
     *    SFMT</a> routine init_by_array, given in
     *    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/cgi-bin/dl/dl.cgi?/MT/SFMT/SFMT-src-1.2.tar.gz
     *    This has better mixing properties that the original MT19937 methods
     *    and the Brown, et al., improvements at a modest (~ 40%) cost in
     *    speed.
     **********************************************************************/
    template<typename UIntType, int BitWidth>
    void SeedToState(UIntType state[], unsigned n) const throw();

    /**
     * A mask to give the low 32 bits of a seed_type.
     **********************************************************************/
    static const seed_type SEED_MASK = seed_type(0xffffffffUL);
    /**
     * A mask to give the low 32 bits in a u32.
     **********************************************************************/
    static const u32 U32_MASK = u32(0xffffffffUL);
    /**
     * A mask to give the low 64 bits in a u64.
     **********************************************************************/
    static const u64 U64_MASK = u64(0xffffffffffffffffULL);

    /**
     * Write integer in 32-bit chunks
     **********************************************************************/
    template<typename UIntType, unsigned BitWidth>
    static void Write32(std::ostream& os, bool bin, int& count, UIntType x)
      throw(std::ios::failure);
    /**
     * Read integer in 32-bit chunks
     **********************************************************************/
    template<typename UIntType, unsigned BitWidth>
    static void Read32(std::istream& is, bool bin, UIntType& x)
      throw(std::ios::failure);
    /**
     * Accumulate a checksum of a integer into a 32-bit check.  This implements
     * a very simple checksum and is intended to avoid accidental corruption
     * only.
     **********************************************************************/
    template<typename UIntType, unsigned BitWidth>
    static void CheckSum(UIntType n, u32& check);

    /**
     * Normalize seed words by masking out all but the low 32 bits.
     **********************************************************************/
    static seed_type SeedMask(seed_type x) { return x & SEED_MASK; }
  };

  inline RandomSeed::~RandomSeed() throw() {};

  /// \cond SKIP

  // Accumulate a checksum of a 32-bit quantity into check
  template<>
  inline void RandomSeed::CheckSum<RandomSeed::u32, 32u>(u32 n, u32& check) {
    // Circular shift left by one bit and add new word.
    check = (check << 1 | check >> 31 & u32(1)) + n;
    check &= U32_MASK;
  }

  // Accumulate a checksum of a 64-bit quantity into check
  template<>
  inline void RandomSeed::CheckSum<RandomSeed::u64, 64u>(u64 n, u32& check) {
    CheckSum<u32, 32u>(static_cast<u32>(n >> 32) & U32_MASK, check);
    CheckSum<u32, 32u>(static_cast<u32>(n      ) & U32_MASK, check);
  }

  /// \endcond
} // namespace RandomLib

#endif	// RANDOMSEED_H
