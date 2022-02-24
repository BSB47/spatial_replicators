/**
 * \file SFMT19937.hpp
 * \brief Header for SFMT19937.
 *
 * This provides an implementation of the SIMD-oriented Fast Mersenne Twister
 * random number generator,
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html">
 * SFMT19937</a>.  See\n Mutsuo Saito,\n An Application of Finite Field: Design
 * and Implementation of 128-bit Instruction-Based Fast Pseudorandom Number
 * Generator,\n Master's Thesis, Dept. of Math., Hiroshima University
 * (Feb. 2007),\n
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/M062821.pdf\n
 * Mutsuo Saito and Makoto Matsumoto,\n
 * SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number
 * Generator,\n accepted in the proceedings of MCQMC2006\n
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/sfmt.pdf
 *
 * Written by <a href="http://charles.karney.info/">Charles Karney</a>
 * <charles@karney.com> and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/random/
 **********************************************************************/

#if !defined(SFMT19937_H)
#define SFMT19937_H

#define RCSID_SFMT19937_H "$Id: SFMT19937.hpp 6406 2007-05-23 13:29:52Z ckarney $"

#include "RandomLib/RandomSeed.hpp"
#include <limits>
#if HAVE_SSE2
#include <emmintrin.h>
#endif

#if HAVE_BOOST_SERIALIZATION
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#endif

/**
 * Which version of seeding and saving to use.  0 = scheme implemented in
 * SFMT-src-1.2.tar.gz by M. Saito.
 **********************************************************************/
#define SFMT19937_VERSION 0

namespace RandomLib {
  /**
   * \brief A generator of random bits.
   *
   * This provides an implementation of the SIMD-oriented Fast Mersenne Twister
   * random number generator,
   * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html">
   * SFMT</a>.  See\n Mutsuo Saito,\n An Application of Finite Field: Design
   * and Implementation of 128-bit Instruction-Based Fast Pseudorandom Number
   * Generator,\n Master's Thesis, Dept. of Math., Hiroshima University
   * (Feb. 2007).\n
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/M062821.pdf
   * Mutsuo Saito and Makoto Matsumoto,\n
   * SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number
   * Generator,\n accepted in the proceedings of MCQMC2006\n
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/sfmt.pdf
   *
   * The template arguments give the type \a UIntType of the "natural" result
   * and the number of random bits \a BitWidth in each result.  The same
   * underlying SFMT19937 algorithm is used in both cases and this is defined
   * in terms of 128-bit words.  BitWidth controls the size of the words used
   * to implement the algorithm and the size of the result returned by Ran().
   *
   * The class supplies the method for advancing the state and generating the
   * random output.  The random results are provided as 32-bit quantities via
   * Ran32(), as 64-bit quantities via Ran64(), and in "natural" units of \a
   * width bits via Ran().  For SFMT19937, \a width is either 32 or 64.  It
   * inherits from RandomSeed the routines for seed management and for
   * converting the seed into state.
   *
   * With the same starting seed, SRandom32::Ran64() and SRandom64::Ran64()
   * produces the same sequences.  Similarly SRandom64::Ran32() produces every
   * other member of the sequence produced by SRandom32::Ran32().
   *
   * This class also provides routines for saving and restoring the state of
   * the generator (including its seed), for stepping the generator forwards or
   * backwards, and for testing the equality of two generators.
   *
   * Written by <a href="http://charles.karney.info/"> Charles Karney</a>
   * <charles@karney.com> and licensed under the LGPL.  For more information,
   * see http://charles.karney.info/random/
   **********************************************************************/
  template<typename UIntType, unsigned BitWidth>
  class SFMT19937 : public RandomSeed {
  public:
    enum {
      /**
       * The number of random bits produced by Ran().
       **********************************************************************/
      width = BitWidth
    };
    /**
     * A type large enough to hold \e width bits.  This is used for the
     * internal state of the generator and the result returned by Ran().
     * (Seeds are specified in terms of the more portable unsigned long.)
     **********************************************************************/
    typedef UIntType result_type;
    /**
     * The minimum result returned by Ran().
     **********************************************************************/
    static const result_type min = 0;
    /**
     * The maximum result returned by Ran() = 2<sup><i>w</i></sup> - 1
     **********************************************************************/
    static const result_type max =
      ~result_type(0) >> std::numeric_limits<result_type>::digits - width;

  protected:
    /**
     * A mask to give the low \e width bits in a result_type.
     **********************************************************************/
    static const result_type RESULT_MASK = max;

  private:
    // Constants
    enum {
      /**
       * The Mersenne prime is 2<sup><i>P</i></sup> - 1
       **********************************************************************/
      P = 19937,
      /**
       * The long lag for SFMT19937 in units of 128-bit words
       **********************************************************************/
      N128 = (P + 128 - 1)/128,
      /**
       * How many width words per 128-bit word.
       **********************************************************************/
      R = 128 / width,
      /**
       * The long lag for SFMT19937
       **********************************************************************/
      N = N128 * R,
      /**
       * The short lag for SFMT19937  in units of 128-bit words
       **********************************************************************/
      M128 = 122,
      /**
       * The short lag for SFMT19937 
       **********************************************************************/
      M = M128 * R
    };
    /**
     * "SMTn" or "SMUn" a signature for Save
     **********************************************************************/
    static const u32 VERSIONID =
      u32(width == 32 ? 0x534d5430UL : 0x534d5530UL) + SFMT19937_VERSION;
    /**
     * Magic matrix for SFMT19937.  Only the low 21 (= 32 - 11) bits need to be
     * set.  (11 is the right shift applied to the words before masking.
     **********************************************************************/
    static const u32 MSK0 = u32(0x1fffefUL),
      MSK1 = u32(0x1ecb7fUL),
      MSK2 = u32(0x1affffUL),
      MSK3 = u32(0x1ffff6UL);
    /**
     * 64-bit versions of MSK.
     **********************************************************************/
    static const u64 LMSK0 = u64(MSK0) | u64(MSK1) << 32,
      LMSK1 = u64(MSK2) | u64(MSK3) << 32;
    /**
     * Mask for simulating u32 << 18 with 64-bit words
     **********************************************************************/
    static const u64 LMSK18 = u64(0xfffc0000fffc0000ULL);
    /**
     * Magic constants needed by "period certification"
     **********************************************************************/
    static const u32 PARITY0 = u32(1UL),
      PARITY1 = u32(0UL),
      PARITY2 = u32(0UL),
      PARITY3 = u32(0x13c9e684UL);
    /**
     * Least significant bit of PARITY
     **********************************************************************/
    static const unsigned PARITY_LSB = 0;
    /**
     * Marker for uninitialized object
     **********************************************************************/
    static const unsigned UNINIT = 0xffffffffU;

    /**
     * \brief Union for SFMT19937's state.
     *
     * A union to hold the state as words of width \a width bits, as 128-bit
     * words (for the SSE2 implemenation) and 32-bit words for seeding.
     **********************************************************************/
    union state {
#if HAVE_SSE2
      /**
       * 128-bit representation.
       **********************************************************************/
      __m128i v[N128];
#endif
      /**
       * 32-bit representation.
       **********************************************************************/
      u32 u[N128 * 4];	// Needed for SeedToState
      /**
       * <i>width</i>-bit representation.
       **********************************************************************/
      result_type w[N];
    };
    /**
     * The state vector
     **********************************************************************/
    state _state;
    /**
     * The index for the next random value
     **********************************************************************/
    unsigned _ptr;
    /**
     * How many times has NextBatch() been called
     **********************************************************************/
    long long _rounds;
    /**
     * Stride for leapfrogging
     **********************************************************************/
    unsigned _stride;

  public:

    /** \name Constructors which set the seed
     **********************************************************************/
    ///@{
    /**
     * Initialize from a vector.  Only the low \e 32 bits of each element are
     * used.
     **********************************************************************/
    template<typename IntType>
    explicit SFMT19937(const std::vector<IntType>& v)
      throw(std::bad_alloc) { Reseed(v); }
    /**
     * Initialize from a pair of iterators setting seed to [\a a, \a b).  The
     * iterator must produce results which can be converted into seed_type.
     * Only the low \e 32 bits of each element are used.
     **********************************************************************/
    template<typename InputIterator>
    SFMT19937(InputIterator a, InputIterator b)
    { Reseed(a, b); }
    /**
     * Initialize with seed [\a n].  Only the low \e width bits of \a n are
     * used.
     **********************************************************************/
    explicit SFMT19937(seed_type n) throw(std::bad_alloc)
    { Reseed(n); }
    /**
     * Initialize with seed [SeedWord()]
     **********************************************************************/
    SFMT19937() throw(std::bad_alloc) { Reseed(); }
    /**
     * Initialize from a string.  See Reseed(const std::string& s)
     **********************************************************************/
    explicit SFMT19937(const std::string& s) throw(std::bad_alloc)
    { Reseed(s); }
    ///@}

    /** \name Functions for returning random data
     **********************************************************************/
    ///@{
    /**
     * Return \e width bits of randomness.  Result is in [0,
     * 2<sup><i>w</i></sup>)
     **********************************************************************/
    result_type operator()() throw() { return Ran(); }
    ///@}

    /** \name Comparing Random objects
     **********************************************************************/
    ///@{
    /**
     * Test equality of two Random objects.  This test that the seeds match and
     * that they have produced the same number of random numbers.
     **********************************************************************/
    bool operator==(const SFMT19937& r) const throw()
    // Ensure that the two Random objects behave the same way.  Note however
    // that the internal states may still be different, e.g., the following all
    // result in Random objects which are == (with Count() == 0) but which all
    // have different internal states:
    //
    // Random r(0);                       _ptr == UNINIT
    // r.StepCount( 1); r.StepCount(-1);  _ptr == 0, _rounds ==  0
    // r.StepCount(-1); r.StepCount( 1);  _ptr == N, _rounds == -1
    { return Count() == r.Count() && _seed == r._seed &&
	_stride == r._stride; }
    /**
     * Test inequality of two Random objects.  See Random::operator==
     **********************************************************************/
    bool operator!=(const SFMT19937& r) const throw()
    { return !operator==(r); }
    ///@}

    /** \name Writing to and reading from a stream
     **********************************************************************/
    ///@{
    /**
     * Save the state of the Random object to an output stream.  Format is a
     * sequence of unsigned 32-bit integers written either in decimal (\a bin
     * false, text format) or in network order with most significant byte first
     * (\a bin true, binary format).  Data consists of:
     *
     *  - VERSIONID (1 word)
     *  - _seed.size() (1 word)
     *  - _seed data (_seed.size() words)
     *  - _ptr (1 word)
     *  - _stride (1 word)
     *  - if _ptr != UNINIT, _rounds (2 words)
     *  - if _ptr != UNINIT, _state (N words or 2 N words)
     *  - checksum
     *
     * Shortest possible saved result consists of 5 words:
     *  - VERSIONID = "MTrn"
     *  - _seed.size() = 0
     *  - _ptr = UNINIT
     *  - _stride = 1
     *  - checksum
     *
     * This corresponds to Seed() = [] and Count() = 0.
     **********************************************************************/
    void Save(std::ostream& os, bool bin = true) const
      throw(std::ios::failure);
    /**
     * Restore the state of the Random object from an input stream.  If \a bin,
     * read in binary, else use text format.  See documentation of Random::Save
     * for the format.  Include error checking on date to make sure the input
     * has not been corrupted.  If an error occurs while reading, the Random
     * object is unchanged.
     **********************************************************************/
    void Load(std::istream& is, bool bin = true)
      throw(std::ios::failure, std::out_of_range, std::bad_alloc) {
      // Read state into temporary so as not to change object on error.
      SFMT19937 t(is, bin);
      _seed.reserve(t._seed.size());
      *this = t;
    }
    ///@}

    /** \name Examining and advancing the Random generator
     **********************************************************************/
    ///@{
    /**
     * Return the number of random numbers used.  This needs to return a long
     * long result since it can reasonably exceed 2<sup>31</sup>.  (On a 1GHz
     * machine, it takes about a minute to produce 2<sup>32</sup> random
     * numbers.)  More precisely this is the (zero-based) index of the next
     * random number to be produced.  (This distinction is important when
     * leapfrogging is in effect.)
     **********************************************************************/
    long long Count() const throw()
    { return _ptr == UNINIT ? 0 : _rounds * N + _ptr; }
    /**
     * Step the generator forwards of backwarks so that the value returned
     * by Count() is \a n
     **********************************************************************/
    void SetCount(long long n) throw() { StepCount(n - Count()); }
    /**
     * Step the generator forward \a n steps.  \a n can be negative.
     **********************************************************************/
    void StepCount(long long n) throw();
    /**
     * Resets the sequence.  Equivalent to SetCount(0), but works by
     * reinitializing the Random object from its seed, rather than by stepping
     * the sequence backwards.  In addition, this undoes leapfrogging.
     **********************************************************************/
    void Reset() throw() { _ptr = UNINIT; _stride = 1; }
    ///@}

    /** \name Leapfrogging
     **********************************************************************/
    ///@{
    /**
     * Set leapfrogging stride to a positive number \a n and increment Count()
     * by \a k < \a n.  If the current Count() is \a i, then normally the next
     * 3 random numbers would have indices \a i, \a i + 1, \a i + 2, and the
     * new Count() is \a i + 2.  However, after SetStride(\a n, \a k) the next
     * 3 random numbers have indices \a i + \a k, \a i + \a k + \a n, \a i + \a
     * k + 2\a n, and the new Count() is \a i + \a k + 3\a n.  With
     * leapfrogging in effect, the time to produce raw random numbers is
     * roughly proportional to 1 + (\a n - 1)/2.  Reseed(...) and Reset() both
     * reset the stride back to 1.  See \ref leapfrog "Leapfrogging" for a
     * description of how to use this facility.
     **********************************************************************/
    void SetStride(unsigned n = 1, unsigned k = 0)
      throw(std::invalid_argument) {
      // Limit stride to UNINIT/2.  This catches negative numbers that have
      // been cast into unsigned.  In reality the stride should be no more than
      // 10-100.
      if (n == 0 || n > UNINIT/2)
	throw std::invalid_argument("SFMT19937: Invalid stride");
      if (k >= n)
	throw std::invalid_argument("SFMT19937: Invalid index");
      _stride = n;
      StepCount(k);
    }
    /**
     * Return leapfrogging stride.
     **********************************************************************/
    unsigned GetStride() const throw() { return _stride; }
    ///@}

    /**
     * Return string description of this generator
     **********************************************************************/
    static std::string Name() {
      std::ostringstream s;
      s << "RandomLib::SFMT19937 " << width << "-bit "
	<< "Version " << SFMT19937_VERSION;
      return s.str();
    }

    /**
     * Tests basic engine.  Throws out_of_range errors on bad results.
     **********************************************************************/
    static void SelfTest() {
      SFMT19937 g(std::vector<seed_type>(0));
      g.SetCount(10000-1);
      if (g() !=
	  result_type(width == 32 ? 2695024307ULL : 1464461649847485149ULL))
	throw std::out_of_range("SFMT19937: Incorrect result with seed " +
				g.SeedString());
      seed_type s[] = {0x1234U, 0x5678U, 0x9abcU, 0xdef0U};
      g.Reseed(s, s+4);
      g.SetCount(-10000);
      std::string save;
      {
	std::ostringstream stream;
	stream << g << std::endl;
	save = stream.str();
      }
      g.Reset();
      {
	std::istringstream stream(save);
	stream >> g;
      }
      g.SetCount(10000);
      {
	std::ostringstream stream;
	g.Save(stream, true);
	save = stream.str();
      }
      {
	std::istringstream stream(save);
	SFMT19937 h(std::vector<seed_type>(0));
	h.Load(stream, true);
	h.SetCount(1000000-1);
	if (h() !=
	    result_type(width == 32 ? 782200760ULL : 5050640804923595109ULL))
	throw std::out_of_range("SFMT19937: Incorrect result with seed " +
				h.SeedString());
	g.SetCount(1000000);
	if (h != g)
	throw std::out_of_range("SFMT19937: Comparison failure");
      }
    }
    
  protected:
    /**
     * Return \e width bits of randomness.  This is the natural unit of random
     * data produced by SFMT19937.
     **********************************************************************/
    result_type Ran() throw() {
      if (_ptr >= N)
	Next();
      result_type y = _state.w[_ptr];
      _ptr += _stride;

      return y;
    }

    /**
     * Return 32 bits of randomness.
     **********************************************************************/
    u32 Ran32() throw();
    /**
     * Return 64 bits of randomness.
     **********************************************************************/
    u64 Ran64() throw();

  private:
    /**
     * Compute initial state from seed
     **********************************************************************/
    void Init() throw();
    /**
     * Advance state by \e N steps
     **********************************************************************/
    void NextBatch(unsigned long long count) throw();
    /**
     * Back up state by \e N steps
     **********************************************************************/
    void PreviousBatch(unsigned long long count) throw();
    /**
     * The interface to NextBatch used by Ran().
     **********************************************************************/
    void Next() throw();

    /**
     * Consistency check on state; mainly used for validating I/O.  This throws
     * an out-of-range error if the state is bad.  It also computes a simple
     * checksum to ensure the integrity of restored data.
     **********************************************************************/
    u32 Check(u32 version) const throw(std::out_of_range);

    /**
     * Read from an input stream.  Potentially corrupts object.  This is used
     * by Load so that it can avoid corrupting its state on bad input.
     **********************************************************************/
    explicit SFMT19937(std::istream& is, bool bin)
      throw(std::ios::failure, std::out_of_range, std::bad_alloc);

#if HAVE_BOOST_SERIALIZATION
    friend class boost::serialization::access;
    /**
     * Save to a boost archive.  Boost versioning isn't very robust.  (It
     * allows a RandomGenerator32 to be read back in as a RandomGenerator64.
     * It doesn't interact well with templates.)  So we do our own versioning
     * and supplement this with a checksum.
     **********************************************************************/
    template<class Archive> void save(Archive& ar, const unsigned int) const {
      u32 _version = VERSIONID,
	_checksum = Check(_version);
      ar & boost::serialization::make_nvp("version" , _version )
	&  boost::serialization::make_nvp("seed"    , _seed    )
	&  boost::serialization::make_nvp("ptr"     , _ptr     )
	&  boost::serialization::make_nvp("stride"  , _stride  );
      if (_ptr != UNINIT)
	ar & boost::serialization::make_nvp("rounds", _rounds  )
	  &  boost::serialization::make_nvp("state" , _state.w );
      ar & boost::serialization::make_nvp("checksum", _checksum);
    }
    /**
     * Load from a boost archive.  Do this safely so that the current object is
     * not corrupted if the archive is bogus.
     **********************************************************************/
    template<class Archive> void load(Archive& ar, const unsigned int) {
      u32 _version, _checksum;
      ar & boost::serialization::make_nvp("version" , _version );
      SFMT19937 t(std::vector<seed_type>(0));
      ar & boost::serialization::make_nvp("seed"    , t._seed  )
	&  boost::serialization::make_nvp("ptr"     , t._ptr   );
      if (_version == VERSIONID - 1) {
	if (t._ptr == N + 1)
	  t._ptr = UNINIT;
	t._stride = 1;
      } else
	ar &  boost::serialization::make_nvp("stride"  , t._stride);
      if (t._ptr != UNINIT)
	ar & boost::serialization::make_nvp("rounds", t._rounds )
	  &  boost::serialization::make_nvp("state" , t._state.w);
      ar & boost::serialization::make_nvp("checksum", _checksum);
      if (t.Check(_version) != _checksum)
	throw std::out_of_range("SFMT19937: Checksum failure");
      _seed.reserve(t._seed.size());
      *this = t;
    }
    /**
     * Glue the boost save and load functionality together -- a bit of boost
     * magic.
     **********************************************************************/
    template<class Archive>
    void serialize(Archive &ar, const unsigned int file_version)
    { boost::serialization::split_member(ar, *this, file_version); }
#endif
  };

  template<>
  inline RandomSeed::u32 SFMT19937<RandomSeed::u32, 32>::Ran32() throw()
  { return Ran(); }

  template<>
  inline RandomSeed::u64 SFMT19937<RandomSeed::u32, 32>::Ran64() throw()
  // Note that the ordering convention is the opposite of MT19937::Rand64().
  // This is so that Ran64() produces the same sequences with both the 32- and
  // 64-bit versions of SFMT19937.
  { const u64 x = Ran(); return static_cast<u64>(Ran()) << width | x; }

  template<>
  inline RandomSeed::u32 SFMT19937<RandomSeed::u64, 64>::Ran32() throw()
  { return static_cast<u32>(Ran()) & U32_MASK; }

  template<>
  inline RandomSeed::u64 SFMT19937<RandomSeed::u64, 64>::Ran64() throw()
  { return Ran(); }
} // namespace RandomLib
#endif	// SFMT19937_H
