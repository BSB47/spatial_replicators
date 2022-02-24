/**
 * \file MT19937.hpp
 * \brief Header for MT19937.
 *
 * This provides an interface to Mersenne Twister random number generator,
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">MT19937</a>.
 * See\n Makoto Matsumoto and Takuji Nishimura,\n Mersenne Twister: A
 * 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator,\n
 * ACM TOMACS 8, 3-30 (1998)
 *
 * Interface routines written by <a href="http://charles.karney.info/">Charles
 * Karney</a> <charles@karney.com> and licensed under the LGPL.  For more
 * information, see http://charles.karney.info/random/
 **********************************************************************/

#if !defined(MT19937_H)
#define MT19937_H

#define RCSID_MT19937_H "$Id: MT19937.hpp 6406 2007-05-23 13:29:52Z ckarney $"

#include "RandomLib/RandomSeed.hpp"
#include <limits>

#if HAVE_BOOST_SERIALIZATION
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#endif

/**
 * Which version of seeding and saving to use.  3 = scheme proposed in C++0X
 * random number proposal, version 3 (2006-09) plus append a checksum to data
 * stream (introduced with 2006-11 version).  4 = add leapfrogging (introduced
 * with 2006-12 version).  In this implementation, I/O routines can read
 * versions 3 and 4.
 **********************************************************************/
#define MT19937_VERSION 4

namespace RandomLib {
  /**
   * \brief A generator of random bits.
   *
   * This provides an interface to Mersenne Twister random number generator,
   * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
   * MT19937</a>.  See\n Makoto Matsumoto and Takuji Nishimura,\n Mersenne
   * Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number
   * Generator,\n ACM TOMACS 8, 3-30 (1998)
   *
   * This is adapted from the 32-bit and 64-bit C versions available at
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html and
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
   *
   * The template arguments give the type \a UIntType of the "natural" result
   * and the number of random bits \a BitWidth in each result.  Although the
   * two versions of MT19937 produce different sequences, the implementations
   * here are portable across 32-bit and 64-bit architectures.  As expected the
   * 64-bit version of MT19937 runs somewhat faster on 64-bit machines and
   * considerably slower on 32-bit machines.
   *
   * The class supplies the method for advancing the state and generating the
   * random output.  The random results are provided as 32-bit quantities via
   * Ran32(), as 64-bit quantities via Ran64(), and in "natural" units of \a
   * width bits via Ran().  For MT19937, \a width is either 32 or 64.  It
   * inherits from RandomSeed the routines for seed management and for
   * converting the seed into state.
   *
   * It also provides routines for saving and restoring the state of the
   * generator (including its seed), for stepping the generator forwards or
   * backwards, and for testing the equality of two generators.
   *
   * Interface routines written by <a href="http://charles.karney.info/">
   * Charles Karney</a> <charles@karney.com> and licensed under the LGPL.  For
   * more information, see http://charles.karney.info/random/
   **********************************************************************/
  template<typename UIntType, unsigned BitWidth>
  class MT19937 : public RandomSeed {
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
       * The long lag for MT19937
       **********************************************************************/
      N = (P + width - 1)/width,
      /**
       * The short lag for MT19937
       **********************************************************************/
      M = width == 32 ? 397 : 156,
      /**
       * The number of ignored bits in the first word of the state
       **********************************************************************/
      R = N * width - P
    };
    /**
     * "MTrn" or "MTsn" a signature for Save
     **********************************************************************/
    static const u32 VERSIONID =
      u32(width == 32 ? 0x4d547230UL : 0x4d547330UL) + MT19937_VERSION;
    /**
     * Magic matrix for MT19937
     **********************************************************************/
    static const result_type MATRIX_A =
      result_type(width == 32 ? 0x9908b0dfULL : 0xb5026f5aa96619e9ULL);
    /**
     * Mask for top \a width - \a R bits of a word
     **********************************************************************/
    static const result_type UPPER_MASK = RESULT_MASK << R & RESULT_MASK;
    /**
     * Mask for low \a R bits of a <i>width</i>-bit word
     **********************************************************************/
    static const result_type LOWER_MASK = ~UPPER_MASK & RESULT_MASK;
    /**
     * Marker for uninitialized object
     **********************************************************************/
    static const unsigned UNINIT = 0xffffffffU;

    /**
     * The state vector
     **********************************************************************/
    result_type _state[N];
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
    template<typename IntType> explicit MT19937(const std::vector<IntType>& v)
      throw(std::bad_alloc) { Reseed(v); }
    /**
     * Initialize from a pair of iterators setting seed to [\a a, \a b).  The
     * iterator must produce results which can be converted into seed_type.
     * Only the low \e 32 bits of each element are used.
     **********************************************************************/
    template<typename InputIterator> MT19937(InputIterator a, InputIterator b)
    { Reseed(a, b); }
    /**
     * Initialize with seed [\a n].  Only the low \e width bits of \a n are
     * used.
     **********************************************************************/
    explicit MT19937(seed_type n) throw(std::bad_alloc)
    { Reseed(n); }
    /**
     * Initialize with seed [SeedWord()]
     **********************************************************************/
    MT19937() throw(std::bad_alloc) { Reseed(); }
    /**
     * Initialize from a string.  See Reseed(const std::string& s)
     **********************************************************************/
    explicit MT19937(const std::string& s) throw(std::bad_alloc)
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
    bool operator==(const MT19937& r) const throw()
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
    bool operator!=(const MT19937& r) const throw()
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
      MT19937 t(is, bin);
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
	throw std::invalid_argument("MT19937: Invalid stride");
      if (k >= n)
	throw std::invalid_argument("MT19937: Invalid index");
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
      s << "RandomLib::MT19937 " << width << "-bit "
	<< "Version " << MT19937_VERSION;
      return s.str();
    }

    /**
     * Tests basic engine.  Throws out_of_range errors on bad results.
     **********************************************************************/
    static void SelfTest() {
      MT19937 g(std::vector<seed_type>(0));
      g.SetCount(10000-1);
      if (g() !=
	  result_type(width == 32 ? 4123659995ULL : 9981545732273789042ULL))
	throw std::out_of_range("MT19937: Incorrect result with seed " +
				g.SeedString());
      seed_type s[] = {1, 2, 3, 4};
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
	MT19937 h(std::vector<seed_type>(0));
	h.Load(stream, true);
	h.SetCount(1000000-1);
	if (h() !=
	    result_type(width == 32 ? 4244229310ULL : 16856309709668453709ULL))
	throw std::out_of_range("MT19937: Incorrect result with seed " +
				h.SeedString());
	g.SetCount(1000000);
	if (h != g)
	throw std::out_of_range("MT19937: Comparison failure");
      }
    }

  protected:
    /**
     * Return \e width bits of randomness.  This is the natural unit of random
     * data produced by MT19937.
     **********************************************************************/
    result_type Ran() throw();
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
    explicit MT19937(std::istream& is, bool bin)
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
	  &  boost::serialization::make_nvp("state" , _state   );
      ar & boost::serialization::make_nvp("checksum", _checksum);
    }
    /**
     * Load from a boost archive.  Do this safely so that the current object is
     * not corrupted if the archive is bogus.
     **********************************************************************/
    template<class Archive> void load(Archive& ar, const unsigned int) {
      u32 _version, _checksum;
      ar & boost::serialization::make_nvp("version" , _version );
      MT19937 t(std::vector<seed_type>(0));
      ar & boost::serialization::make_nvp("seed"    , t._seed  )
	&  boost::serialization::make_nvp("ptr"     , t._ptr   );
      if (_version == VERSIONID - 1) {
	if (t._ptr == N + 1)
	  t._ptr = UNINIT;
	t._stride = 1;
      } else
	ar &  boost::serialization::make_nvp("stride"  , t._stride);
      if (t._ptr != UNINIT)
	ar & boost::serialization::make_nvp("rounds", t._rounds)
	  &  boost::serialization::make_nvp("state" , t._state );
      ar & boost::serialization::make_nvp("checksum", _checksum);
      if (t.Check(_version) != _checksum)
	throw std::out_of_range("MT19937: Checksum failure");
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
  inline RandomSeed::u32 MT19937<RandomSeed::u32, 32>::Ran() throw() {
    // On exit we have _stride <= _ptr < N + _stride.
    if (_ptr >= N)
      Next();
    result_type y = _state[_ptr];
    _ptr += _stride;

    // Specific tempering instantiation for width = 32 given in
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
    y ^= y >> 11;
    y ^= y <<  7 & result_type(0x9d2c5680UL);
    y ^= y << 15 & result_type(0xefc60000UL);
    y ^= y >> 18;

    return y;
  }

  template<>
  inline RandomSeed::u32 MT19937<RandomSeed::u32, 32>::Ran32() throw()
  { return Ran(); }

  template<>
  inline RandomSeed::u64 MT19937<RandomSeed::u32, 32>::Ran64() throw()
  { const u64 x = Ran(); return x << width | static_cast<u64>(Ran()); }

  template<>
  inline RandomSeed::u64 MT19937<RandomSeed::u64, 64>::Ran() throw() {
    if (_ptr >= N)
      Next();
    result_type y = _state[_ptr];
    _ptr += _stride;

    // Specific tempering instantiation for width = 64 given in
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
    y ^= y >> 29 & result_type(0x5555555555555555ULL);
    y ^= y << 17 & result_type(0x71d67fffeda60000ULL);
    y ^= y << 37 & result_type(0xfff7eee000000000ULL);
    y ^= y >> 43;

    return y;
  }

  template<>
  inline RandomSeed::u32 MT19937<RandomSeed::u64, 64>::Ran32() throw()
  { return static_cast<u32>(Ran()) & U32_MASK; }

  template<>
  inline RandomSeed::u64 MT19937<RandomSeed::u64, 64>::Ran64() throw()
  { return Ran(); }
} // namespace RandomLib
#endif	// MT19937_H
