/**
 * \file Random.cpp
 * \brief Code for RandomSeed, MT19937, SFMT19937, RandomCanonical, and Power2
 *
 * RandomSeed provides support for seed management, MT19937 implements the
 * Mersenne-Twister random number generator, SFMT19937 implements the
 * SIMD-oriented Fast Mersenne-Twister random number generator, RandomCanonical
 * supplies uniform distributions of integers and reals, and Power2 multiplies
 * reals by powers of two.
 *
 * The Mersenne Twister random number generator,
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html">
 * mt19937</a> , is described in \n Makoto Matsumoto and Takuji Nishimura,\n
 * Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random
 * Number Generator,\n ACM TOMACS 8, 3-30 (1998)
 *
 * The code for NextBatch is adapted from the code for genrand_int32 by Makoto
 * Matsumoto and Takuji Nishimura in
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
 *
 * PreviousBatch is adapted from revand by Katsumi Hagita.  See
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/REVmt19937b.f
 * The code for SeedToState is adapted from init_genrand{,64} and
 * init_by_array{,64} by Makoto Matsumoto and Takuji Nishimura in
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c and
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
 * following the seed_seq class given in
 * http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2079.pdf
 *
 * The
 * <a href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html">
 * SFMT</a> algorithm is described in\n
 * Mutsuo Saito,\n An Application of Finite Field: Design and Implementation of
 * 128-bit Instruction-Based Fast Pseudorandom Number Generator,\n Master's
 * Thesis, Dept. of Math.  Hiroshima University (Feb. 2007),\n
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/M062821.pdf\n
 * Mutsuo Saito and Makoto Matsumoto,\n
 * SIMD-oriented Fast Mersenne Twister: a 128-bit Pseudorandom Number
 * Generator,\n accepted in the proceedings of MCQMC2006\n
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/sfmt.pdf
 *
 * The adaption to the C++ and the rest of the code was written by
 * <a href="http://charles.karney.info/"> Charles Karney</a>
 * <charles@karney.com> and licensed under the LGPL.  For more information, see
 * http://charles.karney.info/random/
 **********************************************************************/

#include "RandomLib/Random.hpp"

#include <limits>
#include <fstream>		// For SeedWord reading /dev/urandom
#include <ctime>		// For SeedWord calling time()
#if !WINDOWS
#include <sys/time.h>		// For SeedWord calling gettimeofday
#include <unistd.h>		// For SeedWord calling getpid(), gethostid()
#include <netinet/in.h>		// For Save/Load
#else
#include <windows.h>		// For SeedWord calling high prec timer
#include <winbase.h>
#include <process.h>		// For SeedWord calling getpid()
#define getpid _getpid
#endif

namespace {
  char rcsid[] = "$Id: Random.cpp 6405 2007-05-21 15:34:24Z ckarney $";
  char* h_rcsid[] = {
    RCSID_RANDOM_H,
    RCSID_RANDOMCANONICAL_H,
    RCSID_POWER2_H,
    RCSID_MT19937_H,
    RCSID_SFMT19937_H,
    RCSID_RANDOMSEED_H,
  };
}

/**
 * \brief Namespace for random number library
 ***********************************************************************/
namespace RandomLib {

  // RandomSeed implementation

  RandomSeed::seed_type RandomSeed::SeedWord() {
    // Check that the assumptions made about the capabilities of the number
    // system are valid.
    STATIC_ASSERT(std::numeric_limits<seed_type>::radix == 2 &&
		  std::numeric_limits<u32>::radix == 2 &&
		  std::numeric_limits<u64>::radix == 2 &&
		  !std::numeric_limits<seed_type>::is_signed &&
		  !std::numeric_limits<u32>::is_signed &&
		  !std::numeric_limits<u64>::is_signed);
    STATIC_ASSERT(std::numeric_limits<seed_type>::digits >= 32 &&
		  std::numeric_limits<u32>::digits >= 32 &&
		  std::numeric_limits<u64>::digits >= 64);
    seed_type t = 0;
    // Linux has /dev/urandom to initialize the seed randomly.  (Use
    // /dev/urandom instead of /dev/random because it does not block.)
    {
      std::ifstream f("/dev/urandom", std::ios::binary | std::ios::in);
      if (f.good()) {
	u32 x;			// Only need 32 bits from /dev/urandom
	f.read(reinterpret_cast<char *>(&x), sizeof(x));
	t = seed_type(x);
      }
    }
    // XOR with PID to prevent correlations between jobs started close together
    // on the same node.  Multiply by a different prime from that used for the
    // time of day to avoid correlations with the time (e.g., if jobs are
    // started at regular intervals so that both the PID and the time are
    // arithmetic progressions).
    t ^= static_cast<seed_type>(getpid()) * seed_type(10000019UL);
#if !WINDOWS
    timeval tv;
    if (gettimeofday(&tv, 0) == 0)
      // XOR with microsec time.  Multiply secs by a prime, in case
      // gettimeofday lamely returns zero for the usecs.
      t ^= static_cast<seed_type>(tv.tv_sec) * seed_type(1000003UL)
	+ static_cast<seed_type>(tv.tv_usec);
#else
    LARGE_INTEGER taux;
    if (QueryPerformanceCounter((LARGE_INTEGER *)&taux))
      t ^= static_cast<seed_type>(taux.HighPart) * seed_type(1000003UL) +
	static_cast<seed_type>(taux.LowPart);
#endif
    else
      // XOR with plain time as a fall back
      t ^= static_cast<seed_type>(std::time(0));

    return SeedMask(t);
  }

  std::vector<RandomSeed::seed_type> RandomSeed::SeedVector() {
    std::vector<seed_type> v;
    {
      // fine-grained timer
#if !WINDOWS
      timeval tv;
      if (gettimeofday(&tv, 0) == 0)
	v.push_back(static_cast<seed_type>(tv.tv_usec));
#else
      LARGE_INTEGER taux;
      if (QueryPerformanceCounter((LARGE_INTEGER *)&taux)) {
	v.push_back(static_cast<seed_type>(taux.LowPart));
	v.push_back(static_cast<seed_type>(taux.HighPart));
      }
#endif
    }
    // seconds
    const time_t tim = std::time(0);
    v.push_back(static_cast<seed_type>(std::time(0)));
    // PID
    v.push_back(static_cast<seed_type>(getpid()));
#if !WINDOWS
    // host ID
    v.push_back(static_cast<seed_type>(gethostid()));
#endif
    {
      // year
      const tm* gt = std::gmtime(&tim);
      v.push_back((seed_type(1900) + static_cast<seed_type>(gt->tm_year)));
    }
    std::transform(v.begin(), v.end(), v.begin(), SeedMask);
    return v;
  }

  std::vector<RandomSeed::seed_type>
  RandomSeed::StringToVector(const std::string& s)
    throw(std::bad_alloc) {
    std::vector<seed_type> v(0);
    const char* c = s.c_str();
    char* q;
    std::string::size_type p = 0;
    while (true) {
      p = s.find_first_of("0123456789", p);
      if (p == std::string::npos)
	break;
      v.push_back(SeedMask(std::strtoul(c + p, &q, 0)));
      p = q - c;
    }
    return v;
  }

  // Save and restore.

  template<>
  void RandomSeed::Write32<RandomSeed::u32, 32u>(std::ostream& os, bool bin,
						 int& cnt, u32 x)
    throw(std::ios::failure) {
    if (bin) {
      unsigned char buf[4];
      // Use network order -- most significant byte first
      buf[3] = static_cast<unsigned char>(x);
      buf[2] = static_cast<unsigned char>(x >>= 8);
      buf[1] = static_cast<unsigned char>(x >>= 8);
      buf[0] = static_cast<unsigned char>(x >>= 8);
      os.write(reinterpret_cast<const char *>(buf), 4);
    } else {
      const int longsperline = 72/11;
      // No spacing before or after
      if (cnt > 0)
	// Newline every longsperline longs
	os << (cnt % longsperline ? " " : "\n");
      os << x;
      ++cnt;
    }
  }

  template<>
  void RandomSeed::Read32<RandomSeed::u32, 32u>(std::istream& is, bool bin,
						u32& x)
    throw(std::ios::failure) {
    if (bin) {
      unsigned char buf[4];
      is.read(reinterpret_cast<char *>(buf), 4);
      // Use network order -- most significant byte first
      x = u32(buf[0]) << 24 | u32(buf[1]) << 16 |
	u32(buf[2]) << 8 | u32(buf[3]);
    } else {
      u32 n;
      is >> std::ws >> n;
      x = n;
    }
    //    x &= U32_MASK;
  }

  template<>
  void RandomSeed::Write32<RandomSeed::u64, 64u>(std::ostream& os, bool bin,
						 int& cnt, u64 x)
    throw(std::ios::failure) {
    Write32<u32, 32u>(os, bin, cnt, static_cast<u32>(x >> 32) & U32_MASK);
    Write32<u32, 32u>(os, bin, cnt, static_cast<u32>(x      ) & U32_MASK);
  }

  template<>
  void RandomSeed::Read32<RandomSeed::u64, 64u>(std::istream& is, bool bin,
						u64& x)
    throw(std::ios::failure) {
    u32 t;
    Read32<u32, 32u>(is, bin, t);
    x = static_cast<u64>(t) << 32;
    Read32<u32, 32u>(is, bin, t);
    x |= static_cast<u64>(t);
  }

  template<>
  void RandomSeed::SeedToState<RandomSeed::u32, 32u>(u32 state[], unsigned n)
    const throw() {
    // This is the algorithm given in the seed_seq class described in
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2079.pdf It is
    // adapted from
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
    const unsigned s = static_cast<unsigned>(_seed.size());
    const unsigned w = 32u;

    u32 r = (u32(5489UL) + s) & U32_MASK;
    state[0] = r;
    for (unsigned k = 1; k < n; ++k) {
      r = u32(1812433253UL) * (r ^ r >> w - 2) + k;
      r &= U32_MASK;
      state[k] = r;
    }
    if (s > 0) {
      unsigned i1 = 0;
      for (unsigned k = (n > s ? n : s), j = 0;
	   k; --k, i1 = i1 == n - 1 ? 0 : i1 + 1, // i1 = i1 + 1 mod n
	     j = j == s - 1 ? 0 : j + 1 ) { // j = j+1 mod s
	r = state[i1] ^ u32(1664525UL) * (r ^ r >> w - 2);
	r += static_cast<u32>(_seed[j]) + j;
	r &= U32_MASK;
	state[i1] = r;
      }
      for (unsigned k = n; k; --k,
	     i1 = i1 == n - 1 ? 0 : i1 + 1) { // i1 = i1 + 1 mod n
	r = state[i1] ^ u32(1566083941UL) * (r ^ r >> w - 2);
	r -= i1;
	r &= U32_MASK;
	state[i1] = r;
      }
    }
  }

  template<>
  void RandomSeed::SeedToState<RandomSeed::u64, 64u>(u64 state[], unsigned n)
    const throw() {
    // This is adapted from
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
    // with some changes to bring it into line with the seed_seq class given in
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2006/n2079.pdf The
    // differences from the 32-bit version are: different (bigger)
    // multiplicative constants; packing the seed array into 64-bit chunks.

    const unsigned s = static_cast<unsigned>(_seed.size());
    const unsigned w = 64u;

    u64 r = (u64(5489ULL) + s) & U64_MASK;
    state[0] = r;
    for (unsigned k = 1; k < n; ++k) {
      r = u64(6364136223846793005ULL) * (r ^ r >> w - 2) + k;
      r &= U64_MASK;
      state[k] = r;
    }
    if (s > 0) {
      const unsigned key2 = (s + 1)/2;
      unsigned i1 = 0;
      for (unsigned k = (n > key2 ? n : key2), j = 0;
	   k; --k, i1 = i1 == n - 1 ? 0 : i1 + 1, // i1 = i1 + 1 mod n
	     j = j == key2 - 1 ? 0 : j + 1 ) { // j = j + 1 mod key2
	r = state[i1] ^ u64(3935559000370003845ULL) * (r ^ r >> w - 2);
	r += j + static_cast<u64>(_seed[2 * j]) +
	  (2 * j + 1 == s ? u64(0) :
	   static_cast<u64>(_seed[2 * j + 1]) << 32);
	r &= U64_MASK;
	state[i1] = r;
      }
      for (unsigned k = n; k; --k,
	     i1 = i1 == n - 1 ? 0 : i1 + 1) { // i1 = i1 + 1 mod n
	r = state[i1] ^ u64(2862933555777941757ULL) * (r ^ r >> w - 2);
	r -= i1;
	r &= U64_MASK;
	state[i1] = r;
      }
    }
  }

  // Note that BitWidth = 128u is a bogus setting.  This is really a 32-bit
  // SeedToState algorithm.  128u is used since the algorithm was introduced
  // with the 128-bit SFMT19937 and to distinguish it from the 32-bit algorithm
  // used by MT19937
  template<>
  void RandomSeed::SeedToState<RandomSeed::u32, 128u>(u32 state[], unsigned n)
    const throw() {
    // This is adapted from the routine init_by_array by Mutsuo Saito given in
    // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/SFMT-src-1.2.tar.gz

    if (n == 0)
      return;			// Nothing to do

    const unsigned s = static_cast<unsigned>(_seed.size()),
      // Add treatment of small n with lag = (n - 1)/2 for n <= 7.  In
      // particular, the first operation (xor or plus) in each for loop
      // involves three distinct indices for n > 2.
      lag = n >= 623 ? 11 : (n >= 68 ? 7 : (n >= 39 ? 5 :
					    (n >= 7 ? 3 : (n - 1)/2))),
      // count = max( s + 1, n )
      count = s + 1 > n ? s + 1 : n;

    std::fill(state, state + n, u32(0x8b8b8b8bUL));
    const unsigned w = 32u;

    unsigned i = 0, k = (n - lag) / 2, l = k + lag;
    u32 r = state[n - 1];
    for (unsigned j = 0; j < count; ++j,
	   i = i == n - 1 ? 0 : i + 1,
	   k = k == n - 1 ? 0 : k + 1,
	   l = l == n - 1 ? 0 : l + 1) {
      // Here r = state[(j - 1) mod n]
      //      i = j mod n
      //      k = (j + (n - lag)/2) mod n
      //      l = (j + (n - lag)/2 + lag) mod n
      r ^= state[i] ^ state[k];
      r &= U32_MASK;
      r = u32(1664525UL) * (r ^ r >> w - 5);
      state[k] += r;
      r += i + (j > s ? 0 :
		(j ? static_cast<u32>(_seed[j - 1]) : s));
      state[l] += r;
      state[i] = r;
    }

    for (unsigned j = n; j; --j,
	   i = i == n - 1 ? 0 : i + 1,
	   k = k == n - 1 ? 0 : k + 1,
	   l = l == n - 1 ? 0 : l + 1) {
      // Here r = state[(i - 1) mod n]
      //      k = (i + (n - lag)/2) mod n
      //      l = (i + (n - lag)/2 + lag) mod n
      r += state[i] + state[k];
      r &= U32_MASK;
      r = u32(1566083941UL) * (r ^ r >> w - 5);
      r &= U32_MASK;
      state[k] ^= r;
      r -= i;
      r &= U32_MASK;
      state[l] ^= r;
      state[i] = r;
    }
  }

  // MT19937 implementation

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::Init() throw() {
    // On exit we have _ptr == N.

    // Convert the seed into state
    SeedToState<UIntType, BitWidth>(_state, N);

    // Perform the MT-specific sanity check on the resulting state ensuring
    // that the significant 19937 bits are not all zero.
    _state[0] &= UPPER_MASK;	// Mask out unused bits
    unsigned i = 0;
    while (i < N && _state[i] == 0)
      ++i;
    if (i >= N)
      _state[0] = result_type(1) << width - 1; // with prob 2^-19937

    // This sets the low R bits of _state[0] consistent with the rest of the
    // state.  Needed to handle SetCount(-N); Ran32(); immediately following
    // reseeding.  This wasn't required in the original code because a
    // NextBatch was always done first.
    result_type q = _state[N - 1] ^ _state[M - 1], s = q >> width - 1;
    q = (q ^ (s ? MATRIX_A : result_type(0))) << 1 | s;
    _state[0] = _state[0] & UPPER_MASK | q & LOWER_MASK;

    _rounds = -1;
    _ptr = N;
  }

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::StepCount(long long n) throw() {
    // On exit we have 0 <= _ptr <= N.
    if (_ptr == UNINIT)
      Init();
    const long long ncount = n + Count(); // new Count()
    long long nrounds = ncount / N;
    int nptr = static_cast<int>(ncount - nrounds * N);
    if (nptr <= 0) {
      // We pick _ptr = N or _ptr = 0 depending on which choice involves the
      // least work.  We thus avoid doing one (potentially unneeded) called to
      // NextBatch/PreviousBatch.
      --nrounds;
      nptr += N;
    }
    if (nrounds > _rounds)
      NextBatch(nrounds - _rounds);
    else if (nrounds < _rounds) {
      if (nptr == N) {
	nptr = 0;
	++nrounds;
      }
      PreviousBatch(_rounds - nrounds);
    }
    _ptr = nptr;
  }

  template<typename UIntType, unsigned BitWidth>
  RandomSeed::u32 MT19937<UIntType, BitWidth>::Check(u32 version) const
    throw(std::out_of_range) {
    if (version != VERSIONID && version != VERSIONID - 1)
      throw std::out_of_range("MT19937: Unknown version");
    u32 check = version;
    CheckSum<u32, 32u>(_seed.size(), check);
    for (std::vector<seed_type>::const_iterator n = _seed.begin();
	 n != _seed.end(); ++n) {
      if (*n != SeedMask(*n))
	throw std::out_of_range("Illegal seed value");
      CheckSum<u32, 32u>(*n, check);
    }
    if (version == VERSIONID - 1 && _ptr == UNINIT)
      CheckSum<u32, 32u>(N + 1, check);
    else
      CheckSum<u32, 32u>(_ptr, check);
    if (_stride == 0 || _stride > UNINIT/2)
      throw std::out_of_range("MT19937: Invalid stride");
    if (version != VERSIONID - 1)
      CheckSum<u32, 32u>(_stride, check);
    if (_ptr != UNINIT) {
      if (_ptr >= N + _stride)
	throw std::out_of_range("MT19937: Invalid pointer");
      CheckSum<u64, 64u>(static_cast<u64>(_rounds), check);
      result_type x = 0;
      for (unsigned i = 0; i < N; ++i) {
	CheckSum<result_type, width>(_state[i], check);
	x |= _state[i];
      }
      if (x == 0)
	throw std::out_of_range("MT19937: All-zero state");

      // There are only width*(N-1) + 1 = 19937 independent bits of state.
      // Thus the low width-1 bits of _state[0] are derivable from the other
      // bits in _state.  Verify that the redundant bits bits are consistent.
      result_type q = _state[N - 1] ^ _state[M - 1], s = q >> width - 1;
      q = (q ^ (s ? MATRIX_A : result_type(0))) << 1 | s;
      if ((q ^ _state[0]) & LOWER_MASK)
	throw std::out_of_range("MT19937: Invalid state");
    }
    return check;
  }

  // Save and restore.

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::Save(std::ostream& os, bool bin) const
    throw(std::ios::failure) {
    u32 check = Check(VERSIONID);
    int c = 0;
    Write32<u32, 32u>(os, bin, c, VERSIONID);
    Write32<u32, 32u>(os, bin, c, _seed.size());
    for (std::vector<seed_type>::const_iterator n = _seed.begin();
	 n != _seed.end(); ++n)
      Write32<u32, 32u>(os, bin, c, *n);
    Write32<u32, 32u>(os, bin, c, _ptr);
    Write32<u32, 32u>(os, bin, c, _stride);
    if (_ptr != UNINIT) {
      Write32<u64, 64u>(os, bin, c, static_cast<u64>(_rounds));
      for (unsigned i = 0; i < N; ++i)
	Write32<result_type, width>(os, bin, c, _state[i]);
    }
    Write32<u32, 32u>(os, bin, c, check);
  }

  template<typename UIntType, unsigned BitWidth>
  MT19937<UIntType, BitWidth>::MT19937(std::istream& is, bool bin)
    throw(std::ios::failure, std::out_of_range, std::bad_alloc) {
    u32 version, t;
    Read32<u32, 32u>(is, bin, version);
    Read32<u32, 32u>(is, bin, t);
    _seed.resize(static_cast<size_t>(t));
    for (std::vector<seed_type>::iterator n = _seed.begin();
	 n != _seed.end(); ++n) {
      Read32<u32, 32u>(is, bin, t);
      *n = static_cast<seed_type>(t);
    }
    Read32<u32, 32u>(is, bin, t);
    // Don't need to worry about sign extension because _ptr is unsigned.
    _ptr = static_cast<unsigned>(t);
    if (version == VERSIONID - 1) {
      if (_ptr == N + 1)
	_ptr = UNINIT;
      _stride = 1;
    } else {
      Read32<u32, 32u>(is, bin, t);
      _stride = static_cast<unsigned>(t);
    }
    if (_ptr != UNINIT) {
      u64 p;
      Read32<u64, 64u>(is, bin, p);
      _rounds = static_cast<long long>(p);
      // Sign extension in case long long is bigger than 64 bits.
      _rounds <<= 63 - std::numeric_limits<long long>::digits;
      _rounds >>= 63 - std::numeric_limits<long long>::digits;
      for (unsigned i = 0; i < N; ++i)
	Read32<result_type, width>(is, bin, _state[i]);
    }
    Read32<u32, 32u>(is, bin, t);
    if (t != Check(version))
      throw std::out_of_range("MT19937: Checksum failure");
  }

  // Code for NextBatch is from
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
  // and http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html.  Here
  // is the copyright notice.

  // A C-program for MT19937, with initialization improved 2002/1/26.
  // Coded by Takuji Nishimura and Makoto Matsumoto.
  //
  // Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
  // All rights reserved.
  //
  // Redistribution and use in source and binary forms, with or without
  // modification, are permitted provided that the following conditions
  // are met:
  //
  //   1. Redistributions of source code must retain the above copyright
  //	notice, this list of conditions and the following disclaimer.
  //
  //   2. Redistributions in binary form must reproduce the above copyright
  //	notice, this list of conditions and the following disclaimer in the
  //	documentation and/or other materials provided with the distribution.
  //
  //   3. The names of its contributors may not be used to endorse or promote
  //	products derived from this software without specific prior written
  //	permission.
  //
  // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  // A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
  // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  //
  // Any feedback is very welcome.
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  // email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

  // Let y[i] = _state[i] & UPPER_MASK | _state[i + 1] & LOWER_MASK, then
  // _state[N] depends on y[0] and _state[M}.
  //
  // The "state" does NOT depends on the low bits of _state[0].  However, after
  // calls to Init(), NextBatch() and PreviousBatch(), _state[0] needs to
  // contain the right whole word result since it's used (after tempering) as
  // one of the outputs of Ran32().
  //
  // In MT19937_STEP, input is I, J = I + 1, K = I + M; output is I = I + N
  // (mod N)

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::Next() throw() {
    // On exit we have 0 <= _ptr < N.
    if (_ptr == UNINIT)		// Init() has not been called
      Init();
    NextBatch(_ptr/N);
    _ptr %= N;
  }

#define MT19937_STEP(I, J, K) _state[I] = _state[K] ^		\
    (_state[J] & result_type(1) ? MATRIX_A : result_type(0)) ^	\
    (_state[I] & UPPER_MASK | _state[J] & LOWER_MASK) >> 1

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::NextBatch(unsigned long long count)
    throw() {
    _rounds += count;
    for (; count; --count) {
      // This ONLY uses high bit of _state[0]
      unsigned i = 0;
      for (; i < N - M; ++i) MT19937_STEP(i, i + 1, i + M    );
      for (; i < N - 1; ++i) MT19937_STEP(i, i + 1, i + M - N);
      MT19937_STEP(N - 1, 0, M - 1); // i = N - 1
    }
  }

#undef MT19937_STEP

  // This undoes MT19937_STEP.  This piece of coding is from revrand() by
  // Katsumi Hagita.  See
  //
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/REVmt19937b.f
  //
  // The code is cleaned up a little from the Fortran version by getting rid of
  // the unnecessary masking by YMASK and by using simpler logic to restore the
  // correct value of _state[0].
  //
  // y[0] depends on _state[N] and _state[M].
  //
  // Here input is J = I + N - 1, K = I + M - 1, and p = y[I] (only the high
  // bits are used); output _state[I] and p = y[I - 1].

#define MT19937_REVSTEP(I, J, K) {				\
    result_type q = _state[J] ^ _state[K], s = q >> width - 1;	\
    q = (q ^ (s ? MATRIX_A : result_type(0))) << 1 | s;		\
    _state[I] = p & UPPER_MASK | q & LOWER_MASK;		\
    p = q;							\
  }

  template<typename UIntType, unsigned BitWidth>
  void MT19937<UIntType, BitWidth>::PreviousBatch(unsigned long long count)
    throw() {
    _rounds -= count;
    for (; count; --count) {
      // This ONLY uses high bit of _state[0]
      result_type p = _state[0];
      // Fix low bits of _state[0] and compute y[-1]
      MT19937_REVSTEP(0, N - 1, M - 1); // i = N
      unsigned i = N - 1;
      for (; i > N - M; --i) MT19937_REVSTEP(i, i - 1, i + M - 1 - N);
      for (; i        ; --i) MT19937_REVSTEP(i, i - 1, i + M - 1    );
      MT19937_REVSTEP(0, N - 1, M - 1); // i = 0
    }
  }

#undef MT19937_REVSTEP

  // The two standard instantiations
  template class MT19937<RandomSeed::u32, 32u>; // MT19937_32 = mt19937ar
  template class MT19937<RandomSeed::u64, 64u>; // MT19937_64 = mt19937-64

  // RandomCanonical instantiations.

  template<> MRandom32 MRandom32::Global(std::vector<unsigned long>(0));
  template<> MRandom64 MRandom64::Global(std::vector<unsigned long>(0));

  // SFMT19937 implementation

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::Init() throw() {
    // On exit we have _ptr == N.

    // Make sure everything packs together properly.
#if HAVE_SSE2
    STATIC_ASSERT(sizeof(_state.w) == sizeof(_state.v));
#endif

    // Convert the seed into state
    SeedToState<u32, 128u>(_state.u, N128 * 4);

    // Carry out the Period Certification for SFMT19937
    u32 inner = _state.u[0] & PARITY0 ^ _state.u[1] & PARITY1 ^
      _state.u[2] & PARITY2 ^ _state.u[3] & PARITY3;
    for (unsigned s = 16; s; s >>= 1)
      inner ^= inner >> s;
    STATIC_ASSERT(PARITY_LSB < 32 && PARITY0 & 1u << PARITY_LSB);
    // Now inner & 1 is the parity of the number of 1 bits in w_0 & p.
    if ((inner & 1u) == 0)
      // Change bit of w_0 corresponding to LSB of PARITY
      _state.u[PARITY_LSB >> 5] ^= 1u << (PARITY_LSB & 31u);

    // Make sure that writing _state.w[0] will overwrite at most _state.u[0]
    // and _state.u[1].
    STATIC_ASSERT(BitWidth <= 32u ||
		  std::numeric_limits<u64>::digits <=
		  2 * std::numeric_limits<u32>::digits);
    if (BitWidth > 32u)
      // Pack into 64-bit words with LSB ordering.  The loop needs to be in
      // this order in case u32 is a 64-bit quantity.
      for (unsigned j = 0; j < N; ++j)
	_state.w[j] = result_type(u64(_state.u[2*j]) |
				  u64(_state.u[2*j+1]) << 32);

    _rounds = -1;
    _ptr = N;
  }

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::StepCount(long long n) throw() {
    // On exit we have 0 <= _ptr <= N.
    if (_ptr == UNINIT)
      Init();
    const long long ncount = n + Count(); // new Count()
    long long nrounds = ncount / N;
    int nptr = static_cast<int>(ncount - nrounds * N);
    if (nptr <= 0) {
      // We pick _ptr = N or _ptr = 0 depending on which choice involves the
      // least work.  We thus avoid doing one (potentially unneeded) called to
      // NextBatch/PreviousBatch.
      --nrounds;
      nptr += N;
    }
    if (nrounds > _rounds)
      NextBatch(nrounds - _rounds);
    else if (nrounds < _rounds) {
      if (nptr == N) {
	nptr = 0;
	++nrounds;
      }
      PreviousBatch(_rounds - nrounds);
    }
    _ptr = nptr;
  }

  template<typename UIntType, unsigned BitWidth>
  RandomSeed::u32 SFMT19937<UIntType, BitWidth>::Check(u32 version) const
    throw(std::out_of_range) {
    if (version != VERSIONID)
      throw std::out_of_range("SFMT19937: Unknown version");
    u32 check = version;
    CheckSum<u32, 32u>(_seed.size(), check);
    for (std::vector<seed_type>::const_iterator n = _seed.begin();
	 n != _seed.end(); ++n) {
      if (*n != SeedMask(*n))
	throw std::out_of_range("Illegal seed value");
      CheckSum<u32, 32u>(*n, check);
    }
    if (version == VERSIONID - 1 && _ptr == UNINIT)
      CheckSum<u32, 32u>(N + 1, check);
    else
      CheckSum<u32, 32u>(_ptr, check);
    if (_stride == 0 || _stride > UNINIT/2)
      throw std::out_of_range("SFMT19937: Invalid stride");
    if (version != VERSIONID - 1)
      CheckSum<u32, 32u>(_stride, check);
    if (_ptr != UNINIT) {
      if (_ptr >= N + _stride)
	throw std::out_of_range("SFMT19937: Invalid pointer");
      CheckSum<u64, 64u>(static_cast<u64>(_rounds), check);
      result_type x = 0;
      for (unsigned i = 0; i < N; ++i) {
	CheckSum<result_type, width>(_state.w[i], check);
	x |= _state.w[i];
      }
      if (x == 0)
	throw std::out_of_range("SFMT19937: All-zero state");
    }
    return check;
  }

  // Save and restore.

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::Save(std::ostream& os, bool bin) const
    throw(std::ios::failure) {
    u32 check = Check(VERSIONID);
    int c = 0;
    Write32<u32, 32u>(os, bin, c, VERSIONID);
    Write32<u32, 32u>(os, bin, c, _seed.size());
    for (std::vector<seed_type>::const_iterator n = _seed.begin();
	 n != _seed.end(); ++n)
      Write32<u32, 32u>(os, bin, c, *n);
    Write32<u32, 32u>(os, bin, c, _ptr);
    Write32<u32, 32u>(os, bin, c, _stride);
    if (_ptr != UNINIT) {
      Write32<u64, 64u>(os, bin, c, static_cast<u64>(_rounds));
      for (unsigned i = 0; i < N; ++i)
	Write32<result_type, width>(os, bin, c, _state.w[i]);
    }
    Write32<u32, 32u>(os, bin, c, check);
  }

  template<typename UIntType, unsigned BitWidth>
  SFMT19937<UIntType, BitWidth>::SFMT19937(std::istream& is, bool bin)
    throw(std::ios::failure, std::out_of_range, std::bad_alloc) {
    u32 version, t;
    Read32<u32, 32u>(is, bin, version);
    Read32<u32, 32u>(is, bin, t);
    _seed.resize(static_cast<size_t>(t));
    for (std::vector<seed_type>::iterator n = _seed.begin();
	 n != _seed.end(); ++n) {
      Read32<u32, 32u>(is, bin, t);
      *n = static_cast<seed_type>(t);
    }
    Read32<u32, 32u>(is, bin, t);
    // Don't need to worry about sign extension because _ptr is unsigned.
    _ptr = static_cast<unsigned>(t);
    if (version == VERSIONID - 1) {
      if (_ptr == N + 1)
	_ptr = UNINIT;
      _stride = 1;
    } else {
      Read32<u32, 32u>(is, bin, t);
      _stride = static_cast<unsigned>(t);
    }
    if (_ptr != UNINIT) {
      u64 p;
      Read32<u64, 64u>(is, bin, p);
      _rounds = static_cast<long long>(p);
      // Sign extension in case long long is bigger than 64 bits.
      _rounds <<= 63 - std::numeric_limits<long long>::digits;
      _rounds >>= 63 - std::numeric_limits<long long>::digits;
      for (unsigned i = 0; i < N; ++i)
	Read32<result_type, width>(is, bin, _state.w[i]);
    }
    Read32<u32, 32u>(is, bin, t);
    if (t != Check(version))
      throw std::out_of_range("SFMT19937: Checksum failure");
  }

  // Algorithm for NextBatch is from Saito's Master's Thesis
  // http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/M062821.pdf
  //
  // This implements
  //
  //     w_{i+N} = w_i A xor w_M B xor w_{i+N-2} C xor w_{i+N-1} D
  //
  // where w_i is a 128-bit word and
  //
  //     w A = (w << 8)_128 xor w
  //     w B = (w >> 11)_32 & MSK
  //     w C = (w >> 8)_128
  //     w D = (w << 18)_32
  //
  // Here the _128 means shift is of whole 128-bit word.  _32 means the shifts
  // are independently done on each 32-bit word.
  //
  // In SFMT19937_STEP32 and SFMT19937_STEP64 input is I, J = I + M and output
  // is I = I + N (mod N).  On input, s and r give state for I + N - 2 and I +
  // N - 1; on output s and r give state for I + N - 1 and I + N.  The
  // implementation of 128-bit operations is open-coded in a portable fashion
  // (with LSB ordering).
  //
  // N.B. Here N and M are the lags in units of BitWidth words and so are 4
  // (for u32 implementation) or 2 (for u64 implementation) times bigger than
  // defined in Saito's thesis.

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::Next() throw() {
    // On exit we have 0 <= _ptr < N.
    if (_ptr == UNINIT)		// Init() has not been called
      Init();
    NextBatch(_ptr/N);
    _ptr %= N;
  }

#if !HAVE_SSE2		       // The standard portable implementation

#define SFMT19937_STEP32(I, J) {				\
    result_type t = _state.w[I] ^ _state.w[I] << 8 ^		\
      _state.w[J] >> 11 & MSK0 ^				\
      (s0 >> 8 | s1 << 24) ^ r0 << 18;				\
    s0 = r0; r0 = t & RESULT_MASK;				\
    t = _state.w[I + 1] ^					\
      (_state.w[I + 1] << 8 | _state.w[I] >> 24) ^		\
      _state.w[J + 1] >> 11 & MSK1 ^				\
      (s1 >> 8 | s2 << 24) ^ r1 << 18;				\
    s1 = r1; r1 = t & RESULT_MASK;				\
    t = _state.w[I + 2] ^					\
      (_state.w[I + 2] << 8 | _state.w[I + 1] >> 24) ^		\
      _state.w[J + 2] >> 11 & MSK2 ^				\
      (s2 >> 8 | s3 << 24) ^ r2 << 18;				\
    s2 = r2; r2 = t & RESULT_MASK;				\
    t = _state.w[I + 3] ^					\
      (_state.w[I + 3] << 8 | _state.w[I + 2] >> 24) ^		\
      _state.w[J + 3] >> 11 & MSK3 ^ s3 >> 8 ^ r3 << 18;	\
    s3 = r3; r3 = t & RESULT_MASK;				\
    _state.w[I    ] = r0; _state.w[I + 1] = r1;			\
    _state.w[I + 2] = r2; _state.w[I + 3] = r3;			\
  }

  template<>
  void SFMT19937<RandomSeed::u32, 32u>::NextBatch(unsigned long long count)
    throw() {
    _rounds += count;
    // x[i+N] = g(x[i], x[i+M], x[i+N-2], x[i,N-1])
    for (; count; --count) {
      result_type
	s0 = _state.w[N - 8], s1 = _state.w[N - 7],
	s2 = _state.w[N - 6], s3 = _state.w[N - 5],
	r0 = _state.w[N - 4], r1 = _state.w[N - 3],
	r2 = _state.w[N - 2], r3 = _state.w[N - 1];
      unsigned i = 0;
      for (; i + M < N; i += R) SFMT19937_STEP32(i, i + M    );
      for (; i < N    ; i += R) SFMT19937_STEP32(i, i + M - N);
    }
  }

#undef SFMT19937_STEP32

#define SFMT19937_STEP64(I, J) {			\
    result_type t = _state.w[I] ^ _state.w[I] << 8 ^	\
      _state.w[J] >> 11 & LMSK0 ^			\
      (s0 >> 8 | s1 << 56) ^ r0 << 18 & LMSK18;		\
    s0 = r0; r0 = t & RESULT_MASK;			\
    t = _state.w[I + 1] ^				\
      (_state.w[I + 1] << 8 | _state.w[I] >> 56) ^	\
      _state.w[J + 1] >> 11 & LMSK1 ^			\
      s1 >> 8 ^ r1 << 18 & LMSK18;			\
    s1 = r1; r1 = t & RESULT_MASK;			\
    _state.w[I] = r0; _state.w[I + 1] = r1;		\
  }

  template<>
  void SFMT19937<RandomSeed::u64, 64u>::NextBatch(unsigned long long count)
    throw() {
    _rounds += count;
    // x[i+N] = g(x[i], x[i+M], x[i+N-2], x[i,N-1])
    for (; count; --count) {
      result_type
	s0 = _state.w[N - 4], s1 = _state.w[N - 3],
	r0 = _state.w[N - 2], r1 = _state.w[N - 1];
      unsigned i = 0;
      for (; i + M < N; i += R) SFMT19937_STEP64(i, i + M    );
      for (; i < N    ; i += R) SFMT19937_STEP64(i, i + M - N);
    }
  }

#undef SFMT19937_STEP64

  // This undoes SFMT19937_STEP.  Trivially, we have
  //
  //     w_i A = w_{i+N} xor w_{i+M} B xor w_{i+N-2} C xor w_{i+N-1} D
  //
  // Given w_i A we can determine w_i from the observation that A^16 =
  // identity, thus
  //
  //     w_i = (w_i A) A^15
  //
  // Because x A^(2^n) = x << (8*2^n) xor x, the operation y = x A^15 can be
  // implemented as
  //
  //     y'   = (x    << 64)_128 xor x    = x    A^8
  //     y''  = (y'   << 32)_128 xor y'   = y'   A^4 = x A^12
  //     y''' = (y''  << 16)_128 xor y''  = y''  A^2 = x A^14
  //     y    = (y''' <<  8)_128 xor y''' = y''' A   = x A^15
  //
  // Here input is I = I + N, J = I + M, K = I + N - 2, L = I + N -1, and
  // output is I = I.
  //
  // This is about 15-35% times slower than SFMT19937_STEPNN, because (1) there
  // doesn't appear to be a straightforward way of saving intermediate results
  // across calls as with SFMT19937_STEPNN and (2) w A^15 is slower to compute
  // than w A.

#define SFMT19937_REVSTEP32(I, J, K, L) {			\
    result_type							\
      t0 = (_state.w[I] ^ _state.w[J] >> 11 & MSK0 ^		\
	    (_state.w[K] >> 8 | _state.w[K + 1] << 24) ^	\
	    _state.w[L] << 18) & RESULT_MASK,			\
      t1 = (_state.w[I + 1] ^					\
	    _state.w[J + 1] >> 11 & MSK1 ^			\
	    (_state.w[K + 1] >> 8 | _state.w[K + 2] << 24) ^	\
	    _state.w[L + 1] << 18) & RESULT_MASK,		\
      t2 = (_state.w[I + 2] ^					\
	    _state.w[J + 2] >> 11 & MSK2 ^			\
	    (_state.w[K + 2] >> 8 | _state.w[K + 3] << 24) ^	\
	    _state.w[L + 2] << 18) & RESULT_MASK,		\
      t3 = (_state.w[I + 3] ^					\
	    _state.w[J + 3] >> 11 & MSK3 ^			\
	    _state.w[K + 3] >> 8 ^				\
	    _state.w[L + 3] << 18) & RESULT_MASK;		\
    t3 ^= t1; t2 ^= t0; t3 ^= t2; t2 ^= t1; t1 ^= t0;		\
    t3 ^= t2 >> 16 | t3 << 16 & RESULT_MASK;			\
    t2 ^= t1 >> 16 | t2 << 16 & RESULT_MASK;			\
    t1 ^= t0 >> 16 | t1 << 16 & RESULT_MASK;			\
    t0 ^=            t0 << 16 & RESULT_MASK;			\
    _state.w[I    ] = t0 ^             t0 << 8 & RESULT_MASK;	\
    _state.w[I + 1] = t1 ^ (t0 >> 24 | t1 << 8 & RESULT_MASK);	\
    _state.w[I + 2] = t2 ^ (t1 >> 24 | t2 << 8 & RESULT_MASK);	\
    _state.w[I + 3] = t3 ^ (t2 >> 24 | t3 << 8 & RESULT_MASK);	\
 }

  template<>
  void SFMT19937<RandomSeed::u32, 32u>::PreviousBatch(unsigned long long count)
    throw() {
    _rounds -= count;
    for (; count; --count) {
      unsigned i = N;
      for (; i + M > N;) {
	i -= R; SFMT19937_REVSTEP32(i, i + M - N, i - 2 * R, i - R);
      }
      for (; i > 2 * R;) {
	i -= R; SFMT19937_REVSTEP32(i, i + M    , i - 2 * R, i - R);
      }
      SFMT19937_REVSTEP32(R, M + R, N -     R, 0    ); // i = R
      SFMT19937_REVSTEP32(0, M    , N - 2 * R, N - R); // i = 0
    }
  }

#undef SFMT19937_REVSTEP32

  // In combining the left and right shifts to simulate a 128-bit shift we
  // usually use or.  However we can equivalently use xor (e.g., t1 << 8 ^ t0
  // >> 56 instead of t1 ^ t1 << 8 | t0 >> 56) and this speeds up the code if
  // used in some places.
#define SFMT19937_REVSTEP64(I, J, K, L) {			\
    result_type							\
      t0 = _state.w[I] ^ _state.w[J] >> 11 & LMSK0 ^		\
      (_state.w[K] >> 8 | _state.w[K + 1] << 56 & RESULT_MASK)	\
      ^ _state.w[L] << 18 & LMSK18,				\
      t1 = _state.w[I + 1] ^ _state.w[J + 1] >> 11 & LMSK1 ^	\
      _state.w[K + 1] >> 8 ^ _state.w[L + 1] << 18 & LMSK18;	\
    t1 ^= t0;							\
    t1 ^= t0 >> 32 ^ t1 << 32 & RESULT_MASK;			\
    t0 ^=            t0 << 32 & RESULT_MASK;			\
    t1 ^= t0 >> 48 ^ t1 << 16 & RESULT_MASK;			\
    t0 ^=            t0 << 16 & RESULT_MASK;			\
    _state.w[I    ] = t0 ^            t0 << 8 & RESULT_MASK;	\
    _state.w[I + 1] = t1 ^ t0 >> 56 ^ t1 << 8 & RESULT_MASK;	\
 }

  template<>
  void SFMT19937<RandomSeed::u64, 64u>::PreviousBatch(unsigned long long count)
    throw() {
    _rounds -= count;
    for (; count; --count) {
      unsigned i = N;
      for (; i + M > N;) {
	i -= R; SFMT19937_REVSTEP64(i, i + M - N, i - 2 * R, i - R);
      }
      for (; i > 2 * R;) {
	i -= R; SFMT19937_REVSTEP64(i, i + M    , i - 2 * R, i - R);
      }
      SFMT19937_REVSTEP64(R, M + R, N -     R, 0    ); // i = R
      SFMT19937_REVSTEP64(0, M    , N - 2 * R, N - R); // i = 0
    }
  }

#undef SFMT19937_REVSTEP64

#else  // The SSE2 implementation

  // This is adapted from SFMT-sse.c in the SFMT 1.2 distribution.
  // The order of instructions has been rearranged to increase the
  // speed slightly
#define SFMT19937_STEP128(I, J) {		\
    __m128i x = _mm_load_si128(_state.v + I),	\
      y = _mm_srli_epi32(_state.v[J], 11),	\
      z = _mm_srli_si128(s, 1);			\
    s = _mm_slli_epi32(r, 18);			\
    z = _mm_xor_si128(z, x);			\
    x = _mm_slli_si128(x, 1);			\
    z = _mm_xor_si128(z, s);			\
    y = _mm_and_si128(y, m);			\
    z = _mm_xor_si128(z, x);			\
    s = r;					\
    r = _mm_xor_si128(z, y);			\
    _mm_store_si128(_state.v + I, r);		\
  }

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::NextBatch(unsigned long long count)
    throw() {
    _rounds += count;
    for (; count; --count) {
      const __m128i m = _mm_set_epi32(MSK3, MSK2, MSK1, MSK0);
      __m128i s = _mm_load_si128(_state.v + N128 - 2),
	r = _mm_load_si128(_state.v + N128 - 1);
      unsigned i = 0;
      for (; i + M128 < N128; ++i) SFMT19937_STEP128(i, i + M128       );
      for (; i < N128       ; ++i) SFMT19937_STEP128(i, i + M128 - N128);
    }
  }

#undef SFMT19937_STEP128

#define SFMT19937_REVSTEP128(I, J, K, L) {			\
    __m128i x = _mm_load_si128(_state.v + I),			\
      y = _mm_srli_epi32(_state.v[J], 11),			\
      z = _mm_slli_epi32(_state.v[L], 18);			\
    y = _mm_and_si128(y, m);					\
    x = _mm_xor_si128(x, _mm_srli_si128(_state.v[K], 1));	\
    x = _mm_xor_si128(x, z);					\
    x = _mm_xor_si128(x, y);					\
    x = _mm_xor_si128(_mm_slli_si128(x, 8), x);			\
    x = _mm_xor_si128(_mm_slli_si128(x, 4), x);			\
    x = _mm_xor_si128(_mm_slli_si128(x, 2), x);			\
    x = _mm_xor_si128(_mm_slli_si128(x, 1), x);			\
    _mm_store_si128(_state.v + I, x);				\
  }

  template<typename UIntType, unsigned BitWidth>
  void SFMT19937<UIntType, BitWidth>::PreviousBatch(unsigned long long count)
    throw() {
    _rounds -= count;
    for (; count; --count) {
      const __m128i m = _mm_set_epi32(MSK3, MSK2, MSK1, MSK0);
      unsigned i = N128;
      for (; i + M128 > N128;) {
	--i; SFMT19937_REVSTEP128(i, i + M128 - N128, i - 2, i - 1);
      }
      for (; i > 2;) {
	--i; SFMT19937_REVSTEP128(i, i + M128, i - 2, i - 1);
      }
      SFMT19937_REVSTEP128(1, M128 + 1, N128 - 1, 0       ); // i = 1
      SFMT19937_REVSTEP128(0, M128    , N128 - 2, N128 - 1); // i = 0
    }
  }

#undef SFMT19937_REVSTEP128

#endif	// HAVE_SSE2

  // The two standard instantiations
  template class SFMT19937<RandomSeed::u32, 32u>;
  template class SFMT19937<RandomSeed::u64, 64u>;

  // RandomCanonical instantiations.

  template<> SRandom32 SRandom32::Global(std::vector<unsigned long>(0));
  template<> SRandom64 SRandom64::Global(std::vector<unsigned long>(0));

  // Power2 implementation

#if RANDOM_POWERTABLE
  // Powers of two.  Just use floats here.  As long as there's no overflow
  // or underflow these are exact.  In particular they can be cast to
  // doubles or long doubles with no error.
  const float Power2::power2[maxpow - minpow + 1] = {
#if RANDOM_LONGDOUBLEPREC > 64
    // It would be nice to be able to use the C99 notation of 0x1.0p-120
    // for 2^-120 here.
    1/1329227995784915872903807060280344576.f, // 2^-120
    1/664613997892457936451903530140172288.f, // 2^-119
    1/332306998946228968225951765070086144.f, // 2^-118
    1/166153499473114484112975882535043072.f, // 2^-117
    1/83076749736557242056487941267521536.f, // 2^-116
    1/41538374868278621028243970633760768.f, // 2^-115
    1/20769187434139310514121985316880384.f, // 2^-114
    1/10384593717069655257060992658440192.f, // 2^-113
    1/5192296858534827628530496329220096.f, // 2^-112
    1/2596148429267413814265248164610048.f, // 2^-111
    1/1298074214633706907132624082305024.f, // 2^-110
    1/649037107316853453566312041152512.f, // 2^-109
    1/324518553658426726783156020576256.f, // 2^-108
    1/162259276829213363391578010288128.f, // 2^-107
    1/81129638414606681695789005144064.f, // 2^-106
    1/40564819207303340847894502572032.f, // 2^-105
    1/20282409603651670423947251286016.f, // 2^-104
    1/10141204801825835211973625643008.f, // 2^-103
    1/5070602400912917605986812821504.f, // 2^-102
    1/2535301200456458802993406410752.f, // 2^-101
    1/1267650600228229401496703205376.f, // 2^-100
    1/633825300114114700748351602688.f, // 2^-99
    1/316912650057057350374175801344.f, // 2^-98
    1/158456325028528675187087900672.f, // 2^-97
    1/79228162514264337593543950336.f, // 2^-96
    1/39614081257132168796771975168.f, // 2^-95
    1/19807040628566084398385987584.f, // 2^-94
    1/9903520314283042199192993792.f, // 2^-93
    1/4951760157141521099596496896.f, // 2^-92
    1/2475880078570760549798248448.f, // 2^-91
    1/1237940039285380274899124224.f, // 2^-90
    1/618970019642690137449562112.f, // 2^-89
    1/309485009821345068724781056.f, // 2^-88
    1/154742504910672534362390528.f, // 2^-87
    1/77371252455336267181195264.f, // 2^-86
    1/38685626227668133590597632.f, // 2^-85
    1/19342813113834066795298816.f, // 2^-84
    1/9671406556917033397649408.f, // 2^-83
    1/4835703278458516698824704.f, // 2^-82
    1/2417851639229258349412352.f, // 2^-81
    1/1208925819614629174706176.f, // 2^-80
    1/604462909807314587353088.f, // 2^-79
    1/302231454903657293676544.f, // 2^-78
    1/151115727451828646838272.f, // 2^-77
    1/75557863725914323419136.f, // 2^-76
    1/37778931862957161709568.f, // 2^-75
    1/18889465931478580854784.f, // 2^-74
    1/9444732965739290427392.f,	// 2^-73
    1/4722366482869645213696.f,	// 2^-72
    1/2361183241434822606848.f,	// 2^-71
    1/1180591620717411303424.f,	// 2^-70
    1/590295810358705651712.f,	// 2^-69
    1/295147905179352825856.f,	// 2^-68
    1/147573952589676412928.f,	// 2^-67
    1/73786976294838206464.f,	// 2^-66
    1/36893488147419103232.f,	// 2^-65
#endif
    1/18446744073709551616.f,	// 2^-64
    1/9223372036854775808.f,	// 2^-63
    1/4611686018427387904.f,	// 2^-62
    1/2305843009213693952.f,	// 2^-61
    1/1152921504606846976.f,	// 2^-60
    1/576460752303423488.f,	// 2^-59
    1/288230376151711744.f,	// 2^-58
    1/144115188075855872.f,	// 2^-57
    1/72057594037927936.f,	// 2^-56
    1/36028797018963968.f,	// 2^-55
    1/18014398509481984.f,	// 2^-54
    1/9007199254740992.f,	// 2^-53
    1/4503599627370496.f,	// 2^-52
    1/2251799813685248.f,	// 2^-51
    1/1125899906842624.f,	// 2^-50
    1/562949953421312.f,	// 2^-49
    1/281474976710656.f,	// 2^-48
    1/140737488355328.f,	// 2^-47
    1/70368744177664.f,		// 2^-46
    1/35184372088832.f,		// 2^-45
    1/17592186044416.f,		// 2^-44
    1/8796093022208.f,		// 2^-43
    1/4398046511104.f,		// 2^-42
    1/2199023255552.f,		// 2^-41
    1/1099511627776.f,		// 2^-40
    1/549755813888.f,		// 2^-39
    1/274877906944.f,		// 2^-38
    1/137438953472.f,		// 2^-37
    1/68719476736.f,		// 2^-36
    1/34359738368.f,		// 2^-35
    1/17179869184.f,		// 2^-34
    1/8589934592.f,		// 2^-33
    1/4294967296.f,		// 2^-32
    1/2147483648.f,		// 2^-31
    1/1073741824.f,		// 2^-30
    1/536870912.f,		// 2^-29
    1/268435456.f,		// 2^-28
    1/134217728.f,		// 2^-27
    1/67108864.f,		// 2^-26
    1/33554432.f,		// 2^-25
    1/16777216.f,		// 2^-24
    1/8388608.f,		// 2^-23
    1/4194304.f,		// 2^-22
    1/2097152.f,		// 2^-21
    1/1048576.f,		// 2^-20
    1/524288.f,			// 2^-19
    1/262144.f,			// 2^-18
    1/131072.f,			// 2^-17
    1/65536.f,			// 2^-16
    1/32768.f,			// 2^-15
    1/16384.f,			// 2^-14
    1/8192.f,			// 2^-13
    1/4096.f,			// 2^-12
    1/2048.f,			// 2^-11
    1/1024.f,			// 2^-10
    1/512.f,			// 2^-9
    1/256.f,			// 2^-8
    1/128.f,			// 2^-7
    1/64.f,			// 2^-6
    1/32.f,			// 2^-5
    1/16.f,			// 2^-4
    1/8.f,			// 2^-3
    1/4.f,			// 2^-2
    1/2.f,			// 2^-1
    1.f,			// 2^0
    2.f,			// 2^1
    4.f,			// 2^2
    8.f,			// 2^3
    16.f,			// 2^4
    32.f,			// 2^5
    64.f,			// 2^6
    128.f,			// 2^7
    256.f,			// 2^8
    512.f,			// 2^9
    1024.f,			// 2^10
    2048.f,			// 2^11
    4096.f,			// 2^12
    8192.f,			// 2^13
    16384.f,			// 2^14
    32768.f,			// 2^15
    65536.f,			// 2^16
    131072.f,			// 2^17
    262144.f,			// 2^18
    524288.f,			// 2^19
    1048576.f,			// 2^20
    2097152.f,			// 2^21
    4194304.f,			// 2^22
    8388608.f,			// 2^23
    16777216.f,			// 2^24
    33554432.f,			// 2^25
    67108864.f,			// 2^26
    134217728.f,		// 2^27
    268435456.f,		// 2^28
    536870912.f,		// 2^29
    1073741824.f,		// 2^30
    2147483648.f,		// 2^31
    4294967296.f,		// 2^32
    8589934592.f,		// 2^33
    17179869184.f,		// 2^34
    34359738368.f,		// 2^35
    68719476736.f,		// 2^36
    137438953472.f,		// 2^37
    274877906944.f,		// 2^38
    549755813888.f,		// 2^39
    1099511627776.f,		// 2^40
    2199023255552.f,		// 2^41
    4398046511104.f,		// 2^42
    8796093022208.f,		// 2^43
    17592186044416.f,		// 2^44
    35184372088832.f,		// 2^45
    70368744177664.f,		// 2^46
    140737488355328.f,		// 2^47
    281474976710656.f,		// 2^48
    562949953421312.f,		// 2^49
    1125899906842624.f,		// 2^50
    2251799813685248.f,		// 2^51
    4503599627370496.f,		// 2^52
    9007199254740992.f,		// 2^53
    18014398509481984.f,	// 2^54
    36028797018963968.f,	// 2^55
    72057594037927936.f,	// 2^56
    144115188075855872.f,	// 2^57
    288230376151711744.f,	// 2^58
    576460752303423488.f,	// 2^59
    1152921504606846976.f,	// 2^60
    2305843009213693952.f,	// 2^61
    4611686018427387904.f,	// 2^62
    9223372036854775808.f,	// 2^63
    18446744073709551616.f,	// 2^64
  };
#endif
} // namespace RandomLib
