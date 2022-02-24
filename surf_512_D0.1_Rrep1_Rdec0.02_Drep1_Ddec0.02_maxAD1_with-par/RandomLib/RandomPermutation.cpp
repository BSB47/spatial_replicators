/**
 * \file RandomPermutation.cpp
 * \brief Prints a random permutation of integers
 *
 * Written by <a href="http://charles.karney.info/">Charles
 * Karney</a> <charles@karney.com> and licensed under the GPL.
 * For more information, see http://charles.karney.info/random/
 **********************************************************************/

#include "RandomLib/Random.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace {
  char rcsid[] = "$Id: RandomPermutation.cpp 6380 2007-04-20 23:54:15Z ckarney $";
}

int main(int argc, char* argv[]) {
  unsigned n = 100;
  unsigned base = 10;
  for (int m = 1; m < argc; ++m) {
    std::string arg = std::string(argv[m]);
    if (arg == "-o")
      base = 8;
    else if (arg == "-d")
      base = 10;
    else if (arg == "-x")
      base = 16;
    else if (arg == "-h") {
      std::cerr
	<< "Usage: RandomPermutation [-o] [-d] [-x] [-h] [num]\n"
	<< "\n"
	<< "Print a random permutation of numbers from 0 thru num-1\n"
	<< "on standard output.  num is supplied on the command line\n"
	<< "as a decimal number (default is 100).  Optional arguments\n"
	<< "-o, -d, and -x selection octal, decimal, and hexadecimal\n"
	<< "output base (default decimal). -h prints this help.\n";
      return 0;
    }
    else {
      std::istringstream str(arg);
      str >> n;
    }
  }

  unsigned k = 0;
  for (unsigned i = n - 1; i; i /= base) k++;

  std::vector<unsigned> a(n);
  for (unsigned i = n; i;) { --i; a[i] = i; }

  RandomLib::Random r(RandomLib::Random::SeedVector());
  std::random_shuffle(a.begin(), a.end(), r);
  
  std::cout << std::setfill('0')
	    << (base == 16 ? std::hex : (base == 8 ? std::oct : std::dec));
  for (unsigned i = n; i > 0;)
    std::cout << std::setw(k) << a[--i] << "\n";

  return 0;
}
