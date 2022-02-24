/*
  PURPOSE

  I made a wrapper for Random library by Karney in order to use Random
  in Intel Threading Building Blocks. I want each thread to have a
  separate copy of Random class with a distinct
  seed. Tbb::enumerable_thread_specific class will be used to achive
  this. Additionally, I include 'int index' as a thread-specific
  variable, which can uniquely identify different thread. 

  ------------------------------------------------------------------
  USAGE

  Before using random number generators, we first have to set the
  seed. For this, random_wrapper.cpp defines RandomWrapper::seed. It is
  std::vector<unsigned long> type and is default constructed, so that it
  is initially empty. Thus, in order to set a seed, do for example

  RandomWrapper::seed.push_back(12345UL);

  After setting a seed, we can use random number generator.
  random_wrapper.cpp defines random_tbb. In order to use a random number
  generator, do

  RandomTBB::reference what_ever = random_tbb.local();

  Now, what_ever.random is the "thread-specific" Random object of the
  library by Karney. For example, one can now do
  
  row = static_cast<unsigned>(what_ever.random.Real()*nrow + 1);

  Moreover, random_tbb.index will give an integer unique to each
  thread. This index will be greater than 0 (see IMPLEMENTATION)
  
  ------------------------------------------------------------------
  IMPLEMENTATION

  According to some experiments, it appeared that
  enumerable_thread_specific<T>, when constructed, constructs and hold
  one instance of T with the default-constructor (unless an exepmpler is
  supplied to enumerable_thread_specific's constructor, in which case
  copy-constructor is used). This instance of T contained is
  thread-unspecific (not "owned" by any thread). Instead, from this,
  enumerable_thread_specific<T>, when requested by a thread constructs,
  if it's hasn't constructed, constructs a unique instance of T for each
  thread by using copy-constructor. Once the unique instance is
  constructed, enumerable_thread_specific will simply return it when
  requested by the same thread again. Thus, a proper definition of
  copy-constructor is essential for our purpose.

  The default constructor of RandomWrapper will construct Random by
  Random(0) and set "index" to 1. This Random is a dummy (i.e. unused)
  that is cheaper than the default constructor of Random. The copy
  constructor of RandomWrapper will construct Random with Random(0) and
  Reseed() it with a thread-specific seed by using a local copy of the
  static member "seed" appended with the "index", which counts the
  number of instances ever constructed (hence, a first copy-constructed
  instance's index will be 2). Hence, the thread-unspecific instance of
  Random will have [0] as its seed, whereas the thread-specific
  instances of Random will have [seed,index] as its seed (unless
  something else than enumerable_thread_specific<Random> is calling the
  constructors of RandomWrapper, which are not supposed to happen).
*/

#include <iostream>
#include <vector>
#include <unistd.h>
#include "RandomLib/Random.hpp"
#include "tbb/tbb.h"


#ifndef RANDOM_WRAPPER
#define RANDOM_WRAPPER

class RandomWrapper {
  static tbb::atomic<int> counter;
public:
  static std::vector<unsigned long> seed;

  /* A unique index of each instance. It runs 1, 2, ... */
  const int index;
  RandomLib::Random random;

  inline RandomWrapper();
  inline RandomWrapper(const RandomWrapper& rhs);
};

typedef tbb::enumerable_thread_specific<RandomWrapper> RandomTBB;
extern RandomTBB random_tbb;

RandomWrapper::RandomWrapper()
  : index(++counter), random(0)
{
  std::cerr << "RandomWrapper::RandomWrapper() Default constructed." << std::endl;
}

RandomWrapper::RandomWrapper(const RandomWrapper& rhs)
  : index(++counter), random(0)
{
  if(seed.empty()){
    std::cerr << "RandomWrapper::RandomWrapper(const RandomWrapper&) Error, RandomWrapper::seed is not initialized." << std::endl;
    exit(-1);
  }
  std::vector<unsigned long> myseed(seed);
  myseed.push_back(index);
  random.Reseed(myseed);

  std::cerr << "RandomWrapper::RandomWrapper(const RandomWrapper&) Random copy-constructed with seed = " << random.SeedString() << std::endl;
  
}

#endif
