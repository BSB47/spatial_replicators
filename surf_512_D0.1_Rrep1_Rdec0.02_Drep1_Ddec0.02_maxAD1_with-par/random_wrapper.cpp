#include "random_wrapper.hpp"

// This will be zeroed by default.
tbb::atomic<int> RandomWrapper::counter;

// This will be empty by default
std::vector<unsigned long> RandomWrapper::seed;

RandomTBB random_tbb;

