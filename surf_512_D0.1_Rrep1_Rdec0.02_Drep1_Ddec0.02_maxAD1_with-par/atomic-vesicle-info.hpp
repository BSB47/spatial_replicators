#include <iostream>
#include <cmath>
#include <deque>

#include "assert.hpp"

/* Intel TBB */
#include "tbb/tbb.h"

/* This is an atomic version of VesicleInfo. This will be used by
   AtomicVesicles. */

/* This class holds information specific about vesicles. The definition
   of this class does not depend on how replicator dynamics is assumed
   in the model. */

#ifndef ATOMIC_VESICLE_INFO
#define ATOMIC_VESICLE_INFO

class AtomicVesicleInfo {
  tbb::atomic<int> volume;

public:
  tbb::atomic<int> target_volume;

  /* if this flag is true, this vesicle is supposed to divide. This flag
     is initilized to false by AtomicVesicles::PCA_reset() and is updated
     by MyCA::ParallelDivideCollectData::operator(). */
  tbb::atomic<bool> is_divide;

  /* the followings are used to calculate principle components of the CA
     area this vesicle occupies */
  tbb::atomic<unsigned> sum_x;
  tbb::atomic<unsigned> sum_y;
  tbb::atomic<unsigned> sum_xx;
  tbb::atomic<unsigned> sum_xy;
  tbb::atomic<unsigned> sum_yy;
  /* the follwoings are the coffeficients a, b & c in f(X,Y)=aX+bY+c,
     where X=row and Y=col. The sign of f(X,Y) determins how to devide a
     vesicle into two. These coefficients are calculated from principle
     components. For the details of calculation, look at myca.[ch]pp. */
  double a;
  double b;
  double c;
  /* the followings are eigen values */
  double eigen1;
  double eigen2;

  /* This variable holds the index of the daughter vesicle. This is used
     (at least) when division takes place. If daughter_ind==-1, it means
     its daughoter hasn't got an index. To signify this, daughter_ind
     should be initilized to -1 for those vesicles that should
     devide. This will be done by AtomicVesicles::PCA_reset(). */
  tbb::atomic<int> daughter_ind;

  // For a moment, I disable this feature

//   /* This tells when this vesicle is born. This and the vesicle index
//      give a unique identifier of individuals. */
//   long birth_time;

  // For a moment, I disable this feature

//   /* This gives not only an unique indentifier of individual vesicles,
//      but also a unique way to denote lineage relationship among
//      vesicles. Because of this, the identifier becomes too lengthy for
//      long simulations. It is used to make time plots of internal
//      dynamics of vesicles. */
// #ifdef PLOT_VESICLE_INTERNAL_DYNAMICS
//   std::deque<bool> v_identity;
// #endif

  AtomicVesicleInfo() {volume=-1;}
  /* this method tells whether the position (row,col) belongs to the
     daughter vesicle for those vesicles that is supposed to divide. The
     name, mother & daughtor, is arbtrary. Before using this method, PCA
     must be calculated by using Vesicles::PCA_collect_data() and
     Vesicles::PCA_calculate(). */

  // is_daughter must be re-implemented for atomic-vesicle-info
  bool is_daughter(unsigned row,unsigned col) const {return (a*col+b*row+c > 0.);}

  // calculate_PCA() must be re-implemented for atomic-vesicle-info
  void calculate_PCA();

  int get_volume() const {return volume;}
  void add_volume(const int incl) {volume += incl;}
  void set_volume(const int vol) {volume = vol;}

  /* This checks if volume==-1. If so, it sets volume=0, and return
     -1. If volume!=-1, then it the variable's value will remain the
     same and is returned by the method. */
  int resurrect() {return volume.compare_and_swap(0,-1);}

  void set_is_divide(const bool state) {is_divide=state;}
  bool get_is_divide(const bool state) {return is_divide;}
};

#endif
