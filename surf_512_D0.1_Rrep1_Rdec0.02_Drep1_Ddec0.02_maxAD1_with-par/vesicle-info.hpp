#include <iostream>
#include <cmath>
#include <deque>
#include <limits>

#include "assert.hpp"

/* This class holds information specific about vesicles. The definition
   of this class does not depend on how replicator dynamics is assumed
   in the model. */

#ifndef VESICLE_INFO
#define VESICLE_INFO

class VesicleInfo {
  int volume;

public:
  int target_volume;

  /* if this flag is true, this vesicle is supposed to divide. This flag
     is initialized by Vesicles::PCA_reset in MyCA::divide_vesicle and
     set by MyCA::divide_vesicle(). */
  bool is_divide;

  /* the followings are used to calculate principle components of the CA
     area this vesicle occupies */
  unsigned sum_x;
  unsigned sum_y;
  unsigned sum_xx;
  unsigned sum_xy;
  unsigned sum_yy;
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
     devide. This will be done by Vesicles::PCA_reset(). */
  int daughter_ind;
  /* This tells when this vesicle is born. This and the vesicle index
     give a unique identifier of individuals. */
  long birth_time;

  /* This gives not only an unique indentifier of individual vesicles,
     but also a unique way to denote lineage relationship among
     vesicles. Because of this, the identifier becomes too lengthy for
     long simulations. It is used to make time plots of internal
     dynamics of vesicles. */
#ifdef PLOT_VESICLE_INTERNAL_DYNAMICS
  std::deque<bool> v_identity;
#endif
  VesicleInfo(): volume(-1) {}
  /* this method tells whether the position (row,col) belongs to the
     daughter vesicle for those vesicles that is supposed to divide. The
     name, mother & daughtor, is arbtrary. Before using this method, PCA
     must be calculated by using Vesicles::PCA_collect_data() and
     Vesicles::PCA_calculate(). */
  bool is_daughter(unsigned row,unsigned col) const {return (a*col+b*row+c > 0.);}
  void calculate_PCA();
  int get_volume() const {return volume;}
  void add_volume(const int incl) {volume += incl;}
  void set_volume(const int vol) {volume = vol;}
};

#endif
