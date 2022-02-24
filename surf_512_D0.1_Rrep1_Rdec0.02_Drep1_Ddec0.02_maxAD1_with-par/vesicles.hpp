#include <iostream>
#include <fstream>
#include <unistd.h>
#include <limits>
#include <cmath>
#include <deque>

#include "assert.hpp"

#include "para.hpp"
#include "my_dist.hpp"
#include "vesicle-info.hpp"

#ifndef VESICLES
#define VESICLES


class SerialVesicles {
private:
  /* The size of array. The number of vesicles is bounded by this. */
  int length;
  /* If volume[i]==-1, then the i-th vesicle is non-existent. volume[]
     is used to flag which indexes are unused  */
  typedef std::vector<VesicleInfo> VesicleInfoV;
  VesicleInfoV vesicle_info_v;

public:
  SerialVesicles(int a_length);

  int get_unused_index() const;
  /* It initilize vesicle info that are not used and returns its
     index. */
  int make_new_vesicle();
  /* It tries to initilize vesicle infor with a given index. If this
     index has been used, it exits with error diagnosis. */
  void make_new_vesicle(const int index);
  void destroy_vesicle(const int index);

  /* This method sets the temporary variables (sum_x & sum_y etc. except
     for a, b & c) for principle component analysis to zero and set
     is_divide to false for all existing vesicles (i.e. volume>-1).*/
  void PCA_reset();
  /* This method calculate sum_x etc. for vesicles whose index is
     greater or equal to the division threshold. */
  void PCA_collect_data(int index,unsigned row,unsigned col);
  /* This method calculate PCA from the data collected by
     PCA_collect_data. It also sets is_divide according to whether the
     volume is smaller than the threshold, and it also initialize
     daughter_ind by setting it to -1. It will return true if there is
     any vesicle supposed to divide. */
  bool PCA_calculate();
  /* This method calculate PCA, using the data collected by
     PCA_collect_data(), for those vesicles that have is_divide=true. */
  void calculate_PCA_for_dividing_vesicles();

  void divide_target_volume();
  void decay();
  void decay(const int);

  int get_volume(const int index) const {return vesicle_info_v[index].get_volume();}
  void add_volume(const int index,const int incl=1) {vesicle_info_v[index].add_volume(incl);}

  int get_target_volume(const int index) const {return vesicle_info_v[index].target_volume;}
  void add_target_volume(const int index,const int incl=1){if(vesicle_info_v[index].target_volume < Para::max_target_volume) vesicle_info_v[index].target_volume += incl;}
  void set_target_volume(const int index,const int val);
  void set_all_target_volume_to_0();
  void neutral_growth();

  bool is_divide(const int index) const;
  void set_is_divide(const int index);
  bool is_daughter(const int index,const unsigned row,const unsigned col) const;
  int get_daughter_ind(const int index) const;
  void set_daughter_ind(const int index,const int daughter_ind);

  long get_birth_time(const int index) const;
  void set_birth_time(const int index,const long Time);

  /* output the distribution of eigen_1 and volume */
  void dump_eigen(const long Time) const;

#ifdef  PLOT_VESICLE_INTERNAL_DYNAMICS
  void set_identity(const int index, const std::deque<bool>& _identity);
  /* This method changes the identity of this vesicle (mother vesicle)
     that has just gave birth to a new vesicle, and it also gives an
     identity to a daughter vesicle. Thus, this method gives an identity
     to a new born vesicle, while updating the identity of the vesicle
     that produced the new vesicle. */
  void make_identity(const int mother_ind,const int daughter_ind);
  const std::deque<bool>& get_identity(int index) const;
#endif //PLOT_VESICLE_INTERNAL_DYNAMICS

private:
  void set_volume(const int index,const int vol) {vesicle_info_v[index].set_volume(vol);}
};

#endif //VESICLES
