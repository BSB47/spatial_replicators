#include <iostream>
#include <fstream>
#include <unistd.h>
#include <limits>
#include <cmath>
#include <deque>

#include "assert.hpp"

#include "para.hpp"
#include "my_dist.hpp"
#include "atomic-vesicle-info.hpp"
#include "tbb/tbb.h"

/*
  This is the version of Vesicles that should be used in the
  parallelized version of CPM. It contains an array of
  AtomicVesicleInfo, the member variables of which are mostly atomic.

*/

#ifndef ATOMIC_VESICLES
#define ATOMIC_VESICLES

class AtomicVesicles {
private:
  /* The size of array. The number of vesicles is bounded by this. */
  int length;
  /* If volume[i]==-1, then the i-th vesicle is non-existent. volume[]
     is used to flag which indexes are unused  */
  typedef std::vector<AtomicVesicleInfo> AtomicVesicleInfoV;
  AtomicVesicleInfoV vesicle_info_v;

public:
  AtomicVesicles(int a_length);

  /*** THE METHODS BELOW ARE THREAD SAFE ***/

  /*** THE METHODS BELOW ARE **NOT** THREAD SAFE ***/

  int get_unused_index() const;

  /* It initilize vesicle info that are not used and returns its
     index. */
  int make_new_vesicle();

  /* It tries to initilize vesicle infor with a given index. If this
     index has been used, it exits with error diagnosis. */
  void make_new_vesicle(const int index);

  /* It sets all volume whose value is 0 to -1. It uses
     parallelization. */
  inline void destroy_all_empty_vesicles();

  /* This method sets the temporary variables (sum_x & sum_y etc. except
     for a, b & c) for principle component analysis to zero and set
     is_divide to false for all existing vesicles (i.e. volume>-1).*/
  inline void PCA_reset();

  /* This method calculate sum_x etc. for vesicles whose index is
     greater or equal to the division threshold. */
  void PCA_collect_data(int index,unsigned row,unsigned col);

  /* This method calculate PCA, using the data collected by
     PCA_collect_data(), for those vesicles that have is_divide=true. */
  inline void calculate_PCA_for_dividing_vesicles();

  /* This method divides the target volume between mother and daughter
     vesicles. */
  inline void divide_target_volume();
  void decay();
  void decay(const int);

  int get_volume(const int index) const {return vesicle_info_v[index].get_volume();}
  void add_volume(const int index,const int incl=1) {vesicle_info_v[index].add_volume(incl);}
  
  int get_target_volume(const int index) const {return vesicle_info_v[index].target_volume;}
  void add_target_volume(const int index,const int incl=1){if(vesicle_info_v[index].target_volume < Para::max_target_volume) vesicle_info_v[index].target_volume += incl;}
  void set_target_volume(const int index,const int val);
  inline void set_all_target_volume_to_0();
  inline void neutral_growth();

  bool is_divide(const int index) const {return vesicle_info_v[index].is_divide;}
  void set_is_divide(const int index);
  bool is_daughter(const int index,const unsigned row,const unsigned col) const {return vesicle_info_v[index].is_daughter(row,col);}
  int get_daughter_ind(const int index) const {return vesicle_info_v[index].daughter_ind;}
  void set_daughter_ind(const int index);

  /* output the distribution of eigen_1 and volume */
  void dump_eigen(const long Time) const;
  
  
  // I disable this feature until I finish the development of the core
  // part of the parallelized CPM..
  //   long get_birth_time(const int index) const;
  //   void set_birth_time(const int index,const long Time);
  
  // #ifdef  PLOT_VESICLE_INTERNAL_DYNAMICS
  //   void set_identity(const int index, const std::deque<bool>& _identity);
  //   /* This method changes the identity of this vesicle (mother vesicle)
  //      that has just gave birth to a new vesicle, and it also gives an
  //      identity to a daughter vesicle. Thus, this method gives an identity
  //      to a new born vesicle, while updating the identity of the vesicle
  //      that produced the new vesicle. */
  //   void make_identity(const int mother_ind,const int daughter_ind);
  //   const std::deque<bool>& get_identity(int index) const;
  // #endif //PLOT_VESICLE_INTERNAL_DYNAMICS

private:
  void set_volume(const int index,const int vol) {vesicle_info_v[index].set_volume(vol);}

  /* This sets the volume of vesicle_info_v[index] to -1. */
  void destroy_vesicle(const int index);

  /* To parallelize destroy_all_empty_vesicles() */
  class DestroyVesicleParallel {
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
    DestroyVesicleParallel(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };

  /* To parallelize PCA_reset() */
  class PCAReset {
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    void operator()(const tbb::blocked_range<size_t>& r) const;
    PCAReset(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };

  /* To parallelize calculate_PCA_for_dividing_vesicles() */
  class CalculatePCAForDividingVesicles {
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
    CalculatePCAForDividingVesicles(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };

  /* To parallelize neutral_growth() */
  class NeutralGrowth{
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    void operator()(const tbb::blocked_range<size_t>& r) const;
    NeutralGrowth(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };

  /* To parallelize divide_target_volume() */
  class DivideTargetVolume{
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    void operator()(const tbb::blocked_range<size_t>& r) const;
    DivideTargetVolume(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };

  /* To parallelize set_all_target_volume_to_0() */
  class SetAllTargetVolumeTo0 {
    AtomicVesicleInfoV& my_vesicle_info_v;
  public:
    void operator()(const tbb::blocked_range<size_t>& r) const;
    SetAllTargetVolumeTo0(AtomicVesicleInfoV& _my_vesicle_info_v)
      : my_vesicle_info_v(_my_vesicle_info_v) {}
  };
};

void AtomicVesicles::DestroyVesicleParallel::operator()(const tbb::blocked_range<size_t>& r) const
{
  AtomicVesicleInfoV& local_my_vesicle_info_v = my_vesicle_info_v;
  for(size_t i=r.begin();i!=r.end();++i){
    if(local_my_vesicle_info_v[i].get_volume()==0)
      local_my_vesicle_info_v[i].set_volume(-1);
  }

  /* By the way, I seem to have hit a mistake in the core standard of
     C++ here (before C++0x). According to Stroustrup 2003, C.11.3, "the
     members of a member class have no special access to members of an
     enclosing class". Thus, the above code is illegal because
     DestroyVesicleParallel is a member class of Vesicles, and
     vesicle_infor_v is a PRIVATE member of Vesicles. However, according
     to Defect Report 45
     (http://www.open-std.org/jtc1/sc22/wg21/docs/cwg_defects.html#45),
     this was undone by saying, "a nested class is a member and as such
     has the same access rights as any other member". In my opinion, the
     latter seems more consistent too. Seemingly, this correction is to
     be included in C++0x
     (http://cboard.cprogramming.com/cplusplus-programming/110207-can-nested-class-access-parent-classs-private-member.html). */
}

void AtomicVesicles::destroy_all_empty_vesicles()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),DestroyVesicleParallel(vesicle_info_v)); 
}

void AtomicVesicles::PCA_reset()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),PCAReset(vesicle_info_v)); 
}

void AtomicVesicles::CalculatePCAForDividingVesicles::operator()(const tbb::blocked_range<size_t>& r) const
{
  AtomicVesicleInfoV& loc_my_vesicle_info_v = my_vesicle_info_v;
  for(size_t i=r.begin();i!=r.end();++i){
    if(loc_my_vesicle_info_v[i].get_volume()!=-1){
      if(loc_my_vesicle_info_v[i].is_divide){
	loc_my_vesicle_info_v[i].calculate_PCA();
      }
    }
  }
}

void AtomicVesicles::calculate_PCA_for_dividing_vesicles()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),CalculatePCAForDividingVesicles(vesicle_info_v)); 
}

/* This is used by MyCA::vesicle_neutral_growth. */
void AtomicVesicles::neutral_growth()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),NeutralGrowth(vesicle_info_v)); 
}

void AtomicVesicles::divide_target_volume()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),DivideTargetVolume(vesicle_info_v));
}

/* This is used by MyCA::vesicle_neutral_growth. */
void AtomicVesicles::set_all_target_volume_to_0()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,vesicle_info_v.size()),SetAllTargetVolumeTo0(vesicle_info_v));
}

#endif //VESICLES
