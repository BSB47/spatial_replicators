#include <iostream>
#include <vector>
#include <algorithm>
#include <unistd.h>

#include "para.hpp"
#include "random_wrapper.hpp"

#ifndef UPD_ORDER
#define UPD_ORDER

/*
  See also the explanation in myca.hpp
*/


namespace _UpdOrderShuffler{
  /* The number of UpdOrder hold by UpdOrderShuffler. When the
     parallelized algorithm is used, this is the number of arrays
     shuffled in parallel. */
  const unsigned upd_order_array_size = 8;
}

/* this is used for update ordre */
struct Point{
  unsigned row;
  unsigned col;
};

typedef std::vector<Point> UpdOrder;

/* 
   The class contains upd_order_array, which is an array of UpdOrder the
   size of which is upd_order_array_size. In the serial mode, it simply
   shuffles upd_order_array[0] and returns a pointer to it. In the
   parallel mode, in the 1st call, it shuffles all elements of
   upd_order_array and returns a const reference to
   upd_order_array[0]. In the 2nd call, it returns a const reference to
   upd_order_array[1]. In the [upd_order_array_size+1]'th call, since it
   has ran out all the shuffled upd_order_array's elements, the class
   shuffles all the elements of upd_order_array again and returns a
   const reference to upd_order_array[0].
*/

class UpdOrderShuffler {
public:
  UpdOrderShuffler(const unsigned _nrow, const unsigned _ncol);
  inline const UpdOrder& get_upd_order_serial();
  const UpdOrder& get_upd_order_parallel();

private:
  /* Members */
  /* the index of upd_order_array which will be returned when the class
     is requested for a new update order next time. */
  unsigned curr_pos_upd_order_array;
  const unsigned nrow;
  const unsigned ncol;
  const unsigned length;
  typedef std::vector<UpdOrder> UpdOrderArray;
  UpdOrderArray upd_order_array;
  
  /* Methods & Others */
  class ParallelShuffle {
    UpdOrderArray& my_upd_order_array;
  public:
    ParallelShuffle(UpdOrderArray& _upd_order_array)
      : my_upd_order_array(_upd_order_array)
    {}
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
  };
};

void UpdOrderShuffler::ParallelShuffle::operator()(const tbb::blocked_range<size_t>& r) const
{
  /* Get the thread-specific random number generator. */
  RandomTBB::reference random_wrapper = random_tbb.local();
  
  for(size_t i=r.begin(), e=r.end();i!=e;++i){
    std::random_shuffle(my_upd_order_array[i].begin(),my_upd_order_array[i].end(),random_wrapper.random);
  }
}

const UpdOrder& UpdOrderShuffler::get_upd_order_serial()
{
  std::random_shuffle(upd_order_array[0].begin(),upd_order_array[0].end(),rand_karney);
  return upd_order_array[0];
}

#endif
