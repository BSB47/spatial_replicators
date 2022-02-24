#include "upd_order.hpp"

UpdOrderShuffler::UpdOrderShuffler(const unsigned _nrow, const unsigned _ncol)
  : curr_pos_upd_order_array(_UpdOrderShuffler::upd_order_array_size),
    nrow(_nrow),
    ncol(_ncol),
    length(nrow*ncol)
{
  /* A simple error check */
  if(_UpdOrderShuffler::upd_order_array_size==0){
    std::cerr << "UpdOrderShuffler::UpdOrderShuffler() upd_order_array_size=0 is not allowed." << std::endl;
    exit(-1);
  }


  /* Initialize upd_order_array. First, we make upd_order_template,
     which contains the coordinate of all points. Then, we insdert n
     copies of it in upd_order_array. */
  {
    UpdOrder upd_order_template(length);
    UpdOrder::iterator iter=upd_order_template.begin();
    for(unsigned r=1,nr=nrow;r<=nr;++r)
      for(unsigned c=1,nc=ncol;c<=nc;++c){
	iter->row = r;
	iter->col = c;
	++iter;
      }
    upd_order_array.insert(upd_order_array.end(),_UpdOrderShuffler::upd_order_array_size,upd_order_template);
  }
}

const UpdOrder& UpdOrderShuffler::get_upd_order_parallel()
{
  /* If there are still upd_oder_array's emelements that have been
     shuffled and haven't been used. */
  if(curr_pos_upd_order_array<_UpdOrderShuffler::upd_order_array_size){
    return upd_order_array[curr_pos_upd_order_array++];
  }
  else{
    tbb::parallel_for(tbb::blocked_range<size_t>(0,_UpdOrderShuffler::upd_order_array_size),ParallelShuffle(upd_order_array));
    curr_pos_upd_order_array=1;
    return upd_order_array[0];
  }
}
