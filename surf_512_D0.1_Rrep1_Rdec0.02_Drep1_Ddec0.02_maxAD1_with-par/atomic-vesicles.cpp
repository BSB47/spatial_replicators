#include "atomic-vesicles.hpp"

AtomicVesicles::AtomicVesicles(int a_length)
  : length(a_length),
    vesicle_info_v(a_length)
{}

/* This function is not used often, so I include error cheking. */
void AtomicVesicles::set_target_volume(const int index,const int val)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  
  catch(GeneralError){
    std::cerr << "AtomicVesicles::set_target_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].target_volume = val;
}

void AtomicVesicles::SetAllTargetVolumeTo0::operator()(const tbb::blocked_range<size_t>& r) const
{
  AtomicVesicleInfoV& local_vesicle_info_v = my_vesicle_info_v;
  for(size_t i=r.begin();i!=r.end();++i){
    if(local_vesicle_info_v[i].get_volume()!=-1){
      local_vesicle_info_v[i].target_volume = 0;
    }
  }
}

void AtomicVesicles::NeutralGrowth::operator()(const tbb::blocked_range<size_t>& r) const
{
  const double factor = Para::neutral_growth_expand_factor;
  AtomicVesicleInfoV& local_vesicle_info_v = my_vesicle_info_v;
  for(size_t i=r.begin();i!=r.end();++i){
    if(local_vesicle_info_v[i].get_volume()!=-1){
      local_vesicle_info_v[i].target_volume = static_cast<int>(factor * local_vesicle_info_v[i].target_volume + 0.5);
    }
  }
}

/* This function is not used often, so I include error cheking. */
void AtomicVesicles::set_is_divide(const int index)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::set_is_divide() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::set_is_divide() Error, try to get set_is_divide of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].is_divide = true;
}

/* Thread-safe version of set_daughter_ind() requires a modification to
   the non-thread-safe version of it (i.e. that in vesicles.cpp). This
   is due to the way AtomicVesicles::set_daughter_ind() is used by
   ParallelDivideDivideVes::operator(), which is in
   MyCA::parallel_divide_vesicle(). Therein, each thread first checks if
   daughter_ind of a dividing vesicle has already received a valid
   value. If not (i.e. daughter_ind==-1), the thread "tries" to put a
   valid value in it. However, the problem is that if multiple threads
   concurrently try to set daughter_ind of an identical vesicle, there
   is a chance that one or more of vesicle indexes are set to alive
   state (i.e. vesicle_info_v[index]>-1), but, in fact, those indexes
   are never used. And, I want to avoid this non-thread-safe nature
   (despite the fact that AtomicVesicle::DestroyVesicleParallel()
   actually solves this problem by seting to -1 the volume of all the
   elements of vesicle_info_v that have volume=0).  
*/
void AtomicVesicles::set_daughter_ind(const int index)
{
  /* It's non-inline function, so I use error check. */
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::set_daughter_ind() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::set_daughter_ind() Error, try to set daughter index of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  /* For the purpose explained above, I first obtain a possible index of
   a daughter vesicle through AtomicVesicles::make_new_vesicle() and use
   swap_and_compare() to try to give this value to the daughter_ind of
   vesicle_infor_v[index]. If it's unsuccessful, I return the index
   obtained from make_new_vesicle() into the pool. */
  const int daughter_ind = make_new_vesicle();

  /* if swap_and_compare() returns a value other than -1, setting a
     value to the daughter_ind is unsuccessful: somebody else has
     already given a valid value to the daughter_ind. */
  if(vesicle_info_v[index].daughter_ind.compare_and_swap(daughter_ind,-1)!=-1){
    destroy_vesicle(daughter_ind);
  }
}

/* We use compare_and_swap() to safely make a new vesicle (i.e. find a
   vesicle_index that isn't used and set the volume of it to 0). This
   makes it thread safe. */
int AtomicVesicles::make_new_vesicle()
{
  int index=-1;

  for(AtomicVesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->resurrect()==-1){
      index = pos-vesicle_info_v.begin();
      break;
    }
  }

  if(index==-1){
    std::cerr << "AtomicVesicles::make_new_vesicle() We run out vesicle index" << std::endl;
    exit(-1);
  }

  /* Setting the volume to 0 is unnecessary here, because resurrect()
     has done it already. To set is_divide=false is crucial in correct
     working of MyCA::parallel_divide_vesicle(). This method goes though
     an entire CA field to change the CPM state of those CA squares that
     belong to new born daughter vesicles. To decide whether the method
     should change the CPM state or not is dependent on the value of
     is_divide of the (mother) vesicles. If a new born daugher vesicle
     has is_divide=true, the method incorrectly consider that this
     daughter vesicle is a dividing mother vesicle. (Strictly speaking,
     is_divide and daughter_ind are redundant. Actually, I can dispence
     only with daughter_ind. But I didn't do that.) */
  vesicle_info_v[index].set_is_divide(false);
  
  return index;
}

/* This is similar to make_new_vesicle(void). It uses compare_and_swap()
   to make a given vesicle_index "alive". */
void AtomicVesicles::make_new_vesicle(const int index)
{
  /* These error checks must be here even without debugging, because
     they can check the error while reading an input file. */
  if(index >= length){
    std::cerr << "AtomicVesicles::make_new_vesicle(int) Error, index is out of range" << std::endl;
    exit(-1);
  }

  if(vesicle_info_v[index].resurrect()!=-1){
    std::cerr << "AtomicVesicles::make_new_vesicle(int) Error, the vesicle with this index has already been made." << std::endl;
    exit(-1);
  }
  
  vesicle_info_v[index].set_is_divide(false);
}

/* This function is used in AtomicVesicles::set_daughter_ind(). */
void AtomicVesicles::destroy_vesicle(const int index)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::destroy_vesicle() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()==0);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::destroy_vesicle() Error, try to destroy a vesicle whose volume is not 0." << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].set_volume(-1);
}

void AtomicVesicles::PCAReset::operator()(const tbb::blocked_range<size_t>& r) const
{
  AtomicVesicleInfoV& local_my_vesicle_info_v = my_vesicle_info_v;

  for(size_t i=r.begin();i!=r.end();++i){
    /* this is just-in-case error check */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||local_my_vesicle_info_v[i].get_volume()>-2);
    }
    catch(GeneralError){
      std::cerr << "AtomicVesicles::PCA_reset::operator() Error, volume has an abnormal value=" << local_my_vesicle_info_v[i].get_volume() << std::endl;
      exit(-1);
    }

    /* If this vesicle exists, initilize PCA variables & set is_divide
       flag */
    if(local_my_vesicle_info_v[i].get_volume() != -1){
      local_my_vesicle_info_v[i].eigen1 = local_my_vesicle_info_v[i].sum_x = local_my_vesicle_info_v[i].sum_y = local_my_vesicle_info_v[i].sum_xx = local_my_vesicle_info_v[i].sum_yy = local_my_vesicle_info_v[i].sum_xy = 0;
      local_my_vesicle_info_v[i].is_divide = false;
      local_my_vesicle_info_v[i].daughter_ind = -1;
    }

  }

  /* See my mutter in AtomicVesicles::DestroyVesicleParallel::operator(). */
}

void AtomicVesicles::PCA_collect_data(int index,unsigned row,unsigned col)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1 && index<length));
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::PCA_collect_data() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].sum_x += col;
  vesicle_info_v[index].sum_y += row;
  vesicle_info_v[index].sum_xx += col*col;
  vesicle_info_v[index].sum_yy += row*row;
  vesicle_info_v[index].sum_xy += row*col;
}

void AtomicVesicles::DivideTargetVolume::operator()(const tbb::blocked_range<size_t>& r) const
{
  int target;
  AtomicVesicleInfoV& local_vesicle_info_v = my_vesicle_info_v;
  for(size_t i=r.begin();i!=r.end();++i){
    if(local_vesicle_info_v[i].get_volume()!=-1){
      if(local_vesicle_info_v[i].is_divide){
	try{
	  Assert<GeneralError>(!ASSERT::ERROR_CHECK||(local_vesicle_info_v[i].get_volume()>0 && local_vesicle_info_v[local_vesicle_info_v[i].daughter_ind].get_volume() > 0));
	}
	catch(GeneralError){
	  std::cerr << "AtomicVesicles::DivideTargetVolume::operator() Error, volume of the mother vesicle=" << local_vesicle_info_v[i].get_volume() << " that of the daughter=" << local_vesicle_info_v[local_vesicle_info_v[i].daughter_ind].get_volume() << std::endl;
	  exit(-1);
	}

	target=local_vesicle_info_v[i].target_volume;
	local_vesicle_info_v[i].target_volume = static_cast<int>(static_cast<double>(target*local_vesicle_info_v[i].get_volume())/(local_vesicle_info_v[i].get_volume()+local_vesicle_info_v[local_vesicle_info_v[i].daughter_ind].get_volume()) + 0.5);
	local_vesicle_info_v[local_vesicle_info_v[i].daughter_ind].target_volume = target - local_vesicle_info_v[i].target_volume;
      }
    }
  }
}

void AtomicVesicles::decay()
{
  const double decay_rate = Para::vesicle_decay_rate;
  for(AtomicVesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      pos->target_volume -= static_cast<int>(MyDist::bnldev(decay_rate,pos->target_volume));
    }
  }
}

void AtomicVesicles::decay(const int index)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "AtomicVesicles::get_birth_time() Error, index is out of range" << std::endl;
    exit(-1);
  }
  vesicle_info_v[index].target_volume -= static_cast<int>(MyDist::bnldev(Para::vesicle_decay_rate,vesicle_info_v[index].target_volume));
}

void AtomicVesicles::dump_eigen(const long Time) const
{
  std::ofstream fout("time_vol_eigen.dat",std::ios_base::app);
  for(AtomicVesicleInfoV::const_iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      fout << Time << ' ' << pos->get_volume() << ' ' << pos->eigen1 << '\n';
    }
  }

}

// I disable this feature for a moment for AtomicVesicles.
// long AtomicVesicles::get_birth_time(const int index) const
// {
//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::get_birth_time() Error, index is out of range" << std::endl;
//     exit(-1);
//   }

//   return vesicle_info_v[index].birth_time;
// }

// void AtomicVesicles::set_birth_time(const int index,const long Time)
// {
//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::get_birth_time() Error, index is out of range" << std::endl;
//     exit(-1);
//   }

//   vesicle_info_v[index].birth_time = Time;
// }

// I disable this feature for a moment for AtomicVesicles.
// #ifdef PLOT_VESICLE_INTERNAL_DYNAMICS
// void AtomicVesicles::set_identity(const int index, const std::deque<bool>& _identity)
// {
//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::set_identity() Error, index is out of range" << std::endl;
//     exit(-1);
//   }

//   vesicle_info_v[index].v_identity = _identity;
// }

// /* This should be called after a vesicle gives birth to a new
//    vesicle. This gives a new identity to a mother vesicle and an
//    identity to a daughter vesicle. */
// void AtomicVesicles::make_identity(const int mother_ind,const int daughter_ind)
// {
//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(mother_ind>-1&&mother_ind<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::give_identity() Error, mother_ind is out of range" << std::endl;
//     exit(-1);
//   }

//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(daughter_ind>-1&&daughter_ind<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::give_identity() Error, daughter_ind is out of range" << std::endl;
//     exit(-1);
//   }

//   vesicle_info_v[daughter_ind].v_identity = vesicle_info_v[mother_ind].v_identity;
//   vesicle_info_v[daughter_ind].v_identity.push_back(true);
//   vesicle_info_v[mother_ind].v_identity.push_back(false);
// }

// const std::deque<bool>& AtomicVesicles::get_identity(const int index) const
// {
//   try{
//     Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
//   }
//   catch(GeneralError){
//     std::cerr << "AtomicVesicles::give_identity() Error, mother_ind is out of range" << std::endl;
//     exit(-1);
//   }

//   return vesicle_info_v[index].v_identity;
// }
// #endif //PLOT_VESICLE_INTERNAL_DYNAMICS

