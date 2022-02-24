#include "vesicles.hpp"

SerialVesicles::SerialVesicles(int a_length)
  : length(a_length),
    vesicle_info_v(a_length)
{}

void SerialVesicles::calculate_PCA_for_dividing_vesicles()
{
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    /* this is just-in-case error check */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||pos->get_volume()>-2);
    }
    catch(GeneralError){
      std::cerr << "SerialVesicles::calculate_PCA_for_dividing_vesicles() Error, volume has an abnormal value=" << pos->get_volume() << std::endl;
      exit(-1);
    }

    if(pos->get_volume()!=-1){
      if(pos->is_divide){
	pos->calculate_PCA();
      }
    }
  }
}

/* This function is not used often, so I include error cheking. */
void SerialVesicles::set_target_volume(const int index,const int val)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_target_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].target_volume = val;
}

/* This is used by MyCA::vesicle_neutral_growth. */
void SerialVesicles::set_all_target_volume_to_0()
{
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      pos->target_volume = 0;
    }
  }
}

/* This is used by MyCA::vesicle_neutral_growth. */
void SerialVesicles::neutral_growth()
{
  const double factor = Para::neutral_growth_expand_factor;
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      pos->target_volume = static_cast<int>(factor * pos->target_volume + 0.5);
    }
  }
}

void SerialVesicles::destroy_vesicle(const int index)
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

bool SerialVesicles::is_divide(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::is_divide() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::is_divide() Error, try to get is_divide of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  return vesicle_info_v[index].is_divide;
}

void SerialVesicles::set_is_divide(const int index)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_is_divide() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_is_divide() Error, try to get set_is_divide of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].is_divide = true;
}

bool SerialVesicles::is_daughter(const int index, const unsigned row, const unsigned col) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::is_daughter() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::is_daughter() Error, try to use is_daughter of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  return vesicle_info_v[index].is_daughter(row,col);
}

int SerialVesicles::get_daughter_ind(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::get_daughter_ind() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::get_daughter_ind() Error, try to get daughter index of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  return vesicle_info_v[index].daughter_ind;
}


void SerialVesicles::set_daughter_ind(const int index,const int daughter_ind)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_daughter_ind() Error, index is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||vesicle_info_v[index].get_volume()>-1);
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_daughter_ind() Error, try to set daughter index of a vesicle whose index is not used." << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(daughter_ind>-1&&daughter_ind<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_daughter_ind() Error, daughter_ind is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].daughter_ind = daughter_ind;
}

int SerialVesicles::get_unused_index() const
{
  for(VesicleInfoV::const_iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()==-1){
      return (pos-vesicle_info_v.begin());
    }
  }
  
  std::cerr << "SerialVesicles::get_unused() We run out vesicle index" << std::endl;
  exit(-1);
}

int SerialVesicles::make_new_vesicle()
{
  const int index = get_unused_index();
  vesicle_info_v[index].set_volume(0);
  vesicle_info_v[index].is_divide=false;
  
  return index;
}

void SerialVesicles::make_new_vesicle(const int index)
{
  if(index >= length){
    std::cerr << "SerialVesicles::make_new_vesicle(int) Error, index is out of range" << std::endl;
    exit(-1);
  }

  if(vesicle_info_v[index].get_volume()!=-1){
    std::cerr << "SerialVesicles::make_new_vesicle(int) Error, the vesicle with this index has already been made." << std::endl;
    exit(-1);
  }
  
  vesicle_info_v[index].set_volume(0);
  vesicle_info_v[index].is_divide=false;
}


void SerialVesicles::PCA_reset()
{
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    /* this is just-in-case error check */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||pos->get_volume()>-2);
    }
    catch(GeneralError){
      std::cerr << "SerialVesicles::PCA_reset() Error, volume has an abnormal value=" << pos->get_volume() << std::endl;
      exit(-1);
    }

    /* If this vesicle exists, initilize PCA variables & set is_divide
       flag */
    if(pos->get_volume() != -1){
      pos->eigen1 = pos->sum_x = pos->sum_y = pos->sum_xx = pos->sum_yy = pos->sum_xy = 0;
      pos->is_divide = false;
      pos->daughter_ind = -1;
    }
  }
}

void SerialVesicles::PCA_collect_data(int index,unsigned row,unsigned col)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1 && index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::PCA_collect_data() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].sum_x += col;
  vesicle_info_v[index].sum_y += row;
  vesicle_info_v[index].sum_xx += col*col;
  vesicle_info_v[index].sum_yy += row*row;
  vesicle_info_v[index].sum_xy += row*col;
}

bool SerialVesicles::PCA_calculate()
{
  bool there_is_no_division=true;
  /* m is mean; v is variance*/
  double m_x,m_y,v_xx,v_yy,v_xy,vol_inver,d;
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    /* this is just-in-case error check */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||pos->get_volume()>-2);
    }
    catch(GeneralError){
      std::cerr << "SerialVesicles::PCA_calculate() Error, volume has an abnormal value=" << pos->get_volume() << std::endl;
      exit(-1);
    }

    /* if volume==1, we cannot calculate PCA. */
    if(pos->get_volume()>1){
      m_x = pos->sum_x/static_cast<double>(pos->get_volume());
      m_y = pos->sum_y/static_cast<double>(pos->get_volume());
      vol_inver = 1./(pos->get_volume()-1);
      v_xx = (pos->sum_xx - m_x*pos->sum_x)*vol_inver;
      v_yy = (pos->sum_yy - m_y*pos->sum_y)*vol_inver;
      v_xy = (pos->sum_xy - m_x*pos->sum_y)*vol_inver;
    
      /* According to my old code for vesicles, which is according to
	 Biometry (by Sokal):
       
	 greatest eigen value: eigen_1 = (v_xx + v_yy + d)/2 smallest
	 eigen value: eigen_2 = (v_xx + v_yy - d)/2
	 d=sqrt((v_xx+v_yy)*(v_xx+v_yy)-4*(v_xx*v_yy-v_xy*v_xy))
       
	 Then, for PC1, Y=beta1 X + alpha1, beta1 = v_xy/(eigen_1-v_yy),
	 alpha1 = m_y-beta1*m_x
       
	 For PC2, beta2 = -1/beta1, alpha2 = m_y - beta2/m_x
       
	 From these, I can obtain, for PC2,
       
	 L*X + v_xy*Y - m_y*v_xy - m_x*L = 0, where L = eigen_1 - v_yy =
	 (v_xx - v_yy + d)/2

	 Note that by chance, I took x for row and y for column, which
	 is different from an intuitive choice of name.
      */
      d = sqrt((v_xx-v_yy)*(v_xx-v_yy)+4*v_xy*v_xy);

      /* if v_xy==0, there is no correlation between x and y. This means
	 that PC are parallel to either x-axis or y-axis; i.e. the slope
	 of PC is either infinite or zero. Thus, we need special
	 care. If v_xy!=0, then the slope of PC is finite and not
	 zero. Thus, we don't need special care */
      if(fabs(v_xy) > std::numeric_limits<double>::epsilon()){
	pos->a = 0.5*(v_xx-v_yy+d);
	pos->b = v_xy;
	pos->c = -m_y*v_xy - m_x*pos->a;
      }else{
	/* if the cell is elongated in x direction (i.e. v_xx > v_yy),
	   then the 2nd PC is along y-axis. If the cell is elongated in
	   y direction (i.e. v_yy > v_xx), then the 2nd PC is along
	   x-axis. In the following, I include v_xx==v_yy case in
	   v_xx>v_yy. This means that a completely symmetric vesicle in
	   any direction (i.e. circle) will be divided along y-axis. */
	if(v_xx >= v_yy){
 	  pos->a = 1.;
 	  pos->b = 0.;
 	  pos->c = -m_x;//m_x is never zero
	}else{
 	  pos->a = 0.;
 	  pos->b = 1.;
 	  pos->c = -m_y;
	}
      }
      //      std::cout << "a=" << pos->a << " b=" << pos->b << " c=" <<
      //      pos->c << std::endl;

      pos->eigen1 = pos->a+v_yy;
      pos->eigen2 = pos->eigen1 - d;

      if(pos->get_volume()>=Para::volume_threshold){
	pos->is_divide = true;
	pos->daughter_ind = -1;
	if(there_is_no_division)
	  there_is_no_division = false;
      }
      /* I also use eigen values to divide vesicles that have an extreme
	 shape. In particular, when vesicles are forced to be split in a
	 few parts, eigen values tend to be a very large value. I set a
	 threshold in the largest eigen value. If the eigen value
	 exceeds the threshold, I let a vesicle divide. The threshold
	 value is arbitrary determined from the emperical distribution
	 of eigen values whem such a threhsold does not exist. */
      /* Later, I realzed that this part of code causes more troubles
	 than the troubles it solves. Thus, I do not use this feature
	 anymore.  else
	 if(pos->eigen1>pos->get_volume()*0.112704+0.2*Para::volume_threshold){
	 pos->is_divide = true; pos->daughter_ind = -1;
	 if(there_is_no_division) there_is_no_division = false; }
      */
      else{
	pos->is_divide = false;
	/* Note that those vesicles that have is_divide==false do not
	   have to initilize daughter_ind */
      }
    }
    /* This is the case if volume is either 0 or 1. */
    else if(pos->get_volume()>-1){
      pos->is_divide = false;
    }
  }

  return !there_is_no_division;
}

void SerialVesicles::divide_target_volume()
{
  int target;
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      if(pos->is_divide){
	try{
	  Assert<GeneralError>(!ASSERT::ERROR_CHECK||(pos->get_volume()>0 && vesicle_info_v[pos->daughter_ind].get_volume() > 0));
	}
	catch(GeneralError){
	  std::cerr << "SerialVesicles::divide_target_volume() Error, volume of the mother vesicle=" << pos->get_volume() << " that of the daughter=" << vesicle_info_v[pos->daughter_ind].get_volume() << std::endl;
	  exit(-1);
	}
	target=pos->target_volume;
	pos->target_volume = static_cast<int>(static_cast<double>(target*pos->get_volume())/(pos->get_volume()+vesicle_info_v[pos->daughter_ind].get_volume()) + 0.5);
	vesicle_info_v[pos->daughter_ind].target_volume = target - pos->target_volume;
      }
    }
  }
}

void SerialVesicles::decay()
{
  const double decay_rate = Para::vesicle_decay_rate;
  for(VesicleInfoV::iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      pos->target_volume -= static_cast<int>(MyDist::bnldev(decay_rate,pos->target_volume));
    }
  }
}

void SerialVesicles::decay(const int index)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::get_birth_time() Error, index is out of range" << std::endl;
    exit(-1);
  }
  vesicle_info_v[index].target_volume -= static_cast<int>(MyDist::bnldev(Para::vesicle_decay_rate,vesicle_info_v[index].target_volume));
}

void SerialVesicles::dump_eigen(const long Time) const
{
  std::ofstream fout("time_vol_eigen.dat",std::ios_base::app);
  for(VesicleInfoV::const_iterator pos=vesicle_info_v.begin();pos!=vesicle_info_v.end();++pos){
    if(pos->get_volume()!=-1){
      fout << Time << ' ' << pos->get_volume() << ' ' << pos->eigen1 << '\n';
    }
  }

}

long SerialVesicles::get_birth_time(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::get_birth_time() Error, index is out of range" << std::endl;
    exit(-1);
  }

  return vesicle_info_v[index].birth_time;
}

void SerialVesicles::set_birth_time(const int index,const long Time)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::get_birth_time() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].birth_time = Time;
}

#ifdef PLOT_VESICLE_INTERNAL_DYNAMICS
void SerialVesicles::set_identity(const int index, const std::deque<bool>& _identity)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::set_identity() Error, index is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[index].v_identity = _identity;
}

/* This should be called after a vesicle gives birth to a new
   vesicle. This gives a new identity to a mother vesicle and an
   identity to a daughter vesicle. */
void SerialVesicles::make_identity(const int mother_ind,const int daughter_ind)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(mother_ind>-1&&mother_ind<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::give_identity() Error, mother_ind is out of range" << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(daughter_ind>-1&&daughter_ind<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::give_identity() Error, daughter_ind is out of range" << std::endl;
    exit(-1);
  }

  vesicle_info_v[daughter_ind].v_identity = vesicle_info_v[mother_ind].v_identity;
  vesicle_info_v[daughter_ind].v_identity.push_back(true);
  vesicle_info_v[mother_ind].v_identity.push_back(false);
}

const std::deque<bool>& SerialVesicles::get_identity(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&index<length));
  }
  catch(GeneralError){
    std::cerr << "SerialVesicles::give_identity() Error, mother_ind is out of range" << std::endl;
    exit(-1);
  }

  return vesicle_info_v[index].v_identity;
}
#endif //PLOT_VESICLE_INTERNAL_DYNAMICS

