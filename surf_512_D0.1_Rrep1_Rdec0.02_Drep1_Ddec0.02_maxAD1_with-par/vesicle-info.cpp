#include "vesicle-info.hpp"

void VesicleInfo::calculate_PCA()
{
  
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||volume>1);
  }
  catch(GeneralError){
    std::cerr << "VesicleInfo::calculate_PCA() Error, PCA cannot be calculated if volume<1" << std::endl;
    exit(-1);
  }
  
  double m_x,m_y,v_xx,v_yy,v_xy,vol_inver,d;
  m_x = sum_x/static_cast<double>(volume);
  m_y = sum_y/static_cast<double>(volume);
  vol_inver = 1./(volume-1);
  v_xx = (sum_xx - m_x*sum_x)*vol_inver;
  v_yy = (sum_yy - m_y*sum_y)*vol_inver;
  v_xy = (sum_xy - m_x*sum_y)*vol_inver;
    
  /* According to my old code for vesicles, which is according to
     Biometry (by Sokal):
       
     greatest eigen value: eigen_1 = (v_xx + v_yy + d)/2
     smallest eigen value: eigen_2 = (v_xx + v_yy - d)/2
     d=sqrt((v_xx+v_yy)*(v_xx+v_yy)-4*(v_xx*v_yy-v_xy*v_xy))
       
     Then, for PC1, Y=beta1 X + alpha1,
     beta1 = v_xy/(eigen_1-v_yy),
     alpha1 = m_y-beta1*m_x
       
     For PC2,
     beta2 = -1/beta1,
     alpha2 = m_y - beta2/m_x
       
     From these, I can obtain, for PC2,
       
     L*X + v_xy*Y - m_y*v_xy - m_x*L = 0,
     where L = eigen_1 - v_yy
     = (v_xx - v_yy + d)/2

     Note that by chance, I took x for row and y for column, which is
     different from an intuitive choice of name.
  */
  d = sqrt((v_xx-v_yy)*(v_xx-v_yy)+4*v_xy*v_xy);

  /* If v_xy==0, there is no correlation between x and y. This means
     that PC are parallel to either x-axis or y-axis; i.e. the slope of
     PC is either infinite or zero. Thus, we need special care. If
     v_xy!=0, then the slope of PC is finite and not zero. Thus, we
     don't need special care */
  if(fabs(v_xy) > std::numeric_limits<double>::epsilon()){
    a = 0.5*(v_xx-v_yy+d);
    b = v_xy;
    c = -m_y*v_xy - m_x*a;
  }else{
    /* If the cell is elongated in x direction (i.e. v_xx > v_yy), then
       the 2nd PC is along y-axis. If the cell is elongated in y
       direction (i.e. v_yy > v_xx), then the 2nd PC is along x-axis. In
       the following, I include v_xx==v_yy case in v_xx>v_yy. This means
       that a completely symmetric vesicle in any direction
       (i.e. circle) will be divided along y-axis. */
    if(v_xx >= v_yy){
      a = 1.;
      b = 0.;
      c = -m_x;//m_x is never zero
    }else{
      a = 0.;
      b = 1.;
      c = -m_y;
    }
  }
  //      std::cout << "a=" << a << " b=" << b << " c=" << c << std::endl;
  eigen1 = a+v_yy;
  eigen2 = eigen1 - d;
}

