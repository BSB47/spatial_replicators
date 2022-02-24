#include "vesicle-stat.hpp"

void VesicleObserver::calculate_mean()
{
  /* Initilization */
  /* vesicle population */
  unsigned local_population_size = 0;
  unsigned cum_volume = 0;
  unsigned cum_target_volume = 0;
  unsigned cum_n_total = 0;

  for(VesicleStat::iterator iter=vesicle_stat.begin();iter!=vesicle_stat.end();iter++){
    /* Population size of vesicle & volume */
    if(iter->is_existing){
      /* We count only those vesicles that contain at least one
	 molecule. */
      if(iter->n_total>0){
	++local_population_size;
	cum_volume += iter->volume;
	cum_target_volume += iter->target_volume;
	cum_n_total += iter->n_total;
      }
    }
  }

  if(local_population_size>0){
    /* Population size of vesicle & volume */
    double inv_population_size = 1./local_population_size;
    population_size = local_population_size;
    m_volume = static_cast<double>(cum_volume)*inv_population_size;
    m_target_volume = static_cast<double>(cum_target_volume)*inv_population_size;
    m_n_total = static_cast<double>(cum_n_total)*inv_population_size;
  }
  else{
    population_size = 0;
    m_volume = 0.;
    m_target_volume = 0.;
    m_n_total = 0.;
  }
}
