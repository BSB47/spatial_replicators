#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

#include "para.hpp"
#include "assert.hpp"

/* This class holds the statistics of replicator population in
   a vesicle. How this class should be defined depends on how
   molecules are represented in the model. */
struct VesicleStatistics {
  /* If at least one square in vesicle plane contains an index
     of this vesicle, then this vesicle is considered as
     existing (even if it's empty).  */
  bool is_existing;
  /* The simulation algorithm uses volume & target_volume
     stored in VesicleInfo (vesicles.hpp) for the actual
     vesicle dynamics. Information contained here are used
     only for the purpose of monitoring.*/
  int volume;
  int target_volume;

  /* Total number of molecules */
  int n_total;
  /* Number of RNA polymerase RNA. */
  int n_RNAp_RNA;
  /* Number of RNA polymerase DNA. */
  int n_RNAp_DNA;
  /* Number of DNA polymerase RNA. */
  int n_DNAp_RNA;
  /* Number of DNA polymerase DNA. */
  int n_DNAp_DNA;
  /* Number of parasites RNA */
  int n_para_RNA;
  /* Number of parasites DNA */
  int n_para_DNA;
  /* Number of Junk RNA */
  int n_junk_RNA;
  /* Number of Junk DNA */
  int n_junk_DNA;

  /* Default constructor is used for resetting all the statistics. */
  VesicleStatistics()
    : is_existing(false),
      volume(0),
      target_volume(0),
      n_total(0),
      n_RNAp_RNA(0),
      n_RNAp_DNA(0),
      n_DNAp_RNA(0),
      n_DNAp_DNA(0),
      n_para_RNA(0),
      n_para_DNA(0),
      n_junk_RNA(0),
      n_junk_DNA(0)
  {}
};

/* This class manages the statistics of vesicles.*/
class VesicleObserver {
  typedef std::vector<VesicleStatistics> VesicleStat;
  VesicleStat vesicle_stat;

  /* Population size of vesicles. Only potentially surviving vesicles
     are considered. What is potentially surviving depends on the
     assumption on the contiribution of molecules to the growth of
     target volume. */
  int population_size;
  /* Mean volume & target_volume. Only potentially surviving vesicles
     are considered. What is potentially surviving depends on the
     assumption on the contiribution of molecules to the growth of
     target volume. */
  double m_volume;
  double m_target_volume;
  double m_n_total;

  /* The mean of the characters of internal replicators of the
     vesicles. */
  /* nothing yet */

public:
  inline VesicleObserver(const int length);
  
  /* reset all the statistics for the new measurement */
  inline void reset();

  /* The method calculates the mean volume and mean target_volume, and
     it sets the population size of vesicles (population_size). */
  void calculate_mean();


  inline unsigned length(){return vesicle_stat.size();}

  inline const VesicleStatistics& get_statistics(const int index) const;

  inline void set_is_existing(const int index);
  inline bool get_is_existing(const int index) const;

  inline void set_volume(const int index,const int vol);
  inline int get_volume(const int index) const;

  inline void set_target_volume(const int index,const int vol);
  inline int get_target_volume(const int index) const;

  inline void add_n_total(const int index,const int incl=1);
  inline int get_n_total(const int index) const;

  inline void add_n_RNAp_RNA(const int index,const int incl=1);
  inline int get_n_RNAp_RNA(const int index) const;

  inline void add_n_RNAp_DNA(const int index,const int incl=1);
  inline int get_n_RNAp_DNA(const int index) const;

  inline void add_n_DNAp_RNA(const int index,const int incl=1);
  inline int get_n_DNAp_RNA(const int index) const;

  inline void add_n_DNAp_DNA(const int index,const int incl=1);
  inline int get_n_DNAp_DNA(const int index) const;

  inline void add_n_para_RNA(const int index,const int incl=1);
  inline int get_n_para_RNA(const int index) const;

  inline void add_n_para_DNA(const int index,const int incl=1);
  inline int get_n_para_DNA(const int index) const;

  inline void add_n_junk_RNA(const int index,const int incl=1);
  inline int get_n_junk_RNA(const int index) const;

  inline void add_n_junk_DNA(const int index,const int incl=1);
  inline int get_n_junk_DNA(const int index) const;

  double get_mean_volume() const {return m_volume;}
  double get_mean_target_volume() const {return m_target_volume;}
  double get_mean_n_total() const {return m_n_total;}
  int get_population_size() const {return population_size;}
};


VesicleObserver::VesicleObserver(const int length)
  : vesicle_stat(length)
{}

void VesicleObserver::reset()
{
  /* Default constructor of VesicleStatistics has all the data members
     reset. */
  const VesicleStatistics zero;
  for(VesicleStat::iterator iter=vesicle_stat.begin();iter!=vesicle_stat.end();iter++){
    (*iter) = zero;
  }
}

const VesicleStatistics& VesicleObserver::get_statistics(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_birth_time() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index];
}

void VesicleObserver::set_is_existing(const int index)
{
  if(ASSERT::ERROR_CHECK){
    try{
      vesicle_stat.at(index).is_existing = true;
    }catch(std::out_of_range){
      std::cerr << "VesicleObserver::set_is_existing() out of range error has occured." << std::endl;
      exit(-1);
    }
  }else{
    vesicle_stat[index].is_existing = true;
  }
}

bool VesicleObserver::get_is_existing(const int index) const
{
  if(ASSERT::ERROR_CHECK){
    bool result;
    try{
      result = vesicle_stat.at(index).is_existing;
    }catch(std::out_of_range){
      std::cerr << "VesicleObserver::set_is_existing() out of range error has occured." << std::endl;
      exit(-1);
    }
    return result;
  }else{
    return vesicle_stat[index].is_existing;
  }
}

void VesicleObserver::set_volume(const int index,const int vol)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::set_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }
  vesicle_stat[index].volume = vol;
}

int VesicleObserver::get_volume(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::set_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }
  return vesicle_stat[index].volume;
}

void VesicleObserver::set_target_volume(const int index,const int vol)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::set_target_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }
  vesicle_stat[index].target_volume = vol;
}

int VesicleObserver::get_target_volume(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::set_target_volume() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].target_volume;
}

void VesicleObserver::add_n_RNAp_RNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_RNAp_RNA += incl;
}

int VesicleObserver::get_n_RNAp_RNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_RNAp_RNA;
}

void VesicleObserver::add_n_RNAp_DNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_RNAp_DNA += incl;
}

int VesicleObserver::get_n_RNAp_DNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_RNAp_DNA;
}

void VesicleObserver::add_n_DNAp_RNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_DNAp_RNA += incl;
}

int VesicleObserver::get_n_DNAp_RNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_DNAp_RNA;
}

void VesicleObserver::add_n_DNAp_DNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_DNAp_DNA += incl;
}

int VesicleObserver::get_n_DNAp_DNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_DNAp_DNA;
}

void VesicleObserver::add_n_para_RNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_para_RNA += incl;
}

int VesicleObserver::get_n_para_RNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_para_RNA;
}

void VesicleObserver::add_n_para_DNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_para_DNA += incl;
}

int VesicleObserver::get_n_para_DNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_para_DNA;
}

void VesicleObserver::add_n_junk_RNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_junk_RNA += incl;
}

int VesicleObserver::get_n_junk_RNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_junk_RNA;
}

void VesicleObserver::add_n_junk_DNA(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_junk_DNA += incl;
}

int VesicleObserver::get_n_junk_DNA(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_simple_replicase() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_junk_DNA;
}

void VesicleObserver::add_n_total(const int index,const int incl)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::add_n_total() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  vesicle_stat[index].n_total += incl;
}

int VesicleObserver::get_n_total(const int index) const
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(index>-1&&static_cast<unsigned>(index)<vesicle_stat.size()));
  }
  catch(GeneralError){
    std::cerr << "VesicleObserver::get_n_total() Error, index is out of range" << std::endl;
    exit(-1);
  }
  
  return vesicle_stat[index].n_total;
}
