/* Library */
#include <vector>
#include <algorithm>
#include <functional>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <sstream>
#include <deque>
#include <unistd.h>
#include <cmath>

#include "tbb/tbb.h"

/* Generic */
#include "cellular-automata.hpp"
#include "assert.hpp"
#include "cash-display.hpp"
#include "random_wrapper.hpp"

/* Specific */
#include "para.hpp"
#include "automaton.hpp"
#include "vesicles.hpp"
#include "vesicle-stat.hpp"
#include "upd_order.hpp"
#include "atomic-vesicles.hpp"

#ifndef MY_CA
#define MY_CA

/*
  The general problem is that myca.cpp has two slightly different
  models, one is serial version, and the other is parallel version. The
  general question is how to define and use the members (especially
  those that are a class) of MyCA for these two versions. 

  One possibility is to use one and the same class for the two versions,
  and use it differently. This is the case of CA2D<Automaton>
  plane. Automaton contains a member which plays a role of lock (namely,
  tbb:atomic<int> is_busy). In the parallel version, every access to a
  CA square must first check if the square is_busy, so that only one
  thread accesses to the square at the same time. In the serial version,
  this process is skiped.

  UpdOrderShuffler falls in the same category as above. From the problem
  of concurrency piont of view, UpdOrderShufler doesn't actually have to
  be specially arranged for the parallel version because there isn't
  (and should not be) any situations where multiple threads try to
  access the same element of the upd_order_array. However, I still
  implement UpdOrderShuffler because shuffling the upd_order_array
  itself can take advantage of parallelization. In the parallel version,
  the simulation algorithm in MyCA uses UpdOrdershuffler through
  get_upd_order_parallel(), and shuffling is done parallely. So, this is
  a case where a care for parallelization is not necessary, but it is
  advantageous to take it into account.

  Another possibility is to use two different classes for the two
  versions. This will be the case for ParallelVesicles and
  AtomicVesicles. The reason for doing this is that I need an almost
  complete set of members of VesicleInfo to be atomic (free of
  concurrency problems). Luckly, locking the element of the VesicleInfor
  array isn't necessary because access to an element involves only
  either reading operation or incrementing operation (except for one
  case), and it doesn't involve multiple operations that depend on the
  result of previous operations (in CA2D<Automaton>, this is mostly the
  case; e.g., if the cell is empty, change the sate of the cell in
  specific way, which has a danger that another thread might make the
  cell non-empty between the current thread's reading and writing
  operations). One exception is when we want to create a new "alive"
  element in the VesicleInfo array (i.e. an element whose index is---or
  will be---the CPM state of any square in the system). This involves
  multiple operations that depends on the result of the result of the
  previous operations; namely, we first have to find the array element
  that is "dead" (i.e. the index of that element is not---and will not
  be---used by any CPM squares), and then convert it to "alive"; but
  another thread might convert it to living before the first thread
  does. So, in this case, we have to lock the elements. But the lucky
  thing is that we can convert these multiple operations into one atomic
  operation, i.e. compare_and_swap(). The value of member variable
  "volume" designates whether the element is alive or dead, and we can
  simply continue to use the same method of designation.


*/


class MyCA{
private:
  /************** Simulation *****************/
  /* System size */
  const unsigned nrow;
  const unsigned ncol;

  /* Total number of cells (nrow*ncol) */
  const unsigned length;
  /* "Meat" of the class. It contains automata */
  CA2D<Automaton> plane;

  /* This takes care of vesicle related attributes. See vesicles.hpp for
     details. */
  SerialVesicles serial_vesicles;
  AtomicVesicles atomic_vesicles;

  /* To update vesicles, one can use random sequencial updating. */
  UpdOrderShuffler upd_order_shuffler;

  /* An array that holds probabilities as a function of energy
     differencial. */
  double probabilities[Para::PROBABILITY_ARRAY_SIZE];

  /* This holds information about protocells (e.g. number of molecules
     in a vesicle). */
  VesicleObserver vesicle_observer;

  /************** Visualization *****************/
  /* for visualization */
  enum DisplayPanel {FUNCTION,PANEL_FOLD,HIST_RNAp_RNArecog,HIST_RNAp_DNArecog,HIST_RNAp_fold,HIST_DNAp_RNArecog,HIST_DNAp_DNArecog,HIST_DNAp_fold};
  CashDisplay* display_p;
  std::vector<unsigned char> color_of_gradient;
  enum MyCAColor {COLOR_BLACK,COLOR_WHITE,COLOR_CELL_INSIDE,COLOR_CELL_WALL,COLOR_RNAp,COLOR_DNAp,COLOR_PARAS,COLOR_JUNK,COLOR_HIST_BARS,COLOR_HIST_BACK,COLOR_MARGIN};

  /*************** Output ****************/
  std::ofstream plot_fout;

public:
  MyCA(const unsigned a_nrow,const unsigned a_ncol);
  ~MyCA();
  
  /*** INITIALIZATION ***/
  void initialize(const long Time);
  
  /* Some miscerenious manipulation of initial conditions */
  void get_rid_of_some_species();

  /*** VISUALIZATION ***/
  void visualize(const long t);
  void mouse_event(void (*button2)(int x,int y),void (*button3)(int x,int y)){display_p->handle_mouse_event(button2,button3);}
  void button2(int x,int y);
  void button3(int x,int y);
  /* This is used when the initial time != 0. */
  void reset_movie_frame(const long init_time);

  /*** SIMULATION OF VESICLES ***/

  /* Serial versions */
  inline void update_vesicle_whole();
  void divide_vesicle(const long Time);
  /* Vesicle, neutral case: target vol is a*n_mol. */
  void vesicle_neutral_growth();
  /* Vesicle, lipid synthetase case: target vol is inclemented when
     lipid synthesizing reaction happens, and it linearly decays. */
  void vesicle_decay(){serial_vesicles.decay();}

  /* Parallelized version */
  inline void parallel_update_vesicle_whole();
  /* Vesicle, neutral case: target vol is a*n_mol. */
  void parallel_divide_vesicle(const long Time);
  inline void parallel_vesicle_neutral_growth();
  /* The target volume of vesicles is set to Para::fixed_target_volume
     if a vesicle contains at least one replicator. Otherwise, it sets
     the target volume to zero. */
  inline void parallel_vesicle_fixed_growth();

  /*** SIMULATION OF REPLICATORS ***/
  void update_replicator_whole_EQA();
  void update_replicator_whole_no_complex_EQA();

  /* Parallelized version */
  inline void parallel_update_replicator_whole_no_complex_EQA();
  inline void parallel_update_replicator_whole_EQA();

  /*** OUTPUT ***/
  bool output(const long Time);
  /* The usage of .sav files. I simply dump the values of the parameters
     of Molecule in these files, so that the effect of rate_scale has
     already been taken into account. This means that one cannot simply
     use different rate_scale from that having used used in the
     simulation that produced file.sav. If the modification of
     rate_scale is necessary, we have to re-scale the rates when reading
     in the file, which hasn't been implemented yet. */
  void save(const char* fname) const;
  void save_wrapper(const long Time) const;

  /*** INPUT ***/
  void read(const std::string& file_name);
  void read_only_vesicle(const std::string &fname);
  void read_half_half(const std::string &fname1, const std::string &fname2);

private:
  /*** VISUALIZATION ***/
  void scan_plane_to_display();
  unsigned get_bin(const unsigned n_bins,const double min,const double max,const double value); 

  /*** SIMULATION OF VESICLES ***/
  /* This is for serial version */
  void update_vesicle_1square(const unsigned row,const unsigned col);

  /* This is for parallel version */
  void parallel_update_vesicle_1square(const unsigned row,const unsigned col);

  /* This method can be used by both serial and parallel versions */
  inline void energy_J(int& energy,const int index1,const int index2) const;

  /* This method is for serial version. */
  inline bool is_copy_state(const int energy) const;

  /* This method is for parallel version. It is distinguished by the
     function parameters. Sorry for the confusing way of coding. */
  inline bool is_copy_state(const int energy,RandomTBB::reference random_wrapper) const;

  /* This nested class is used by parallel_update_vesicle_whole() in
     tbb::parallel_for(). It simply calls
     parallel_update_vesicle_1square() with a coordinate read from
     upd_order. parallel_update_vesicle_1square() does the actual
     updating in such a way that it is concurrently safe. */
  class ParallelUpdateVesicle1Square {
    MyCA* const my_myca;
    const UpdOrder& upd_order;
  public:
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
    ParallelUpdateVesicle1Square(MyCA* _my_myca, const UpdOrder& _upd_order)
      : my_myca(_my_myca), upd_order(_upd_order){}
  };

  /* This nested class is used by parallel_divide_vesicle(), which calls
     tbb::parallel_for(). It looks at one CA square and records its
     coordinate information in atomic_vesicles for the later use if the
     vesicle in that CA square has a volume that is greater than the
     division threshold. */
  class ParallelDivideCollectData {
    MyCA& my_myca;
    tbb::atomic<bool>& my_is_divide;
  public:
    void operator()(const tbb::blocked_range2d<size_t>& r) const;
    ParallelDivideCollectData(MyCA& _my_myca,tbb::atomic<bool>& _my_is_divide)
      : my_myca(_my_myca), my_is_divide(_my_is_divide){}
  };

  /* This nested class is used by parallel_divide_vesicle(), which calls
     tbb::parallel_for(). It looks at one CA square and, if it belongs
     to a new daughter vesicle, it assigns a new CPM state. It also
     dissociate complex molecules that are spaning between a mother and
     daughter vesicle. */
  class ParallelDivideDivideVes {
    MyCA& my_myca;
  public:
    void operator()(const tbb::blocked_range2d<size_t>& r) const;
    ParallelDivideDivideVes(MyCA& _my_myca)
      : my_myca(_my_myca){}
  };

  /* This nested class is used by parallel_vesicle_neutral_growth().  It
     counts the number of replicators in each vesicle and its the target
     volume to the number of replicators times some factor. */
  class ParallelVesicleNeutralGrowthCount {
    MyCA& my_myca;
  public:
    void operator()(const tbb::blocked_range2d<size_t>& r) const;
    ParallelVesicleNeutralGrowthCount(MyCA& _my_myca)
      : my_myca(_my_myca){}
  };

  /* This nested class is used by parallel_vesicle_fixed_growth(). It
     counts the number of replicators in each vesicle and, if the
     number of replicators is more than zero, it sets the vesicle's
     target volume to Para::fixed_target_volume. */
  class ParallelVesicleFixedGrowth {
    MyCA& my_myca;
  public:
    void operator()(const tbb::blocked_range2d<size_t>& r) const;
    ParallelVesicleFixedGrowth(MyCA& _my_myca)
      : my_myca(_my_myca){}
  };

  /*** SIMULATION OF REPLICATORS ***/
  /* Common methods with complex formation */
  unsigned neighbor_converter(const unsigned nei) const {const unsigned result[9] = {0,4,3,2,1,8,7,6,5};  return result[nei];}
  void update_replicator_1square_empty(const unsigned row,const unsigned col);
  void update_replicator_1square_complex(const unsigned row,const unsigned col);
  /* Under equilibrium assumption with complex formation */
  inline void update_replicator_1square_EQA(const unsigned row,const unsigned col);
  void update_replicator_1square_simple_EQA(const unsigned row,const unsigned col);

  /* Without complex formation under equilibrium assumption. The rate of
     replication is cauclated in Molecule::rate_repli_no_complex() */
  void update_replicator_1square_empty_no_complex_EQA(const unsigned row,const unsigned col);
  void update_replicator_1square_simple_no_complex_EQA(const unsigned row,const unsigned col);

  /* PARALLEL_SURFACE_EQUIL_NO_COMPLEX */
  void parallel_update_replicator_1square_no_complex_EQA();
  bool parallel_update_replicator_1square_empty_no_complex_EQA(const unsigned row,const unsigned col);
  bool parallel_update_replicator_1square_simple_no_complex_EQA(const unsigned row,const unsigned col);
  class ParallelUpdateReplicator1SquareNoComplexEQA {
    MyCA* const  my_myca;
  public:
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
    ParallelUpdateReplicator1SquareNoComplexEQA(MyCA* myca_p) : my_myca(myca_p) {}
  };

  /* PARALLEL_SURFACE_EQUIL */
  void parallel_update_replicator_1square_EQA();
  bool parallel_update_replicator_1square_empty(const unsigned row,const unsigned col);
  bool parallel_update_replicator_1square_simple_EQA_ugly(const unsigned row,const unsigned col);
  bool parallel_update_replicator_1square_simple_EQA(const unsigned row,const unsigned col);
  bool parallel_update_replicator_1square_complex(const unsigned row,const unsigned col);
  class ParallelUpdateReplicator1SquareEQA {
    MyCA* const  my_myca;
  public:
    inline void operator()(const tbb::blocked_range<size_t>& r) const;
    ParallelUpdateReplicator1SquareEQA(MyCA* myca_p) : my_myca(myca_p) {}
  };

  /*** OUTPUT ***/
  void call_vesicle_observer(const unsigned row,const unsigned col);

  /*** INPUT ***/
  void read_plane_one_line(std::ifstream& fin,const int row, const int col,const int vesicle_index_offset=0);
  void read_only_vesicle_one_line(std::ifstream& fin,const int row,const int col,const int vesicle_index_offset=0);
  void set_color_half_half();
  void set_color_RNA_DNA();

  /*** MISCELLANEOUS ***/
  /* All methods that are not the part of the simulation algorithm such
     as those for initialization, input, output and observation do not
     need different implementation whether the simulation algorithm is
     parallelized or not (because they are always serial
     anyway). However, they must access different array of vesicle
     information, viz. either serial_vesicles or atomic_vesicles,
     depending on which simulation algorithm is used. The following
     methods offer an automatic choice. The member functions that use
     the following methods are contained in myca.cpp.

     In the simulation algorithm, I don't use these methods, but instead
     directly use serial_vesicles or atomic_vesicles. I thought this
     would make the difference in the algorithms more explicit. The
     member functions that implement the simulation algorithms are
     either in serial-myca.cpp or in parallel-myca.cpp.
  */
  inline int vesicle_get_volume(const int index) const;
  inline int vesicle_get_target_volume(const int index) const;
  inline void vesicle_set_target_volume(const int index,const int val);
  inline void vesicle_add_volume(const int index);
  inline void vesicle_make_new_vesicle(const int index);
  inline int vesicle_make_new_vesicle();

#ifdef _COUNT_RNA_DNA
public:
  void output_count_rna_dna(const long Time);
#endif //_COUNT_RNA_DNA

};


/* This function add energy associated to surfaces to the first argument
   "energy" depending on the second and third arguments that take the
   index of vesicles. We only use two vesicle types: vesicle and
   media. Media has index=-1; vesicles have index that is greater than
   or equal to zero. */
void MyCA::energy_J(int& energy,const int index1,const int index2) const
{
  if(index1 != index2){
    if(index1==-1 || index2 == -1){
      energy += Para::J_media_cell;
      return;
    }else{
      energy += Para::J_cell_cell;
      return;
    }
  }
}

bool MyCA::is_copy_state(const int energy) const
{
  if(energy<=0){
    return true;
  }
  else if(energy<Para::PROBABILITY_ARRAY_SIZE){
    return (rand_karney.Real() < probabilities[energy]);
  }
  else{
    return false;
  }
}

void MyCA::update_vesicle_whole()
{
  const UpdOrder& upd_order = upd_order_shuffler.get_upd_order_serial();
  for(UpdOrder::const_iterator pos=upd_order.begin();pos!=upd_order.end();++pos){
    update_vesicle_1square(pos->row,pos->col);
  }
}


////////////////////////////PARALLEL////////////////////////////////
bool MyCA::is_copy_state(const int energy,RandomTBB::reference random_wrapper) const
{
  if(energy<=0){
    return true;
  }
  else if(energy<Para::PROBABILITY_ARRAY_SIZE){
    return (random_wrapper.random.Real() < probabilities[energy]);
  }
  else{
    return false;
  }
}

void MyCA::ParallelUpdateVesicle1Square::operator()(const tbb::blocked_range<size_t>& r) const 
{
  MyCA* const myca_p = my_myca;
  for(size_t i=r.begin(), e=r.end();i!=e;++i){
    myca_p->parallel_update_vesicle_1square(upd_order[i].row,upd_order[i].col);
  }
}

void MyCA::parallel_update_vesicle_whole()
{
  const UpdOrder& upd_order = upd_order_shuffler.get_upd_order_parallel();
  tbb::parallel_for(tbb::blocked_range<size_t>(0,upd_order.size()),ParallelUpdateVesicle1Square(this,upd_order));
  atomic_vesicles.destroy_all_empty_vesicles();
}

void MyCA::ParallelUpdateReplicator1SquareNoComplexEQA::operator()(const tbb::blocked_range<size_t>& r) const
{
  MyCA* myca_p = my_myca;
  for(size_t i=r.begin();i!=r.end();++i){
    myca_p->parallel_update_replicator_1square_no_complex_EQA();
  }
}

void MyCA::parallel_update_replicator_whole_no_complex_EQA()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,Para::sys_nrow*Para::sys_ncol),ParallelUpdateReplicator1SquareNoComplexEQA(this));
}

void MyCA::ParallelUpdateReplicator1SquareEQA::operator()(const tbb::blocked_range<size_t>& r) const
{
  MyCA* myca_p = my_myca;
  for(size_t i=r.begin();i!=r.end();++i){
    myca_p->parallel_update_replicator_1square_EQA();
  }
}

void MyCA::parallel_update_replicator_whole_EQA()
{
  tbb::parallel_for(tbb::blocked_range<size_t>(0,Para::sys_nrow*Para::sys_ncol),ParallelUpdateReplicator1SquareEQA(this));
}

int MyCA::vesicle_get_volume(const int index) const
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM ||
     Para::model_type==Para::INORG_EQUIL){

    return serial_vesicles.get_volume(index);

  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    return atomic_vesicles.get_volume(index);

  }else{
    std::cerr << "MyCA::vesicle_get_volume() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

int MyCA::vesicle_get_target_volume(const int index) const
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM ||
     Para::model_type==Para::INORG_EQUIL){

    return serial_vesicles.get_target_volume(index);

  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    return atomic_vesicles.get_target_volume(index);

  }else{
    std::cerr << "MyCA::vesicle_get_target_volume() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

void MyCA::vesicle_set_target_volume(const int index,const int val)
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM ||
     Para::model_type==Para::INORG_EQUIL){

    serial_vesicles.set_target_volume(index,val);

  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    atomic_vesicles.set_target_volume(index,val);

  }else{
    std::cerr << "MyCA::vesicle_set_target_volume() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

void MyCA::vesicle_make_new_vesicle(const int index)
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM ||
     Para::model_type==Para::INORG_EQUIL){

     serial_vesicles.make_new_vesicle(index);

  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL || 
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    atomic_vesicles.make_new_vesicle(index);

  }else{
    std::cerr << "MyCA::vesicle_make_new_vesicle() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

int MyCA::vesicle_make_new_vesicle()
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM ||
     Para::model_type==Para::INORG_EQUIL){

    return serial_vesicles.make_new_vesicle();
    
  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL || 
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){
    
    return atomic_vesicles.make_new_vesicle();

  }else{
    std::cerr << "MyCA::vesicle_make_new_vesicle() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

void MyCA::vesicle_add_volume(const int index)
{
  if(Para::model_type == Para::VESIC_NEUT_EQUIL ||
     Para::model_type == Para::VESIC_LIPID_EQUIL ||
     Para::model_type == Para::SERIAL_ONLY_CPM||
     Para::model_type==Para::INORG_EQUIL){

    serial_vesicles.add_volume(index);

  }else if(Para::model_type == Para::PARALLEL_ONLY_CPM ||
	   Para::model_type == Para::PARALLEL_VESIC_NEUT_EQUIL ||
	   Para::model_type==Para::PARALLEL_INORG_EQUIL ||
	   Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    atomic_vesicles.add_volume(index);

  }else{
    std::cerr << "MyCA::vesicle_ad_volume() Error, tring to access vesicle info array although the model doesn't involve vesicles?" << std::endl;
    exit(-1);
  }
}

void MyCA::parallel_vesicle_neutral_growth()
{
  atomic_vesicles.set_all_target_volume_to_0();
  tbb::parallel_for(tbb::blocked_range2d<size_t>(1,nrow+1,1,ncol+1),ParallelVesicleNeutralGrowthCount(*this));
  atomic_vesicles.neutral_growth();
}

void MyCA::parallel_vesicle_fixed_growth()
{
  atomic_vesicles.set_all_target_volume_to_0();
  tbb::parallel_for(tbb::blocked_range2d<size_t>(1,nrow+1,1,ncol+1),ParallelVesicleFixedGrowth(*this));
}

#endif //MY_CA
