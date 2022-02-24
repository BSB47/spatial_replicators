#include <iostream>
#include <string>
#include <boost/test/floating_point_comparison.hpp>

/* Intel TBB */
#include "tbb/tbb.h"

#include "assert.hpp"

#include "para.hpp"

#ifndef AUTOMATON
#define AUTOMATON

class Molecule {
  /*
    SIMPL_TEMPL: a simple molecule that can be a template for
                 replication. It can form a complex with a
                 replicase. Under the equilibrium assumption, only
                 SIMPLE_TEMPL is used (SIMPL_FOLDE is not used).

    SIMPL_FOLDE: a simple molecule that can catalyze some reaction and
                 cannot be a template for replication. This is not used
                 if the system is under the equilibrium assumption.

    COMPL_ACTIV: a complexed molecule that is a potential polymerase.

    COMPL_PASSI: a complexed molecule that is a potential template.

    BOUND: boundary.
  */
  enum MolType {EMPTY,SIMPL_TEMPL,SIMPL_FOLDE,COMPL_ACTIV,COMPL_PASSI,BOUND};
  /* The following variable will be deep in the algorithm of
     reaction. */
  /* This variable represents the type of molecules. The algorithm of
     reaction will not directly access to this. */
  MolType mol_type;
public:
  /* If bon_nei==0, it's not making a complexed molecule. This variable
     will be directly accessed by the reaction algorith in myca.cpp */
  unsigned bon_nei;
  /* the follwoing is used only for well-mixed simulation */
  unsigned bon_row; unsigned bon_col;

private:
  /* The following variables will be used as a kind of temporary
     variables. These are used for a complexed molecule. The value of
     these depends on the partner molecule. */
  double disso_rate;
  /* This does not specify if the replication produces RNA or DNA. It is
     a sum of both rates. */
  double repli_rate;

  /* The following variables will specify the property of molecules. The
     meaning & existence of these variables depend on how we represent
     molecules. */

  /* Activity is used to determine the way in which the parameters are
     mutated. I assume JUNK molecules don't bind. */
  enum Activity {REPLI_RNA,REPLI_DNA,PARAS_RNA,PARAS_DNA,JUNK_RNA,JUNK_DNA};
  Activity activity;

  double decay_rate;

  /* This specifies the association rate of this molecule to replicase
     as template. Currently, it is set to max_assoc_diss_rate. For
     details, see Molecule::rate_assoc_disso_to(). */
  double bind_param;

  /* The next two are abstract parameters to determin
     association-dissociation rate of this replicase towards RNA or DNA
     template. They take value in [0,1] and multiplied to
     bind_param. For details, see Molecule::rate_assoc_disso_to(). */
  double rna_recog_param;
  double dna_recog_param;

  /* The next two determine the speed of RNA or DNA polymerization of
     this replicase. */
  double rna_repli_rate;
  double dna_repli_rate;

  /* These are the rate of folding/unfolding rates. If the system is
     under the equilibrium assumption (which is most likely the case in
     this study), fold_rate+unfold_rate=1, i.e. they represents the
     (averaged) fraction of time molecules are in the corresponding
     states. Note that I assume DNA is always available as template. */
  double fold_rate;
  double unfold_rate;

  /* This specifies the rate of "lipid synthesis", i.e. the rate of
     reaction that increases the target volume of a vesicle. */
  double lipid_synth_rate;

  /* This is used in several ways. See MyCA::initialize()

     1. When doing half & half competition experiments, I color
     molecules according to which file molecules came from. This will be
     activated if Para::read_fname2 is not empty. See
     MyCA::set_color_half_half().

     2. It is also used to do a simple ancestor tracing. I color
     molecules according to whether they are RNA or DNA templates. This
     will be activated if Para::read_fname2 is empty. See
     MyCA::set_color_RNA_DNA().
  */
  int color;

public:
  /* Generic methods that should exist irrespective of how molecules are
     represented (the detailes of the code do depend on the
     representation ) */
  /* Produce empty molecule */
  Molecule(): mol_type(EMPTY){} 
  /* produce boundary */
  Molecule(const char dummy): mol_type(BOUND){} 

  bool is_empty() const {return (mol_type==EMPTY);}
  bool is_simple() const {return (mol_type==SIMPL_TEMPL || mol_type==SIMPL_FOLDE);}
  bool is_complex() const {return (mol_type==COMPL_ACTIV || mol_type==COMPL_PASSI);}
  bool is_boundary() const {return (mol_type==BOUND);}
  bool is_folded() const {return (mol_type==SIMPL_FOLDE);}
  bool is_simple_templ() const {return (mol_type==SIMPL_TEMPL);}
  bool is_empty_or_boundary() const {return (mol_type==BOUND||mol_type==EMPTY);}

  /* Save & read the simulation */
  friend std::ostream& operator<<(std::ostream& s,const Molecule& mol);
  friend std::istream& operator>>(std::istream& s,Molecule& mol);

  /* Specific methods that depend on how molecules are represented */
  /* Return status */
  bool is_replicase_rna() const {return (activity==REPLI_RNA);}
  bool is_replicase_dna() const {return (activity==REPLI_DNA);}
  bool is_replicase() const {return (activity==REPLI_RNA||activity==REPLI_DNA);}
  bool is_parasite_rna() const {return (activity==PARAS_RNA);}
  bool is_parasite_dna() const {return (activity==PARAS_DNA);}
  bool is_parasite() const {return (activity==PARAS_RNA||activity==PARAS_DNA);}
  bool is_junk_rna() const {return (activity==JUNK_RNA);}
  bool is_junk_dna() const {return (activity==JUNK_DNA);}
  bool is_junk() const {return (activity==JUNK_RNA||activity==JUNK_DNA);}
  bool is_rna() const {return (activity==REPLI_RNA||activity==PARAS_RNA||activity==JUNK_RNA);}
  bool is_dna() const {return (activity==REPLI_DNA||activity==PARAS_DNA||activity==JUNK_DNA);}

  /*** With complex formation ***/
  /* Return a rate */
  double rate_repli() const {return repli_rate;}
  double rate_rna_repli() const {return rna_repli_rate;}
  double rate_dna_repli() const {return dna_repli_rate;}
  double rate_decay() const {/*std::cout << decay_rate << std::endl;*/return decay_rate;}
  void rate_assoc_disso_to(const Molecule& partner,double& assoc,double& disso) const;
  double rate_lipid_synth() const;
  double rate_disso() const {return disso_rate;}

  double get_rna_recog_param() const {return rna_recog_param;}
  void set_rna_recog_param(double x) {rna_recog_param=x;}

  double get_dna_recog_param() const {return dna_recog_param;}
  void set_dna_recog_param(double x) {dna_recog_param=x;}

  double rate_fold() const {return fold_rate;}
  void set_fold_rate(double x) {fold_rate = x; unfold_rate = Para::max_fold_unfold_rate-fold_rate;}
  
  /* Do some action (change state) */
  void initialize(const std::string& _activity,const double _fold_rate,const double _rna_recog_param=0.,const double _dna_recog_param=0.);
  void decay() {mol_type = EMPTY;}
  void replicate(Molecule& mol1,Molecule& mol2);
  void copy_with_mutation(const Molecule& original);
  void copy_with_mutation_new(const Molecule& original);
  inline void convert_to_rna();
  inline void convert_to_dna();
  inline void swap_with(Molecule& rhs);
  inline void rotate(Molecule& a,Molecule& b);
  inline void push(Molecule& a,Molecule& b);
  void binds_to(Molecule& partner,const double disso);
  void dissociate() {mol_type = SIMPL_TEMPL;}

  /*** Without complex formation ***/
  /* Return a rate */
  double rate_repli_no_complex(const Molecule& templ) const;
  
  /* Do some action (change state) */
  void replicate_no_complex(Molecule& templ,Molecule& repli);

  /* For debugging */
  int  get_mol_type() const {return static_cast<int>(mol_type);}
  bool is_complex_passive() const {return (mol_type==COMPL_PASSI);}

  /*** Visualization ***/
  void set_color(int _color) {color=_color;}
  int get_color() {return color;}

private:
  bool validate_moltype(const int _i);
  bool validate_activity(const int _i);
  /* Do some action */
  inline void mutate_parameter(double& parameter,const double min,const double max,const double step);
  inline void mutate_parameter_step_wrapped(double& parameter,const double min,const double max,const double step);
  static bool relative_close(double x,double y,double tolerance=1e-4);

#ifdef _COUNT_RNA_DNA
  /* These counts are inclemented according to whether this molecule
     is RNA or DNA template. */
  unsigned long _rna_count;
  unsigned long _dna_count;
  /* This is called in replicate(). */
  inline void _inclement_rna_dna_count();
  
public:
  void output_rna_dna_count(std::ostream& s);
#endif //_COUNT_RNA_DNA
};

/* About locking & is_busy, see also the explanation in myca.hpp */

struct Automaton {
  /* If is_busy=1, this automaton is being updated by a thread */
  tbb::atomic<int> is_busy;
  int vesicle_index;
  Molecule molecule; 
  /* According to TBB tutorial, I shouldn't do is_busy(0). */
  Automaton() : is_busy(), vesicle_index(-1){}

  /* Simple locking. It won't remember which thread locks a square. */
  bool lock(){return (is_busy.compare_and_swap(1,0)==0);}

  /* Simple releasing. It always releases a square (independent of
     which thread has locked the square). */
  void release(){is_busy = 0;}
};

/* Swap this and rhs */
void Molecule::swap_with(Molecule& rhs)
{
  const Molecule tmp(rhs);
  rhs = *this;
  *this = tmp;
}

/* This method is used in the diffusion between a simple (non-empty!)
   and complex molecule.

   Rotate (this,a,b) in the following order: this=a, a=b, b=this
*/
void Molecule::rotate(Molecule& a,Molecule& b)
{
  const Molecule tmp(*this);
  
  *this = a;
  a = b;
  b = tmp;
}
/* This method is used in the diffusion between an emty square and a
   complex molecule. "This" should be the empty square.

   Push (this,a,b) in the following way: this=a, a=b, b.decay(). */
void Molecule::push(Molecule& a,Molecule& b)
{
  *this = a;
  a = b;
  b.decay();
}

void Molecule::mutate_parameter(double& parameter,const double min,const double max,const double step)
{
  if(Para::switch_mutation_parameter_method==0){
    mutate_parameter_step_wrapped(parameter,min,max,step);
  }
  else{
    std::cerr << "Molecule::mutate_parameter() Unknown case" << std::endl;
    exit(-1);
  }
}

void Molecule::mutate_parameter_step_wrapped(double& parameter,const double min,const double max,const double step)
{
  parameter +=step*rand_karney.FixedS();
  if(parameter<min){
    parameter = 2.*min - parameter;
  }else if(parameter>max){
    parameter = 2.*max - parameter;
  }
}


void Molecule::convert_to_dna()
{
  /* If this is RNA, convert it to DNA. */
  switch(activity){
  case REPLI_RNA:
    activity = REPLI_DNA;
  case REPLI_DNA:
    decay_rate = Para::dna_decay_rate;
    break;
  case PARAS_RNA:
    activity = PARAS_DNA;
  case PARAS_DNA:
    decay_rate = Para::dna_decay_rate;
    break;
  case JUNK_RNA:
    activity = JUNK_DNA;
  case JUNK_DNA:
    decay_rate = Para::dna_decay_rate;
    break;
  default:
    break;
  }
}

void Molecule::convert_to_rna()
{
  /* If this is DNA, convert it to RNA. */
  switch(activity){
  case REPLI_DNA:
    activity = REPLI_RNA;
  case REPLI_RNA:
    decay_rate = Para::rna_decay_rate;
    break;
  case PARAS_DNA:
    activity = PARAS_RNA;
  case PARAS_RNA:
    decay_rate = Para::rna_decay_rate;
    break;
  case JUNK_DNA:
    activity = JUNK_RNA;
  case JUNK_RNA:
    decay_rate = Para::rna_decay_rate;
    break;
  default:
    break;
    }
}

#ifdef _COUNT_RNA_DNA
void Molecule::_inclement_rna_dna_count()
{
  if(is_rna()){
    ++_rna_count;
  }else{
    ++_dna_count;
  }
}
#endif //_COUNT_RNA_DNA

#endif
