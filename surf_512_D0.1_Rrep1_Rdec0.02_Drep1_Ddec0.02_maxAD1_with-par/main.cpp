#include <iostream>
#include <functional>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <limits>
#include <cstring>
#include <csignal>
#include <string>
#include <cstdlib>

#include "tbb/tbb.h"

#include "para.hpp"
#include "myca.hpp"
#include "random_wrapper.hpp"

/* Core of any simulation */
void args(int argc, char **argv);
void scale_rates();
void core_simulation(MyCA& myca,const long Time);

/* Types of the simulation */
void single_run();

/* This allows continuation of a simulation. It changes the initial
   value of Time in single_run() and properly set the name of
   Para::read_fname. */
long setup_continuation(); 

/* This is an ugly way of using class methods via function pointer in C
   codes. */
void button2(int x,int y);
void button3(int x,int y);

/* Signal handlers */
void handler_usr1(int);
void handler_usr2(int);

MyCA* myca_p;

/* Flag to stop a simulation */
/* Signal USR1 causes the program to output .sav file. Singnal USR2
   tells the program to output .sav file and terminate.  */
bool is_end=false;
bool is_usr1=false;

int main(int argc, char **argv)
{
  args(argc, argv);

  // Seed the random number generator
  rand_karney.Reseed(Para::seed);
  RandomWrapper::seed.push_back(Para::seed);

  scale_rates();

  /* set signals */
  if(SIG_ERR==signal(SIGUSR1,handler_usr1)){
    std::cerr << "main() Error, signal handler cannot be set." << std::endl;
  }
  if(SIG_ERR==signal(SIGUSR2,handler_usr2)){
    std::cerr << "main() Error, signal handler cannot be set." << std::endl;
  }

  tbb::task_scheduler_init init(5);

  tbb::tick_count t0 = tbb::tick_count::now();

  switch(Para::simulation_type){
  case 0:
    single_run(); 
    break;
  default: std::cerr << "main(): Unknown simulation_type" << std::endl; break;
  }
  
  tbb::tick_count t1 = tbb::tick_count::now();

  double t = (t1-t0).seconds();
  std::cout << " time = " << t << "\n";

  return (0);
}

void core_simulation(MyCA& myca,const long Time)
{
  /****** core of simulations ******/
  if(Para::model_type==Para::VESIC_NEUT_EQUIL){
    myca.update_replicator_whole_EQA();

    if(Time%Para::vesicle_update_interval==0){
      myca.vesicle_neutral_growth();
      myca.update_vesicle_whole();
      myca.divide_vesicle(Time);
    }
  }
  else if(Para::model_type==Para::VESIC_LIPID_EQUIL){
    myca.update_replicator_whole_EQA();
    myca.vesicle_decay();
    if(Time%Para::vesicle_update_interval==0){
      myca.update_vesicle_whole();
      myca.divide_vesicle(Time);
    }
  }
  else if(Para::model_type==Para::SURFACE_EQUIL){
    myca.update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::WELL_MIXED_EQUIL){
    myca.update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX){
    myca.update_replicator_whole_no_complex_EQA();
  }
  else if(Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX){
    myca.update_replicator_whole_no_complex_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX){
    myca.parallel_update_replicator_whole_no_complex_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_SURFACE_EQUIL){
    myca.parallel_update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){
    std::cerr << "main.c core_simulation() Do not use PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX" << std::endl;
    exit(-1);
    myca.parallel_update_replicator_whole_no_complex_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL){
    std::cerr << "main.c core_simulation() Do not use PARALLEL_WELL_MIXED_EQUIL" << std::endl;
    exit(-1);
    myca.parallel_update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL){
    myca.parallel_update_replicator_whole_EQA();
    if(Time%Para::vesicle_update_interval==0){
      myca.parallel_vesicle_neutral_growth();
      myca.parallel_update_vesicle_whole();
      myca.parallel_divide_vesicle(Time);
    }
  }
  else if(Para::model_type==Para::SERIAL_ONLY_CPM){
    myca.update_vesicle_whole();
  }
  else if(Para::model_type==Para::PARALLEL_ONLY_CPM){
    myca.parallel_update_vesicle_whole();
  }
  else if(Para::model_type==Para::PARALLEL_INORG_EQUIL){
    myca.parallel_update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::INORG_EQUIL){
    myca.update_replicator_whole_EQA();
  }
  else if(Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){
    myca.parallel_update_replicator_whole_EQA();

    if(Time%Para::vesicle_update_interval==0){
      myca.parallel_vesicle_fixed_growth();
      myca.parallel_update_vesicle_whole();
      myca.parallel_divide_vesicle(Time);
    }
  }
  else{
    std::cerr << "main() Error, unknown model_type=" << Para::model_type << std::endl;
    exit(-1);
  }
}


void single_run()
{
  /* create a world */
  long init_time = 0;
  if(Para::is_continue_last_save){
    init_time = setup_continuation();
  }

  long Time=init_time-1;
  MyCA myca(Para::sys_nrow,Para::sys_ncol);
  myca_p = &myca;
  myca.initialize(Time);
  myca.reset_movie_frame(init_time);

  /* begin the dynamics */
  for(Time=init_time;Time<=Para::max_time;++Time){
    myca.visualize(Time);

    if(Time%Para::plot_interval==0){
      is_end = (!myca.output(Time)) && !(Para::model_type==Para::PARALLEL_ONLY_CPM);
    }

#ifdef _COUNT_RNA_DNA   
    if(Time%Para::count_RNA_DNA_interval==0 && Time!=0){
      myca.output_count_rna_dna(Time);
    }
#endif // _COUNT_RNA_DNA

    if(Para::is_display)
      myca.mouse_event(button2,button3);
    
    if(is_end || is_usr1 || (Time%Para::save_interval==0 && Time!=init_time)){
      myca.save_wrapper(Time);

      if(is_end){
	break;
      }
      else if(is_usr1){
	is_usr1 = false;
      }
    }
    
    core_simulation(myca,Time);

    if(is_end){
      myca.save_wrapper(Time+1);
      break;
    }
  }
}



void scale_rates()
{
  /*** Scale reaction rate constants ***/
  if(Para::model_type==Para::VESIC_NEUT_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
     Para::model_type==Para::VESIC_LIPID_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL ||
     Para::model_type==Para::WELL_MIXED_EQUIL ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL ||
     Para::model_type==Para::PARALLEL_INORG_EQUIL ||
     Para::model_type==Para::INORG_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){
    /* Para::max_rna_recog_param, Para::max_dna_recog_param and
       Para::std_bind_param are not rates. Thus, we shouldn't scale
       them. */
    Para::diffusion_rate *= Para::rate_scale;
    Para::rna_decay_rate *= Para::rate_scale;
    Para::dna_decay_rate *= Para::rate_scale;
    Para::max_asso_diss_rate *= Para::rate_scale;
    Para::max_rna_repli_rate *= Para::rate_scale;
    Para::max_dna_repli_rate *= Para::rate_scale;
    Para::max_lipid_synth_rate *= Para::rate_scale;
    Para::vesicle_decay_rate *= Para::rate_scale;

    /* We want this simulation model to behave in qualitatively the same
       way as the ODE model when we use the same parameters. For this
       sake, we have to scale rate constants. Let me first explain how to
       scale the complex association reaction. Assuming that the system is
       well-mixed, we pick up two molecules. The chance that we pick up X
       and Y is (X/N)*(Y/N-1) + (Y/N)*(X/N-1) = 2*X*Y/N/N-1, whereas the
       chance of picking up X and X is (X/N)*(X-1/N-1) = X*X/N/N-1, where
       N is the total number of squares in CA. Therefore, if the chance of
       reaction upon the meeting of two molecules is equal for (X,X) and
       (X,Y), then the rate of reaction is twice greater for (X,Y) than
       for (X,X). Let me call this Advantage 1. This is the case in the
       Gillespie model. In the ODE model, if we write X'=a*X*X and
       Y'=b*X*Y, this would mean that X has twice chance of going through
       reaction upon the meeting of two X's. Biologically, we can
       interpret this as meaning that X can be either catalyst or
       template, so that there are two possible ways of reaction occuring
       for two X's. On the other hand, X and Y have only one possible
       way. Let me call this Advantage 2. Now, Advantage 1 gives a factor
       of two to the reaction rate of X+Y, whereas Advantage 2 also gives
       a factor of two to the reaction rate of X+X. Thus, effectively the
       advantages cancel each other, so we get the situation in the
       ODE. This simulation model takes into account Advantage 1
       implicitly (for it does not matter whether molecule is chosen
       firstly or chosen seconly as a neighbor), whereas it also takes
       into account Advantage 2 explicitly in the reaction
       algorithm. Then, assuming that the system is well-mixed and
       denoting the reaction rate constant by alpha, the net reaction rate
       of X+X is 2*alpha*X*X/N/N-1, whereas the reaction rate of X+Y is
       2*alpha*X*Y/N/N-1. Hence, we want to halven alpha to make the
       system behave in the way similar to the behavior of the ODE. */
    /* Next, let me explain how to scale the complex dissociation
       reaction. A complex molecule occupies two squares, so that it has
       the chance of being picked up twice as much as a simple
       molecule. Therefore, the rate constant of complex dissociation
       must be halved. */
    /* Scale the diffusion rate of complex molecules. Since a complex
       molecule occupies two squares, it has the chance of being picked up
       twice as much as a simple molecule. Thus, we simply have to halve
       the dissociation rate. We now see that both association rate
       constant and dissociation rate constant have to be halved to
       properly scale them. This boils down to halving
       max_asso_diss_rate. */
    Para::max_asso_diss_rate *=0.5;

    /* For the replication reaction, we have to scale them by dividing it
       by 4. This can be seen as follows. The reaction is C+E->3X. Here, C
       occupies two squares, so it has twice chance of being picked
       up. Moreover, it is a second-order reaction. Thus, the argument we
       made for the association reaction between X and Y also applies
       here. Thus, we have quadrupling effect, in that the rate of
       replication is 4*alpha*C*E, where alpha is the rate constant. Thus,
       we have to divide it by 4. By the way, note that we do some special
       thing about chosing a neighbor if we choose a complex molecule in
       the first place. Namely, we exclude the partner molecule of the
       chosen molecule as a possible neighbor. In other words, I chose one
       square from seven neighboring squares which does not contain the
       square where the partner molecule lies. If I didn't do this special
       trick, something strange would happen. Assuming that the system is
       well-mixed and we still stick to the algorithm of finite diffusion
       (i.e.\ we randomly choose a neighbor), the chance of choosing an
       empty square is (7/8)*E/N-1. This factor of 7/8 is the strange
       thing.  */
    Para::max_rna_repli_rate *= 0.25;
    Para::max_dna_repli_rate *= 0.25;

    /* As I previously explained for the complex dissociation reaction,
       diffusion of complex molecules must be halved too. By the way, we
       don't do anything special for the decay reaction of complex
       molecules. This is because the decay of a complex molecule is not a
       separate chemical reaction, but consists of "independent" decay of
       constituent molecules, which has the same rate constant as the
       decay of simple molecules. */
    Para::diffusion_rate_complex = 0.5*Para::diffusion_rate;

  }
  else if(Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){
    /* Para::max_rna_recog_param, Para::max_dna_recog_param and
       Para::std_bind_param are not rates. Thus, we shouldn't scale
       them. */
    Para::diffusion_rate *= Para::rate_scale;
    Para::rna_decay_rate *= Para::rate_scale;
    Para::dna_decay_rate *= Para::rate_scale;
    /* I do not scale max_asso_diss_rate. The rate of replication is
       written as m*fp*(1-ft)*rp*k (see below for notation). */
    Para::max_rna_repli_rate *= Para::rate_scale;
    Para::max_dna_repli_rate *= Para::rate_scale;
    Para::max_lipid_synth_rate *= Para::rate_scale;
    Para::vesicle_decay_rate *= Para::rate_scale;

    /* If we don't take complex formation into account, replication is
       tri-molecular reaction. This means that there are six (3!)
       different ways of choosing three molecules for one reaction. Note
       that I can factor in this six-fold speeding of reaction in my
       algorithm by deviding either max_asso_diss_rate or
       max_rna(dna)_repli_rate. I arbitrary chose repli_rate because
       it's more intuitive to do scaling there. */
    Para::max_rna_repli_rate /= 6.;
    Para::max_dna_repli_rate /= 6.;
  }
  else if(Para::model_type==Para::PARALLEL_ONLY_CPM ||
	  Para::model_type==Para::SERIAL_ONLY_CPM){
    /* do nothing */
  }
  else{
    std::cerr << "scale_rates() Unkown model_type." << std::endl;
    exit(-1);
  }


  /* Scale vesicle_update_inverval */
  Para::vesicle_update_interval = static_cast<int>(0.5+ Para::vesicle_update_interval_d/Para::rate_scale);
  if(Para::vesicle_update_interval<1){
    Para::vesicle_update_interval = 1;
    std::cerr << "scale_rates() Warning, rate_scale is too small so that vesicle_update_interval had to be set to 1." << std::endl;
  }

  /*** Error check ***/
  /* Under the equilibrium assumption, we do not consider fold_rate and
     unfold_rate as a rate, but as something related to equilibrium
     constant. Namely, we assume that fold_rate+unfold_rate=1. Thus, we
     should set Para::max_fold_unfold_rate to 1. */
  if(Para::model_type==Para::VESIC_NEUT_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
     Para::model_type==Para::VESIC_LIPID_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL ||
     Para::model_type==Para::WELL_MIXED_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX || 
     Para::model_type==Para::PARALLEL_INORG_EQUIL ||
     Para::model_type==Para::INORG_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){
    if(Para::max_fold_unfold_rate != 1.){
      std::cerr << "scale_rates() Although the system is under the equilibrium condition, max_fold_unfold_rate is not set to 1. Abort." << std::endl;
      exit(-1);
    }
  }
  else if(Para::model_type==Para::PARALLEL_ONLY_CPM ||
	  Para::model_type==Para::SERIAL_ONLY_CPM){
    /* do nothing */
  }
  else{
    std::cerr << "scale_rates() Unkown model_type. Do I have to scale init_fold_rate_repli and init_fold_rate_para?" << std::endl;
    exit(-1);
  }

   /* Mutations are mutually exclusive. The sum of mutation rates must be
     less than one. */
  if(1.<Para::mut_rate_rna_dna_repl_rate + 2.*Para::mut_rate_recog_param + Para::mut_rate_to_junk + Para::mut_rate_to_para 
     || 1.<Para::mut_rate_lipid_synth_rate + Para::mut_rate_fold_unfold_rate + Para::mut_rate_to_junk){
    std::cerr << "scale_rates() The sum of mutation rates is not less than 1!" << std::endl;
    exit(-1);
  }
  
  /* If the system is well-mixed, we have to set diffusion to zero. */
  if((Para::model_type==Para::WELL_MIXED_EQUIL ||
      Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
      Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX ||
      Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL)
     && Para::diffusion_rate!=0.){
    std::cerr << "scale_rates() If the system is well-mixed, diffusion_rate must be 0." << std::endl;
    exit(-1);
  }

  /* If the model doesn't assume complex formation, we don't use
     Automaton::Molecule::bind_param. However, despite this fact,
     bind_param is still set to max_asso_diss_rate. Then, setting
     max_asso_diss_rate to 0 is a good practice because this will reduce
     the size of save files and also makes it easy to detect any mistake
     such as mixing up the two different models. So, I force myself to
     do this. */
  if(Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){
    if(Para::max_asso_diss_rate != 0.){
      std::cerr << "scale_rates() max_asso_diss_rate != 0 even though the model assumes complex formation." << std::endl;
      exit(-1);
    }

    /* Currently, I disable mutations to parasites and junks [see
       Molecule::copy_with_mutation()]. So, I check if I am intending
       to do simulations with those mutations. */
    if(Para::mut_rate_to_junk!=0 || Para::mut_rate_to_para!=0){
      std::cerr << "scale_rates() mut_rate_junk or mut_rate_para is not zero. If you are intending to do simulations that require mutations to junk molecules or parasites, comment out appropriate regiouns in Molecule::copy_with_mutation()." << std::endl;
      exit(-1);
    }
  }
  
  /*** Report on the rates ***/
  /* Maximum replication rate depends on rna_dna_repl_constrain. See
     para.hpp */
  double max_repli_rate;
  if(Para::rna_dna_repl_constrain==0 ||
     Para::rna_dna_repl_constrain==1){
    max_repli_rate = (Para::max_rna_repli_rate>Para::max_dna_repli_rate) ? Para::max_rna_repli_rate : Para::max_dna_repli_rate;
  }
  else{
    std::cerr << "scale_rates() Unkown rna_dna_repl_constrain." << std::endl;
    exit(-1);
  }
  if(Para::model_type==Para::VESIC_NEUT_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
     Para::model_type==Para::VESIC_LIPID_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL ||
     Para::model_type==Para::WELL_MIXED_EQUIL ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL ||
     Para::model_type==Para::PARALLEL_INORG_EQUIL ||
     Para::model_type==Para::INORG_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    if(Para::fold_unfold_constrain==0){
      std::cerr << "EMPTY " << max_repli_rate + Para::diffusion_rate << std::endl;
      /* As I wrote before, the complex association between X and X has
	 two fold advantage in its rate, and we here want to know the
	 maximum association rate . In the current model, we consider the
	 folding/unfolding reaction not only for parasites, but also for
	 replicases (in my previous study, I considered only the folding
	 of parasites). It is impossible for two replicases to associate
	 with each other if both are always folded. Of course, if both are
	 always unfolded, they cannot associate either. The sum of the two
	 association rates can be written as A[X(1-Y)+Y(1-X)] where X is
	 the fold_rate of the one molecule, and Y is that of the other
	 molecule, and A is max_asso_diss_rate. It can be shown that this
	 function is maximal either at (1,0) or (0,1) and the maximum is
	 A. Hence, as the maximal association rate, I can use
	 max_asso_diss_rate. */
      std::cerr << "REPLI_RNA " << Para::max_asso_diss_rate + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
      /* Similaly, parasites cannot be a perfect template and a perfect
	 lipid synthetase at the same time because of folding. The sum of
	 these two rates are A(1-X)+BX where A is the maximal rate of
	 association and B is the maximal rate of lipid synthesis. If A>B,
	 then the maximum is A at X=0; if B>A, then the maximum is B at
	 X=1. I calculate this maximum rate as max_paras_reaction. */
      double max_paras_reaction = (Para::max_asso_diss_rate > Para::max_lipid_synth_rate) ? Para::max_asso_diss_rate : Para::max_lipid_synth_rate;
      std::cerr << "PARAS_RNA " << max_paras_reaction + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
      std::cerr << "DNA " << Para::max_asso_diss_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;
      std::cerr << "COMPLEX (RNA/DNA) " << max_repli_rate + Para::max_asso_diss_rate + Para::diffusion_rate_complex + ((Para::rna_decay_rate>Para::dna_decay_rate) ? Para::rna_decay_rate:Para::dna_decay_rate) << std::endl;
    }
    else if(Para::fold_unfold_constrain==1){
      std::cerr << "EMPTY " << max_repli_rate + Para::diffusion_rate << std::endl;
      /* As I wrote before, the complex association between X and X has
	 two fold advantage in its rate, and we here want to know the
	 maximum association rate. If
	 Para::switch_mol_recognition_assumption==1, we don't consider
	 the folding/unfolding. For REPLI_RNA, the sum of the two
	 association rates is as follows:

	 REPLI_RNA_1 + REPLI_RNA_2: bind_param*(rna_recog_param_1+rna_recog_param_2)
	 REPLI_RNA_1 + REPLI_DNA_2: bind_param*dna_recog_param_1
	 REPLI_RNA_1 + PARAS_RNA_2: parasite_advantage*bind_param*rna_recog_param_1
	 REPLI_RNA_1 + PARAS_DNA_2: parasite_advantage*bind_param*dna_recog_param_1
	 
	 Bind_param is set to Para::max_asso_diss_rate. Hence, the
	 maximum of the two numbers, 2*max_asso_diss_rate and
	 paras_advantage*max_asso_diss_rate is the maximum association
	 rate. */
      std::cerr << "REPLI_RNA " << (2.>Para::parasite_advantage ? 2.*Para::max_asso_diss_rate : Para::parasite_advantage*Para::max_asso_diss_rate) + Para::diffusion_rate + Para::rna_decay_rate << std::endl;

      std::cerr << "REPLI_DNA " << Para::max_asso_diss_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;      

      /* Since we don't have folding, parasites can freely do lipid
	 synthesis or complex formation. */
      std::cerr << "PARAS_RNA " << Para::parasite_advantage*Para::max_asso_diss_rate + Para::max_lipid_synth_rate + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
      std::cerr << "PARAS_DNA " << Para::parasite_advantage*Para::max_asso_diss_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;

      std::cerr << "COMPLEX (RNA/DNA) " << max_repli_rate + Para::max_asso_diss_rate + Para::diffusion_rate_complex + ((Para::rna_decay_rate>Para::dna_decay_rate) ? Para::rna_decay_rate:Para::dna_decay_rate) << std::endl;
    }
    else if(Para::fold_unfold_constrain==2){
      std::cerr << "EMPTY " << max_repli_rate + Para::diffusion_rate << std::endl;
      /* The calculation of the maximum association rate for REPLI_RNA
	 is almost the same as when fold_unfold_constrain==0. However,
	 the critical difference is that we here have explicit
	 parasites. When fold_unfold_constrain=1, we had to check if
	 parasite_advantage>2 to determine the max association
	 rate. Here, we instead have to check if
	 parasite_advantage>1. */
      std::cerr << "REPLI_RNA " << (Para::parasite_advantage > 1. ? Para::parasite_advantage*Para::max_asso_diss_rate : Para::max_asso_diss_rate) + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
      /* Similaly, parasites cannot be a perfect template and a perfect
	 lipid synthetase at the same time because of folding. The sum of
	 these two rates are A(1-X)+BX where A is the maximal rate of
	 association and B is the maximal rate of lipid synthesis. If A>B,
	 then the maximum is A at X=0; if B>A, then the maximum is B at
	 X=1. I calculate this maximum rate as max_paras_reaction. */
      double max_paras_reaction = (Para::parasite_advantage*Para::max_asso_diss_rate > Para::max_lipid_synth_rate) ? Para::parasite_advantage*Para::max_asso_diss_rate : Para::max_lipid_synth_rate;
      std::cerr << "PARAS_RNA " << max_paras_reaction + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
      std::cerr << "REPLI_DNA " << Para::max_asso_diss_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;
      std::cerr << "PARAS_DNA " << Para::parasite_advantage*Para::max_asso_diss_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;
      std::cerr << "COMPLEX (RNA/DNA) " << max_repli_rate + Para::max_asso_diss_rate + Para::diffusion_rate_complex + ((Para::rna_decay_rate>Para::dna_decay_rate) ? Para::rna_decay_rate:Para::dna_decay_rate) << std::endl;
    }
    else{
      std::cerr << "scale_rates(): Unknown value in switch_mol_recognition_assumption." << std::endl;
      exit(-1);
    }
    
  }
  else if(Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX ||
	  Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){ 
    /* For a system without complex formation, the rate of replication
       reaction is m*fp*(1-ft)*rp*k where m is max_asso_diss_rate; fp is
       the folding rate of polymerase; ft is the folding rate of
       template; rp is the recognition rate of polymerase; and k is the
       bare rate of replication. As I wrote above, the maximum rate of
       replication for a replicase (i.e. non-lipid synthesizing
       molecule) is obtained when one replicase is always folded, and
       the other is unfolded. Therefore, the maximum of fp*(1-ft) is
       1. Since rp is 1 at maximum, the maximum rate of replication is
       m*k. */
    std::cerr << "EMPTY " << Para::max_asso_diss_rate*max_repli_rate + Para::diffusion_rate << std::endl;
    std::cerr << "REPLI_RNA " << Para::max_asso_diss_rate*max_repli_rate + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
    /* As I wrote above, parasites have the form of A(1-X) + BX kind of
       reaction rate, where X is the folding rate. The maximum of A is
       m*k, whereas that of B is max_lipid_synth_rate. In the next line,
       I am comparing A and B, and I take the greater one. */
    double max_paras_reaction = (Para::max_asso_diss_rate*max_repli_rate > Para::max_lipid_synth_rate) ? Para::max_asso_diss_rate*max_repli_rate : Para::max_lipid_synth_rate;
    std::cerr << "PARAS_RNA " << max_paras_reaction  + Para::diffusion_rate + Para::rna_decay_rate << std::endl;
    std::cerr << "DNA_RNA " << Para::max_asso_diss_rate*max_repli_rate + Para::diffusion_rate + Para::dna_decay_rate << std::endl;
  }
  else if(Para::model_type==Para::PARALLEL_ONLY_CPM ||
	  Para::model_type==Para::SERIAL_ONLY_CPM){
    std::cerr << "Only vesicles: no need to output rates." << std::endl;
  }
  else{
    std::cerr << "scale_rates() Unkown model_type." << std::endl;
    exit(-1);
  }



  return;
}

void args(int argc, char **argv)
{
  std::string argument;
  for (int i = 1; i < argc; i++) {
    argument = argv[i];
    if(argument=="-Display"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it==1)
	Para::is_display=true;
      else if (it==0)
	Para::is_display=false;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-Scale"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it>0){
	Para::scale = it;
      }
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-Movie"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it==1)
	Para::is_movie=true;
      else if (it==0)
	Para::is_movie=false;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-DisplayInt"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::display_interval=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-SaveInt"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::save_interval=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MovieInt"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::movie_interval=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MovieDir"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      strncpy(Para::movie_directory_name,argv[i],255);
    }
    else if(argument=="-Seed"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      Para::seed = boost::lexical_cast<unsigned long>(argv[i]);
    }
    else if(argument=="-MaxTime"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      unsigned ui = boost::lexical_cast<unsigned long>(argv[i]);
      if(ui==0)
	Para::max_time = std::numeric_limits<long>::max();
      else
	Para::max_time = ui;
    }
    else if(argument=="-RateScale"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::rate_scale=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-DiffRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::diffusion_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-CrossDiff"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::cross_membrane=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-RNADecayRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::rna_decay_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-DNADecayRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::dna_decay_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-VesDecayRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::vesicle_decay_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxAssoDissRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::max_asso_diss_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxRNARepliRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::max_rna_repli_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxDNARepliRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::max_dna_repli_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxFoldUnfoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::max_fold_unfold_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxLipidSynthRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::max_lipid_synth_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApFoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_fold_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApFoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_fold_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaFoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_para_fold_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitLipidSynthRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_lipid_synth_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateRecogParam"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.){
	Para::mut_rate_recog_param=db;
	Para::mut_rate_RNA_recog_param_Rp=db;
	Para::mut_rate_DNA_recog_param_Rp=db;
	Para::mut_rate_RNA_recog_param_Dp=db;
	Para::mut_rate_DNA_recog_param_Dp=db;
      }
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutStepRecogParam"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_step_recog_param=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateRNADNAReplRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_rna_dna_repl_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateFoldUnfoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if(db>=0.){
	Para::mut_rate_fold_unfold_rate=db;
	Para::mut_rate_fold_unfold_rate_Rp=db;
	Para::mut_rate_fold_unfold_rate_Dp=db;
      }
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutStepFoldUnfoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_step_fold_unfold_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateLipidSynthRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_lipid_synth_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutStepLipidSynthRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_step_lipid_synth_rate=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateToPara"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_to_para=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateToJunk"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_to_junk=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApRNArecog"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_rna_recog_param=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApDNArecog"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_dna_recog_param=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApRNArecog"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_rna_recog_param=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApDNArecog"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_dna_recog_param=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApRNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNApRNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApDNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNApDNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApRNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNApRNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApDNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNApDNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaRNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_paraRNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaDNAdens"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_paraDNA_dens=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-VolThr"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it>0){
	Para::volume_threshold = it;
      }
      else{
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MaxTarVol"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it>0){
	Para::max_target_volume = it;
      }
      else{
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    } 
    else if(argument=="-FixedTarVol"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it>0){
	Para::fixed_target_volume = it;
      }
      else{
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
   else if(argument=="-MinRNApFoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::min_fold_rate_RNAp=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-PlotInt"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::plot_interval=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-PlotFile"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      Para::plot_fname=argv[i];
    }
    else if(argument=="-ReadFile"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      Para::read_fname=argv[i];
    }
    else if(argument=="-ReadFile2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      Para::read_fname2=argv[i];
    }
    else if(argument=="-OldRateScale"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::old_rate_scale=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-HistMaxY"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::hist_max_y=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-HistFoldMinX"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::hist_fold_min_x=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-HistFoldMaxX"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::hist_fold_max_x=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-HistRecoMinX"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::hist_reco_min_x=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-HistRecoMaxX"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::hist_reco_max_x=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-CompetitionExp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it==0){
	Para::init_competition_experiment = 0;
      }
      else{
	Para::init_competition_experiment = 1;
      }
    }
    else if(argument=="-InitRNApRNArecog2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_rna_recog_param_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApDNArecog2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_dna_recog_param_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApRNArecog2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_rna_recog_param_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApDNArecog2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_dna_recog_param_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApRNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNApRNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApDNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNApDNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApRNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNApRNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApDNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNApDNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaRNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_paraRNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaDNAdens2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_paraDNA_dens_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitRNApFoldRate2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_RNAp_fold_rate_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitDNApFoldRate2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_DNAp_fold_rate_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-InitParaFoldRate2"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::init_para_fold_rate_2=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-NeutralGrowthFactor"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::neutral_growth_expand_factor=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-FileNrow"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::file_nrow=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-FileNcol"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::file_ncol=it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-VesicleFile"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      Para::read_vesicle_fname=argv[i];
    }
    /* If Continue=1, the program looks for the last ".sav" or
       ".sav.bz2" file in ./sav. If it's bzipped, it extracts it through
       bzcat. Then, the program sets Para::read_fname to the name of the
       last ".sav" file and also sets the initial time to the number in
       the ".sav" file name. */
    else if(argument=="-Continue"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it==1)
	Para::is_continue_last_save=true;
      else if (it==0)
	Para::is_continue_last_save=false;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-RemoveSpecies"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      int it = boost::lexical_cast<int>(argv[i]);
      if (it==1)
	Para::is_removing_some_species = true;
      else if (it==0)
	Para::is_removing_some_species = false;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-CountRNADNAInt"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      long it = boost::lexical_cast<long>(argv[i]);
      if (it>0)
	Para::count_RNA_DNA_interval = it;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-ParasiteAdv"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::parasite_advantage=db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutInitRNARecogRp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::mut_init_rna_recog_param_Rp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutInitDNARecogRp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::mut_init_dna_recog_param_Rp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutInitRNARecogDp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::mut_init_rna_recog_param_Dp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutInitDNARecogDp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0. && db<=1.)
	Para::mut_init_dna_recog_param_Dp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutInitFoldRate"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_init_fold_rate = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateRNARecogParamRp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_RNA_recog_param_Rp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateDNARecogParamRp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_DNA_recog_param_Rp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateFoldUnfoldRateRp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_fold_unfold_rate_Rp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateRNARecogParamDp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_RNA_recog_param_Dp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateDNARecogParamDp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_DNA_recog_param_Dp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else if(argument=="-MutRateFoldUnfoldRateDp"){
      i++; if (i == argc) {std::cerr << "No parameter after " << argument << std::endl; exit(-1);}
      double db = boost::lexical_cast<double>(argv[i]);
      if (db>=0.)
	Para::mut_rate_fold_unfold_rate_Dp = db;
      else {
	std::cerr << "Args() wrong value in " << argument << std::endl;
        exit(-1);
      }
    }
    else{
      std::cerr << "Unknown option: \"" << argument << '"' << std::endl;
      std::cerr << "-SaveInt -Scale -Display -DisplayInt -Movie -MovieInt -MovieDir -Seed -MaxTime -RateScale -DiffRate -CrossDiff -RNADecayRate -DNADecayRate -MaxRNARepliRate -MaxDNARepliRate -MaxFoldUnfoldRate -MaxLipidSynthRate -InitRNApFoldRate -InitDNApFoldRate -InitParaFoldRate -InitLipidSynthRate -MutRateRecogParam -MutStepRecogParam -MutRateRNADNAReplRate -MutRateFoldUnfoldRate -MutStepFoldUnfoldRate -MutRateLipidSynthRate -MutStepLipidSynthRate -MutRateToPara -MutRateToJunk -InitRNApRNArecog -InitRNApDNArecog -InitDNApRNArecog -InitDNApDNArecog -InitRNApRNAdens -InitRNApDNAdens -InitDNApRNAdens -InitDNApDNAdens -InitParaRNAdens -InitParaDNAdens -VolThr -NeutralGrowthFactor -FixedTarVol -VesDecayRate -MinRNApFoldRate -PlotInt -PlotFile -ReadFile -ReadFile2 -OldRateScale -HistMaxY -HistRecoMinX -HistRecoMaxX -HistFoldMinX -HistFoldMaxX -CompetitionExp -FileNrow -FileNcol -VesicleFile -Continue -RemoveSpecies -CountRNADNAInt -ParasiteAdv -MutInitRNARecogRp -MutInitDNARecogRp -MutInitRNARecogDp -MutInitDNARecogDp -MutInitFoldRate -MutRateRNARecogParamRp -MutRateDNARecogParamRp -MutRateFoldUnfoldRateRp -MutRateRNARecogParamDp -MutRateDNARecogParamDp -MutRateFoldUnfoldRateDp" << std::endl;
      exit(-1);
    }
  }
}

void button2(int x,int y)
{
  myca_p->button2(x,y);
}

void button3(int x,int y)
{
  myca_p->button3(x,y);
}

void handler_usr1(int)
{
  signal(SIGUSR1,handler_usr1);
  is_usr1 = true;
  std::cerr << "hanlder_usr1() is called." << std::endl;
}

void handler_usr2(int)
{
  signal(SIGUSR2,handler_usr2);
  is_end = true;
  std::cerr << "hanlder_usr2() is called." << std::endl;
}

long setup_continuation()
{
  if(system("$HOME/study/dna-simple/program/find_last_save/find_last_save.sh")){
    std::cerr << "setup_continuation() Error occured while calling find_last_save.sh" << std::endl;
    exit(-1);
  }

  std::ifstream fin("last_save");
  if(!fin.is_open()){
    std::cerr << "MyCA::setup_continuation() Error. Cannot open last_save." << std::endl;
    exit(-1);
  }

  getline(fin,Para::read_fname);

  std::string line;
  getline(fin,line);
  return boost::lexical_cast<long>(line);
}
