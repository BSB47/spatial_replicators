#include "RandomLib/Random.hpp"
#include <limits>

#ifndef PARAMETERS
#define PARAMETERS


extern RandomLib::Random rand_karney;

namespace Para {
  const unsigned sys_nrow=512;
  const unsigned sys_ncol=512;
  /* 
     The types of the model: 

     VESIC_NEUT_EQUIL: 

        Compartmentalized system (equilibrium assumption). The growth of
     target volume depends only on the total number of molecules.

        Equilibrium assumption is made.

        The replication rate of replicases is either rna_replication>0 &
     dna_replication=0, or rna_replication=0 & dna_replication>0.

     ------------------------------------------------------------------

     VESIC_LIPID_EQUIL: 

        Compartmentalized system (equilibrium assumption). The growth of
     target volume depends on the lipid synthesis reaction. Target
     volume decays linearly.

        The replication rate of replicases is either rna_replication>0 &
     dna_replication=0, or rna_replication=0 & dna_replication>0.

     ------------------------------------------------------------------

     SURFACE_EQUIL: 

        Surface-bound system (equilibrium assumption).

        The replication rate of replicases is either rna_replication>0 &
     dna_replication=0, or rna_replication=0 & dna_replication>0.

     ------------------------------------------------------------------

     SURFACE_EQUIL_NO_COMPLEX:

        Surface-bound system with equilibrium assumption without complex
        formation.
  */
  enum ModelType {VESIC_NEUT_EQUIL,VESIC_LIPID_EQUIL,SURFACE_EQUIL,WELL_MIXED_EQUIL,SURFACE_EQUIL_NO_COMPLEX,WELL_MIXED_EQUIL_NO_COMPLEX,PARALLEL_SURFACE_EQUIL_NO_COMPLEX,PARALLEL_SURFACE_EQUIL,PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX,PARALLEL_WELL_MIXED_EQUIL,PARALLEL_VESIC_NEUT_EQUIL,PARALLEL_ONLY_CPM,SERIAL_ONLY_CPM,PARALLEL_INORG_EQUIL,INORG_EQUIL,PARALLEL_VESIC_FIXED_EQUIL};
  const ModelType model_type = SURFACE_EQUIL;//PARALLEL_VESIC_NEUT_EQUIL;//WELL_MIXED_EQUIL;//VESIC_NEUT_EQUIL;//SURFACE_EQUIL;//PARALLEL_INORG_EQUIL;//PARALLEL_VESIC_FIXED_EQUIL;//INORG_EQUIL;//PARALLEL_ONLY_CPM;//SERIAL_ONLY_CPM;//VESIC_LIPID_EQUIL;//PARALLEL_SURFACE_EQUIL_NO_COMPLEX;//SURFACE_EQUIL_NO_COMPLEX;//WELL_MIXED_EQUIL_NO_COMPLEX;//PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX;//

  /* the type of boundary */
#define NEIGH_X neigh_fix
#define XY_NEIGH_X xy_neigh_fix

  /* 0: RNA recognition and DNA recognition is independent. That is,
     rna_recog_param is in [0,1] and dna_recog_param is in [0,1]. 

     1: Almost the same as 0, but it disables DNA->DNA replication (Dp's
     dna_recog=0). This option requires Para::rna_dna_repl_constrain=0
     because we need to make a sharp distinction between RNAp and
     DNAp. Note that the initial state of the system is not modified by
     this option. So, we must additionaly do something with the
     initialization so that there is no Dp that does DNA->DNA
     replication. See MyCA::get_rid_of_some_species().

     2: Almost the same as 0, but it disables RNA->DNA reverse
     transcription (Dp's rna_recog=0). Note that this option requires
     Para::rna_dna_repl_constrain=0 because we need to make a sharp
     distinction between RNAp and DNAp. Note that the initial state of
     the system is not modified by this option. So, we must additionaly
     do something with the initialization so that there is no Dp that
     does RNA->DNA reverse_transcription. See
     MyCA::get_rid_of_some_species().

     3: It disables transcription: DNA->RNA. So, Rp's dna_recog=0. This
     option requires Para::rna_dna_repl_constrain=0 because we need to
     make a sharp distinction between RNAp and DNAp. Note that the
     initial state of the system is not modified by this option. So, we
     must additionaly do something with the initialization so that there
     is no Rp that does DNA->RNA replication. See
     MyCA::get_rid_of_some_species().
  */
  const int rna_dna_recog_constrain = 0;

  /* 0: Binary decision (i.e. perfect specificity is assumed);
     i.e. rna_repli_rate=max_rna_repli_rate & dna_repli_rate=0 or
     rna_repli_rate=0 & dna_repli_rate=max_dna_repli_rate. 

     1: Binary decision, plus non-copying initialization of folding &
     recognition parameters. Upon Rp<=>Dp mutation, parameters are
     re-set to pre-specified values.
  */
  const int rna_dna_repl_constrain = 0;

   /* This switch controles how the recognition between molecules are
     modeled. 

     0: fold_param + unfold_param = Para::max_fold_unfold_rate.

     Thus, if the system is under the equilibrium assumption, set
     max_fold_unfold_rate to 1. [For a historical reason, I allowed the
     possibility of fold_param + unfold_param = 1. But this feature will
     likeliy to be not used because we always assume equilibrium
     assumption.]

     We assume that everybody recognizes everybody with an equal
     affinity; i.e. we don't assume any specificity in the recognition
     between molecules. Moreover, we assume that ribozymes must fold in
     order to execute its catalytic activity. While a molecule is
     folded, it cannot be replicated. Thus, the difference between
     replicase and parasite can be made by the difference in folding. We
     also assume that DNA is always unfolded.

     1: We do not assume folding. Ribozymes can catalyze replication
     without folding and can also always be replicated (unless it's
     making a complex). We assume that the parasite has an advantage in
     being recognized by replicases.

     So, fold_param + unfold_param is undefined. They are not used. I
     set fold_param=0, unfold_param=Para::max_fold_unfold_rate.

     2: We assume folding. The difference from fold_unfold_constrain=1
     is that we also take folding into account in the dissociation
     rate. If fold_unfold_constrain=1, we have

     dissoc=bind_param*(1-rna_recog_param)

     If fold_unfold_constrain=2, we have

     dissoc=unfold_rate*bind_param*(1-rna_recog_param)*partner.fold_rate

     This will favor complex formation more than before, and it is hoped
     that this will make the model with folding more similar to the
     model without folding with respect to the timescale of evolutionary
     dynamics.

  */

  const int fold_unfold_constrain = 1;

  /* Tye type of simulation run.
     0: Normal run. (Single evolutionary simulation.) */
  const int simulation_type = 0;

  /* The switch for type of mutation.
     0: Independent mutation on rna_dna_recog_param,
     rna_dna_repli_rate and fold_unfold_rate. */
  const int mutation_type = 0;

  /* The switch for method parameter mutation. Either by step
     or uniform distribution in entire parameter space: 
     0: by step with wrapping boundary condition
  */
  const int switch_mutation_parameter_method = 0;

  /* 0: colored by fold parameter.  

     1: clored according to which file they come from. This is used in
     competition experiment using MyCA::read_half_half().

     2: colored by whether the original template is RNA or DNA
     template. If this option is chosen, the program also outputs the
     color of molecules in ".sav" files for the porpse of later analysis
     (see operator<<(std::ostream& s,const Molecule& mol) in
     automaton.cpp). However, note that at a moment of writing (June
     11), ".sav" file containing the information abot the color cannot
     be read by MyCA::read().

     3: colored by whether the molecule is RNA or DNA template.
   */
  const int switch_color_code_1st_panel = 3;

  /**** Conditional compilation ****/
  /* This option turns on the counting of the cumulative number of RNA
 and DNA templates along the ancestral lineage. */
#undef _COUNT_RNA_DNA

  /*******************************************************
   ************** From here, look at para.cpp ************
   *******************************************************/
  extern unsigned long seed;

  extern long max_time;

  /* These parameters are never modified in my model. To speed up, I
     make them const. */
  const int J_media_cell = 1;
  const int J_cell_cell = 3;
  const double temperature = 1.;
  const int lambda = 1;

  extern int volume_threshold;
  extern double neutral_growth_expand_factor;
  extern int max_target_volume;
  extern int fixed_target_volume;
  extern double vesicle_decay_rate;
  extern double cross_membrane;
  extern double vesicle_update_interval_d;
  extern int vesicle_update_interval;
  extern int only_cpm_init_tar_vol;
  extern unsigned only_cpm_init_n_vesicle;

  extern double rate_scale;
  extern double diffusion_rate;
  extern double diffusion_rate_complex;
  extern double rna_decay_rate;
  extern double dna_decay_rate;
  extern double max_asso_diss_rate;
  extern double max_rna_repli_rate;
  extern double max_dna_repli_rate;
  extern double max_fold_unfold_rate;
  extern double min_fold_rate_RNAp;
  extern double max_lipid_synth_rate;
  extern double parasite_advantage;

  extern double init_RNAp_rna_recog_param;
  extern double init_RNAp_dna_recog_param;
  extern double init_DNAp_rna_recog_param;
  extern double init_DNAp_dna_recog_param;
  extern double init_RNApRNA_dens;
  extern double init_RNApDNA_dens;
  extern double init_DNApRNA_dens;
  extern double init_DNApDNA_dens;
  extern double init_paraRNA_dens;
  extern double init_paraDNA_dens;
  extern double init_RNAp_fold_rate;
  extern double init_DNAp_fold_rate;
  extern double init_para_fold_rate;
  extern double init_lipid_synth_rate;

  extern int   init_competition_experiment;
  extern double init_RNAp_rna_recog_param_2;
  extern double init_RNAp_dna_recog_param_2;
  extern double init_DNAp_rna_recog_param_2;
  extern double init_DNAp_dna_recog_param_2;
  extern double init_RNApRNA_dens_2;
  extern double init_RNApDNA_dens_2;
  extern double init_DNApRNA_dens_2;
  extern double init_DNApDNA_dens_2;
  extern double init_paraRNA_dens_2;
  extern double init_paraDNA_dens_2;
  extern double init_RNAp_fold_rate_2;
  extern double init_DNAp_fold_rate_2;
  extern double init_para_fold_rate_2;

  extern bool is_removing_some_species;

  extern double mut_rate_recog_param;
  extern double mut_step_recog_param;
  extern double mut_rate_rna_dna_repl_rate;
  extern double mut_rate_fold_unfold_rate;
  extern double mut_step_fold_unfold_rate;
  extern double mut_rate_lipid_synth_rate;
  extern double mut_step_lipid_synth_rate;
  extern double mut_rate_to_para;
  extern double mut_rate_to_junk;

  extern double mut_init_rna_recog_param_Rp;
  extern double mut_init_dna_recog_param_Rp;
  extern double mut_init_rna_recog_param_Dp;
  extern double mut_init_dna_recog_param_Dp;
  extern double mut_init_fold_rate;

  extern double mut_rate_RNA_recog_param_Rp;
  extern double mut_rate_DNA_recog_param_Rp;
  extern double mut_rate_fold_unfold_rate_Rp;
  extern double mut_rate_RNA_recog_param_Dp;
  extern double mut_rate_DNA_recog_param_Dp;
  extern double mut_rate_fold_unfold_rate_Dp;

  extern bool is_display;
  extern bool is_movie;
  extern char movie_directory_name[255];
  extern int margin;
  extern int scale;
  extern long display_interval;
  extern long movie_interval;

  extern double hist_max_y;
  extern unsigned hist_n_bins;
  extern double hist_fold_min_x;
  extern double hist_fold_max_x;
  extern double hist_reco_min_x;
  extern double hist_reco_max_x;

  extern std::string plot_fname;
  extern long plot_interval;
  extern long save_interval;
  extern char save_directory_name[255];
  extern long count_RNA_DNA_interval;
  extern char count_rna_dna_directory_name[255];

  extern std::string read_fname;
  extern unsigned file_nrow;
  extern unsigned file_ncol;
  extern std::string read_vesicle_fname;
  extern bool is_continue_last_save;
  extern std::string read_fname2;
  extern double old_rate_scale;

  /*** The followings are not a parameter of the model ***/
  const int PROBABILITY_ARRAY_SIZE=25;
  const int COLOR_OFFSET = 156;
}
#endif 
