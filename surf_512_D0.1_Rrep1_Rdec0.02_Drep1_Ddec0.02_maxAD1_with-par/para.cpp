#include "para.hpp"

RandomLib::Random rand_karney(0);

namespace Para {
  unsigned long seed = 125ul;

  long max_time = std::numeric_limits<long>::max();

  /*** CPM parameters ***/
  /* CPM parameters are now const */
  /* lambda is coefficient for volume constrain term */
  /* note that temperature cannot be modified once MyCA is
     instantiated. See MyCA::MyCA() */

  /*** Vesicle parameters ***/
  int volume_threshold = 1000;
  /* neutral_growth_expand_factor is used when Para::model_type =
     VESIC_NEUT_EQUIL. The target volume of vesicles will be set to
     num_mol*neutral_growth_expand_factor. */
  double neutral_growth_expand_factor = 1.1;
  int max_target_volume = std::numeric_limits<int>::max();
  int fixed_target_volume = 1100;
  double vesicle_decay_rate = 0.06;
  double cross_membrane = 0.;
  /* This is the times of updating replicase plane between vesicle
     updating events. This value is scaled through multiplying
     rate_scale. Thus, the actual number is equal to the times of
     updating replicase plane betwen vesicle updating events when
     rate_scale is set to 1. If rate_scale = 0.331 and
     vesicle_update_interval_d=6.62, vesicle_update_interval becomes
     20. */
  double vesicle_update_interval_d=6.62; 
  int vesicle_update_interval;
  /* if model==(PARALLEL_)ONLY_CPM */
  int only_cpm_init_tar_vol=500;
  unsigned only_cpm_init_n_vesicle=550;

  /*** Boundaries for rate constants ***/
  double rate_scale = 1.;
  double diffusion_rate = 0.1;
  /* This should be a half of diffusion_rate.*/
  double diffusion_rate_complex;
  double rna_decay_rate = 0.02;
  double dna_decay_rate = 0.002;
  double max_asso_diss_rate = 1.;
  double max_rna_repli_rate = 1.;
  double max_dna_repli_rate = 0.09;
  /* If the system is under the equilibrium assumption, this should be
     1. See also para.hpp fold_unfold_constrain and
     Automaton::initialize() */
  double max_fold_unfold_rate = 1.;
  double min_fold_rate_RNAp = 0.;
  double max_lipid_synth_rate = 0.;
  double parasite_advantage = 1.1;

  /*** Initial values of rate constants ***/
  /* This specifies the rate of association such that association_rate =
     max_asso_diss_rate * RNA_or_DNA_recog_para * fold_rate_repli *
     fold_rate_template. dissociation_rate = max_asso_diss_rate *
     (1-RNA_or_DNA_recog_para)*/
  double init_RNAp_rna_recog_param = 0.6; //[0,1]
  double init_RNAp_dna_recog_param = 0.6; //[0,1]
  double init_DNAp_rna_recog_param = 0.6; //[0,1]
  double init_DNAp_dna_recog_param = 0.6; //[0,1]
  double init_RNApRNA_dens = 1.;
  double init_RNApDNA_dens = 0.;
  double init_DNApRNA_dens = 0.;
  double init_DNApDNA_dens = 0.;
  double init_paraRNA_dens = 0.;
  double init_paraDNA_dens = 0.;
  double init_RNAp_fold_rate = 0.3;
  double init_DNAp_fold_rate = 0.3;
  double init_para_fold_rate = 0.1;
  double init_lipid_synth_rate = 0.;

  bool is_removing_some_species = false;

  int   init_competition_experiment = 0;
  double init_RNAp_rna_recog_param_2;
  double init_RNAp_dna_recog_param_2;
  double init_DNAp_rna_recog_param_2;
  double init_DNAp_dna_recog_param_2;
  double init_RNApRNA_dens_2;
  double init_RNApDNA_dens_2;
  double init_DNApRNA_dens_2;
  double init_DNApDNA_dens_2;
  double init_paraRNA_dens_2;
  double init_paraDNA_dens_2;
  double init_RNAp_fold_rate_2;
  double init_DNAp_fold_rate_2;
  double init_para_fold_rate_2;

  /*** Mutation rate of rate constants ***/
  double mut_rate_recog_param = 0.01;
  /* mut_step is the total width of the uniform distribution. */
  double mut_step_recog_param = 0.1;
  double mut_rate_rna_dna_repl_rate = 0.00001;
  double mut_rate_fold_unfold_rate = 0.01;
  double mut_step_fold_unfold_rate = 0.1;
  double mut_rate_lipid_synth_rate = 0.01;
  double mut_step_lipid_synth_rate = 0.1;
  double mut_rate_to_para = 0.0000;
  double mut_rate_to_junk = 0.;

  double mut_init_rna_recog_param_Rp=0.;
  double mut_init_dna_recog_param_Rp=0.;
  double mut_init_rna_recog_param_Dp=0.;
  double mut_init_dna_recog_param_Dp=0.;
  double mut_init_fold_rate=0.;

  double mut_rate_RNA_recog_param_Rp=0.01;
  double mut_rate_DNA_recog_param_Rp=0.01;
  double mut_rate_fold_unfold_rate_Rp=0.01;
  double mut_rate_RNA_recog_param_Dp=0.01;
  double mut_rate_DNA_recog_param_Dp=0.01;
  double mut_rate_fold_unfold_rate_Dp=0.01;

  /*** Visualization ***/
  bool is_display = true;
  bool is_movie = false;
  char movie_directory_name[255] = "movie";
  int margin=10;
  int scale=1;
  long display_interval = 500;
  long movie_interval = 100000;

  /*** Histogram ***/
  double hist_max_y = 0.3;
  unsigned hist_n_bins = 100;
  double hist_fold_min_x = 0.;
  double hist_fold_max_x = 1.;
  double hist_reco_min_x = 0.;
  double hist_reco_max_x = 1.;

  /*** Output ***/
  std::string plot_fname = "plot.dat";
  long plot_interval = 1000;
  long save_interval = 1000000;
  char save_directory_name[255] = "save";
  long count_RNA_DNA_interval = 100000;
  char count_rna_dna_directory_name[255] = "count";

  /*** Input ***/
  std::string read_fname;
  unsigned file_nrow = Para::sys_nrow;
  unsigned file_ncol = Para::sys_ncol;
  std::string read_vesicle_fname;
  bool is_continue_last_save = false;
  std::string read_fname2;
  double old_rate_scale = 0.;
}
