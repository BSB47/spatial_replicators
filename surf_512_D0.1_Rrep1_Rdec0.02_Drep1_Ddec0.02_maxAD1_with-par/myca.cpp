#include "myca.hpp"

MyCA::MyCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow(a_nrow), ncol(a_ncol), length(a_nrow * a_ncol),
      plane(a_nrow, a_ncol), serial_vesicles(a_nrow * a_ncol),
      atomic_vesicles(a_nrow * a_ncol), upd_order_shuffler(a_nrow, a_ncol),
      vesicle_observer(a_nrow * a_ncol) {
  /* Initilize probabilities[] */
  for (int i = 0; i != Para::PROBABILITY_ARRAY_SIZE; ++i) {
    probabilities[i] = exp(-i / Para::temperature);
  }

  /* Set boundaries (note that vesicle_index is set to -1 by default,
     which will produce the desired behavior in the boundary in vesicle
     updating step) */
  for (unsigned row = 0; row <= nrow + 1; ++row) {
    plane.cell(row, 0).molecule = Molecule('b');
    plane.cell(row, ncol + 1).molecule = Molecule('b');
  }
  for (unsigned col = 1; col <= ncol; ++col) {
    plane.cell(0, col).molecule = Molecule('b');
    plane.cell(nrow + 1, col).molecule = Molecule('b');
  }

  /* Designing the display */
  /* 0th panel: molecules by function (RNAp, DNAp, parasite & junk)
     1st panel: molecules by RNA or DNA
  */
  std::vector<CashPanelInfo> panel_info(8);
  int window_row = Para::margin + nrow + Para::margin + 100 + Para::margin;
  int window_col =
      (3 * Para::margin + 2 * (int)ncol > 6 * 100 + 7 * Para::margin)
          ? 3 * Para::margin + 2 * ncol
          : 6 * 100 + 7 * Para::margin;

  panel_info[FUNCTION].n_row = nrow;
  panel_info[FUNCTION].n_col = ncol;
  panel_info[FUNCTION].o_row = Para::margin;
  panel_info[FUNCTION].o_col = Para::margin;

  panel_info[PANEL_FOLD].n_row = nrow;
  panel_info[PANEL_FOLD].n_col = ncol;
  panel_info[PANEL_FOLD].o_row = Para::margin;
  panel_info[PANEL_FOLD].o_col = Para::margin + nrow + Para::margin;

  panel_info[HIST_RNAp_RNArecog].n_row = 100;
  panel_info[HIST_RNAp_RNArecog].n_col = 100;
  panel_info[HIST_RNAp_RNArecog].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_RNAp_RNArecog].o_col = Para::margin;

  panel_info[HIST_RNAp_DNArecog].n_row = 100;
  panel_info[HIST_RNAp_DNArecog].n_col = 100;
  panel_info[HIST_RNAp_DNArecog].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_RNAp_DNArecog].o_col = Para::margin + 100 + Para::margin;

  panel_info[HIST_RNAp_fold].n_row = 100;
  panel_info[HIST_RNAp_fold].n_col = 100;
  panel_info[HIST_RNAp_fold].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_RNAp_fold].o_col =
      Para::margin + 100 + Para::margin + 100 + Para::margin;

  panel_info[HIST_DNAp_RNArecog].n_row = 100;
  panel_info[HIST_DNAp_RNArecog].n_col = 100;
  panel_info[HIST_DNAp_RNArecog].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_DNAp_RNArecog].o_col = Para::margin + 100 + Para::margin +
                                         100 + Para::margin + 100 +
                                         Para::margin;

  panel_info[HIST_DNAp_DNArecog].n_row = 100;
  panel_info[HIST_DNAp_DNArecog].n_col = 100;
  panel_info[HIST_DNAp_DNArecog].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_DNAp_DNArecog].o_col = Para::margin + 100 + Para::margin +
                                         100 + Para::margin + 100 +
                                         Para::margin + 100 + Para::margin;

  panel_info[HIST_DNAp_fold].n_row = 100;
  panel_info[HIST_DNAp_fold].n_col = 100;
  panel_info[HIST_DNAp_fold].o_row = Para::margin + nrow + Para::margin;
  panel_info[HIST_DNAp_fold].o_col = Para::margin + 100 + Para::margin + 100 +
                                     Para::margin + 100 + Para::margin + 100 +
                                     Para::margin + 100 + Para::margin;

  /* Setting up the display */
  try {
    display_p = new CashDisplay(window_row, window_col, panel_info, Para::scale,
                                COLOR_MARGIN);
  } catch (std::bad_alloc) {
    std::cerr << "MyCA::MyCA(): Error, memory exhaustion" << std::endl;
    exit(-1);
  }

  /* Initilize display */
  /* COLOR_OFFSET means that the index of gradient color starts at
     COLOR_OFFSET and ends at COLOR_OFFSET+hist_n_bins-1. */
  display_p->color_yellow2red(Para::COLOR_OFFSET, Para::hist_n_bins);
  for (unsigned i = Para::COLOR_OFFSET;
       i < Para::COLOR_OFFSET + Para::hist_n_bins; ++i)
    color_of_gradient.push_back(static_cast<unsigned char>(i));

  display_p->color_rgb(COLOR_BLACK, 0, 0, 0);
  display_p->color_rgb(COLOR_WHITE, 255, 255, 255);
  display_p->color_rgb(COLOR_CELL_INSIDE, 20, 20, 20);
  display_p->color_rgb(COLOR_CELL_WALL, 255, 255, 255);
  display_p->color_rgb(COLOR_RNAp, 0, 0, 255);     // blue
  display_p->color_rgb(COLOR_DNAp, 0, 200, 0);     // dark green
  display_p->color_rgb(COLOR_PARAS, 255, 0, 0);    // red
  display_p->color_rgb(COLOR_JUNK, 188, 143, 143); // brown
  display_p->color_rgb(COLOR_HIST_BARS, 255, 255, 255);
  display_p->color_rgb(COLOR_HIST_BACK, 0, 0, 0);
  display_p->color_rgb(COLOR_MARGIN, 150, 150, 150);

  if (Para::is_display)
    display_p->open_window();
  if (Para::is_movie)
    display_p->open_png(Para::movie_directory_name);
}

MyCA::~MyCA() {
  if (display_p)
    delete display_p;
}

/* The method takes the initial time as arguments. If
   init_time%movie_interval==0 (i.e. a png file is created when ".sav"
   file is generated), the program will write a png file at
   Time=init_time as (init_time/movie_interval)-th png file. This will
   overwrite the existing file (I want to make sure that I have a png
   slide depicting the initial condition). If
   init_time%movie_interval!=0, the program won't output a png file when
   Time=init_time, so no overwriting.  */
void MyCA::reset_movie_frame(const long init_time) {
  int nframes = (init_time % Para::movie_interval)
                    ? (init_time / Para::movie_interval + 1)
                    : init_time / Para::movie_interval;
  display_p->reset_movie_frame(nframes);
}

void MyCA::initialize(const long Time) {
  if (Para::read_fname.empty()) {
    /* Note that at this point every square contains empty-molecule and
       vesicle_index=-1, except for boundaries.  */
    if (Para::model_type == Para::SURFACE_EQUIL ||
        Para::model_type == Para::SURFACE_EQUIL_NO_COMPLEX ||) {
      for (unsigned row = 1; row <= nrow; ++row)
        for (unsigned col = 1; col <= ncol; ++col) {
          if ((row - nrow / 2.) * (row - nrow / 2.) +
                  (col - ncol / 2.) * (col - ncol / 2.) <
              100 * 100) {
            if (!Para::init_competition_experiment || col >= ncol / 2) {
              double p = rand_karney.Real();
              if (p < Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "RNAp_RNA", Para::init_RNAp_fold_rate,
                    Para::init_RNAp_rna_recog_param,
                    Para::init_RNAp_dna_recog_param);
              } else if (p <
                         Para::init_RNApDNA_dens + Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "RNAp_RNA", Para::init_RNAp_fold_rate,
                    Para::init_RNAp_rna_recog_param,
                    Para::init_RNAp_dna_recog_param);
                plane.cell(row, col).molecule.convert_to_dna();
              } else if (p < Para::init_DNApRNA_dens + Para::init_RNApDNA_dens +
                                 Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "DNAp_RNA", Para::init_DNAp_fold_rate,
                    Para::init_DNAp_rna_recog_param,
                    Para::init_DNAp_dna_recog_param);
              } else if (p < Para::init_DNApDNA_dens + Para::init_DNApRNA_dens +
                                 Para::init_RNApDNA_dens +
                                 Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "DNAp_RNA", Para::init_DNAp_fold_rate,
                    Para::init_DNAp_rna_recog_param,
                    Para::init_DNAp_dna_recog_param);
                plane.cell(row, col).molecule.convert_to_dna();
              } else if (p < Para::init_paraRNA_dens + Para::init_DNApDNA_dens +
                                 Para::init_DNApRNA_dens +
                                 Para::init_RNApDNA_dens +
                                 Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "PARAS_RNA", Para::init_para_fold_rate);
              } else if (p < Para::init_paraDNA_dens + Para::init_paraRNA_dens +
                                 Para::init_DNApDNA_dens +
                                 Para::init_DNApRNA_dens +
                                 Para::init_RNApDNA_dens +
                                 Para::init_RNApRNA_dens) {
                plane.cell(row, col).molecule.initialize(
                    "PARAS_RNA", Para::init_para_fold_rate);
                plane.cell(row, col).molecule.convert_to_dna();
              }
            } else if (col < ncol / 2 && Para::init_competition_experiment) {
              double p = rand_karney.Real();
              if (p < Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "RNAp_RNA", Para::init_RNAp_fold_rate_2,
                    Para::init_RNAp_rna_recog_param_2,
                    Para::init_RNAp_dna_recog_param_2);
              } else if (p < Para::init_RNApDNA_dens_2 +
                                 Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "RNAp_RNA", Para::init_RNAp_fold_rate_2,
                    Para::init_RNAp_rna_recog_param_2,
                    Para::init_RNAp_dna_recog_param_2);
                plane.cell(row, col).molecule.convert_to_dna();
              } else if (p < Para::init_DNApRNA_dens_2 +
                                 Para::init_RNApDNA_dens_2 +
                                 Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "DNAp_RNA", Para::init_DNAp_fold_rate_2,
                    Para::init_DNAp_rna_recog_param_2,
                    Para::init_DNAp_dna_recog_param_2);
              } else if (p < Para::init_DNApDNA_dens_2 +
                                 Para::init_DNApRNA_dens_2 +
                                 Para::init_RNApDNA_dens_2 +
                                 Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "DNAp_RNA", Para::init_DNAp_fold_rate_2,
                    Para::init_DNAp_rna_recog_param_2,
                    Para::init_DNAp_dna_recog_param_2);
                plane.cell(row, col).molecule.convert_to_dna();
              } else if (p < Para::init_paraRNA_dens_2 +
                                 Para::init_DNApDNA_dens_2 +
                                 Para::init_DNApRNA_dens_2 +
                                 Para::init_RNApDNA_dens_2 +
                                 Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "PARAS_RNA", Para::init_para_fold_rate_2);
              } else if (p < Para::init_paraDNA_dens_2 +
                                 Para::init_paraRNA_dens_2 +
                                 Para::init_DNApDNA_dens_2 +
                                 Para::init_DNApRNA_dens_2 +
                                 Para::init_RNApDNA_dens_2 +
                                 Para::init_RNApRNA_dens_2) {
                plane.cell(row, col).molecule.initialize(
                    "PARAS_RNA", Para::init_para_fold_rate_2);
                plane.cell(row, col).molecule.convert_to_dna();
              }
            }
          }
        }

      /* If this is inorganic compartment smulation, we need to read in
         the vesicle only file. */
      if (Para::model_type == Para::PARALLEL_INORG_EQUIL ||
          Para::model_type == Para::INORG_EQUIL) {
        if (!Para::read_vesicle_fname.empty()) {
          read_only_vesicle(Para::read_vesicle_fname);
        } else {
          std::cerr << "MyCA::initialize() Error, although "
                       "model_type==PARALLEL_INORG_EQUIL, read_vesicle_fname "
                       "is not set."
                    << std::endl;
          exit(-1);
        }
      }
    } else if (Para::model_type == Para::WELL_MIXED_EQUIL ||
               Para::model_type == Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
               Para::model_type == Para::PARALLEL_WELL_MIXED_EQUIL ||
               Para::model_type == Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX) {
      for (unsigned row = 1; row <= nrow; ++row)
        for (unsigned col = 1; col <= ncol; ++col) {
          double p = rand_karney.Real();
          if (p < Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize(
                "RNAp_RNA", Para::init_RNAp_fold_rate,
                Para::init_RNAp_rna_recog_param,
                Para::init_RNAp_dna_recog_param);
          } else if (p < Para::init_RNApDNA_dens + Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize(
                "RNAp_RNA", Para::init_RNAp_fold_rate,
                Para::init_RNAp_rna_recog_param,
                Para::init_RNAp_dna_recog_param);
            plane.cell(row, col).molecule.convert_to_dna();
          } else if (p < Para::init_DNApRNA_dens + Para::init_RNApDNA_dens +
                             Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize(
                "DNAp_RNA", Para::init_DNAp_fold_rate,
                Para::init_DNAp_rna_recog_param,
                Para::init_DNAp_dna_recog_param);
          } else if (p < Para::init_DNApDNA_dens + Para::init_DNApRNA_dens +
                             Para::init_RNApDNA_dens +
                             Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize(
                "DNAp_RNA", Para::init_DNAp_fold_rate,
                Para::init_DNAp_rna_recog_param,
                Para::init_DNAp_dna_recog_param);
            plane.cell(row, col).molecule.convert_to_dna();
          } else if (p < Para::init_paraRNA_dens + Para::init_DNApDNA_dens +
                             Para::init_DNApRNA_dens + Para::init_RNApDNA_dens +
                             Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize("PARAS_RNA",
                                                     Para::init_para_fold_rate);
          } else if (p < Para::init_paraDNA_dens + Para::init_paraRNA_dens +
                             Para::init_DNApDNA_dens + Para::init_DNApRNA_dens +
                             Para::init_RNApDNA_dens +
                             Para::init_RNApRNA_dens) {
            plane.cell(row, col).molecule.initialize("PARAS_RNA",
                                                     Para::init_para_fold_rate);
            plane.cell(row, col).molecule.convert_to_dna();
          }
        }
    } else if (Para::model_type == Para::SERIAL_ONLY_CPM ||
               Para::model_type == Para::PARALLEL_ONLY_CPM) {
      /* We put a square vesicle of a small size in randomly chosen
         points for  Para::only_cpm_init_n_vesicle times. */

      /* Set the width of squared vesicles */
      unsigned half_width =
          static_cast<unsigned>(sqrt(Para::only_cpm_init_tar_vol * 0.5) * 0.5);
      if (half_width == 0)
        half_width = 1;

      /* Randomize the order of choosing a point */
      const UpdOrder &upd_order = upd_order_shuffler.get_upd_order_serial();
      UpdOrder::const_iterator pos = upd_order.begin();

      int vesicle_ind = -1;
      bool vesicle_has_not_been_produced;
      for (unsigned n_ves = 0; n_ves < Para::only_cpm_init_n_vesicle; ++n_ves) {
        vesicle_has_not_been_produced = true;
        do {
          for (unsigned row = pos->row - half_width;
               row < pos->row + half_width; ++row)
            for (unsigned col = pos->col - half_width;
                 col < pos->col + half_width; ++col) {
              /* If this coordinate is legal & this square is empty  */
              if (row >= 1 && row <= nrow && col >= 1 && col <= ncol &&
                  plane.cell(row, col).vesicle_index == -1) {
                /* Make a vesicle if not made yet */
                if (vesicle_has_not_been_produced) {
                  vesicle_has_not_been_produced = false;
                  vesicle_ind = vesicle_make_new_vesicle();
                  vesicle_set_target_volume(vesicle_ind,
                                            Para::only_cpm_init_tar_vol);
                }

                /* Place a vesicle */
                plane.cell(row, col).vesicle_index = vesicle_ind;
                vesicle_add_volume(vesicle_ind);
              }
            }
          /* Go to the next random position */
          ++pos;
          if (pos == upd_order.end()) {
            std::cerr << "MyCA::initialize() The initial number of vesicle "
                         "seems too large. Currently, n_ves="
                      << n_ves << std::endl;
            exit(-1);
          }
        } while (vesicle_has_not_been_produced);
      }
    } else {
      std::cerr << "MyCA::initialize() Unknown model_type" << std::endl;
      exit(-1);
    }
  }
  /* Read the file */
  else if (Para::read_fname2.empty()) {
    read(Para::read_fname);

    if (!Para::read_vesicle_fname.empty()) {
      read_only_vesicle(Para::read_vesicle_fname);
    }

    if (Para::is_removing_some_species) {
      get_rid_of_some_species();
    }
  } else {
    read_half_half(Para::read_fname, Para::read_fname2);

    if (Para::is_removing_some_species) {
      get_rid_of_some_species();
    }
  }

  /* If this switch is 1, we assume that we are using read_half_half()
     to initialize the system. In this case, we color molecules
     according to which file molecles come from.

     If this switch is 2, we color them according to whether they are
     RNA or DNA templates, which allows a simple ancestor tracing. */
  if (Para::switch_color_code_1st_panel == 0) {
    /* Do nothing */
    ;
  } else if (Para::switch_color_code_1st_panel == 1) {
    set_color_half_half();
  } else if (Para::switch_color_code_1st_panel == 2) {
    set_color_RNA_DNA();
  } else if (Para::switch_color_code_1st_panel == 3) {
    /* Do nothing */
    ;
  } else {
    std::cerr << "MyCA::initialize() Unknown value in "
                 "Para::switch_color_code_1st_panel."
              << std::endl;
    exit(-1);
  }
}

void MyCA::get_rid_of_some_species() {
  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (!plane.cell(row, col).molecule.is_empty()) {
        if (plane.cell(row, col).molecule.is_replicase() &&
            plane.cell(row, col).molecule.rate_rna_repli() > 0. &&
            plane.cell(row, col).molecule.get_dna_recog_param() < 0.5) {
          if (plane.cell(row, col).molecule.is_simple()) {
            plane.cell(row, col).molecule.decay();
          } else {
            plane.NEIGH_X(row, col, plane.cell(row, col).molecule.bon_nei)
                .molecule.dissociate();
            plane.cell(row, col).molecule.decay();
          }
        }
      }
    }
}

void MyCA::visualize(const long t) {
  if ((Para::is_display && t % Para::display_interval == 0) ||
      (Para::is_movie && t % Para::movie_interval == 0)) {
    scan_plane_to_display();

    if (Para::is_display && t % Para::display_interval == 0) {
      display_p->draw_window();
    }

    if (Para::is_movie && t % Para::movie_interval == 0) {
      display_p->draw_png();
    }
  }
}

void MyCA::scan_plane_to_display() {
  /* This is used for 0th & 1st panel */
  bool is_boundary;

  /* This is used to draw histograms */
  unsigned n_RNAp = 0, n_DNAp = 0;
  const unsigned n_bins = Para::hist_n_bins;
  std::vector<unsigned> occu_RNAp_RNArecog(n_bins, 0);
  std::vector<unsigned> occu_RNAp_DNArecog(n_bins, 0);
  std::vector<unsigned> occu_RNAp_fold(n_bins, 0);
  std::vector<unsigned> occu_DNAp_RNArecog(n_bins, 0);
  std::vector<unsigned> occu_DNAp_DNArecog(n_bins, 0);
  std::vector<unsigned> occu_DNAp_fold(n_bins, 0);

  const double hist_fold_min_x = Para::hist_fold_min_x;
  const double hist_fold_max_x = Para::hist_fold_max_x;
  const double hist_reco_min_x = Para::hist_reco_min_x;
  const double hist_reco_max_x = Para::hist_reco_max_x;

  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      plane.cell(row, col).vesicle_index = plane.cell(row, col).vesicle_index;

      /* We will directly write down a picture in the 0th & 1st
         panel. For histograms, we gather the data. */
      if (plane.cell(row, col).molecule.is_empty()) {
        display_p->put_pixel(FUNCTION, row, col, COLOR_BLACK);
        display_p->put_pixel(PANEL_FOLD, row, col, COLOR_BLACK);
      } else {
        /* Draw the 0th panel by function, and gather the data for
           histograms. */
        if (plane.cell(row, col).molecule.is_replicase()) {

          /* If this is RNAp */
          if (plane.cell(row, col).molecule.rate_rna_repli() > 0) {

            /* 0th panel */
            display_p->put_pixel(FUNCTION, row, col, COLOR_RNAp);

            /* Histogram */
            ++n_RNAp;
            ++occu_RNAp_RNArecog[get_bin(
                n_bins, hist_reco_min_x, hist_reco_max_x,
                plane.cell(row, col).molecule.get_rna_recog_param())];
            ++occu_RNAp_DNArecog[get_bin(
                n_bins, hist_reco_min_x, hist_reco_max_x,
                plane.cell(row, col).molecule.get_dna_recog_param())];
            ++occu_RNAp_fold[get_bin(
                n_bins, hist_fold_min_x, hist_fold_max_x,
                plane.cell(row, col).molecule.rate_fold())];
          }
          /* If this is DNAp */
          else {

            /* 0th panel */
            display_p->put_pixel(FUNCTION, row, col, COLOR_DNAp);

            /* Histogram */
            ++n_DNAp;
            ++occu_DNAp_RNArecog[get_bin(
                n_bins, hist_reco_min_x, hist_reco_max_x,
                plane.cell(row, col).molecule.get_rna_recog_param())];
            ++occu_DNAp_DNArecog[get_bin(
                n_bins, hist_reco_min_x, hist_reco_max_x,
                plane.cell(row, col).molecule.get_dna_recog_param())];
            ++occu_DNAp_fold[get_bin(
                n_bins, hist_fold_min_x, hist_fold_max_x,
                plane.cell(row, col).molecule.rate_fold())];
          }
        }
        /* If this is a parasite */
        else if (plane.cell(row, col).molecule.is_parasite()) {
          display_p->put_pixel(FUNCTION, row, col, COLOR_PARAS);
        } else if (plane.cell(row, col).molecule.is_junk()) {
          display_p->put_pixel(FUNCTION, row, col, COLOR_JUNK);
        } else {
          std::cerr << "MyCA::scan_plane_to_display() Unknown activity type!"
                    << std::endl;
          exit(-1);
        }

        /* 1st by folding */
        if (Para::switch_color_code_1st_panel == 0) {
          display_p->put_pixel(PANEL_FOLD, row, col,
                               color_of_gradient[get_bin(
                                   n_bins, hist_fold_min_x, hist_fold_max_x,
                                   plane.cell(row, col).molecule.rate_fold())]);
        } else if (Para::switch_color_code_1st_panel == 1 ||
                   Para::switch_color_code_1st_panel == 2) {
          display_p->put_pixel(PANEL_FOLD, row, col,
                               plane.cell(row, col).molecule.get_color());
        } else if (Para::switch_color_code_1st_panel == 3) {
          if (plane.cell(row, col).molecule.is_rna())
            display_p->put_pixel(PANEL_FOLD, row, col, COLOR_RNAp);
          else
            display_p->put_pixel(PANEL_FOLD, row, col, COLOR_DNAp);
        } else {
          std::cerr << "MyCA::scan_plane_to_display() Unknown value in "
                       "Para::switch_color_code_1st_panel."
                    << std::endl;
          exit(-1);
        }
      }

      /* Inside of vesicle (only for 0th and 1st panels ) */
      if (plane.cell(row, col).vesicle_index != -1) {

        /* Check if this pixel is a vesicle boundary */
        is_boundary = false;
        for (unsigned i = 1; i != 5; ++i) {
          if (plane.cell(row, col).vesicle_index !=
              plane.NEIGH_X(row, col, i).vesicle_index) {
            is_boundary = true;
            break;
          }
        }

        /* Draw vesicle boundaries */
        if (is_boundary) {
          display_p->put_pixel(FUNCTION, row, col, COLOR_CELL_WALL);
          display_p->put_pixel(PANEL_FOLD, row, col, COLOR_CELL_WALL);
        } else if (plane.cell(row, col).molecule.is_empty()) {
          display_p->put_pixel(FUNCTION, row, col, COLOR_CELL_INSIDE);
          display_p->put_pixel(PANEL_FOLD, row, col, COLOR_CELL_INSIDE);
        }
      }
    }

  /* Convert occurence to frequencies */
  std::vector<double> freq_RNAp_RNArecog(n_bins, 0.);
  std::vector<double> freq_RNAp_DNArecog(n_bins, 0.);
  std::vector<double> freq_RNAp_fold(n_bins, 0.);
  if (n_RNAp > 0) {
    for (unsigned i = 0; i < n_bins; ++i) {
      freq_RNAp_RNArecog[i] =
          static_cast<double>(occu_RNAp_RNArecog[i]) / n_RNAp;
      freq_RNAp_DNArecog[i] =
          static_cast<double>(occu_RNAp_DNArecog[i]) / n_RNAp;
      freq_RNAp_fold[i] = static_cast<double>(occu_RNAp_fold[i]) / n_RNAp;
    }
  }
  std::vector<double> freq_DNAp_RNArecog(n_bins, 0.);
  std::vector<double> freq_DNAp_DNArecog(n_bins, 0.);
  std::vector<double> freq_DNAp_fold(n_bins, 0.);
  if (n_DNAp > 0) {
    for (unsigned i = 0; i < n_bins; ++i) {
      freq_DNAp_RNArecog[i] =
          static_cast<double>(occu_DNAp_RNArecog[i]) / n_DNAp;
      freq_DNAp_DNArecog[i] =
          static_cast<double>(occu_DNAp_DNArecog[i]) / n_DNAp;
      freq_DNAp_fold[i] = static_cast<double>(occu_DNAp_fold[i]) / n_DNAp;
    }
  }

  display_p->put_histogram(HIST_RNAp_RNArecog, freq_RNAp_RNArecog,
                           Para::hist_max_y, COLOR_HIST_BACK, COLOR_BLACK,
                           COLOR_HIST_BARS);
  display_p->put_histogram(HIST_RNAp_DNArecog, freq_RNAp_DNArecog,
                           Para::hist_max_y, COLOR_HIST_BACK, COLOR_BLACK,
                           COLOR_HIST_BARS);
  display_p->put_histogram(HIST_RNAp_fold, freq_RNAp_fold, Para::hist_max_y,
                           COLOR_HIST_BACK, COLOR_BLACK, color_of_gradient);
  display_p->put_histogram(HIST_DNAp_RNArecog, freq_DNAp_RNArecog,
                           Para::hist_max_y, COLOR_HIST_BACK, COLOR_BLACK,
                           COLOR_HIST_BARS);
  display_p->put_histogram(HIST_DNAp_DNArecog, freq_DNAp_DNArecog,
                           Para::hist_max_y, COLOR_HIST_BACK, COLOR_BLACK,
                           COLOR_HIST_BARS);
  display_p->put_histogram(HIST_DNAp_fold, freq_DNAp_fold, Para::hist_max_y,
                           COLOR_HIST_BACK, COLOR_BLACK, color_of_gradient);
}

bool MyCA::output(const long Time) {
  if (!plot_fout.is_open()) {
    plot_fout.open(Para::plot_fname.c_str(), std::ios_base::app);
    if (!plot_fout.is_open()) {
      std::cerr << "MyCA::output() Error in opening file" << std::endl;
      exit(-1);
    }
    plot_fout << "# time(1) RNAp_RNA(2) RNAp_DNA(3) DNAp_RNA(4) DNAp_DNA(5) "
                 "para_RNA(6) para_DNA(7) RNAp_RNArecog(8) RNAp_DNArecog(9) "
                 "RNAp_fold(10) DNAp_RNArecog(11) DNAp_DNArecog(12) "
                 "DNAp_fold(13) para_fold(14) vesicle_pop(15) mean_vol(16) "
                 "mean_target_vol(17) mean_repli_density(18)"
              << std::endl;
  }

  unsigned rna_poly_rna = 0, rna_poly_dna = 0, dna_poly_rna = 0,
           dna_poly_dna = 0, paras_rna = 0, paras_dna = 0;
  double rna_poly_rna_recog_param = 0., rna_poly_dna_recog_param = 0.,
         rna_poly_fold_rate = 0.;
  double dna_poly_rna_recog_param = 0., dna_poly_dna_recog_param = 0.,
         dna_poly_fold_rate = 0.;
  double paras_fold_rate = 0.;

  vesicle_observer.reset();

  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (!plane.cell(row, col).molecule.is_empty()) {
        /* Statistics of vesicles (If this part becomes complicated,
           consider to marge this part and the next part that processes
           the statistics of replicators for seep up.) */
        call_vesicle_observer(row, col);

        /* Statistics of replicators */
        if (plane.cell(row, col).molecule.is_replicase_rna()) {
          if (plane.cell(row, col).molecule.rate_dna_repli() == 0.) {
            ++rna_poly_rna;
            rna_poly_rna_recog_param +=
                plane.cell(row, col).molecule.get_rna_recog_param();
            rna_poly_dna_recog_param +=
                plane.cell(row, col).molecule.get_dna_recog_param();
            rna_poly_fold_rate += plane.cell(row, col).molecule.rate_fold();
          } else {
            ++dna_poly_rna;
            dna_poly_rna_recog_param +=
                plane.cell(row, col).molecule.get_rna_recog_param();
            dna_poly_dna_recog_param +=
                plane.cell(row, col).molecule.get_dna_recog_param();
            dna_poly_fold_rate += plane.cell(row, col).molecule.rate_fold();
          }
        } else if (plane.cell(row, col).molecule.is_replicase_dna()) {
          if (plane.cell(row, col).molecule.rate_dna_repli() == 0.) {
            ++rna_poly_dna;
            rna_poly_rna_recog_param +=
                plane.cell(row, col).molecule.get_rna_recog_param();
            rna_poly_dna_recog_param +=
                plane.cell(row, col).molecule.get_dna_recog_param();
            rna_poly_fold_rate += plane.cell(row, col).molecule.rate_fold();
          } else {
            ++dna_poly_dna;
            dna_poly_rna_recog_param +=
                plane.cell(row, col).molecule.get_rna_recog_param();
            dna_poly_dna_recog_param +=
                plane.cell(row, col).molecule.get_dna_recog_param();
            dna_poly_fold_rate += plane.cell(row, col).molecule.rate_fold();
          }
        } else if (plane.cell(row, col).molecule.is_parasite_rna()) {
          ++paras_rna;
          paras_fold_rate += plane.cell(row, col).molecule.rate_fold();
        } else if (plane.cell(row, col).molecule.is_parasite_dna()) {
          ++paras_dna;
          paras_fold_rate += plane.cell(row, col).molecule.rate_fold();
        }
      }
    }

  /* Calculate mean values for replicators */
  if (rna_poly_rna + rna_poly_dna > 0) {
    rna_poly_rna_recog_param /= (rna_poly_rna + rna_poly_dna);
    rna_poly_dna_recog_param /= (rna_poly_rna + rna_poly_dna);
    rna_poly_fold_rate /= (rna_poly_rna + rna_poly_dna);
  }

  if (dna_poly_rna + dna_poly_dna > 0) {
    dna_poly_rna_recog_param /= (dna_poly_rna + dna_poly_dna);
    dna_poly_dna_recog_param /= (dna_poly_rna + dna_poly_dna);
    dna_poly_fold_rate /= (dna_poly_rna + dna_poly_dna);
  }

  if (paras_rna + paras_dna > 0) {
    paras_fold_rate /= (paras_rna + paras_dna);
  }

  /* Calculate mean for vesicles */
  vesicle_observer.calculate_mean();
  double mean_repli_density = 0.;
  if (vesicle_observer.get_population_size() > 0) {
    mean_repli_density =
        static_cast<double>(vesicle_observer.get_mean_n_total()) /
        vesicle_observer.get_mean_volume();
  }

  plot_fout << Time << ' ' << rna_poly_rna << ' ' << rna_poly_dna << ' '
            << dna_poly_rna << ' ' << dna_poly_dna << ' ' << paras_rna << ' '
            << paras_dna << ' ' << rna_poly_rna_recog_param << ' '
            << rna_poly_dna_recog_param << ' ' << rna_poly_fold_rate << ' '
            << dna_poly_rna_recog_param << ' ' << dna_poly_dna_recog_param
            << ' ' << dna_poly_fold_rate << ' ' << paras_fold_rate << ' '
            << vesicle_observer.get_population_size() << ' '
            << vesicle_observer.get_mean_volume() << ' '
            << vesicle_observer.get_mean_target_volume() << ' '
            << mean_repli_density << ' ' << std::endl;

  if (rna_poly_rna + dna_poly_rna == 0) {
    return false;
  } else {
    return true;
  }
}

/* The method collect data of vesicles. Note that this method must be
   called for CA squares that contain a molecule. is_empty() is checked
   in MyCA::output() */
void MyCA::call_vesicle_observer(const unsigned row, const unsigned col) {
  const int v_ind = plane.cell(row, col).vesicle_index;
  if (v_ind != -1) {
    /* Wether this vesicle has a molecule or not, I set this
       vesicle as is_existing=true in vesicle_observer as long
       as it exists in vesicle planne. */
    if (!vesicle_observer.get_is_existing(v_ind)) {
      vesicle_observer.set_is_existing(v_ind);
      vesicle_observer.set_volume(v_ind, vesicle_get_volume(v_ind));
      vesicle_observer.set_target_volume(v_ind,
                                         vesicle_get_target_volume(v_ind));
    }

    vesicle_observer.add_n_total(v_ind);

    /* Collect data from this vesicle */
    //   if(plane.cell(row,col).molecule.is_replicase_rna()){
    //     if(Para::rna_dna_repl_constrain == 0){
    //       if(plane.cell(row,col).molecule.rate_dna_repli() == 0.){
    //         vesicle_observer.add_n_RNAp_RNA(v_ind);
    //       }
    //       else{
    //         vesicle_observer.add_n_DNAp_RNA(v_ind);
    //       }
    //     }
    //     else{
    //       std::cerr << "MyCA::call_vesicle_observer() Error.
    //       Para::rna_dna_repl_constrain is not 0, for which this method
    //       failss." << std::endl; exit(-1);
    //     }
    //   }
    //   else if(plane.cell(row,col).molecule.is_replicase_dna()){
    //     if(Para::rna_dna_repl_constrain == 0){
    //       if(plane.cell(row,col).molecule.rate_dna_repli() == 0.){
    //         vesicle_observer.add_n_RNAp_DNA(v_ind);
    //       }
    //       else{
    //         vesicle_observer.add_n_DNAp_DNA(v_ind);
    //       }
    //     }
    //   }
    //   else if(plane.cell(row,col).molecule.is_parasite_rna()){
    //     vesicle_observer.add_n_para_RNA(v_ind);
    //   }
    //   else if(plane.cell(row,col).molecule.is_parasite_dna()){
    //     vesicle_observer.add_n_para_DNA(v_ind);
    //   }
    // }
  }
}

void MyCA::button2(int x, int y) {
  /* convert window's xy-coordinate into a simulation's
     (row,col) */
  unsigned panel_ind;
  int row, col;
  const bool is_inside_panel =
      display_p->xy_window_to_rc_panel(x, y, panel_ind, row, col);

  if (is_inside_panel) {
    if (panel_ind == 0) {
      if (plane.cell(row, col).vesicle_index != -1) {
        std::cout
            << "ind=" << plane.cell(row, col).vesicle_index << " vol="
            << vesicle_observer.get_volume(plane.cell(row, col).vesicle_index)
            << " target="
            << vesicle_observer.get_target_volume(
                   plane.cell(row, col).vesicle_index)
            << " tot_mol="
            << vesicle_observer.get_n_total(plane.cell(row, col).vesicle_index)
            << std::endl;
      } else {
        std::cout << "row=" << row << " col=" << col << " no-vesicle"
                  << std::endl;
      }
    } else {
      std::cout << "MyCA::button2() this is panel_ind=" << panel_ind
                << std::endl;
    }
  } else {
    std::cout << "MyCA::button2() this coordinate is out of panels"
              << std::endl;
  }
}

void MyCA::button3(int x, int y) {
  std::cout << "MyCA::button3() is not functional" << std::endl;
}

unsigned MyCA::get_bin(const unsigned n_bins, const double min,
                       const double max, const double value) {
  const double delta = (max - min) / n_bins;
  if (value >= min && value < max) {
    return static_cast<unsigned>((value - min) / delta);
  } else if (value >= max) {
    return n_bins - 1;
  } else {
    std::cout << value << std::endl;
    return 0;
  }
}

void MyCA::save(const char *fname) const {
  if (Para::model_type == Para::WELL_MIXED_EQUIL ||
      Para::model_type == Para::WELL_MIXED_EQUIL_NO_COMPLEX ||
      Para::model_type == Para::PARALLEL_WELL_MIXED_EQUIL ||
      Para::model_type == Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX) {
    std::cerr << "MyCA::save() For  model_type=WELL_MIXED_EQUIL, saved data "
                 "doesn't work because we don't output bon_row & bon_col."
              << std::endl;
    return;
  }

  std::ofstream fout(fname);
  if (!fout.is_open()) {
    std::cerr << "MyCA::save() Error. Cannot open " << fname << std::endl;
    exit(-1);
  }

  /* volume target_volume molecule */
  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (plane.cell(row, col).vesicle_index != -1) {
        fout << plane.cell(row, col).vesicle_index << ' '
             << vesicle_get_volume(plane.cell(row, col).vesicle_index) << ' '
             << vesicle_get_target_volume(plane.cell(row, col).vesicle_index);
      } else {
        fout << "-1 0 0";
      }
      fout << ' ' << plane.cell(row, col).molecule << '\n';
    }
}

void MyCA::save_wrapper(const long Time) const {
  std::string ofname;
  std::string command("mkdir -p ");
  command += Para::save_directory_name;
  if (system(command.c_str()) == -1) {
    std::cerr << "MyCA::save_wrapper() Failed to mkdir "
              << Para::save_directory_name
              << ". So, save in the current directory" << std::endl;
    ofname = "/t" + boost::lexical_cast<std::string>(Time) + ".sav";

  } else {
    ofname = Para::save_directory_name;
    ofname += "/t";
    ofname += boost::lexical_cast<std::string>(Time);
    ofname += ".sav";
  }

  save(ofname.c_str());
}

/* This method reads ".sav" file to initilize a simulation. The size of
   CA in a file can be different from the size of the CA of the model.

   HOWEVER, NOTE THAT this function assumes that the boundary of CA is
   fixed.

   I cannot completely check the consistency of a save file in
   read(). Call check_consistency() after input */
void MyCA::read(const std::string &fname) {
  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    std::cerr << "MyCA::read() Error. Cannot open " << fname << std::endl;
    exit(-1);
  }

  unsigned f_nrow = Para::file_nrow;
  unsigned f_ncol = Para::file_ncol;

  unsigned max_nrow = (f_nrow > nrow) ? f_nrow : nrow;
  unsigned max_ncol = (f_ncol > ncol) ? f_ncol : ncol;

  for (unsigned row = 1; row <= max_nrow; ++row) {
    for (unsigned col = 1; col <= max_ncol; ++col) {
      if (row <= nrow && col <= ncol && row <= f_nrow && col <= f_ncol) {
        read_plane_one_line(fin, row, col);

        if (row == nrow && plane.cell(row, col).molecule.is_complex()) {
          switch (plane.cell(row, col).molecule.bon_nei) {
          case 4:
          case 7:
          case 8:
            plane.cell(row, col).molecule.dissociate();
            break;
          default:
            break;
          }
        }

        if (col == ncol && plane.cell(row, col).molecule.is_complex()) {
          switch (plane.cell(row, col).molecule.bon_nei) {
          case 3:
          case 6:
          case 8:
            plane.cell(row, col).molecule.dissociate();
            break;
          default:
            break;
          }
        }

        /* If this is a well-mixed model, we have to convert bon_nei to
           bon_row & bon_col. */
        if (Para::model_type == Para::WELL_MIXED_EQUIL) {
          unsigned bon_row = 0, bon_col = 0;
          plane.XY_NEIGH_X(row, col, plane.cell(row, col).molecule.bon_nei,
                           bon_row, bon_col);
          plane.cell(row, col).molecule.bon_row = bon_row;
          plane.cell(row, col).molecule.bon_col = bon_col;
        }

      } else if ((row > nrow || col > ncol) && row <= f_nrow && col <= f_ncol) {
        while (fin.get() != '\n')
          ;
      } else if (row <= nrow && col <= ncol && (row > f_nrow || col > f_ncol)) {
        /* Do nothing */
        ;
      } else {
        std::cerr << "MyCA::read() There is a bug." << std::endl;
        exit(-1);
      }
    }
  }

  if (fin.peek() != std::char_traits<char>::eof()) {
    std::cerr << "MyCA::read() Error: there is more to read in " << fname
              << ". Did you set -FileNrow/Ncol properly?" << std::endl;
    exit(-1);
  }

  if (fin.fail()) {
    std::cerr << "MyCA::read() Error in reading " << fname << std::endl;
    exit(-1);
  }
}

/* The state of the input file stream is not checked
   below. This should be checked by the method that calls this
   method. */
void MyCA::read_plane_one_line(std::ifstream &fin, const int row, const int col,
                               const int vesicle_index_offset) {
  /*** VESICLE ***/
  read_only_vesicle_one_line(fin, row, col, vesicle_index_offset);

  /*** MOLECULE ****/
  fin >> plane.cell(row, col).molecule;

  /* throw away \n */
  while (fin.get() != '\n')
    ;
}

void MyCA::read_only_vesicle_one_line(std::ifstream &fin, const int row,
                                      const int col,
                                      const int vesicle_index_offset) {
  int vesicle_index, volume, target_volume;

  /*** VESICLE_INDEX ***/
  fin >> vesicle_index;

  /* error check */
  if (vesicle_index >= static_cast<int>(nrow * ncol) || vesicle_index < -2) {
    std::cerr << "MyCA::read_only_vesicle_plane_one_line() Error, invalid "
                 "vesicle_index="
              << vesicle_index << " nrow=" << nrow << " ncol=" << ncol
              << std::endl;
    exit(-1);
  }

  if (vesicle_index != -1 && vesicle_index_offset != 0)
    vesicle_index += vesicle_index_offset;

  /* input */
  plane.cell(row, col).vesicle_index = vesicle_index;

  /*** VOLUME & TARGET_VOLUME ***/
  fin >> volume >> target_volume;

  if (vesicle_index != -1) {
    /* error check */
    if (volume < 0 || volume > static_cast<int>(nrow * ncol) ||
        target_volume < 0) {
      std::cerr
          << "MyCA::read_only_vesicle_plane_one_line() Error, vesicle_index="
          << vesicle_index << " has invalid volume=" << volume
          << " or invalid target_volume=" << target_volume << std::endl;
      exit(-1);
    }

    /* if this is the first time we see this vesicle index,
       make this vesicle. */
    if (vesicle_get_volume(vesicle_index) == -1) {
      vesicle_make_new_vesicle(vesicle_index);
    }
    /* inclement the volume of this vesicle
       (make_new_vesicle() will set the volume to zero.) */
    vesicle_add_volume(vesicle_index);

    /* division by zero has been checked */
    if (vesicle_get_volume(vesicle_index) < volume) {
      /* It might appear that the following (inclemental increase of
         target_volume) is useless, it is actually important when we
         want read a file the size of which is greater than the system
         size (for some vesicles, only part of their whole volume will
         be inside of the system). */
      vesicle_set_target_volume(
          vesicle_index,
          static_cast<int>(
              static_cast<double>(vesicle_get_volume(vesicle_index) *
                                  target_volume) /
                  volume +
              0.5));
    } else if (vesicle_get_volume(vesicle_index) == volume) {
      vesicle_set_target_volume(vesicle_index, target_volume);
    } else {
      std::cerr << "MyCA::read_only_vesicle_plane_one_line() Error, the actual "
                   "volume [="
                << vesicle_get_volume(vesicle_index)
                << "] is greather than the volume in the saved file [="
                << volume << "]" << std::endl;
      exit(-1);
    }
  }
}

void MyCA::read_only_vesicle(const std::string &fname) {
  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    std::cerr << "MyCA::read() Error. Cannot open " << fname << std::endl;
    exit(-1);
  }

  unsigned f_nrow = Para::file_nrow;
  unsigned f_ncol = Para::file_ncol;

  for (unsigned row = 1; row <= f_nrow; ++row) {
    for (unsigned col = 1; col <= f_ncol; ++col) {

      if (row <= nrow && col <= ncol) {
        read_only_vesicle_one_line(fin, row, col);
      }

      /* Discard the rest of the file. */
      while (fin.get() != '\n')
        ;
    }
  }

  if (fin.fail()) {
    std::cerr << "MyCA::read() Error in reading " << fname << std::endl;
    exit(-1);
  }
}

/* This method reads two ".sav" files. The right half and left half of
   the CA will be initialized by each file. It will also color molecules
   according to which which file they come from. This can be displayed
   by using Para::switch_color_code_1st_panel.

   NOTE THAT this function assumes that the boundary of CA is fixed. */
void MyCA::read_half_half(const std::string &fname1,
                          const std::string &fname2) {
  std::ifstream fin1(fname1.c_str());
  if (!fin1.is_open()) {
    std::cerr << "MyCA::read_half_half() Error. Cannot open " << fname1
              << std::endl;
    exit(-1);
  }

  std::ifstream fin2(fname2.c_str());
  if (!fin2.is_open()) {
    std::cerr << "MyCA::read_half_half() Error. Cannot open " << fname2
              << std::endl;
    exit(-1);
  }

  /* volume target_volume molecule */
  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (col <= ncol / 2) {
        read_plane_one_line(fin1, row, col);
        while (fin2.get() != '\n')
          ;
      } else {
        while (fin1.get() != '\n')
          ;
        read_plane_one_line(fin2, row, col, nrow * ncol / 2);
      }

      /* complex molecules in the middle must be uncomplexed */
      if (col == ncol / 2 && plane.cell(row, col).molecule.is_complex()) {
        /* if the partner is in the right side, it would not
           exist. So, uncomplex this molecule. */
        switch (plane.cell(row, col).molecule.bon_nei) {
        case 3:
        case 6:
        case 8:
          plane.cell(row, col).molecule.dissociate();
        default:
          break;
        }
      } else if (col == ncol / 2 + 1 &&
                 plane.cell(row, col).molecule.is_complex()) {
        switch (plane.cell(row, col).molecule.bon_nei) {
        case 2:
        case 5:
        case 7:
          plane.cell(row, col).molecule.dissociate();
        default:
          break;
        }
      }
    }
  if (fin1.fail()) {
    std::cerr << "MyCA::read_half_half() Error in reading " << fname1
              << std::endl;
    exit(-1);
  }

  if (fin2.fail()) {
    std::cerr << "MyCA::read_half_half() Error in reading " << fname2
              << std::endl;
    exit(-1);
  }
}

void MyCA::set_color_half_half() {
  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (col <= ncol / 2) {
        /* color of the left half */
        if (!plane.cell(row, col).molecule.is_empty()) {
          plane.cell(row, col).molecule.set_color(COLOR_RNAp);
        }
      } else {
        /* color of the right half */
        if (!plane.cell(row, col).molecule.is_empty()) {
          plane.cell(row, col).molecule.set_color(COLOR_PARAS);
        }
      }
    }
}

void MyCA::set_color_RNA_DNA() {
  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      if (!plane.cell(row, col).molecule.is_empty()) {
        if (plane.cell(row, col).molecule.is_rna()) {
          plane.cell(row, col).molecule.set_color(COLOR_RNAp);
        } else if (plane.cell(row, col).molecule.is_dna()) {
          plane.cell(row, col).molecule.set_color(COLOR_PARAS);
        } else {
          std::cerr << "set_color_RNA_DNA() Error, there is a molecule that is "
                       "neither RNA nor DNA."
                    << std::endl;
          exit(-1);
        }
      }
    }
}

#ifdef _COUNT_RNA_DNA
void MyCA::output_count_rna_dna(const long Time) {
  std::string ofname;
  std::string command("mkdir -p ");
  command += Para::count_rna_dna_directory_name;
  if (system(command.c_str()) == -1) {
    std::cerr << "MyCA::output_count_rna_dna() Failed to mkdir "
              << Para::count_rna_dna_directory_name
              << ". So, save in the current directory" << std::endl;
    ofname = "/t" + boost::lexical_cast<std::string>(Time) + ".dat";

  } else {
    ofname = Para::count_rna_dna_directory_name;
    ofname += "/t";
    ofname += boost::lexical_cast<std::string>(Time);
    ofname += ".dat";
  }

  std::ofstream fout(ofname.c_str());
  if (!fout.is_open()) {
    std::cerr << "MyCA::output_count_rna_dna() Error. Cannot open " << ofname
              << std::endl;
    exit(-1);
  }

  /* Output the header */
  fout << "#polymerase_type(1) template_type(2) RNA_recog(3) DNA_recog(4) "
          "fold(5) RNA_count(6) DNA_count(7)"
       << std::endl;

  for (unsigned row = 1; row <= nrow; ++row)
    for (unsigned col = 1; col <= ncol; ++col) {
      plane.cell(row, col).molecule.output_rna_dna_count(fout);
    }
}

#endif // _COUNT_RNA_DNA
