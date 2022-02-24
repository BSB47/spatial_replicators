#include "automaton.hpp"

bool Molecule::validate_moltype(const int i)
{
  MolType _i = static_cast<MolType>(i);
  switch(_i)
    {
    case EMPTY: 
    case SIMPL_TEMPL:
    case SIMPL_FOLDE:
    case COMPL_ACTIV:
    case COMPL_PASSI:
    case BOUND:
      return true;
    default: return false;
    }
}

bool Molecule::validate_activity(const int i)
{
  const Activity _i = static_cast<Activity>(i);
  switch(_i)
    {
    case REPLI_RNA:
    case REPLI_DNA: 
    case PARAS_RNA:
    case PARAS_DNA:
    case JUNK_RNA:
    case JUNK_DNA:
      return true;
    default: return false;
    }
}

std::ostream& operator<<(std::ostream& s,const Molecule& mol)
{
  /* EMPTY, SIMPL_TEMPL etc. */
  std::streamsize size = s.precision();
  s.precision(10);
  s << mol.mol_type;

  /* If not empty, */
  if(!mol.is_empty()){
    /* Only when the model assumes complex formation, there can be
       complex molecules. */
    if(mol.is_complex()){
      s  << ' ' << mol.bon_nei 
	/* dissociation rate is bind_param * (1-rna/dna_regoc_param). In
	   the model with complex formation, bind_param is set to
	   max_asso_diss_rate, which is scaled by rate_scale in
	   main(). Thus, I renormalize disso_rate. */
	 << ' ' << mol.disso_rate/Para::rate_scale
	/* replication arte is also normalized by rate_scale. So, I
	   remormalize it too. */
	 << ' ' << mol.repli_rate/Para::rate_scale;
    }else{
      s << " 0 0 0";
    }

    s << ' ' << mol.activity
      << ' ' << mol.decay_rate/Para::rate_scale;

    if(mol.is_replicase()){
      /* In the model with complex formaton, bind_param is set to
	 max_asso_diss, which is scaled by rate_scale in main(). In the
	 model without complex formation, bind_param is simply not
	 used. So, it doesn't matter what I do with it. So, I always
	 renormalize it with rate_scale whether the model assumes
	 complex formation or not.  */
      s << ' ' << mol.bind_param/Para::rate_scale
	<< ' ' << mol.rna_recog_param
	<< ' ' << mol.dna_recog_param
	<< ' ' << mol.rna_repli_rate/Para::rate_scale
	<< ' ' << mol.dna_repli_rate/Para::rate_scale;
    }
    else{
      s << " 0 0 0 0 0";
    }

    s << ' ' << mol.fold_rate
      << ' ' << mol.unfold_rate;

    if(mol.is_parasite()){
      s << ' ' << mol.lipid_synth_rate/Para::rate_scale;
    }
    else{
      s << " 0";
    }
  }
  else{
    s << " 0 0 0 0 0 0 0 0 0 0 0 0 0";
  }

  /* If we are doing the simple ancestor tracing, we also output the
     color of molecules for the porpose of later analysis. */
  if(Para::switch_color_code_1st_panel==2){
      if(!mol.is_empty()){
	s << ' ' << mol.color;
      }
      else{
	s << " 0";
      }
  }

  s.precision(size);
  return s;
}

/* This is the new imput method. When reading a new save file, comment
   out the older one when using it; and vice versa. */
std::istream& operator>>(std::istream& s,Molecule& mol)
{
  /*** MOL_TYPE ***/
  int i;
  s >> i;
  if(!mol.validate_moltype(i)){   /* Error check */
    std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, unknown mol_type=" << i << std::endl;
    exit(-1);
  }
  mol.mol_type = static_cast<Molecule::MolType>(i);

  /*** BON_NEI ***/
  s >> mol.bon_nei;

  /*** DISSO_RATE ***/
  s >> mol.disso_rate;
  mol.disso_rate *= Para::rate_scale;
  
  /*** REPLI_RATE ***/
  s >> mol.repli_rate;
  mol.repli_rate *= Para::rate_scale;

  /*** ACTIVITY_TYPE ***/
  s >> i;
  if(!mol.validate_activity(i)){
    std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, unknown activity_type=" << i << std::endl;
    exit(-1);
  }
  mol.activity = static_cast<Molecule::Activity>(i);

  /*** DECAY_RATE ***/
  s >> mol.decay_rate;
  mol.decay_rate *= Para::rate_scale;

  /*** BIND_PARAM ***/
  s >> mol.bind_param;
  mol.bind_param *= Para::rate_scale;

  /*** RNA_RECOG_PARAM ***/
  s >> mol.rna_recog_param;

  /*** DNA_RECOG_PARAM ***/
  s >> mol.dna_recog_param;

  /*** RNA_REPLI_RATE ***/
  s >> mol.rna_repli_rate;
  mol.rna_repli_rate *= Para::rate_scale;

  /*** DNA_REPLI_RATE ***/
  s >> mol.dna_repli_rate;
  mol.dna_repli_rate *= Para::rate_scale;

  /*** FOLD_RATE ***/
  s >> mol.fold_rate;

  /*** UNFOLD_RATE ***/
  s >>  mol.unfold_rate;

  /*** LIPID_SYNTH_RATE ***/
  s >> mol.lipid_synth_rate;
  mol.lipid_synth_rate *= Para::rate_scale;

  /*** Error check ***/
  if(mol.is_complex()){
    /*** BON_NEI ***/
    if(mol.bon_nei<0 || mol.bon_nei>8){
      std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range bon_nei=" << mol.bon_nei << std::endl;
      exit(-1);
    }
    
    /*** DISSO_RATE ***/
    /* I cannot completely check disso_rate here because it depends on the
       partner molecule  */
    if(mol.disso_rate<0. || mol.disso_rate>Para::max_asso_diss_rate){
      std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range disso_rate=" << mol.disso_rate << std::endl;
      exit(-1);
    }
    
    /*** REPLI_RATE ***/
    /* Currently, I assume that replicase is either RNAp or DNAp (not
       both), and the replication rate does not evolve. Thus, repli_rate
       should be either of replication rates. */
    if(Para::rna_dna_repl_constrain==0 ||
       Para::rna_dna_repl_constrain==1){
      if(mol.repli_rate != 0. && !Molecule::relative_close(mol.repli_rate,Para::max_rna_repli_rate) && !Molecule::relative_close(mol.repli_rate,Para::max_dna_repli_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range repli_rate=" << mol.repli_rate << std::endl;
	exit(-1);
      }
    }
    else{
      std::cerr <<  "operator>>(std::istream& s,Molecule& mol) Error, unknown Para::rna_dna_repl_constrain" << std::endl;
    }
  }

  if(!mol.is_empty() && !mol.is_boundary()){
    /*** DECAY_RATE ***/
    if(mol.is_dna()){
      if(!Molecule::relative_close(mol.decay_rate,Para::dna_decay_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, decay_rate [=" << mol.decay_rate << "] is not the same as dna_decay_rate=" << Para::dna_decay_rate << std::endl;
	exit(-1);
      }
    }
    else if(mol.is_rna()){
      if(!Molecule::relative_close(mol.decay_rate,Para::rna_decay_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, decay_rate [=" << mol.decay_rate << "] is not the same as rna_decay_rate=" << Para::rna_decay_rate << std::endl;
	exit(-1);
      }
    }

    /*** BIND_PARAM ***/
    /* Currently, I always set bind_param==Para::max_asso_diss_rate. If
       this is not the case, rate_scale may be set wrongly. */
    if(mol.is_replicase()){
      if(!Molecule::relative_close(mol.bind_param,Para::max_asso_diss_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, bind_param != Para::max_asso_diss_rate even though this molecule is replicase." << std::endl;
	exit(-1);
      }
    }
    else{
      if(mol.bind_param != 0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, bind_param != 0 even though this molecule is not replicase." << std::endl;
	exit(-1);
      }
    }

    /*** RNA_RECOG_PARAM ***/
    if(mol.is_replicase()){
      if(mol.rna_recog_param < 0. || mol.rna_recog_param > 1.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range rna_recog_param=" << mol.rna_recog_param << std::endl;
	exit(-1);
      }
    }
    else{
      if(mol.rna_recog_param != 0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna_recog_param [=" << mol.rna_recog_param << "] is not zero even though this molecule is not replicase." << std::endl;
	exit(-1);
      }
    }

    /*** DNA_RECOG_PARAM ***/
    if(mol.is_replicase()){
      if(mol.dna_recog_param < 0. || mol.dna_recog_param > 1.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range dna_recog_param=" << mol.dna_recog_param << std::endl;
	exit(-1);
      }
    }
    else{
      if(mol.dna_recog_param != 0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_recog_param [=" << mol.dna_recog_param << "] is not zero even though this molecule is not replicase." << std::endl;
	exit(-1);
      }
    }

    /*** RNA_REPLI_RATE & DNA_REPLI_RATE ***/
    if(mol.is_replicase()){
      if(mol.dna_repli_rate == 0. && !Molecule::relative_close(mol.rna_repli_rate,Para::max_rna_repli_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna_repli_rate is not Para::max_rna_repli_rate although dna_repli_rate=0" << std::endl;
	exit(-1);
      }
      else if(mol.rna_repli_rate == 0. && !Molecule::relative_close(mol.dna_repli_rate,Para::max_dna_repli_rate)){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_repli_rate [=" << mol.dna_repli_rate << "] is not Para::max_dna_repli_rate [="<< Para::max_dna_repli_rate <<"] although rna_repli_rate=0" << std::endl;
	std::cout.precision(15);
	std::cout << mol.dna_repli_rate << " " << Para::max_dna_repli_rate << std::endl;
	exit(-1);
      }
      else if(mol.rna_repli_rate > 0. && mol.dna_repli_rate > 0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_repli_rate>0 and  rna_repli_rate>0" << std::endl;
	exit(-1);
      }
    }
    else{
      if(mol.rna_repli_rate!=0. || mol.dna_repli_rate!=0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna/dna_repli_rate is not 0 although this molecule is not replicase" << std::endl;
	exit(-1);
      }
    }

    /*** FOLD_RATE & UNFOLD_RATE ***/
    if(!Molecule::relative_close(mol.fold_rate+mol.unfold_rate,Para::max_fold_unfold_rate)){
      std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, fold_rate+unfold_rate!=Para::max_fold_unfold_rate." << std::endl;
      exit(-1);
    }

    /*** LIPID_SYNTH_RATE ***/
    if(mol.is_parasite()){
      if(mol.lipid_synth_rate < 0. || mol.lipid_synth_rate > Para::max_lipid_synth_rate){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range lipid_synth_rate=" << mol.lipid_synth_rate << std::endl;
	exit(-1);
      }
    }
    else{
      if(mol.lipid_synth_rate != 0.){
	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error,  lipid_synth_rate [=" << mol.lipid_synth_rate << "] is not zero although this is not parasite." << std::endl;
	exit(-1);
      }
    }
  } 

#ifdef _COUNT_RNA_DNA
  /* This molecule generates a new lineage of replicators (i.e. it
     doesn't have an original template from which it is produced). So,
     we zero _rna_coutn & _dna_count and inclement them according to
     whether this molecule is RNA or DNA template.  */
  mol._rna_count = 0;
  mol._dna_count = 0;
  mol._inclement_rna_dna_count();
#endif //_COUNT_RNA_DNA

  return s;
}

/*This is the old method. I need it when I read an old save
   file. Comment out the newer one when using the older one; and vice
   versa. */
// std::istream& operator>>(std::istream& s,Molecule& mol)
// {
//   if(Para::old_rate_scale==0.){
//     std::cerr << "operator>>(std::istream& s, Molecule &mol) Error, -OldRateScale is not set." << std::endl;
//     exit(-1);
//   }

//   /*** MOL_TYPE ***/
//   int i;
//   s >> i;
//   if(!mol.validate_moltype(i)){   /* Error check */
//     std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, unknown mol_type=" << i << std::endl;
//     exit(-1);
//   }
//   mol.mol_type = static_cast<Molecule::MolType>(i);

//   /*** BON_NEI ***/
//   s >> mol.bon_nei;

//   /*** DISSO_RATE ***/
//   s >> mol.disso_rate;
//   mol.disso_rate *= Para::rate_scale/Para::old_rate_scale;
  
//   /*** REPLI_RATE ***/
//   s >> mol.repli_rate;
//   mol.repli_rate *= Para::rate_scale/Para::old_rate_scale;

//   /*** ACTIVITY_TYPE ***/
//   s >> i;
//   if(!mol.validate_activity(i)){
//     std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, unknown activity_type=" << i << std::endl;
//     exit(-1);
//   }
//   mol.activity = static_cast<Molecule::Activity>(i);
  
//   /*** DECAY_RATE ***/
//   s >> mol.decay_rate;
//   mol.decay_rate *= Para::rate_scale/Para::old_rate_scale;

//   /*** BIND_PARAM ***/
//   s >> mol.bind_param;
//   mol.bind_param *= Para::rate_scale/Para::old_rate_scale;

//   /*** RNA_RECOG_PARAM ***/
//   s >> mol.rna_recog_param;

//   /*** DNA_RECOG_PARAM ***/
//   s >> mol.dna_recog_param;

//   /*** RNA_REPLI_RATE ***/
//   s >> mol.rna_repli_rate;
//   mol.rna_repli_rate *= Para::rate_scale/Para::old_rate_scale;

//   /*** DNA_REPLI_RATE ***/
//   s >> mol.dna_repli_rate;
//   mol.dna_repli_rate *= Para::rate_scale/Para::old_rate_scale;

//   /*** FOLD_RATE ***/
//   s >> mol.fold_rate;

//   /*** UNFOLD_RATE ***/
//   s >>  mol.unfold_rate;

//   /*** LIPID_SYNTH_RATE ***/
//   s >> mol.lipid_synth_rate;
//   mol.lipid_synth_rate *= Para::rate_scale/Para::old_rate_scale;

//   /*** Overwrite parameters  ***/
//   if(Para::is_overwriting_parameters){
//     mol.overwrite_parameters();
//   }

//   /* we might want to skip an error check. */
//   // if(Para::is_skipping_error_check_in_read_file){
//   //   return s;
//   // }

//   /*******************
//    *** Error check ***
//    *******************/
//   if(mol.is_complex()){
//     /*** BON_NEI ***/
//     if(mol.bon_nei<0 || mol.bon_nei>8){
//       std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range bon_nei=" << mol.bon_nei << std::endl;
//       exit(-1);
//     }
    
//     /*** DISSO_RATE ***/
//     /* I cannot completely check disso_rate here because it depends on the
//        partner molecule  */
//     if(mol.disso_rate<0. || mol.disso_rate>Para::max_asso_diss_rate){
//       std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range disso_rate=" << mol.disso_rate << std::endl;
//       exit(-1);
//     }
    
//     /*** REPLI_RATE ***/
//     /* Currently, I assume that replicase is either RNAp or DNAp (not
//        both), and the replication rate does not evolve. Thus, repli_rate
//        should be either of replication rates. */
//     if(Para::rna_dna_repl_constrain==0){
//       if(mol.repli_rate != 0. && !Molecule::relative_close(mol.repli_rate,Para::max_rna_repli_rate) && !Molecule::relative_close(mol.repli_rate,Para::max_dna_repli_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range repli_rate=" << mol.repli_rate << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       std::cerr <<  "operator>>(std::istream& s,Molecule& mol) Error, unknown Para::rna_dna_repl_constrain" << std::endl;
//     }
//   }

//   if(!mol.is_empty() && !mol.is_boundary()){
//     /*** DECAY_RATE ***/
//     if(mol.is_dna()){
//       if(!Molecule::relative_close(mol.decay_rate,Para::dna_decay_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, decay_rate [=" << mol.decay_rate << "] is not the same as dna_decay_rate=" << Para::dna_decay_rate << std::endl;
// 	exit(-1);
//       }
//     }
//     else if(mol.is_rna()){
//       if(!Molecule::relative_close(mol.decay_rate,Para::rna_decay_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, decay_rate [=" << mol.decay_rate << "] is not the same as rna_decay_rate=" << Para::rna_decay_rate << std::endl;
// 	exit(-1);
//       }
//     }

//     /*** BIND_PARAM ***/
//     /* Currently, I always set bind_param==Para::max_asso_diss_rate. If
//        this is not the case, rate_scale may be set wrongly. */
//     if(mol.is_replicase()){
//       if(!Molecule::relative_close(mol.bind_param,Para::max_asso_diss_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, bind_param != Para::max_asso_diss_rate even though this molecule is replicase." << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       if(mol.bind_param != 0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, bind_param != 0 even though this molecule is not replicase." << std::endl;
// 	exit(-1);
//       }
//     }

//     /*** RNA_RECOG_PARAM ***/
//     if(mol.is_replicase()){
//       if(mol.rna_recog_param < 0. || mol.rna_recog_param > 1.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range rna_recog_param=" << mol.rna_recog_param << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       if(mol.rna_recog_param != 0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna_recog_param [=" << mol.rna_recog_param << "] is not zero even though this molecule is not replicase." << std::endl;
// 	exit(-1);
//       }
//     }

//     /*** DNA_RECOG_PARAM ***/
//     if(mol.is_replicase()){
//       if(mol.dna_recog_param < 0. || mol.dna_recog_param > 1.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range dna_recog_param=" << mol.dna_recog_param << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       if(mol.dna_recog_param != 0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_recog_param [=" << mol.dna_recog_param << "] is not zero even though this molecule is not replicase." << std::endl;
// 	exit(-1);
//       }
//     }

//     /*** RNA_REPLI_RATE & DNA_REPLI_RATE ***/
//     if(mol.is_replicase()){
//       if(mol.dna_repli_rate == 0. && !Molecule::relative_close(mol.rna_repli_rate,Para::max_rna_repli_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna_repli_rate is not Para::max_rna_repli_rate although dna_repli_rate=0" << std::endl;
// 	exit(-1);
//       }
//       else if(mol.rna_repli_rate == 0. && !Molecule::relative_close(mol.dna_repli_rate,Para::max_dna_repli_rate)){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_repli_rate is not Para::max_dna_repli_rate although rna_repli_rate=0" << std::endl;
// 	std::cout.precision(15);
// 	std::cout << mol.dna_repli_rate << " " << Para::max_dna_repli_rate << std::endl;
// 	exit(-1);
//       }
//       else if(mol.rna_repli_rate > 0. && mol.dna_repli_rate > 0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, dna_repli_rate>0 and  rna_repli_rate>0" << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       if(mol.rna_repli_rate!=0. || mol.dna_repli_rate!=0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, rna/dna_repli_rate is not 0 although this molecule is not replicase" << std::endl;
// 	exit(-1);
//       }
//     }

//     /*** FOLD_RATE & UNFOLD_RATE ***/
//     if(!Molecule::relative_close(mol.fold_rate+mol.unfold_rate,Para::max_fold_unfold_rate)){
//       std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, fold_rate+unfold_rate!=Para::max_fold_unfold_rate." << std::endl;
//       exit(-1);
//     }

//     /*** LIPID_SYNTH_RATE ***/
//     if(mol.is_parasite()){
//       if(mol.lipid_synth_rate < 0. || mol.lipid_synth_rate > Para::max_lipid_synth_rate){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error, out of range lipid_synth_rate=" << mol.lipid_synth_rate << std::endl;
// 	exit(-1);
//       }
//     }
//     else{
//       if(mol.lipid_synth_rate != 0.){
// 	std::cerr << "operator>>(std::istream& s,Molecule& mol) Error,  lipid_synth_rate [=" << mol.lipid_synth_rate << "] is not zero although this is not parasite." << std::endl;
// 	exit(-1);
//       }
//     }
//   } 
// #ifdef _COUNT_RNA_DNA
//   /* When a ".sav" file is fed to a program, we consider that each
//      molecule in the file generates a new lineage of replicators (we
//      ignore the original templates from which they are produced because
//      the original templates cannot be recovered from the file). So, we
//      zero _rna_coutn & _dna_count and inclement them according to
//      whether this molecule is RNA or DNA template. */
//   mol._rna_count = 0;
//   mol._dna_count = 0;
//   mol._inclement_rna_dna_count();
// #endif //_COUNT_RNA_DNA
//   return s;
// }

void Molecule::initialize(const std::string& _activity,const double _fold_rate,const double _rna_recog_param,const double _dna_recog_param)
{
  mol_type = SIMPL_TEMPL;

  if(_activity == "RNAp_RNA"){
    activity = REPLI_RNA;
    decay_rate = Para::rna_decay_rate;
    bind_param = Para::max_asso_diss_rate;

    /* rna/dna_recog_param */
    if(Para::rna_dna_recog_constrain==0 ||
       Para::rna_dna_recog_constrain==1 ||
       Para::rna_dna_recog_constrain==2){
      rna_recog_param = _rna_recog_param;
      dna_recog_param = _dna_recog_param;
    }
    else if(Para::rna_dna_recog_constrain==3){
      rna_recog_param = _rna_recog_param;
      dna_recog_param = 0.;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::rna_dna_recog_constrain." << std::endl;
      exit(-1);
    }
    
    /* rna/dna_repl_rate */
    if(Para::rna_dna_repl_constrain==0 ||
       Para::rna_dna_repl_constrain==1){
      rna_repli_rate = Para::max_rna_repli_rate;
      dna_repli_rate = 0.;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::rna_dna_repl_constrain." << std::endl;
      exit(-1);
    }

    /* fold/unforld_rate */
    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      /* Note that the meaning of fold_rate's depend on whether the
	 system is under the equilibrium assumption. If the system is
	 under the equilibrium assumption, fold_rate and unfold_rate
	 represent the fraction of time rather than rates. Moreover,
	 max_fold_unfold_rate is set to 1. */
      fold_rate = _fold_rate;
      unfold_rate = Para::max_fold_unfold_rate - fold_rate;
    }
    else if(Para::fold_unfold_constrain==1){
      fold_rate = 0.;
      unfold_rate = Para::max_fold_unfold_rate;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* lipid synthesis rate not needed lipid_synth_rate = 0.; */
  }
  else if(_activity == "DNAp_RNA"){
    activity = REPLI_RNA;
    decay_rate = Para::rna_decay_rate;
    bind_param = Para::max_asso_diss_rate;

    /* rna/dna_recog_param */
    if(Para::rna_dna_recog_constrain==0 ||
       Para::rna_dna_recog_constrain==3){
      rna_recog_param = _rna_recog_param;
      dna_recog_param = _dna_recog_param;
    }
    else if(Para::rna_dna_recog_constrain==1){
      if(_dna_recog_param!=0.){
	std::cerr << "Molecule::initialize() Error, although Para::rna_dna_recog_constrain = 1, Para::init_DNAp_rna_recog_param != 0." << std::endl;
	exit(-1);
      }
      rna_recog_param = _rna_recog_param;
      dna_recog_param = 0.;
    }
    else if(Para::rna_dna_recog_constrain==2){
      if(_rna_recog_param!=0.){
	std::cerr << "Molecule::initialize() Error, although Para::rna_dna_recog_constrain = 2, Para::init_DNAp_rna_recog_param != 0." << std::endl;
	exit(-1);
      }
      rna_recog_param = 0.;
      dna_recog_param = _dna_recog_param;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::rna_dna_recog_constrain." << std::endl;
      exit(-1);
    }
    
    /* rna/dna_repl_rate */
    if(Para::rna_dna_repl_constrain==0 ||
       Para::rna_dna_repl_constrain==1){
      rna_repli_rate = 0.;
      dna_repli_rate = Para::max_dna_repli_rate;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::rna_dna_repl_constrain." << std::endl;
      exit(-1);
    }

    /* fold/unforld_rate */
    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      /* Note that the meaning of fold_rate's depend on whether the
	 system is under the equilibrium assumption. If the system is
	 under the equilibrium assumption, fold_rate and unfold_rate
	 represent the fraction of time rather than rates. Moreover,
	 max_fold_unfold_rate is set to 1. */
      fold_rate = _fold_rate;
      unfold_rate = Para::max_fold_unfold_rate - fold_rate;
    }
    else if(Para::fold_unfold_constrain==1){
      fold_rate = 0.;
      unfold_rate = Para::max_fold_unfold_rate;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* lipid synthesis rate not needed lipid_synth_rate = 0.; */
  }
  else if(_activity == "PARAS_RNA"){
    activity = PARAS_RNA;
    decay_rate = Para::rna_decay_rate;
    /* Not needed bind_param = 0.; */

    /* rna/dna_recog_param not needed */
    /* rna/dna_repl_rate not needed */
    
    /* fold/unforld_rate */
    if(Para::fold_unfold_constrain==0||
       Para::fold_unfold_constrain==2){
      fold_rate = _fold_rate;
      unfold_rate = Para::max_fold_unfold_rate - fold_rate;
    }
    else if(Para::fold_unfold_constrain==1){
      fold_rate = 0.;
      unfold_rate = Para::max_fold_unfold_rate;
    }
    else{
      std::cerr << "Molecule::initialize() Error, unknown Para::fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* lipid synthesis rate */
    lipid_synth_rate = Para::init_lipid_synth_rate;
  }
  else{
    std::cerr << "Molecule::initialize() Error, unknown activity." << std::endl;
    exit(-1);
  }

#ifdef _COUNT_RNA_DNA
  /* When a ".sav" file is fed to a program, we consider that each
     molecule in the file generates a new lineage of replicators (we
     ignore the original templates from which they are produced because
     the original templates cannot be recovered from the file). So, we
     zero _rna_coutn & _dna_count and inclement them according to
     whether this molecule is RNA or DNA template. */
  _rna_count = 0;
  _dna_count = 0;
  _inclement_rna_dna_count();
#endif //_COUNT_RNA_DNA
}

/* By this method, (*this) will get a copy of either mol1 or mol2 that
   has mol_type==COMPL_PASSI. The values of bon_nei of mol1 and mol2 are
   invariant in this method.  */
void Molecule::replicate(Molecule& mol1,Molecule& mol2)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(mol1.mol_type==COMPL_PASSI || mol2.mol_type==COMPL_PASSI));
  }
  catch(GeneralError){
    std::cerr << "Molecule::replicate() Error, neither mol1 nor mol2 has mol_type==COMPL_PASSI" << std::endl;
    exit(-1);
  }

  /* When I assume that replicases are either RNA polymerase or DNA
     polymerase (no hybrid). Thus, either dna_repli_rate or
     rna_repli_rate must be 0. */
  if(Para::rna_dna_repl_constrain==0 ||
     Para::rna_dna_repl_constrain==1){
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||
			   ((mol2.mol_type==COMPL_ACTIV && (mol2.rna_repli_rate == 0. || mol2.dna_repli_rate == 0.))||
			    (mol1.mol_type==COMPL_ACTIV && (mol1.rna_repli_rate == 0. || mol1.dna_repli_rate == 0.))));
    }
    catch(GeneralError){
      std::cerr << "Molecule::replicate() Error, even though rna_dna_repl_constrain=0, neither rna_repli_rate nor dna_repli_rate is 0." << std::endl;
      exit(-1);
    }

    /* Mol1 is the template. */
    if(mol1.mol_type==COMPL_PASSI){
      copy_with_mutation_new(mol1);
      if(mol2.rna_repli_rate > 0.){
	convert_to_rna();
      }
      else{
	convert_to_dna();
      }
    }
    /* Mol2 is the template. */
    else{
      copy_with_mutation_new(mol2);
      if(mol1.rna_repli_rate > 0.){
	convert_to_rna();
      }
      else{
	convert_to_dna();
      }
    }

    /* The complex must dissociate. */
    mol1.mol_type = mol2.mol_type = SIMPL_TEMPL;

#ifdef _COUNT_RNA_DNA
    /* Every time a new molecule is produced, _rna_count & _dna_count
       are inclemented according to whether the molecule is DNA or RNA
       template. */
    _inclement_rna_dna_count();
#endif //_COUNT_RNA_DNA
  }
  else{
    std::cerr << "Molecule::replicate() Error, unknown Para::rna_dna_repl_constrain." << std::endl;
    exit(-1);
  }
}

/* It replicates a molecule by using "original" as template. The
   produced molecule will be in SIMPL_TEMPL state. */
void Molecule::copy_with_mutation(const Molecule& original)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||mol_type==EMPTY);
  }
  catch(GeneralError){
    std::cerr << "Molecule::copy_with_mutation() Error, I am asked to copy something into a non-empty square" << std::endl;
    std::cerr << "mol_type=" << original.mol_type << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||
			 ((original.mol_type==COMPL_PASSI && (Para::model_type==Para::VESIC_NEUT_EQUIL || Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL || Para::model_type==Para::VESIC_LIPID_EQUIL || Para::model_type==Para::SURFACE_EQUIL || Para::model_type==Para::WELL_MIXED_EQUIL || Para::model_type==Para::PARALLEL_SURFACE_EQUIL || Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL || Para::model_type==Para::PARALLEL_INORG_EQUIL || Para::model_type==Para::INORG_EQUIL || Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL)) ||
			  (original.mol_type==SIMPL_TEMPL && (Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX || Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX || Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX || Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX))));
  }
  catch(GeneralError){
    std::cerr << "Molecule::copy_with_mutation() Error, I am asked to replicate non-CMPL_PASSI molecule" << std::endl;
    std::cerr << "mol_type=" << original.mol_type << std::endl;
    exit(-1);
  }

  /* This is raw copy. */
  *this = original;
  mol_type = SIMPL_TEMPL;
  
  if(activity==REPLI_RNA || activity==REPLI_DNA){
    double p = rand_karney.Real();

    /* Mutation in fold_rate */
    double cum_rate=0.;
    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      cum_rate += Para::mut_rate_fold_unfold_rate;
      if(p < cum_rate){
	/* If this is RNA polymerase */
	if(rna_repli_rate>0.){
	  mutate_parameter(fold_rate,Para::min_fold_rate_RNAp,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	}else{
	  mutate_parameter(fold_rate,0.,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	}
	unfold_rate = Para::max_fold_unfold_rate - fold_rate;
	return;
      }
    }
    else if(Para::fold_unfold_constrain==1){
      /* We don't mutate fold_rate because when
	 Para::fold_unfold_constrain == 1, we don't use fold_rate. */
      ;
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* Mutation to parasites */
    cum_rate += Para::mut_rate_to_para;
    if(p < cum_rate){
      /* Even if activity is originally REPLI_DNA, we don't have to
    	 make it into PARAS_DNA because Molecule::repliate() will take
    	 care of it. Moreover, we don't have to set
    	 rna/dna_recog_param to 0 because rate_assoc_diss() will take
    	 care of it. */
      activity = PARAS_RNA;
      lipid_synth_rate = 0.;
      return;
    }

    /* Mutation to junk molecules */
    // FOR NOW, I TURN THIS OFF
    // cum_rate += Para::mut_rate_to_junk;
    // if(p < cum_rate){
    //   /* Why not JUNK_DNA? See above. */
    //   activity = JUNK_RNA;
    //   return;
    // }

    /* Mutation on RNA/DNA recognition parameters */
    if(Para::rna_dna_recog_constrain == 0){
      cum_rate += Para::mut_rate_recog_param;
      if(p < cum_rate){
	mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	return;
      }
      cum_rate += Para::mut_rate_recog_param;
      if(p < cum_rate){
	mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	return;
      }
    }
    else if(Para::rna_dna_recog_constrain == 1){
      if(Para::rna_dna_repl_constrain == 0 ||
	 Para::rna_dna_repl_constrain == 1){
	/* If this is RNAp, both rna/dna_recog_param are allowed to evolve. */
	if(dna_repli_rate == 0.){
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	/* If this is DNAp, only rna_recog_param is allowed to evolve. */
	else{
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
      }
      else{
	std::cerr << "Molecule::copy_with_mutation() Error: Although Para::rna_dna_recog_constrain==1, Para::rna_dna_repl_constrain is not set to 1. We need this in order to distinguish RNAp and DNAp." << std::endl;
	exit(-1);
      }
    }
    else if(Para::rna_dna_recog_constrain == 2){
      if(Para::rna_dna_repl_constrain == 0 ||
	 Para::rna_dna_repl_constrain == 1){
	/* If this is RNAp, both rna/dna_recog_param are allowed to evolve. */
	if(dna_repli_rate == 0.){
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	/* If this is DNAp, only dna_recog_param is allowed to evolve. */
	else{
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
      }
      else{
	std::cerr << "Molecule::copy_with_mutation() Error: Although Para::rna_dna_recog_constrain==2, Para::rna_dna_repl_constrain is not set to 1. We need this in order to distinguish RNAp and DNAp." << std::endl;
	exit(-1);
      }
    }
    else if(Para::rna_dna_recog_constrain == 3){
      if(Para::rna_dna_repl_constrain == 0 ||
	 Para::rna_dna_repl_constrain == 1){
	/* If this is RNAp, only rna_recog_param is allowed to evolve. */
	if(dna_repli_rate == 0.){
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	/* If this is DNAp, both rna/dna_recog_param and allowed to evolve. */
	else{
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	  cum_rate += Para::mut_rate_recog_param;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
      }
      else{
	std::cerr << "Molecule::copy_with_mutation() Error: Although Para::rna_dna_recog_constrain==3, Para::rna_dna_repl_constrain is not set to 1. We need this in order to distinguish RNAp and DNAp." << std::endl;
	exit(-1);
      }
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::rna_dna_recog_constrain." << std::endl;
      exit(-1);
    }

    /* Mutation on RNA/DNA replication parameters. */
    if(Para::rna_dna_repl_constrain == 0 ||
       Para::rna_dna_repl_constrain == 1){

      cum_rate += Para::mut_rate_rna_dna_repl_rate;
      if(p < cum_rate){
	if(dna_repli_rate == 0.){
	  dna_repli_rate = Para::max_dna_repli_rate;
	  rna_repli_rate = 0.;

	  /* If rna_dna_repl_constrain==1, we don't copy parameters of the
	     original template. */
	  if(Para::rna_dna_repl_constrain == 1){
	    /* Recognition parameters */
	    rna_recog_param = Para::mut_init_rna_recog_param_Dp;
	    dna_recog_param = Para::mut_init_dna_recog_param_Dp;

	    /* Folding/unfolding */
	    if(Para::fold_unfold_constrain==0 ||
	       Para::fold_unfold_constrain==2){
	      fold_rate = Para::mut_init_fold_rate;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	    else if(Para::fold_unfold_constrain==1){
	      /* Don't do anything */
	      ;
	    }
	    else{
	      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	      exit(-1);
	    }
	  }

	  /* If we are disabling DNA->DNA or RNA->DNA replication, we
	     should change recog_para when RNAp mutates into DNAp */
	  if(Para::rna_dna_recog_constrain == 1){
	    dna_recog_param = 0.;	    
	  }
	  else if(Para::rna_dna_recog_constrain == 2){
	    rna_recog_param = 0.;	    
	  }
	}
	else{
	  rna_repli_rate = Para::max_rna_repli_rate;
	  dna_repli_rate = 0.;
	    
	  /* If rna_dna_repl_constrain==1, we don't copy parameters of the
	     original template. */
	  if(Para::rna_dna_repl_constrain == 1){
	    /* Recognition parameters */
	    rna_recog_param = Para::mut_init_rna_recog_param_Rp;
	    dna_recog_param = Para::mut_init_dna_recog_param_Rp;

	    /* Folding/unfolding */
	    if(Para::fold_unfold_constrain==0 ||
	       Para::fold_unfold_constrain==2){
	      fold_rate = Para::mut_init_fold_rate;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	    else if(Para::fold_unfold_constrain==1){
	      /* Don't do anything */
	      ;
	    }
	    else{
	      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	      exit(-1);
	    }
	  }

	  /* If we are disabling DNA->RNA  replication, we
	     should change recog_para when RNAp mutates into DNAp */
	  if(Para::rna_dna_recog_constrain == 3){
	    dna_recog_param = 0.;	    
	  }

	  /* Check if the folding rate is smaller than allowed */
	  if(Para::fold_unfold_constrain==0 ||
	     Para::fold_unfold_constrain==2){
	    if(fold_rate < Para::min_fold_rate_RNAp){
	      fold_rate = Para::min_fold_rate_RNAp;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	  }
	  else if(Para::fold_unfold_constrain==1){
	    /* Don't do anything */
	    ;
	  }
	  else{
	    std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	    exit(-1);
	  }
	}
	return;
      }
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::rna_dna_repl_constrain." << std::endl;
      exit(-1);
    }
  }
  else if(activity == PARAS_RNA || activity == PARAS_DNA){
    double p = rand_karney.Real();
    double cum_rate = 0;

    /* Mutation of fold/unfold rate */
    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      cum_rate += Para::mut_rate_fold_unfold_rate;
      if(p < cum_rate){
	mutate_parameter(fold_rate,0.,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	return;
      }
    }
    else if(Para::fold_unfold_constrain==1){
      /* We don't mutate fold_rate because when
	 Para::fold_unfold_constrain == 1, we don't use fold_rate. */
      ;
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* Mutation of lipid synthesis rate */
    cum_rate += Para::mut_rate_lipid_synth_rate;
    if(p < cum_rate){
      mutate_parameter(lipid_synth_rate,0.,Para::max_lipid_synth_rate,Para::mut_step_lipid_synth_rate);
      return;
    }

    /* Mutation to junk molecules */
    // FOR NOW, I TURN IT OFF.
    // cum_rate += Para::mut_rate_to_junk;
    // if(p < cum_rate){
    // 	/* Why not JUNK_DNA? See above. */
    // 	activity = JUNK_RNA;
    // 	return;
    // }
  }
  /* JUNK do not mutate */
}

/*
  This is an updated version of copy_with_mutation(). It uses separate
  parameters for the mutation of RNA recognition parameter, DNA
  recognition parameter and folding for Rp and Dp. This allows me to run
  a simulation where only a subset of parameters is allowed to evolve.
*/
void Molecule::copy_with_mutation_new(const Molecule& original)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||mol_type==EMPTY);
  }
  catch(GeneralError){
    std::cerr << "Molecule::copy_with_mutation() Error, I am asked to copy something into a non-empty square" << std::endl;
    std::cerr << "mol_type=" << original.mol_type << std::endl;
    exit(-1);
  }

  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||
			 ((original.mol_type==COMPL_PASSI && (Para::model_type==Para::VESIC_NEUT_EQUIL || Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL || Para::model_type==Para::VESIC_LIPID_EQUIL || Para::model_type==Para::SURFACE_EQUIL || Para::model_type==Para::WELL_MIXED_EQUIL || Para::model_type==Para::PARALLEL_SURFACE_EQUIL || Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL || Para::model_type==Para::PARALLEL_INORG_EQUIL || Para::model_type==Para::INORG_EQUIL || Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL)) ||
			  (original.mol_type==SIMPL_TEMPL && (Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX || Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX || Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX || Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX))));
  }
  catch(GeneralError){
    std::cerr << "Molecule::copy_with_mutation() Error, I am asked to replicate non-CMPL_PASSI molecule" << std::endl;
    std::cerr << "mol_type=" << original.mol_type << std::endl;
    exit(-1);
  }

  /* This is raw copy. */
  *this = original;
  mol_type = SIMPL_TEMPL;

  double p = rand_karney.Real();
  double cum_rate=0.;
  if(activity==REPLI_RNA || activity==REPLI_DNA){
    if(Para::rna_dna_repl_constrain == 0 ||
       Para::rna_dna_repl_constrain == 1) {

      /*****************
       * If this is Rp *
       *****************/
      if(rna_repli_rate>0.){
	/*************************
	 * Mutation in fold_rate *
	 *************************/
	if(Para::fold_unfold_constrain==0 ||
	   Para::fold_unfold_constrain==2){
	  cum_rate += Para::mut_rate_fold_unfold_rate_Rp;
	  if(p < cum_rate){
	    mutate_parameter(fold_rate,Para::min_fold_rate_RNAp,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	    unfold_rate = Para::max_fold_unfold_rate - fold_rate;
	    return;
	  }
	}
	else if(Para::fold_unfold_constrain==1){
	  /* We don't mutate fold_rate because when
	     Para::fold_unfold_constrain == 1, we don't use fold_rate. */
	  ;
	}
	else{
	  std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain." << std::endl;
	  exit(-1);
	}

	/*************************
	 * Mutation to parasites *
	 *************************/
	cum_rate += Para::mut_rate_to_para;
	if(p < cum_rate){
	  /* Even if activity is originally REPLI_DNA, we don't have to
	     make it into PARAS_DNA because Molecule::repliate() will take
	     care of it. Moreover, we don't have to set
	     rna/dna_recog_param to 0 because rate_assoc_diss() will take
	     care of it. */
	  activity = PARAS_RNA;
	  lipid_synth_rate = 0.;
	  return;
	}
	
	/**********************************************
	 * Mutation on RNA/DNA recognition parameters *
	 **********************************************/
	if(Para::rna_dna_recog_constrain == 0 ||
	   Para::rna_dna_recog_constrain == 1 ||
	   Para::rna_dna_recog_constrain == 2){
	  cum_rate += Para::mut_rate_RNA_recog_param_Rp;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	  cum_rate += Para::mut_rate_DNA_recog_param_Rp;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	else if(Para::rna_dna_recog_constrain == 3){
	  /* Only rna_recog_param is allowed to evolve. */
	  cum_rate += Para::mut_rate_RNA_recog_param_Rp;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	else{
	  std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::rna_dna_recog_constrain." << std::endl;
	  exit(-1);
	}

	/***********************************************
	 * Mutation on RNA/DNA replication parameters. *
	 ***********************************************/
	cum_rate += Para::mut_rate_rna_dna_repl_rate;
	if(p < cum_rate){
	  dna_repli_rate = Para::max_dna_repli_rate;
	  rna_repli_rate = 0.;
	  
	  /* If rna_dna_repl_constrain==1, we don't copy parameters of the
	     original template. */
	  if(Para::rna_dna_repl_constrain == 1){
	    /* Recognition parameters */
	    rna_recog_param = Para::mut_init_rna_recog_param_Dp;
	    dna_recog_param = Para::mut_init_dna_recog_param_Dp;

	    /* Folding/unfolding */
	    if(Para::fold_unfold_constrain==0 ||
	       Para::fold_unfold_constrain==2){
	      fold_rate = Para::mut_init_fold_rate;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	    else if(Para::fold_unfold_constrain==1){
	      /* Don't do anything */
	      ;
	    }
	    else{
	      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	      exit(-1);
	    }
	  }

	  /* If we are disabling DNA->DNA or RNA->DNA replication, we
	     should change recog_para when RNAp mutates into DNAp */
	  if(Para::rna_dna_recog_constrain == 1){
	    dna_recog_param = 0.;	    
	  }
	  else if(Para::rna_dna_recog_constrain == 2){
	    rna_recog_param = 0.;	    
	  }

	  return;
	}

      }
      /*****************
       * If this is Dp *
       *****************/
      else{
	/*************************
	 * Mutation in fold_rate *
	 **************************/
	if(Para::fold_unfold_constrain==0 ||
	   Para::fold_unfold_constrain==2){
	  cum_rate += Para::mut_rate_fold_unfold_rate_Dp;
	  if(p < cum_rate){
	    mutate_parameter(fold_rate,0.,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	    unfold_rate = Para::max_fold_unfold_rate - fold_rate;
	    return;
	  }
	}
	else if(Para::fold_unfold_constrain==1){
	  /* We don't mutate fold_rate because when
	     Para::fold_unfold_constrain == 1, we don't use fold_rate. */
	  ;
	}
	else{
	  std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain." << std::endl;
	  exit(-1);
	}
	
	/*************************
	 * Mutation to parasites *
	 *************************/
	cum_rate += Para::mut_rate_to_para;
	if(p < cum_rate){
	  activity = PARAS_RNA;
	  lipid_synth_rate = 0.;
	  return;
	}

	/**********************************************
	 * Mutation on RNA/DNA recognition parameters *
	 **********************************************/
	if(Para::rna_dna_recog_constrain == 0 ||
	   Para::rna_dna_recog_constrain == 3 ){
	  cum_rate += Para::mut_rate_RNA_recog_param_Dp;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	  cum_rate += Para::mut_rate_DNA_recog_param_Dp;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	else if(Para::rna_dna_recog_constrain == 1){
	  cum_rate += Para::mut_rate_RNA_recog_param_Dp;
	  if(p < cum_rate){
	    mutate_parameter(rna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	else if(Para::rna_dna_recog_constrain == 2){
	  cum_rate += Para::mut_rate_DNA_recog_param_Dp;
	  if(p < cum_rate){
	    mutate_parameter(dna_recog_param,0.,1.,Para::mut_step_recog_param);
	    return;
	  }
	}
	else{
	  std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::rna_dna_recog_constrain." << std::endl;
	  exit(-1);
	}


	/***********************************************
	 * Mutation on RNA/DNA replication parameters. *
	 ***********************************************/
	cum_rate += Para::mut_rate_rna_dna_repl_rate;
	if(p < cum_rate){
	  rna_repli_rate = Para::max_rna_repli_rate;
	  dna_repli_rate = 0.;
	    
	  /* If rna_dna_repl_constrain==1, we don't copy parameters of the
	     original template. */
	  if(Para::rna_dna_repl_constrain == 1){
	    /* Recognition parameters */
	    rna_recog_param = Para::mut_init_rna_recog_param_Rp;
	    dna_recog_param = Para::mut_init_dna_recog_param_Rp;

	    /* Folding/unfolding */
	    if(Para::fold_unfold_constrain==0 ||
	       Para::fold_unfold_constrain==2){
	      fold_rate = Para::mut_init_fold_rate;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	    else if(Para::fold_unfold_constrain==1){
	      /* Don't do anything */
	      ;
	    }
	    else{
	      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	      exit(-1);
	    }
	  }

	  /* If we are disabling DNA->RNA  replication, we
	     should change recog_para when RNAp mutates into DNAp */
	  if(Para::rna_dna_recog_constrain == 3){
	    dna_recog_param = 0.;	    
	  }

	  /* Check if the folding rate is smaller than allowed */
	  if(Para::fold_unfold_constrain==0 ||
	     Para::fold_unfold_constrain==2){
	    if(fold_rate < Para::min_fold_rate_RNAp){
	      fold_rate = Para::min_fold_rate_RNAp;
	      unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	    }
	  }
	  else if(Para::fold_unfold_constrain==1){
	    /* Don't do anything */
	    ;
	  }
	  else{
	    std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::fold_unfold_constrain" << std::endl;
	    exit(-1);
	  }
	  return;
	}
      }
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown value in Para::rna_dna_repl_constrain." << std::endl;
      exit(-1);
    }
  }
  else if(activity == PARAS_RNA || activity == PARAS_DNA){
    /* Mutation of fold/unfold rate */
    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      cum_rate += Para::mut_rate_fold_unfold_rate;
      if(p < cum_rate){
	mutate_parameter(fold_rate,0.,Para::max_fold_unfold_rate,Para::mut_step_fold_unfold_rate);
	unfold_rate = Para::max_fold_unfold_rate-fold_rate;
	return;
      }
    }
    else if(Para::fold_unfold_constrain==1){
      /* We don't mutate fold_rate because when
	 Para::fold_unfold_constrain == 1, we don't use fold_rate. */
      ;
    }
    else{
      std::cerr << "Molecule::copy_with_mutation() Unknown fold_unfold_constrain." << std::endl;
      exit(-1);
    }

    /* Mutation of lipid synthesis rate */
    cum_rate += Para::mut_rate_lipid_synth_rate;
    if(p < cum_rate){
      mutate_parameter(lipid_synth_rate,0.,Para::max_lipid_synth_rate,Para::mut_step_lipid_synth_rate);
      return;
    }
  }
  /* JUNK do not mutate */
}


/* This calculates the rate of association and dissociation between two
   simple molecules. Note that this method considers *this to be the
   replicase that takes "partner" as template. */
void Molecule::rate_assoc_disso_to(const Molecule& partner,double& assoc,double& disso) const
{
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
    /* Under the equilibrium assumption, we don't use SIMPL_FOLDE. If
       molecules are not in complex, they are always SIMPLE_TEMPL. */
    if(mol_type != SIMPL_TEMPL || partner.mol_type != SIMPL_TEMPL){
      assoc = 0.;
      /* Since association does not happen, I do not set disso and
	 repli. */
      return;
    }

    /* I assume that binding happens only if the catalyst (i.e. *this)
        is REPLI_RNA. This part depends on how we represent molecules */
    if(activity != REPLI_RNA){
      assoc = 0.;
      return;
    }

    if(Para::fold_unfold_constrain==0){
      switch(partner.activity){
	/* If the template is RNA. */
      case REPLI_RNA: case PARAS_RNA:
	assoc = fold_rate * bind_param * rna_recog_param * partner.unfold_rate;
	disso = bind_param * (1. - rna_recog_param);
	return;
	/* If the template is DNA */
      case REPLI_DNA: case PARAS_DNA:
	/* I assume that DNA molecule does not fold. */
	assoc = fold_rate * bind_param * dna_recog_param;
	disso = bind_param * (1. - dna_recog_param);
	return;
      case JUNK_RNA: case JUNK_DNA:
	/* I assume that JUNK molecules don't bind. */
	assoc = 0.;
	return;
      default:
	std::cerr << "Molecule::rate_assoc_disso_to() Error, unknown activity=" << partner.activity << std::endl;
      }
    }
    else if(Para::fold_unfold_constrain==1){
      switch(partner.activity){
	/* If the template is RNA. */
      case REPLI_RNA:
	assoc = bind_param * rna_recog_param;
	disso = bind_param * (1. - rna_recog_param);
	return;
	/* If the template is DNA */
      case REPLI_DNA:
	/* I assume that DNA molecule does not fold. */
	assoc = bind_param * dna_recog_param;
	disso = bind_param * (1. - dna_recog_param);
	return;
	/* If the template is Parasite RNA. */
      case PARAS_RNA:
	assoc = Para::parasite_advantage * bind_param * rna_recog_param;
	disso = bind_param * (1. - rna_recog_param);
	return;
	/* If the template is Parasite DNA */
      case PARAS_DNA:
	/* I assume that DNA molecule does not fold. */
	assoc = Para::parasite_advantage * bind_param * dna_recog_param;
	disso = bind_param * (1. - dna_recog_param);
	return;
      case JUNK_RNA: case JUNK_DNA:
	/* I assume that JUNK molecules don't bind. */
	assoc = 0.;
	return;
      default:
	std::cerr << "Molecule::rate_assoc_disso_to() Error, unknown activity=" << partner.activity << std::endl;
      }
    }
    else if(Para::fold_unfold_constrain==2){
      switch(partner.activity){
	/* If the template is RNA. */
      case REPLI_RNA:
	assoc = fold_rate * bind_param * rna_recog_param * partner.unfold_rate;
	disso = unfold_rate * bind_param * (1. - rna_recog_param) * partner.fold_rate;
	return;
	/* If the template is DNA */
      case REPLI_DNA:
	/* I assume that DNA molecule does not fold. */
	assoc = fold_rate * bind_param * dna_recog_param;
	disso = 0.;
	return;
      case PARAS_RNA:
	assoc = Para::parasite_advantage * fold_rate * bind_param * rna_recog_param * partner.unfold_rate;
	disso = unfold_rate * bind_param * (1. - rna_recog_param) * partner.fold_rate;
	return;
	/* If the template is Parasite DNA */
      case PARAS_DNA:
	/* I assume that DNA molecule does not fold. */
	assoc = Para::parasite_advantage * fold_rate * bind_param * dna_recog_param;
	disso = 0.;
	return;
      case JUNK_RNA: case JUNK_DNA:
	/* I assume that JUNK molecules don't bind. */
	assoc = 0.;
	return;
      default:
	std::cerr << "Molecule::rate_assoc_disso_to() Error, unknown activity=" << partner.activity << std::endl;
      }
    }
    else{
      std::cerr << "Molecule::rate_assoc_disso_to() Unknown value in fold_unfold_constrain." << std::endl;
      exit(-1);
    }
  }
  else{
    std::cerr << "Molecule::rate_assoc_disso_to() Unkown model_type" << std::endl;
    exit(-1);
  }
}

/* This method makes *this molecule to bind to the "partner"
   molecule. The method considers *this to be the replicase which takes
   "partner" as template. */
void Molecule::binds_to(Molecule& partner,const double disso)
{
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
    /* Under the equilibrium assumption, we don't use SIMPL_FOLDE. If
       molecules are not in complex, they are always SIMPLE_TEMPL. */

    /* Check if the types of molecules are correct */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||
			   (mol_type==SIMPL_TEMPL && partner.mol_type==SIMPL_TEMPL));
    }
    catch(GeneralError){
      std::cerr << "Molecule::binds_to() Error, mol_type is going wrong" << std::endl;
      exit(-1);
    }

    /* Check if the activity types of molecules are correct */
    try{
      Assert<GeneralError>(!ASSERT::ERROR_CHECK||
			   (activity==REPLI_RNA && partner.activity!=JUNK_RNA && partner.activity!=JUNK_DNA));
    }
    catch(GeneralError){
      std::cerr << "Molecule::binds_to() Error, activity_type is going wrong" << std::endl;
      exit(-1);
    }

    /* I set repli_rate to the sum of the two kinds of replication
       rates. Thus, when replication does occur, we must determine which
       molecule (RNA/DNA) is produced.  See Molecule::replicate(). */
    mol_type = COMPL_ACTIV;
    disso_rate = disso;
    repli_rate = rna_repli_rate + dna_repli_rate;

    partner.mol_type = COMPL_PASSI;
    partner.disso_rate = disso;
    partner.repli_rate = rna_repli_rate + dna_repli_rate;
  }
  else{
    std::cerr << "Molecule::binds_to() Unkown model_type" << std::endl;
    exit(-1);
  }
}

double Molecule::rate_lipid_synth() const
{
  if(Para::model_type==Para::VESIC_NEUT_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_NEUT_EQUIL ||
     Para::model_type==Para::VESIC_LIPID_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL ||
     Para::model_type==Para::WELL_MIXED_EQUIL ||
     Para::model_type==Para::SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::WELL_MIXED_EQUIL_NO_COMPLEX || 
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_SURFACE_EQUIL ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL ||
     Para::model_type==Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX ||
     Para::model_type==Para::PARALLEL_INORG_EQUIL ||
     Para::model_type==Para::INORG_EQUIL ||
     Para::model_type==Para::PARALLEL_VESIC_FIXED_EQUIL){

    if(Para::fold_unfold_constrain==0 ||
       Para::fold_unfold_constrain==2){
      /* Explcit checking of activity makes things easier because we also
	 have DNA version of PARAS. */
      return ((activity==PARAS_RNA) ? lipid_synth_rate*fold_rate:0.);
    }
    else if(Para::fold_unfold_constrain==1){
      return ((activity==PARAS_RNA) ? lipid_synth_rate:0.);
    }
    else{
      std::cerr << "Molecule::rate_lipid_synth() Unknown value in Para::fold_unfold_constrain." << std::endl;
      exit(-1);
    }
  }
  else{
    std::cerr << "Molecule::rate_lipid_synth() Unkown model_type" << std::endl;
    exit(-1);
  }
}

/* This calculates the rate of replication reaction between two simple
   molecules "without complex formation". Note that this method
   considers *this to be the replicase that takes "partner" as
   template. The method returns an agglomerated replication rate of DNA
   and RNA (i.e. rna_repli_rate + dna_repli_rate). */
double Molecule::rate_repli_no_complex(const Molecule& templ) const
{
  /* I assume that replication happens only if the catalyst (i.e. *this)
     is REPLI_RNA. In the model without complex formation, I always use
     SIMPL_TEMPL for mol_type of a non-empty square. */
  if(activity != REPLI_RNA || templ.mol_type != SIMPL_TEMPL){
    return 0.;
  }

  switch(templ.activity){
    /* If the template is RNA. */
  case REPLI_RNA: case PARAS_RNA:
    return (fold_rate * rna_recog_param * templ.unfold_rate * (rna_repli_rate + dna_repli_rate));
    /* If the template is DNA */
  case REPLI_DNA: case PARAS_DNA:
    /* I assume that DNA molecule does not fold. */
    return (fold_rate * dna_recog_param * (rna_repli_rate + dna_repli_rate));
  case JUNK_RNA: case JUNK_DNA:
    /* I assume that JUNK molecules don't bind. */
    return 0.;
  default:
    std::cerr << "Molecule::rate_repli_no_complex() Error, unknown activity=" << templ.activity << std::endl;
    exit(-1);
  }
}

/* This replicates a molecule without taking complex formation into
   account. Unlike Molecule::replicate(), this function requires a user
   to specify which molecule is the template and which molecule is the
   replicate. This is because we know which molecule will be
   replicated. (*this) should be the empty square as in
   Molecule()::replicate(). */
void Molecule::replicate_no_complex(Molecule& templ,Molecule& repli)
{
    /* When I assume that replicases are either RNA polymerase or DNA
     polymerase (no hybrid). Thus, either dna_repli_rate or
     rna_repli_rate must be 0. */
  if(Para::rna_dna_repl_constrain==0 ||
     Para::rna_dna_repl_constrain==1){
    copy_with_mutation(templ);
    
    if(repli.rna_repli_rate>0.){
      convert_to_rna();
    }
    else{
      convert_to_dna();
    }
  }
  else{
    std::cerr << "Molecule::replicate_no_complex() Error, unknown Para::rna_dna_repl_constrain." << std::endl;
    exit(-1);
  }

}

bool Molecule::relative_close(double A,double B,double tolerance)
{
  if(A!=B){
    if(fabs((A-B)/B) <= tolerance && fabs((A-B)/A) <= tolerance)
      return true;
    else
      return false;
  }
  else{
    return true;
  }
}

#ifdef _COUNT_RNA_DNA
/* This is a function to output the information collected for RNA count
   and DNA count. It only outputs replicase molecules (ignore parasites
   and junk molecules). */
void Molecule::output_rna_dna_count(std::ostream& s)
{
  if(!is_empty())
    if(is_replicase()){
      /* Rp or Dp */
      if(rate_rna_repli()>0.){
	s << "Rp";
      }
      else{
	s << "Dp";
      }
      
      /* RNA template or DNA template */
      if(is_replicase_rna()){
	s << " R";
      }
      else{
	s << " D";
      }

      /* Recognition parameters & fold */
      s << ' ' << get_rna_recog_param()
	<< ' ' << get_dna_recog_param()
	<< ' ' << rate_fold();

      /* RNA & DNA count */
      s << ' ' << _rna_count
	<< ' ' << _dna_count
	<< std::endl;
    }
}
#endif //_COUNT_RNA_DNA

