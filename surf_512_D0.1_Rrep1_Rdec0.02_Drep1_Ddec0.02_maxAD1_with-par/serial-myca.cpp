#include "myca.hpp"

void MyCA::update_vesicle_1square(const unsigned row,const unsigned col)
{
  /* The state of vesicle of the current square may be modified if this
     square doesn't contain a molecule */
  if(!plane.cell(row,col).molecule.is_empty())
    return;

  /* Get a neighbor from Moore neighborhood */
  unsigned nei = rand_karney.Integer(8u) +1;
  
  /* We need to modify the state only if the first square and its
     neighbor are from different vesicle. The first square can get the
     state of the neighbor square. */
   int index_this,index_nei;
  if((index_this = plane.cell(row,col).vesicle_index)
     != (index_nei = plane.NEIGH_X(row,col,nei).vesicle_index)){

    /* To introduce an empty space between vesicles once in a while */
    if(index_this!=-1 && index_nei!=-1)
      if(rand_karney.Real()<0.001){
	index_nei=-1;
      }

    /* First, calculate the surface free energy differential */
     int energy_before=0, energy_after=0, index_tmp;

    /* Here, it's faster to use a for-loop, though this sounds
       crazy. Note that boundaries have vesicle_index=-1.  Thus, when
       boundaries are chosen as a neighbor, the neighbor will be
       considered as media (empty). Since copying process only modifies
       the vesicle state of the current cell, boundaries' state will be
       never modified. */
    for(unsigned i=1;i!=9;++i){
      index_tmp = plane.NEIGH_X(row,col,i).vesicle_index;
      energy_J(energy_before,index_this,index_tmp);
      energy_J(energy_after,index_nei,index_tmp);
    }

    /* Second, calculate the volume free energy differential */
    int energy_diff;
    if(index_this != -1 && index_nei != -1){
      energy_diff = energy_after - energy_before 
	+ 2*Para::lambda*(serial_vesicles.get_target_volume(index_this)-serial_vesicles.get_volume(index_this)+serial_vesicles.get_volume(index_nei)-serial_vesicles.get_target_volume(index_nei)+1);
    }else if(index_this != -1){
      energy_diff = energy_after - energy_before 
	+ 1 + 2*Para::lambda*(serial_vesicles.get_target_volume(index_this)-serial_vesicles.get_volume(index_this));
    }else{
      energy_diff = energy_after - energy_before 
	+ 1 + 2*Para::lambda*(serial_vesicles.get_volume(index_nei)-serial_vesicles.get_target_volume(index_nei));
    }

    /* Third, determine if we change the state */
    /* In contrast to using energy_J(), using get_probability() doesn't
       slow down the program */
    if(is_copy_state(energy_diff)){

      /* copy the state */
      plane.cell(row,col).vesicle_index = index_nei;

      /* update volume */
      if(index_this != -1){
	serial_vesicles.add_volume(index_this,-1);
	/* if volume becomes zero, delete this vesicle by setting the
	   volue to -1. */
	if(serial_vesicles.get_volume(index_this)==0){
	  serial_vesicles.destroy_vesicle(index_this);
	}
      }
      if(index_nei != -1)
	serial_vesicles.add_volume(index_nei);
    }
  }
}

void MyCA::update_replicator_whole_EQA()
{
  unsigned row,col;
  for(unsigned pos=0;pos!=Para::sys_nrow*Para::sys_ncol;++pos){
    /* Randomly choose a square to update */
    row = rand_karney.Integer(Para::sys_nrow) + 1;
    col = rand_karney.Integer(Para::sys_ncol) + 1;
//     row = (unsigned)(rand_karney.Real()*nrow) + 1;
//     col = (unsigned)(rand_karney.Real()*ncol) + 1;

    if(plane.cell(row,col).molecule.is_empty()){
      update_replicator_1square_empty(row,col);
    }
    /* the square is not empty && the molecule is simple */
    else if(plane.cell(row,col).molecule.is_simple_templ()){
      update_replicator_1square_simple_EQA(row,col);
    }
    /* the square is not empty && the molecule is complexed */
    else if(plane.cell(row,col).molecule.is_complex()){
      update_replicator_1square_complex(row,col);
    }
    else{
      std::cerr << "MyCA::update_replicator_1square_EQA() Unkown activity!" << std::endl;
      exit(-1);
    }
  }
}

void MyCA::update_replicator_1square_empty(const unsigned row,const unsigned col)
{
  const double p = rand_karney.Real();
  double cum_rate = Para::diffusion_rate;
  /*********************
   * DIFFUSION - EMPTY *
   *********************/
  if(p< cum_rate){

    unsigned nei = rand_karney.Integer(8u) + 1;
//     unsigned nei = (unsigned)(rand_karney.Real()*8.) +1;

    /* If neighbor is empty, do nothing */
    if(plane.NEIGH_X(row,col,nei).molecule.is_empty_or_boundary()){
      return;
    }
    /* The method of swapping differs according to whether the neighbor
       molecule is complexed or not */
    else if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
      /* If swapping is across membranes, we check if this swapping is
	 allowed to cross membranes. If not, don't swap*/
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	/* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	   cannot use cum_rate here!*/
	if(p>=Para::cross_membrane*cum_rate)
	  return;

      /* Swap the molecules (empty-simple) */
      plane.cell(row,col).molecule = plane.NEIGH_X(row,col,nei).molecule;
      plane.NEIGH_X(row,col,nei).molecule.decay();
    }
    else if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){
      /* If the neighbor is complexed, and swapping is across membranes,
	 we don't swap (this is to avoid a complex molecule spanning
	 across a membrane ) */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	return;

      /* Swap the molecules (empty-complex) */
      plane.cell(row,col).molecule.push(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).molecule);

      /* Resetting bon_nei */
      plane.cell(row,col).molecule.bon_nei = nei;
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
    }
    else{
      std::cerr << "MyCA::update_replicator_1square_empty() Unkown activity!" << std::endl;
      exit(-1);
    }

    /* The end of DIFFUSION. The next return is necessary because the
       next replication reaction does not start with else-if. */
    return;
  }
  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type!=Para::WELL_MIXED_EQUIL){
    unsigned nei = rand_karney.Integer(8u) + 1;
//     unsigned nei = (unsigned)(rand_karney.Real()*8.) +1;

    /* We should check if neighbor is complexed before calling
       rate_repli(). This is a little short cut to the algorithm and
       should not be used if there is more reactions to follow after
       this, in which case we must have the rate of replication to test
       whether the reaction that follows happen. */
    if(plane.NEIGH_X(row,col,nei).molecule.is_complex() && 
       plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

      cum_rate += plane.NEIGH_X(row,col,nei).molecule.rate_repli();

      if(p<cum_rate){
	plane.cell(row,col).molecule.replicate(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).molecule);
	return;
      }
    }
  }
  /* If the system is well-mixed. */
  else{
    unsigned n_row,n_col;
    do{
      n_row = rand_karney.Integer(Para::sys_nrow) + 1;
      n_col = rand_karney.Integer(Para::sys_ncol) + 1;
    }while(n_row==row && n_col==col);
    /* We should check if neighbor is complexed before calling
       rate_repli() */
    if(plane.cell(n_row,n_col).molecule.is_complex() &&
       plane.cell(row,col).vesicle_index == plane.cell(n_row,n_col).vesicle_index){

      cum_rate += plane.cell(n_row,n_col).molecule.rate_repli();

      if(p<cum_rate){
	plane.cell(row,col).molecule.replicate(plane.cell(n_row,n_col).molecule,plane.cell(plane.cell(n_row,n_col).molecule.bon_row,plane.cell(n_row,n_col).molecule.bon_col).molecule);
	return;
     }
    }
  }
}

/* 
   Under the equilibrium assumption, molecules are neither folded nor
   unfolded, but are folded and unfolded at the same time. When they are
   folded, they (except for junk molecules) can catalyze some chemical
   reaction (the fraction of time folded is multiplied to the rate
   constants of the reaction). Replicases can catalyze replication,
   whereas parasites can catalyze "lipid synthesis", i.e. the growth of
   vesicle target volume.
*/
void MyCA::update_replicator_1square_simple_EQA(const unsigned row,const unsigned col)
{
  const double p = rand_karney.Real();
  double cum_rate = Para::diffusion_rate;
  /**********************
   * DIFFUSION - SIMPLE *
   **********************/
  if(p< cum_rate){
    unsigned nei = rand_karney.Integer(8u) + 1;
//     unsigned nei = (unsigned)(rand_karney.Real()*8.) +1;

    /* The method of swapping differs according to whether the neighbor
       molecule is complexed or not */
    if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
      /* If diffusion is not allowed to cross membrane && if neighbor is
	 accross membranes, then do nothing. */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	/* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	   cannot use cum_rate here!*/
	if(p>=Para::cross_membrane*cum_rate)
	  return;

      /* Swap the molecules (simple-simple) */
      plane.NEIGH_X(row,col,nei).molecule = plane.cell(row,col).molecule;
      plane.cell(row,col).molecule.decay();
    }
    else if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
      /* If diffusion is not allowed to cross membrane && if neighbor is
	 accross membranes, then do nothing. */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	/* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	   cannot use cum_rate here!*/
	if(p>=Para::cross_membrane*cum_rate)
	  return;

      /* Swap the molecules (simple-simple) */
      plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei).molecule);
    }
    else if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){
      /* If the neighbor is complexed, and swapping is across membranes,
	 we don't swap */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	return;

      /* Swap the molecules (simple-complex) */
      plane.cell(row,col).molecule.rotate(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).molecule);

      /* Resetting bon_nei */
      plane.cell(row,col).molecule.bon_nei = nei;
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
    }
    /* This return is necessary because the follwing COMPLEX FORMATION
       does not begien with else-if. */
    return;
  }
  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p < cum_rate){
    plane.cell(row,col).molecule.decay();
    /* This return is necessary. */
    return;
  }
  /*******************
   * LIPID SYNTHESIS *
   *******************/
  cum_rate += plane.cell(row,col).molecule.rate_lipid_synth();
  if(p<cum_rate){
    if(plane.cell(row,col).vesicle_index>-1){
      serial_vesicles.add_target_volume(plane.cell(row,col).vesicle_index);
    }
   /* This return is necessary. */
   return;
  }
  /*********************
   * COMPLEX FORMATION *
   *********************/
  if(Para::model_type!=Para::WELL_MIXED_EQUIL){
    unsigned nei = rand_karney.Integer(8u) + 1;
//     unsigned nei = (unsigned)(rand_karney.Real()*8.) +1;

    /* Molecules cannot form complex across membranes. Checking first if
       we are looking at a neighbor that has the same vesicle_index is a
       little short cut to the algorithm. */
    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

      /* Bond in the direction of this->neighbor */
      double assoc=0.,disso=0.;
      plane.cell(row,col).molecule.rate_assoc_disso_to(plane.NEIGH_X(row,col,nei).molecule,assoc,disso);
      cum_rate += assoc;
      if(p<cum_rate){
	plane.cell(row,col).molecule.bon_nei = nei;
	plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	plane.cell(row,col).molecule.binds_to(plane.NEIGH_X(row,col,nei).molecule,disso);
	/* This return is necessary because the next reaction does not
	   begin with else-if. */
	return;
      }
      /* Bond in the direction of neighbot->this */
      plane.NEIGH_X(row,col,nei).molecule.rate_assoc_disso_to(plane.cell(row,col).molecule,assoc,disso);
      cum_rate += assoc;
      if(p<cum_rate){
	plane.cell(row,col).molecule.bon_nei = nei;
	plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	plane.NEIGH_X(row,col,nei).molecule.binds_to(plane.cell(row,col).molecule,disso);
	/* This return is only here to make code symmetric */
	return;
      }
    }
  }
  /* If the system is well-mixed */ 
  else{
    unsigned n_row,n_col;
    do{
      n_row = rand_karney.Integer(Para::sys_nrow) + 1;
      n_col = rand_karney.Integer(Para::sys_ncol) + 1;
    }while(n_row==row && n_col==col);

    /* Molecules cannot form complex across membranes. Checking first if
       we are looking at a neighbor that has the same vesicle_index is a
       little short cut to the algorithm. */
    if(plane.cell(row,col).vesicle_index == plane.cell(n_row,n_col).vesicle_index){
      double assoc=0.,disso=0.;
      plane.cell(row,col).molecule.rate_assoc_disso_to(plane.cell(n_row,n_col).molecule,assoc,disso);
      cum_rate += assoc;
      if(p<cum_rate){
	plane.cell(row,col).molecule.bon_row = n_row;
	plane.cell(row,col).molecule.bon_col = n_col;
	plane.cell(n_row,n_col).molecule.bon_row = row;
	plane.cell(n_row,n_col).molecule.bon_col = col;
	plane.cell(row,col).molecule.binds_to(plane.cell(n_row,n_col).molecule,disso);
	/* This regurn is necessary. */
	return;
      }

      plane.cell(n_row,n_col).molecule.rate_assoc_disso_to(plane.cell(row,col).molecule,assoc,disso);
      cum_rate += assoc;
      if(p<cum_rate){
	plane.cell(row,col).molecule.bon_row = n_row;
	plane.cell(row,col).molecule.bon_col = n_col;
	plane.cell(n_row,n_col).molecule.bon_row = row;
	plane.cell(n_row,n_col).molecule.bon_col = col;
	plane.cell(n_row,n_col).molecule.binds_to(plane.cell(row,col).molecule,disso);
	return;
      }
    }
  }
}

void MyCA::update_replicator_1square_complex(const unsigned row,const unsigned col)
{
  const double p = rand_karney.Real();
  double cum_rate = Para::diffusion_rate_complex;
  /***********************
   * DIFFUSION - COMPLEX *
   ***********************/
  /* Diffusion of complexed molecules are handled as follows.  The
     diffusion rate is halved for complex molecules because the chance
     that a complex molecule is chosen is twice as much as that of a
     simple molecule. When the chosen neighbor of a chosen complex turns
     out to be the partner molecule, then we "rotate" a complex
     molecule. Note that when the partner molecule is chosen and diffuse
     (move), the non-chosen molecule also diffuse (move) with the same
     distance. Imagine there is only one complex molecule in the
     system. In the above method, the distance one of the two molecules
     comprising a complex travels should be equal to that traveled by a
     simple molecule. */
  if(p< cum_rate){
    unsigned nei = rand_karney.Integer(8u) + 1;
//     unsigned nei = (unsigned)(rand_karney.Real()*8.) +1;

    /* The method of swapping differs according to whether the neighbor
       molecule is complexed or not. If the neighbor molecule is
       complexed, it can be the partner molecule of this molecule or a
       molecule of another complex.  Thus, there are three
       possibilities  */
    /* case 1: neighbor is the partner molecule of itself */
    if(plane.cell(row,col).molecule.bon_nei == nei){
      /* Swapping */
      plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei).molecule);
      /* Reset bon_nei */
      plane.cell(row,col).molecule.bon_nei = nei;
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
    }
    /* Case 2: neighbor is a molecule of another complex */
    else if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){
      /* If diffusion is not allowed to cross membrane && if neighbor
	 is accross membranes, then do nothing. */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	/* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	   cannot use cum_rate here!*/
	if(p>=Para::cross_membrane*cum_rate)
	  return;

      /* Swap the molecules (complex-complex) */
      /* I must remember bon_nei of this and bon_nei of nei because I
	 need them (resetting bon_nei) after overwriting this and
	 nei. */
      const unsigned bon_nei = plane.cell(row,col).molecule.bon_nei;
      const unsigned nei_bon_nei = plane.NEIGH_X(row,col,nei).molecule.bon_nei;
      /* Actual swapping */
      plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei,nei_bon_nei).molecule);
      plane.NEIGH_X(row,col,nei).molecule.swap_with(plane.NEIGH_X(row,col,bon_nei).molecule);

      /* Reset bon_nei */
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = nei_bon_nei;
      plane.NEIGH_X(row,col,nei,nei_bon_nei).molecule.bon_nei = neighbor_converter(nei_bon_nei);
      plane.cell(row,col).molecule.bon_nei = bon_nei;
      plane.NEIGH_X(row,col,bon_nei).molecule.bon_nei = neighbor_converter(bon_nei);

    }
    /* Case 3: neighbor is a simple molecule */
    else if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
      /* If only one of the two molecule is complexed, and swapping is
	 across membranes, we don't swap (or cannot due to the lack of
	 a good algorithm.) */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	return;

      /* Swaping */
      plane.NEIGH_X(row,col,nei).molecule.rotate(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule);
	
      /* Resetting bon_nei */
      plane.cell(row,col).molecule.bon_nei = nei;
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
    }
    /* Case 4: neighbor is empty */
    else if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
      /* If only one of the two molecule is complexed, and swapping is
	 across membranes, we don't swap (or cannot due to the lack of
	 a good algorithm.) */
      if(plane.cell(row,col).vesicle_index != plane.NEIGH_X(row,col,nei).vesicle_index)
	return;

      /* Actual swaping */
      plane.NEIGH_X(row,col,nei).molecule.push(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule);
	
	/* Resetting bon_nei */
      plane.cell(row,col).molecule.bon_nei = nei;
      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
    }

    /* Wether diffusion actually happened or not, we must return because
       p was smaller than cum_rate.*/
    return;
  }
  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p<cum_rate){
    if(Para::model_type!=Para::WELL_MIXED_EQUIL){
      plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule.dissociate();
    }else{
      plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_col).molecule.dissociate();
    }
    plane.cell(row,col).molecule.decay();
    /* The next return is necessary. */
    return;
  }

  /****************
   * DISSOCIATION *
   ****************/
  cum_rate += plane.cell(row,col).molecule.rate_disso();
  if(p<cum_rate){
    if(Para::model_type!=Para::WELL_MIXED_EQUIL){
      plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule.dissociate();
    }else{
      plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_col).molecule.dissociate();
    }
    plane.cell(row,col).molecule.dissociate();
    /* The next return is necessary. */    
    return;
  }

  /***************
   * REPLICATION *
   ***************/
  cum_rate += plane.cell(row,col).molecule.rate_repli();
  if(p<cum_rate){
    if(Para::model_type!=Para::WELL_MIXED_EQUIL){
      /* See scale_rates() in main.cpp for the explanation on why we use
	 neigh_7_select() here. */
      unsigned nei = rand_karney.Integer(7u) + 1;
//       unsigned nei = (unsigned)(rand_karney.Real()*7.) +1;

      nei = plane.neigh_7_select(nei,plane.cell(row,col).molecule.bon_nei);
      
      if(plane.NEIGH_X(row,col,nei).molecule.is_empty()&&
	 plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){
	plane.NEIGH_X(row,col,nei).molecule.replicate(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule);
      }
    }
    /* If the system is well-mixed */
    else{
      unsigned n_row,n_col;
      do{
	n_row = rand_karney.Integer(Para::sys_nrow) + 1;
	n_col = rand_karney.Integer(Para::sys_ncol) + 1;
      } while(n_row==row && n_col==col);
    
      if(plane.cell(n_row,n_col).molecule.is_empty()){
	plane.cell(n_row,n_col).molecule.replicate(plane.cell(row,col).molecule,plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_col).molecule);
      }
    }
  }
}

void MyCA::update_replicator_whole_no_complex_EQA()
{
  unsigned row,col;
  for(unsigned pos=0;pos!=Para::sys_nrow*Para::sys_ncol;++pos){
    /* Randomly choose a square to update */
    row = rand_karney.Integer(Para::sys_nrow) + 1;
    col = rand_karney.Integer(Para::sys_ncol) + 1;

    if(plane.cell(row,col).molecule.is_empty()){
      update_replicator_1square_empty_no_complex_EQA(row,col);
    }
    /* the square is not empty && the molecule is simple */
    else if(plane.cell(row,col).molecule.is_simple_templ()){
      update_replicator_1square_simple_no_complex_EQA(row,col);
    }
    else{
      std::cerr << "MyCA::update_replicator_1square_no_complex_EQA() Unkown activity!" << std::endl;
      exit(-1);
    }
  }
}

void MyCA::update_replicator_1square_empty_no_complex_EQA(const unsigned row,const unsigned col)
{
  const double p = rand_karney.Real();
  double cum_rate = Para::diffusion_rate;
  /*********************
   * DIFFUSION - EMPTY *
   *********************/
  if(p< cum_rate){

    const unsigned nei = rand_karney.Integer(8u) + 1;

    /* The method of swapping differs according to whether the neighbor
       molecule is complexed or not */
    if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
      /* If swapping is across membranes, we check if this swapping is
	 allowed to cross membranes. If not, don't swap*/
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	 p<Para::cross_membrane*cum_rate){

	/* Swap the molecules (empty-simple) */
	plane.cell(row,col).molecule = plane.NEIGH_X(row,col,nei).molecule;
	plane.NEIGH_X(row,col,nei).molecule.decay();
      }
    }
    /* This return is necessary because the next replication reaction
       does not start with else-if. */
    return;
  }
  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type != Para::WELL_MIXED_EQUIL_NO_COMPLEX){
    /* Since replication reaction without complex formation is
       tri-molecular reaction, we randomly choose two neighbors. */
    const unsigned nei1 = rand_karney.Integer(8u) + 1;

    /* We should check if the neighbor is simple_template (i.e. not
       empty) and in the same vesicle or outside before calling
       rate_repli_no_complex(). This is a little short cut to the
       algorithm and should not be used if there is more reactions to
       follow after this and if those reactions do not require this
       condition to be fulfilled in order to occur. Generally speaking,
       we must first check the rate of replication to test whether the
       reaction can happen and then check the other conditions because
       those other conditions are not necessarly shared by differnet
       kinds of reaction. Also, note that although replication reaction
       can happen in two different ways below, it is ok to do the
       checking because they both share the same condition of
       happening. */
    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei1).vesicle_index &&
       plane.NEIGH_X(row,col,nei1).molecule.is_simple_templ()){

      const unsigned nei2 = plane.neigh_7_select(rand_karney.Integer(7u) + 1,nei1);

      /* We check if the second neighbor is simple_template and in the
	 same vesicle or outside before calling
	 rate_repli_no_complex(). For the same reason as explained
	 above, we should be careful with this way of doing things
	 because it cannot be generally used.  */
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei2).vesicle_index &&
	 plane.NEIGH_X(row,col,nei2).molecule.is_simple_templ()){
	
	/* We finally check at the rate of reaction. */
	cum_rate += plane.NEIGH_X(row,col,nei1).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei2).molecule);
	if(p < cum_rate){
	  plane.cell(row,col).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei2).molecule,plane.NEIGH_X(row,col,nei1).molecule);
	  /* This return is necessary because the next reaction does not
	     start with else-if. */
	  return;
	}

	cum_rate += plane.NEIGH_X(row,col,nei2).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei1).molecule);
	if(p < cum_rate){
	  plane.cell(row,col).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei1).molecule,plane.NEIGH_X(row,col,nei2).molecule);
	  return;
	}
      }
    }
  }
  /* If the system is well-mixed. */
  else{
    unsigned n1_row,n1_col;
    do{
      n1_row = rand_karney.Integer(Para::sys_nrow) + 1;
      n1_col = rand_karney.Integer(Para::sys_ncol) + 1;
    }while(n1_row==row && n1_col==col);

    if(plane.cell(n1_row,n1_col).molecule.is_simple_templ()){

      unsigned n2_row,n2_col;
      do{
	n2_row = rand_karney.Integer(Para::sys_nrow) + 1;
	n2_col = rand_karney.Integer(Para::sys_ncol) + 1;
      }while((n2_row==row && n2_col==col) || (n2_row==n1_row && n2_col==n1_col));

      if(plane.cell(n2_row,n2_col).molecule.is_simple_templ()){
	
	cum_rate += plane.cell(n1_row,n1_col).molecule.rate_repli_no_complex(plane.cell(n2_row,n2_col).molecule);

	if(p < cum_rate){
	  plane.cell(row,col).molecule.replicate_no_complex(plane.cell(n2_row,n2_col).molecule,plane.cell(n1_row,n1_col).molecule);
	  /* This return is necessary because the next reaction does not
	     start with else-if. */
	  return;
	}

	cum_rate += plane.cell(n2_row,n2_col).molecule.rate_repli_no_complex(plane.cell(n1_row,n1_col).molecule);
	if(p < cum_rate){
	  plane.cell(row,col).molecule.replicate_no_complex(plane.cell(n1_row,n1_col).molecule,plane.cell(n2_row,n2_col).molecule);
	  return;
	}
      }
    }
  }
}

void MyCA::update_replicator_1square_simple_no_complex_EQA(const unsigned row,const unsigned col)
{
  const double p = rand_karney.Real();
  double cum_rate = Para::diffusion_rate;
  /**********************
   * DIFFUSION - SIMPLE *
   **********************/
  if(p < cum_rate){
    const unsigned nei = rand_karney.Integer(8u) + 1;

    /* The method of swapping differs according to whether the neighbor
       molecule is complexed or not */
    if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
      /* If diffusion is not allowed to cross membrane && if neighbor is
	 accross membranes, then do nothing. */
      if(p < Para::cross_membrane*cum_rate ||
	 plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){
	
	/* Swap the molecules (simple-simple) */
	plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei).molecule);
      }
    }
    else if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
      if(p < Para::cross_membrane*cum_rate ||
	 plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	/* Swap the molecules (simple-simple) */
	plane.NEIGH_X(row,col,nei).molecule = plane.cell(row,col).molecule;
	plane.cell(row,col).molecule.decay();
      }
    }
    /* This return is necessary because the follwing COMPLEX FORMATION
       does not begien with else-if. */
    return;
  }
  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p < cum_rate){
    plane.cell(row,col).molecule.decay();
    /* This return is necessary. */
    return;
  }
  /*******************
   * LIPID SYNTHESIS *
   *******************/
  cum_rate += plane.cell(row,col).molecule.rate_lipid_synth();
  if(p < cum_rate){
    if(plane.cell(row,col).vesicle_index>-1){
      serial_vesicles.add_target_volume(plane.cell(row,col).vesicle_index);
    }
   /* This return is necessary. */
   return;
 }
  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type != Para::WELL_MIXED_EQUIL_NO_COMPLEX){
    const unsigned nei1 = rand_karney.Integer(8u) + 1;

    /* We first check if the neighbors have the same vesicle_index. This
       is a little short cut to the algorithm, and should not be used if
       there is more reactions to follow after this reaction and if
       those reactions do not share the same condition in order to
       happen. In general, we must first check the rate of association
       to test whether the reaction that follows happen. But it would be
       a shame if we didn't exploit the specificity of the reaction to
       speed up the algorithm. */
    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei1).vesicle_index){

      unsigned const nei2 = plane.neigh_7_select(rand_karney.Integer(7u) + 1,nei1);
    
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei2).vesicle_index){
	/* Depending on which square is empty, which molecules we should
	   apply rate_repli_no_complex() differs. So, we check it;
	   otherwise, there are six (3!) ways in which we have to call
	   rate_repli_no_complex(), which would be inefficient.
	   Moreover, we also discard those combinations of mol_type's
	   where no replication occurs. This is again a shortcut to the
	   algorithm and should not be used if there are more reactions
	   to follow and if those reactions do not share the same
	   condition (on the combination of mol_type) in order to
	   happen. */
	if(plane.NEIGH_X(row,col,nei1).molecule.is_simple_templ() && plane.NEIGH_X(row,col,nei2).molecule.is_empty()){
	  cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei1).molecule);
	  if(p < cum_rate){
	    plane.NEIGH_X(row,col,nei2).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei1).molecule,plane.cell(row,col).molecule);
	    /* This return is necessary because the next reaction does
	       not start with else-if. */
	    return;
	  }

	  cum_rate += plane.NEIGH_X(row,col,nei1).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
	  if(p < cum_rate){
	    plane.NEIGH_X(row,col,nei2).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,nei1).molecule);
	    return;
	  }
	}
	else if(plane.NEIGH_X(row,col,nei1).molecule.is_empty() && plane.NEIGH_X(row,col,nei2).molecule.is_simple_templ()){
	  cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei2).molecule);
	  if(p < cum_rate){
	    plane.NEIGH_X(row,col,nei1).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei2).molecule,plane.cell(row,col).molecule);
	    /* This return is necessary because the next reaction does
	       not start with else-if. */
	    return;
	  }

	  cum_rate += plane.NEIGH_X(row,col,nei2).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
	  if(p < cum_rate){
	    plane.NEIGH_X(row,col,nei1).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,nei2).molecule);
	    return;
	  }
	}
      }
    }
  }
  /* If the system is well-mixed */ 
  else{
    unsigned n1_row,n1_col;
    do{
      n1_row = rand_karney.Integer(Para::sys_nrow) + 1;
      n1_col = rand_karney.Integer(Para::sys_ncol) + 1;
    }while(n1_row==row && n1_col==col);

    /* Molecules cannot form complex across membranes. Checking first if
       we are looking at a neighbor that has the same vesicle_index is a
       little short cut to the algorithm, and should not be used if
       there is more reactions to follow after this reaction. */
    unsigned n2_row,n2_col;
    do{
      n2_row = rand_karney.Integer(Para::sys_nrow) + 1;
      n2_col = rand_karney.Integer(Para::sys_ncol) + 1;
    }while((n2_row==row && n2_col==col) || (n2_row==n1_row && n2_col==n1_col));

    if(plane.cell(n1_row,n1_col).molecule.is_simple_templ() && plane.cell(n2_row,n2_col).molecule.is_empty()){

      cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.cell(n1_row,n1_col).molecule);
      if(p < cum_rate){
	plane.cell(n2_row,n2_col).molecule.replicate_no_complex(plane.cell(n1_row,n1_col).molecule,plane.cell(row,col).molecule);
	/* This return is necessary because the next reaction does
	   not start with else-if. */
	return;
      }

      cum_rate += plane.cell(n1_row,n1_col).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
      if(p < cum_rate){
	plane.cell(n2_row,n2_col).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.cell(n1_row,n1_col).molecule);
	return;
      }
    }
    else if(plane.cell(n1_row,n1_col).molecule.is_empty() && plane.cell(n2_row,n2_col).molecule.is_simple_templ()){

      cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.cell(n2_row,n2_col).molecule);
      if(p < cum_rate){
	plane.cell(n1_row,n1_col).molecule.replicate_no_complex(plane.cell(n2_row,n2_col).molecule,plane.cell(row,col).molecule);
	/* This return is necessary because the next reaction does
	   not start with else-if. */
	return;
      }
      cum_rate += plane.cell(n2_row,n2_col).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
      if(p < cum_rate){
	plane.cell(n1_row,n1_col).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.cell(n2_row,n2_col).molecule);
	return;
      }
    }
  }
}

void MyCA::vesicle_neutral_growth()
{
  serial_vesicles.set_all_target_volume_to_0();
  /* We count the number of molecules in vesicles. */
  for(unsigned row=1;row<=nrow;++row)
    for(unsigned col=1;col<=ncol;++col){
      if(plane.cell(row,col).vesicle_index!=-1 
	 && !plane.cell(row,col).molecule.is_empty()){
	serial_vesicles.add_target_volume(plane.cell(row,col).vesicle_index);
      }
    }
  serial_vesicles.neutral_growth();
}

void MyCA::divide_vesicle(const long Time)
{
  /* Collect information to calculate principle component of vesicles
     that must be divided. */
  serial_vesicles.PCA_reset();
  bool is_divide = false;
  for(unsigned row=1;row<=nrow;++row)
    for(unsigned col=1;col<=ncol;++col){
      if(plane.cell(row,col).vesicle_index>-1){
	/* Set vesicles with volume>threshold to divide, and collect
	   information to calculate PCA. */
	if(serial_vesicles.get_volume(plane.cell(row,col).vesicle_index)>=Para::volume_threshold){
	  /* Tell that this vesicle is going to divide. */
	  if(!serial_vesicles.is_divide(plane.cell(row,col).vesicle_index)){
	    serial_vesicles.set_is_divide(plane.cell(row,col).vesicle_index);
	  }
	  /* This is just to know that there is a division event at
	     all. */
	  if(!is_divide){
	    is_divide = true;
	  }
	  serial_vesicles.PCA_collect_data(plane.cell(row,col).vesicle_index,row,col);
	}
      }
    }

  /* Calculate principle component of vesicles that must be divide */
  serial_vesicles.calculate_PCA_for_dividing_vesicles();

  /* If there are vesicles that have to divide, is_divide is true. */
  if(is_divide){
    /* Go through a whole plane to divide vesicles. We assign some
       squares to new vesicles */
    int mother_ind;
    for(unsigned row=1;row<=nrow;++row)
      for(unsigned col=1;col<=ncol;++col){
	/* Here, I use 'mother_ind' because I cannot use
	   plane.cell(row,col).vesicle_index directly, for this might be
	   changed during the division. For the same reason, don't use
	   const reference either! */
	mother_ind = plane.cell(row,col).vesicle_index;
	if(mother_ind>-1){
	  /* Check if this vesicle has to be divided (if so, is_divide()
	     return true). */
	  if(serial_vesicles.is_divide(mother_ind)){
	    /* is_daughter() tells whether this square belongs to a
	       daughter vesicle. If so, we change the state of this
	       square */
	    /* Note that it may be possible that the next 'is_daughter'
	       never returns true. This happens, e.g., when the volume
	       of the mother vesicle is very small. In that case, even
	       if the mother vesicle's is_divide flag is true, it will
	       not divide. */
	    if(serial_vesicles.is_daughter(mother_ind,row,col)){
	      /* If a daughter hasn't been assigned to this mother
		 vesicle, get_daughter_ind() returns -1. In that case,
		 we assign one. */
	      if(serial_vesicles.get_daughter_ind(mother_ind)==-1){
		serial_vesicles.set_daughter_ind(mother_ind,serial_vesicles.make_new_vesicle());
	      }
	      /* A daughter index has already been assined */
	      /* This square is now belonging to the daughter */
	      plane.cell(row,col).vesicle_index = serial_vesicles.get_daughter_ind(mother_ind);
	      
	      /* update volume */
	      serial_vesicles.add_volume(mother_ind,-1);
	      serial_vesicles.add_volume(serial_vesicles.get_daughter_ind(mother_ind));
	      
	      /* a complex molecule that spans between mother and
		 daughter vesicles must dissociate. */
	      if(plane.cell(row,col).molecule.is_complex()){
		unsigned nei_row=0,nei_col=0;
		plane.XY_NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei,nei_row,nei_col);
		if(!serial_vesicles.is_daughter(mother_ind,nei_row,nei_col)){
		  plane.cell(row,col).molecule.dissociate();
		  plane.cell(nei_row,nei_col).molecule.dissociate();
		}
	      }
	    }
	  }
	}
      }
    /* set target_volume (the mother and daughter get the fraction of
       the original target_volume of the mother proportional to their
       volume */
    serial_vesicles.divide_target_volume();
  }
}
