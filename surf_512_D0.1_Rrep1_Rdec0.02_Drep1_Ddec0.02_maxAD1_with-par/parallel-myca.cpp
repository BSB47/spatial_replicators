#include "myca.hpp"

void MyCA::parallel_update_replicator_1square_no_complex_EQA()
{
  /* First, get the thread-specific random number generator. */
  RandomTBB::reference random_wrapper = random_tbb.local();

  /* Second, we reserve all squares that are potentially updated. In the
     model without complex formation, we need to reserve two distinct
     neighbors for replication reaction and one neighbor for diffusion
     that can be identical to one of the two neighbors for
     replication. */
  unsigned row,col;

  while(true){
    row = random_wrapper.random.Integer(Para::sys_nrow) + 1;
    col = random_wrapper.random.Integer(Para::sys_ncol) + 1;

    if(plane.cell(row,col).lock()){
      if(plane.cell(row,col).molecule.is_empty()){
	if(parallel_update_replicator_1square_empty_no_complex_EQA(row,col)){
	  plane.cell(row,col).release();
	  break;
	}
      }
      else if(plane.cell(row,col).molecule.is_simple_templ()){
	if(parallel_update_replicator_1square_simple_no_complex_EQA(row,col)){
	  plane.cell(row,col).release();
	  break;
	}
      }
      else{
	std::cerr << "MyCA::parallel_update_replicator_1square_no_complex_EQA() Unkown activity!" << std::endl;
	exit(-1);
      }
      plane.cell(row,col).release();
    }
  }
}


/*
  The following principles are adopted to write the program in the
  parallelized algorithm. 

  Principle I: Always return if p<cum_rate is true.

  if(p<cum_rate){
     ...
     ...
     return true;
  }



  Principle II: Always release X in the same clause if(X.lock()) unless
  it contradicts with Principle I.

  if(X.lock()){
     ...
     ...
     if(p<cum_rate){
        ...
        ...
        X.release();    <= This is an exception.
        return true;
     }
     ...
     ...
     X.release();       <= This is the general rule.
  }
  else{
     ...
     ...
     return false;
  }
*/

bool MyCA::parallel_update_replicator_1square_empty_no_complex_EQA(const unsigned row,const unsigned col)
{
  /* Get the thread-specific random number generator. */
  RandomTBB::reference wrapper = random_tbb.local();
  const double p = wrapper.random.Real();
  double cum_rate = Para::diffusion_rate;
  /*********************
   * DIFFUSION - EMPTY *
   *********************/
  if(p < cum_rate){
    const unsigned nei = wrapper.random.Integer(8u) + 1;
    if(plane.NEIGH_X(row,col,nei).lock()){
      /* The method of swapping differs according to whether the neighbor
	 molecule is complexed or not */
      if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){

	/* If swapping is across membranes, we check if this swapping is
	   allowed to cross membranes. If not, don't swap*/
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   /* At this point, cum_rate equals to diffusion_rate. */
	   p<Para::cross_membrane*cum_rate){
	  
	  /* Swap the molecules (empty-simple) */
	  plane.cell(row,col).molecule = plane.NEIGH_X(row,col,nei).molecule;
	  plane.NEIGH_X(row,col,nei).molecule.decay();
	}
      }
      plane.NEIGH_X(row,col,nei).release();
    }
    else{
      /* Locking failed. Tell the caller that updating is unfinished. */
      return false;
    }
    return true;
  }

  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type != Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){
    /* Since replication reaction without complex formation is
       tri-molecular reaction, we randomly choose two neighbors. */
    const unsigned nei1 = wrapper.random.Integer(8u) + 1;
    if(plane.NEIGH_X(row,col,nei1).lock()){

      /* We should check if the neighbor is simple_template (i.e. not
	 empty) and in the same vesicle or outside before calling
	 rate_repli_no_complex(). This is a little shortcut to the
	 algorithm. */
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei1).vesicle_index &&
	 plane.NEIGH_X(row,col,nei1).molecule.is_simple_templ()){

	const unsigned nei2 = plane.neigh_7_select(wrapper.random.Integer(7u) + 1,nei1);
	if(plane.NEIGH_X(row,col,nei2).lock()){

	  /* We check if the second neighbor is simple_template and in the
	     same vesicle or outside before calling
	     rate_repli_no_complex(). */
	  if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei2).vesicle_index &&
	     plane.NEIGH_X(row,col,nei2).molecule.is_simple_templ()){
	
	    /* We finally check at the rate of reaction. */
	    cum_rate += plane.NEIGH_X(row,col,nei1).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei2).molecule);
	    if(p<cum_rate){
	      plane.cell(row,col).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei2).molecule,plane.NEIGH_X(row,col,nei1).molecule);

	      /* This return is necessary; */
	      plane.NEIGH_X(row,col,nei2).release();
	      plane.NEIGH_X(row,col,nei1).release();
	      return true;
	    }

	    cum_rate += plane.NEIGH_X(row,col,nei2).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei1).molecule);
	    if(p<cum_rate){
	      plane.cell(row,col).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei1).molecule,plane.NEIGH_X(row,col,nei2).molecule);

	      /* This return is necessary; */
	      plane.NEIGH_X(row,col,nei2).release();
	      plane.NEIGH_X(row,col,nei1).release();
	      return true;
	    }
	  }
	  plane.NEIGH_X(row,col,nei2).release();
	}
	else{
	  /* Locking failed. Tell the caller that updating isn't finished. */
	  plane.NEIGH_X(row,col,nei1).release();
	  return false;
	}
      }
      plane.NEIGH_X(row,col,nei1).release();
    }
    else{
      /* Locking failed. Tell the caller that updating isn't finished. */
      return false;
    }
  }
  /* If the system is well-mixed. */
  else{
    unsigned n1_row,n1_col;
    do{
      n1_row = wrapper.random.Integer(Para::sys_nrow) + 1;
      n1_col = wrapper.random.Integer(Para::sys_ncol) + 1;
    }while(n1_row==row && n1_col==col);

    if(plane.cell(n1_row,n1_col).lock()){

      if(plane.cell(row,col).vesicle_index == plane.cell(n1_row,n1_col).vesicle_index &&
	 plane.cell(n1_row,n1_col).molecule.is_simple_templ()){

	unsigned n2_row,n2_col;
	do{
	  n2_row = wrapper.random.Integer(Para::sys_nrow) + 1;
	  n2_col = wrapper.random.Integer(Para::sys_ncol) + 1;
	}while((n2_row==row && n2_col==col) || (n2_row==n1_row && n2_col==n1_col));

	if(plane.cell(n2_row,n2_col).lock()){

	  if(plane.cell(n2_row,n2_col).molecule.is_simple_templ()){
	
	    cum_rate += plane.cell(n1_row,n1_col).molecule.rate_repli_no_complex(plane.cell(n2_row,n2_col).molecule);
	    if(p<cum_rate){
	      plane.cell(row,col).molecule.replicate_no_complex(plane.cell(n2_row,n2_col).molecule,plane.cell(n1_row,n1_col).molecule);

	      /* This return is necessary. */
	      plane.cell(n2_row,n2_col).release();
	      plane.cell(n1_row,n1_col).release();
	      return true;
	    }

	    cum_rate += plane.cell(n2_row,n2_col).molecule.rate_repli_no_complex(plane.cell(n1_row,n1_col).molecule);
	    if(p<cum_rate){
	      plane.cell(row,col).molecule.replicate_no_complex(plane.cell(n1_row,n1_col).molecule,plane.cell(n2_row,n2_col).molecule);

	      /* This return is necessary. */
	      plane.cell(n2_row,n2_col).release();
	      plane.cell(n1_row,n1_col).release();
	      return true;
	    }
	  }
	  plane.cell(n2_row,n2_col).release();
	}
	else{
	  plane.cell(n1_row,n1_col).release();
	  return false;
	}
      }
      plane.cell(n1_row,n1_col).release();
    }
    else{
      return false;
    }
  }
  
  /* If we reach here, nothing has happened. Updating is finished. More
     reactions can follow here. */
  return true;
}

bool MyCA::parallel_update_replicator_1square_simple_no_complex_EQA(const unsigned row,const unsigned col)
{
  /* Get the thread-specific random number generator. */
  RandomTBB::reference wrapper = random_tbb.local();
  const double p = wrapper.random.Real();
  double cum_rate = Para::diffusion_rate;

  /**********************
   * DIFFUSION - SIMPLE *
   **********************/
  if(p < cum_rate){
    const unsigned nei = wrapper.random.Integer(8u) + 1;
    if(plane.NEIGH_X(row,col,nei).lock()){
      /* The method of swapping differs according to whether the neighbor
	 molecule is complexed or not */
      if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
	/* If diffusion is not allowed to cross membrane && if neighbor is
	   accross membranes, then do nothing. */
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   /* Be careful with cum_rate here. Officially, it has to be
	      diffusion_rate. */
	   p<Para::cross_membrane*cum_rate ){

	  /* Swap the molecules (simple-simple) */
	    plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei).molecule);
	}
      }
      else if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   p<Para::cross_membrane*cum_rate ){

	  plane.NEIGH_X(row,col,nei).molecule = plane.cell(row,col).molecule;
	  plane.cell(row,col).molecule.decay();
	}
      }
      plane.NEIGH_X(row,col,nei).release();
    }
    else{
      /* Locking failed. Tell the caller that updating is unfinished. */
      return false;
    }
    return true;
  }

  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p < cum_rate){
    /* Let molecule decay. */
    plane.cell(row,col).molecule.decay();
    /* This return is necessary because the following Replication
       doesn't start with else-if. */
    return true;
  }

  /*******************
   * LIPID SYNTHESIS *
   *******************/
  cum_rate += plane.cell(row,col).molecule.rate_lipid_synth();
  if(p < cum_rate){
    /* Let molecule synthesize lipid. */
    if(plane.cell(row,col).vesicle_index>-1){
      atomic_vesicles.add_target_volume(plane.cell(row,col).vesicle_index);
    }
    return true;
  }

  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type != Para::PARALLEL_WELL_MIXED_EQUIL_NO_COMPLEX){

    const unsigned nei1 = wrapper.random.Integer(8u) + 1;
    if(plane.NEIGH_X(row,col,nei1).lock()){

      /* We first check if nei1 have the same vesicle_index. */
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei1).vesicle_index){

	const unsigned nei2 = plane.neigh_7_select(wrapper.random.Integer(7u) + 1,nei1);
	if(plane.NEIGH_X(row,col,nei2).lock()){

	  /* We first check if nei1 have the same vesicle_index. */
	  if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei2).vesicle_index){

	    /* Depending on which square is empty, which molecules we
	       should apply rate_repli_no_complex() differs. So, we check
	       it; otherwise, there are six (3!) ways in which we have to
	       call rate_repli_no_complex(), which would be inefficient.
	       Moreover, we also discard those combinations of mol_type's
	       where no replication occurs. */
	    if(plane.NEIGH_X(row,col,nei1).molecule.is_simple_templ() && plane.NEIGH_X(row,col,nei2).molecule.is_empty()){

	      /* Nei1 is replicated by this into nei2. */
	      cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei1).molecule);
	      if(p < cum_rate){
		/* Replicate molecule. */
		plane.NEIGH_X(row,col,nei2).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei1).molecule,plane.cell(row,col).molecule);

		/* Updating is done. Release all the locked squares. */
		plane.NEIGH_X(row,col,nei2).release();
		plane.NEIGH_X(row,col,nei1).release();

		/* The next return is necessary */
		return true;
	      }
      
	      /* This is replicated by nei1 into nei2. */
	      cum_rate += plane.NEIGH_X(row,col,nei1).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);

	      if(p < cum_rate){
		/* Replicate molecule. */
		plane.NEIGH_X(row,col,nei2).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,nei1).molecule);

		plane.NEIGH_X(row,col,nei2).release();
		plane.NEIGH_X(row,col,nei1).release();
		return true;
	      }
	    }
	    else if(plane.NEIGH_X(row,col,nei1).molecule.is_empty() && plane.NEIGH_X(row,col,nei2).molecule.is_simple_templ()){
	      /* Nei2 is replicated by this into nei1. */
	      cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.NEIGH_X(row,col,nei2).molecule);
	      if(p < cum_rate){
		/* Replicate molecule. */
		plane.NEIGH_X(row,col,nei1).molecule.replicate_no_complex(plane.NEIGH_X(row,col,nei2).molecule,plane.cell(row,col).molecule);

		plane.NEIGH_X(row,col,nei2).release();
		plane.NEIGH_X(row,col,nei1).release();
		return true;
	      }
      
	      /* This is replicated by nei2 into nei1. */
	      cum_rate += plane.NEIGH_X(row,col,nei2).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
	      if(p < cum_rate){
		/* Replicate molecule. */
		plane.NEIGH_X(row,col,nei1).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,nei2).molecule);

		plane.NEIGH_X(row,col,nei2).release();
		plane.NEIGH_X(row,col,nei1).release();
		return true;
	      }
	    }
	    /* Updating is done including the case where replication didn't
	       happen. Release all the locked squares. */
	    plane.NEIGH_X(row,col,nei1).release();
	    plane.NEIGH_X(row,col,nei2).release();
	    plane.cell(row,col).release();
	  }
	  plane.NEIGH_X(row,col,nei2).release();
	}
	else{
	  /* Couldn't lock nei2. */
	  plane.NEIGH_X(row,col,nei1).release();
	  return false;
	}
      }
      plane.NEIGH_X(row,col,nei1).release();
    }
    else{
      /* Couldn't lock nei1. */
      return false;
    }
  }
  else{
    unsigned n1_row,n1_col;
    do{
      n1_row = wrapper.random.Integer(Para::sys_nrow) + 1;
      n1_col = wrapper.random.Integer(Para::sys_ncol) + 1;
    }while(n1_row==row && n1_col==col);

    if(plane.cell(n1_row,n1_col).lock()){

      unsigned n2_row,n2_col;
      do{
	n2_row = wrapper.random.Integer(Para::sys_nrow) + 1;
	n2_col = wrapper.random.Integer(Para::sys_ncol) + 1;
      }while((n2_row==row && n2_col==col) || (n2_row==n1_row && n2_col==n1_col));

      if(plane.cell(n2_row,n2_col).lock()){

	if(plane.cell(n1_row,n1_col).molecule.is_simple_templ() && plane.cell(n2_row,n2_col).molecule.is_empty()){
	  /* Nei1 is replicated by this into nei2. */
	  cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.cell(n1_row,n1_col).molecule);
	  if(p < cum_rate){
	    plane.cell(n2_row,n2_col).molecule.replicate_no_complex(plane.cell(n1_row,n1_col).molecule,plane.cell(row,col).molecule);
	    
	    plane.cell(n2_row,n2_col).release();
	    plane.cell(n1_row,n1_col).release();
	    return true;
	  }

	  /* This is replicated by nei1 into nei2. */
	  cum_rate += plane.cell(n1_row,n1_col).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
	  if(p < cum_rate){
	    plane.cell(n2_row,n2_col).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.cell(n1_row,n1_col).molecule);

	    plane.cell(n2_row,n2_col).release();
	    plane.cell(n1_row,n1_col).release();
	    return true;
	  }
	}
	else if(plane.cell(n1_row,n1_col).molecule.is_empty() && plane.cell(n2_row,n2_col).molecule.is_simple_templ()){

	  /* Nei2 is replicated by this into nei1. */
	  cum_rate += plane.cell(row,col).molecule.rate_repli_no_complex(plane.cell(n2_row,n2_col).molecule);
	  if(p < cum_rate){
	    plane.cell(n1_row,n1_col).molecule.replicate_no_complex(plane.cell(n2_row,n2_col).molecule,plane.cell(row,col).molecule);

	    plane.cell(n2_row,n2_col).release();
	    plane.cell(n1_row,n1_col).release();
	    return true;
	  }

	  /* This is replicated by nei2 into nei1. */
	  cum_rate += plane.cell(n2_row,n2_col).molecule.rate_repli_no_complex(plane.cell(row,col).molecule);
	  if(p < cum_rate){
	    plane.cell(n1_row,n1_col).molecule.replicate_no_complex(plane.cell(row,col).molecule,plane.cell(n2_row,n2_col).molecule);

	    plane.cell(n2_row,n2_col).release();
	    plane.cell(n1_row,n1_col).release();
	    return true;
	  }
	}
	plane.cell(n2_row,n2_col).release();
      }
      else{
	plane.cell(n1_row,n1_col).release();
	return false;
      }
      plane.cell(n1_row,n1_col).release();
    }
    else{
      return false;
    }
  }

  /* If we reach here, nothing has happened. Updating is finished. More
     reactions can follow here. */
  return true;
}

/*
  General idea is the following. We try to lock every square this thread
  might access (either read or write) until we succeed. Within this
  method, it is possible to release some of the chosen neighbors
  according to its content; however, we do not release them because of a
  potential overlap among neighbors---that is, two chosen neighbors can
  actually be the same square. When the two neighbors are the same
  square and if we release one of them, we will actually release both
  squares, which is illegal. Thus, if we are to release it, we must
  check the neighbors are not overlapping. I feel will be an overhead,
  so releasing will be done at the end of updating.
*/
void MyCA::parallel_update_replicator_1square_EQA()
{
  /* First, get the thread-specific random number generator. */
  RandomTBB::reference random_wrapper = random_tbb.local();

  /* Second, we reserve all squares that are potentially updated. In the
     model with complex formation, if updated squares contain complex
     molecules, we need to reserve a square where the partner molecules
     reside. In generall, we have to choose two neighbors: one for
     diffusion and one for second-order reaction, which are either
     replication or complex association. */

  unsigned row,col;
  bool not_finished=true;
  do{
    row = random_wrapper.random.Integer(Para::sys_nrow) + 1;
    col = random_wrapper.random.Integer(Para::sys_ncol) + 1;

    if(plane.cell(row,col).lock()){
      if(plane.cell(row,col).molecule.is_empty()){
	not_finished=parallel_update_replicator_1square_empty(row,col);
      }
      /* the square is not empty && the molecule is simple */
      else if(plane.cell(row,col).molecule.is_simple_templ()){
	not_finished=parallel_update_replicator_1square_simple_EQA(row,col);
      }
      /* the square is not empty && the molecule is complexed */
      else if(plane.cell(row,col).molecule.is_complex()){
	not_finished=parallel_update_replicator_1square_complex(row,col);
      }
      else{
	std::cerr << "MyCA::update_replicator_1square_EQA() Unkown activity!" << std::endl;
	exit(-1);
      }
      plane.cell(row,col).release();
    }
  }while(not_finished);
}

/*
  This is generally applicable to
  parallel_update_replicator_1square_empty(),
  parallel_update_replicator_1square_simple(),
  parallel_update_replicator_1square_complex(). The idea is the
  following. We first proceed updating until we know which neighbor we
  need access. We then try to lock it. If we cannot lock it, we abort
  the updating and return with true, which signifies that the updating
  is not_finished. If we succesfully update, we return false, which
  signifies the updating is not not_finished (i.e. finished). This
  algorithm does introduce some bias such that 1st-order reaction
  happens more often than 2nd-order reaction in comparison with the
  canonical algorithm. This is simply because we abort only 2nd-order
  reactions due to neighbor retrieval. But, given that "collision" is
  infrequent, I hope this do not introduce too much bias.

  Basic philosophy of the following code is this. We lock and realease
  on the same depth of claudes. We return on two occasions. 

  (1) p<cum_rate was true. In this case, return will be placed at the
  end of the clause which starts with if(p<cum_rates).

  (2) Neighbors weren't successfully retrieved. In this case, return
  will be such that

  if(X.lock()){
  ...do_reaction...
  }
  else{
  ...release_other_neighbors... 
  return true;
  }
  
  (3) Nothing happend. Return will be placed at the end of the method.
*/

bool MyCA::parallel_update_replicator_1square_empty(const unsigned row,const unsigned col)
{
  RandomTBB::reference wrapper = random_tbb.local();
  const double p = wrapper.random.Real();
  double cum_rate = Para::diffusion_rate;
  /*********************
   * DIFFUSION - EMPTY *
   *********************/
  if(p< cum_rate){

    /* Try to lock a neighbor for diffusion */
    unsigned nei = wrapper.random.Integer(8u) + 1;

    if(plane.NEIGH_X(row,col,nei).lock()){

      /* The method of swapping differs according to whether the neighbor
	 molecule is complexed or not */
      if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
	/* If swapping is across membranes, we check if this swapping is
	   allowed to cross membranes. If not, don't swap*/
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   /* Note that we can do the next line because
	      cum_rate==Para::diffusion_rate, which is not generally true,
	      in which case we cannot use cum_rate here!*/
	   p<Para::cross_membrane*cum_rate){
	  /* Updating is done. Release nei and return not
	     not_finished. */

	  /* Swap the molecules (empty-simple) */
	  plane.cell(row,col).molecule = plane.NEIGH_X(row,col,nei).molecule;
	  plane.NEIGH_X(row,col,nei).molecule.decay();
	}
      }
      else if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){

	/* If the neighbor is complexed, and diffusion is prohibited
	   across membranes. This is to avoid a complex molecule
	   spanning across a membrane */
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	  /* In order to swap a complex molecule, we must lock the partner
	     molecule of nei. Be careful: we must store the value of bon_nei
	     of nei here, because we will lose it through swapping. */
	  const unsigned neipar = plane.NEIGH_X(row,col,nei).molecule.bon_nei;
	  if(plane.NEIGH_X(row,col,nei,neipar).lock()){

	    /* Swap the molecules (empty-complex) */
	    plane.cell(row,col).molecule.push(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,neipar).molecule);
	    /* Resetting bon_nei */
	    plane.cell(row,col).molecule.bon_nei = nei;
	    plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);

	    /* Updating is done. Release the bon_nei of nei. I could also
	       release nei and and return not not_finished here, but I do it
	       later. */
	    plane.NEIGH_X(row,col,nei,neipar).release();
	  }
	  else{
	    /* We couldn't lock it, so we abort this updating. */
	    plane.NEIGH_X(row,col,nei).release();
	    return true;
	  }
	}
      }
      plane.NEIGH_X(row,col,nei).release();
    }
    /* We couldn't lock the neighbor, so we abort and go find another
       squares to update. We return true, which means it is
       not_finished. */
    else{
      return true;
    }
    /* Updating is done. Return not not_finished.*/
    return false;
  }
  /***************
   * REPLICATION *
   ***************/
  if(Para::model_type!=Para::PARALLEL_WELL_MIXED_EQUIL){
    unsigned nei = wrapper.random.Integer(8u) + 1;

    /* We must first lock nei for replication. */
    if(plane.NEIGH_X(row,col,nei).lock()){
      /* We check if neighbor is complexed before calling
	 rate_repli(). This is a little shortcut to the algorithm. We
	 should be careful with this kind of shortcut because it's not
	 very clear when the updating is finished (see below). */
      if(plane.NEIGH_X(row,col,nei).molecule.is_complex() && 
	 plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	/* The condition for replication is now almost complete. We must
	   now also lock the partner of nei. */
	if(plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).lock()){
      
	  cum_rate += plane.NEIGH_X(row,col,nei).molecule.rate_repli();

	  if(p<cum_rate){
	    plane.cell(row,col).molecule.replicate(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).molecule);
	    /* Updating is done. Release the neighbors and return not
	       not_finished. Note that the next three lines are actually
	       redundant because we release neighbors and return false
	       anyway after this clause. However, if there were more
	       reactions to follow, these lines must be placed here. So, for
	       generality, I did write them down here. */
	    plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).release();
	    plane.NEIGH_X(row,col,nei).release();
	    return false;
	  }
	  /* Note that updating can be considered as finished here
	     because there is no reaction to follow. But, if there were
	     more reactions, we would not be able to consider this point
	     as the end of updating. So, for generality, I don't return
	     here, but only release bon_nei of nei.  */
	  plane.NEIGH_X(row,col,nei,plane.NEIGH_X(row,col,nei).molecule.bon_nei).release();
	}
	else{
	  /* We couldn't lock it, abort, release nei and report that
	     updating is not_finished. */
	  plane.NEIGH_X(row,col,nei).release();

	  return true;
	}
      }
      plane.NEIGH_X(row,col,nei).release();
    }
    else{
      /* We couldn't lock it, abort and report that updating is
	 not_finished. */
      return true;
    }
  }
  /* If the system is well-mixed. */
  else{
    unsigned n_row,n_col;
    do{
      n_row = wrapper.random.Integer(Para::sys_nrow) + 1;
      n_col = wrapper.random.Integer(Para::sys_ncol) + 1;
    }while(n_row==row && n_col==col);

    /* We must first lock nei for replication. */
    if(plane.cell(n_row,n_col).lock()){

      /* We should check if neighbor is complexed before calling
	 rate_repli(). This is a shortcut. If the condition is not
	 fulfiled, we don't have to increase cum_rate, so we can
	 immediately go to next reaction (but currently there is no more
	 reaction to follow). */
      if(plane.cell(n_row,n_col).molecule.is_complex()){

	/* We must now also lock the partner of nei. */
	if(plane.cell(plane.cell(n_row,n_col).molecule.bon_row,plane.cell(n_row,n_col).molecule.bon_col).lock()){

	  cum_rate += plane.cell(n_row,n_col).molecule.rate_repli();

	  if(p<cum_rate){
	    plane.cell(row,col).molecule.replicate(plane.cell(n_row,n_col).molecule,plane.cell(plane.cell(n_row,n_col).molecule.bon_row,plane.cell(n_row,n_col).molecule.bon_col).molecule);

	    /* Updating is done. Release the neighbors and return not
	       not_finished. */
	    plane.cell(plane.cell(n_row,n_col).molecule.bon_row,plane.cell(n_row,n_col).molecule.bon_col).release();
	    plane.cell(n_row,n_col).release();
	    return false;
	  }
	  plane.cell(plane.cell(n_row,n_col).molecule.bon_row,plane.cell(n_row,n_col).molecule.bon_col).release();
	}
	else{
	  /* We couldn't lock it, so we abort. Release nei and report that
	     updating is not_finished. */
	  plane.cell(n_row,n_col).release();
	  return true;
	}
      }
      plane.cell(n_row,n_col).release();
    }
    else{      
      /* We couldn't lock it, abort and report that updating is
	 not_finished. */
      return true;
    }
  }
  /* Nothing happend up to here. HERE IS THE POINT WHERE MORE REACTIONS
     CAN FOLLOW. Note that cum_rate is correctly updated at this point:
     in Replicatoin, if nei was not complexed, cum_rate was not
     increased; if nei was complexed, cum_rate is most likely to have
     been increased. This is consistent whether there are more reactions
     to follow or not. */
  return false;
}

bool MyCA::parallel_update_replicator_1square_simple_EQA(const unsigned row,const unsigned col)
{
  RandomTBB::reference wrapper = random_tbb.local();
  const double p = wrapper.random.Real();
  double cum_rate = Para::diffusion_rate;
  /**********************
   * DIFFUSION - SIMPLE *
   **********************/
  if(p< cum_rate){
    unsigned nei = wrapper.random.Integer(8u) + 1;

    if(plane.NEIGH_X(row,col,nei).lock()){

      /* The method of swapping differs according to whether the neighbor
	 molecule is complexed or not */
      if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
	/* Diffusion always happens within a vesicle. It can happen
	   across membrane with a small chance. */
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   /* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	      cannot use cum_rate here!*/
	   p<Para::cross_membrane*cum_rate){
	
	  /* Swap the molecules (simple-simple) */
	  plane.NEIGH_X(row,col,nei).molecule = plane.cell(row,col).molecule;
	  plane.cell(row,col).molecule.decay();
	}
      }
      else if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
	/* Diffusion always happens within a vesicle. It can happen
	   across membrane with a small chance. */
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	   /* Note that cum_rate==Para::diffusion_rate; othwerwise, we
	      cannot use cum_rate here!*/
	   p<Para::cross_membrane*cum_rate){

	  /* Swap the molecules (simple-simple) */
	  plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei).molecule);
	}
      }
      else if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){
	/* If the neighbor is complexed, swapping happens only within
	   the same vesicle because this is simple molecule. (We must
	   not place a complex half-inside half-outside of a
	   vesicle.) */
	if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	  /* In order to swap a complex molecule, we must lock the partner
	     molecule of nei. Caution: we must remember the value of bon_nei
	     of nei, because we later swap molecules, which will change the
	     value of bon_nei. */
	  const unsigned neipar = plane.NEIGH_X(row,col,nei).molecule.bon_nei;
	  if(plane.NEIGH_X(row,col,nei,neipar).lock()){

	    /* Swap the molecules (simple-complex) */
	    plane.cell(row,col).molecule.rotate(plane.NEIGH_X(row,col,nei).molecule,plane.NEIGH_X(row,col,nei,neipar).molecule);

	    /* Resetting bon_nei */
	    plane.cell(row,col).molecule.bon_nei = nei;
	    plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	  
	    /* Release the bon_nei of nei. */
	    plane.NEIGH_X(row,col,nei,neipar).release();
	  }
	  else{
	    /* We couldn't lock neipar, so we abort this updating. */
	    plane.NEIGH_X(row,col,nei).release();
	    return true;
	  }
	}
      }
      /* This return is necessary because the follwing COMPLEX FORMATION
	 does not begien with else-if. */
      plane.NEIGH_X(row,col,nei).release();
    }
    else{
      /* We couldn't lock nei. */
      return true;
    }
    return false;
  }

  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p < cum_rate){
    plane.cell(row,col).molecule.decay();
    /* This return is necessary. */
    return false;
  }
  /*******************
   * LIPID SYNTHESIS *
   *******************/
  cum_rate += plane.cell(row,col).molecule.rate_lipid_synth();
  if(p < cum_rate){
    if(plane.cell(row,col).vesicle_index>-1){
      atomic_vesicles.add_target_volume(plane.cell(row,col).vesicle_index);
    }
   /* This return is necessary. */
   return false;
  }
  /*********************
   * COMPLEX FORMATION *
   *********************/
  if(Para::model_type!=Para::PARALLEL_WELL_MIXED_EQUIL){
    unsigned nei = wrapper.random.Integer(8u) + 1;

    /* We must first lock nei for replication. */
    if(plane.NEIGH_X(row,col,nei).lock()){

      /* Molecules cannot form complex across membranes. Checking first if
	 nei has the same vesicle_index. Also, we check if nei is simple
	 molecule. These kinds of checking are shortcut to the reaction
	 algorithm. We must be careful with the shortcut because this
	 makes unclear where is the end of the reaction. */
      if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index && 
	 plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){

	/* Bond in the direction of this->neighbor */
	double assoc=0.,disso=0.;
	plane.cell(row,col).molecule.rate_assoc_disso_to(plane.NEIGH_X(row,col,nei).molecule,assoc,disso);
	cum_rate += assoc;
	if(p<cum_rate){
	  plane.cell(row,col).molecule.bon_nei = nei;
	  plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	  plane.cell(row,col).molecule.binds_to(plane.NEIGH_X(row,col,nei).molecule,disso);

	  /* Release nei */
	  plane.NEIGH_X(row,col,nei).release();
	  /* This return is necessary because the next reaction does not
	     begin with else-if. */
	  return false;
	}

	/* Bond in the direction of neighbot->this */
	plane.NEIGH_X(row,col,nei).molecule.rate_assoc_disso_to(plane.cell(row,col).molecule,assoc,disso);
	cum_rate += assoc;
	if(p<cum_rate){
	  plane.cell(row,col).molecule.bon_nei = nei;
	  plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	  plane.NEIGH_X(row,col,nei).molecule.binds_to(plane.cell(row,col).molecule,disso);

	  /* Release nei */
	  plane.NEIGH_X(row,col,nei).release();
	  /* This return is only here to make code generic */

	  return false;
	}
      }
      plane.NEIGH_X(row,col,nei).release();
    }
    else{
      /* We couldn't lock it, abort and report that updating is
	 not_finished. */
      return true;
    }
    /* Caution: Arriving here means that p<cum_rate was never true,
       rather than p<cum_rate was true, but the other condition for
       complex association wasn't fulfilled. So, we must not return if
       there are more reactions to follow. */
  }
  /* If the system is well-mixed */ 
  else{
    unsigned n_row,n_col;
    do{
      n_row = wrapper.random.Integer(Para::sys_nrow) + 1;
      n_col = wrapper.random.Integer(Para::sys_ncol) + 1;
    }while(n_row==row && n_col==col);

    /* We must first lock nei for replication. */
    if(plane.cell(n_row,n_col).lock()){

      /* Molecules cannot form complex across membranes. Checking first if
	 nei has the same vesicle_index. Also, we check if nei is simple
	 molecule. These kinds of checking are shortcut to the reaction
	 algorithm. */
      if(plane.cell(n_row,n_col).molecule.is_simple_templ()){

	double assoc=0.,disso=0.;
	plane.cell(row,col).molecule.rate_assoc_disso_to(plane.cell(n_row,n_col).molecule,assoc,disso);
	cum_rate += assoc;
	if(p<cum_rate){
	  plane.cell(row,col).molecule.bon_row = n_row;
	  plane.cell(row,col).molecule.bon_col = n_col;
	  plane.cell(n_row,n_col).molecule.bon_row = row;
	  plane.cell(n_row,n_col).molecule.bon_col = col;
	  plane.cell(row,col).molecule.binds_to(plane.cell(n_row,n_col).molecule,disso);

	  /* Release nei */
	  plane.cell(n_row,n_col).release();
	  /* This return is necessary. */
	  return false;
	}

	plane.cell(n_row,n_col).molecule.rate_assoc_disso_to(plane.cell(row,col).molecule,assoc,disso);
	cum_rate += assoc;
	if(p<cum_rate){
	  plane.cell(row,col).molecule.bon_row = n_row;
	  plane.cell(row,col).molecule.bon_col = n_col;
	  plane.cell(n_row,n_col).molecule.bon_row = row;
	  plane.cell(n_row,n_col).molecule.bon_col = col;
	  plane.cell(n_row,n_col).molecule.binds_to(plane.cell(row,col).molecule,disso);

	  /* Release nei */
	  plane.cell(n_row,n_col).release();
	  /* This return is necessary. */
	  return false;
	}
      }
      /* Release nei */
      plane.cell(n_row,n_col).release();
    }
    else{
      return true;
    }
  }
  /* Nothing happend up to here. HERE IS THE POINT WHERE MORE
     REACTIONS CAN FOLLOW. Note that cum_rate is correctly updated at
     this point: if nei was not complexed, cum_rate was not increased;
     if nei was complexed, cum_rate is most likely to have been
     increased. This is consistent whether there are more reactions to
     follow or not. */
  return false;
}

bool MyCA::parallel_update_replicator_1square_complex(const unsigned row,const unsigned col)
{
  RandomTBB::reference wrapper = random_tbb.local();
  const double p = wrapper.random.Real();
  /***********************
   * DIFFUSION - COMPLEX *
   ***********************/
  /* Diffusion of complexed molecules are handled as follows.  The
     diffusion rate is halved for complex molecules because the chance
     that a complex molecule is chosen is twice as much as that of a
     simple molecule. When the chosen neighbor of a chosen complex
     turns out to be the partner molecule, then we "rotate" a complex
     molecule. Note that when the partner molecule is chosen and
     diffuse (move), the non-chosen molecule also diffuse (move) with
     the same distance. Imagine there is only one complex molecule in
     the system. In the above method, the distance one of the two
     molecules comprising a complex travels should be equal to that
     traveled by a simple molecule. */
  double cum_rate = Para::diffusion_rate_complex;
  if(p< cum_rate){

    /* Lock the partner of this. We must remember the value of bon_nei
       of this, because the original value will be modified by
       diffusion. */
    const unsigned par=plane.cell(row,col).molecule.bon_nei;

    if(plane.NEIGH_X(row,col,par).lock()){
      unsigned nei = wrapper.random.Integer(8u) + 1;

      /* The method of swapping differs according to whether the neighbor
	 molecule is complexed or not. If the neighbor molecule is
	 complexed, it can be the partner molecule of this molecule or a
	 molecule of another complex.  Thus, there are three
	 possibilities  */

      /* If nei!=par (i.e. nei is not the partner molecule of this), we
	 must next lock nei. */
      if(nei!=par){
	if(plane.NEIGH_X(row,col,nei).lock()){

	  /* Case 1: neighbor is a molecule of another complex */
	  if(plane.NEIGH_X(row,col,nei).molecule.is_complex()){
	    /* Diffusion always happen within a vesicle. It can happen
	       across membrane with a small chance. */
	    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index ||
	       /* Note that cum_rate==Para::diffusion_rate; othwerwise, we
		  cannot use cum_rate here!*/
	       p<Para::cross_membrane*cum_rate){

	      /* Lock the partner molecule of nei. I must remember
		 bon_nei of nei because I need them (when resetting
		 bon_nei) after overwriting nei. */
	      const unsigned neipar = plane.NEIGH_X(row,col,nei).molecule.bon_nei;
	      if(plane.NEIGH_X(row,col,nei,neipar).lock()){

		/* Swap the molecules (complex-complex) */
		plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,nei,neipar).molecule);
		plane.NEIGH_X(row,col,nei).molecule.swap_with(plane.NEIGH_X(row,col,par).molecule);

		/* Reset bon_nei */
		plane.NEIGH_X(row,col,nei).molecule.bon_nei = neipar;
		plane.NEIGH_X(row,col,nei,neipar).molecule.bon_nei = neighbor_converter(neipar);
		plane.cell(row,col).molecule.bon_nei = par;
		plane.NEIGH_X(row,col,par).molecule.bon_nei = neighbor_converter(par);

		/* Release the partner of nei */
		plane.NEIGH_X(row,col,nei,neipar).release();
	      }
	      else{
		/* Updating is failed. */
		plane.NEIGH_X(row,col,nei).release();
		plane.NEIGH_X(row,col,par).release();
		return true;
	      }
	    }
	  }
	  /* Case 2: neighbor is a simple molecule */
	  else if(plane.NEIGH_X(row,col,nei).molecule.is_simple_templ()){
	    /* If only one of the two molecule is complexed,  swapping
	       happens only if diffusion is within a vesicle. */
	    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	      /* Swaping */
	      plane.NEIGH_X(row,col,nei).molecule.rotate(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,par).molecule);
	
	      /* Resetting bon_nei */
	      plane.cell(row,col).molecule.bon_nei = nei;
	      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	    }
	  }
	  /* Case 3: neighbor is empty */
	  else if(plane.NEIGH_X(row,col,nei).molecule.is_empty()){
	    /* If only one of the two molecule is complexed,  swapping
	       happens only if diffusion is within a vesicle. */
	    if(plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){

	      /* Swaping */
	      plane.NEIGH_X(row,col,nei).molecule.push(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule);
	
	      /* Resetting bon_nei */
	      plane.cell(row,col).molecule.bon_nei = nei;
	      plane.NEIGH_X(row,col,nei).molecule.bon_nei = neighbor_converter(nei);
	    }
	  }
	  plane.NEIGH_X(row,col,nei).release();
	}
	else{
	  /* We couldn't lock nei. We abort. */
	  plane.NEIGH_X(row,col,par).release();
	  return true;
	}
      }
      /* Case 4: neighbor is the partner molecule of this. Rotation
	 of the complex happens. */
      else{
	/* Swapping */
	plane.cell(row,col).molecule.swap_with(plane.NEIGH_X(row,col,par).molecule);
	/* Reset bon_nei */
	plane.cell(row,col).molecule.bon_nei = par;
	plane.NEIGH_X(row,col,par).molecule.bon_nei = neighbor_converter(par);
      }
      plane.NEIGH_X(row,col,par).release();
    }
    else{
      /* We couldn't lock the partner molecule of this. Abort. */
      return true;
    }

    /* Updating is done.*/
    return false;
  }

  /*********
   * DECAY *
   *********/
  cum_rate += plane.cell(row,col).molecule.rate_decay();
  if(p<cum_rate){
    if(Para::model_type!=Para::PARALLEL_WELL_MIXED_EQUIL){
      
      /* Lock the partner molecule of this. */
      if(plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).lock()){
	plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule.dissociate();
	plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).release();
      }
      else{
	/* We couldn't lock the partner molecule. Abort. */
	return true;
      }
    }
    else{
      if(plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_row).lock()){
	plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_col).molecule.dissociate();
	plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_row).release();
      }
      else{
	return true;
      }
    }
    /* The next line should not be placed before locking the partner
       molecule of this is done. We do the reaction only after it is
       assured that we can do the reaction. */
    plane.cell(row,col).molecule.decay();
    /* The next return is necessary. */
    return false;
  }

  /****************
   * DISSOCIATION *
   ****************/
  cum_rate += plane.cell(row,col).molecule.rate_disso();
  if(p<cum_rate){
    if(Para::model_type!=Para::PARALLEL_WELL_MIXED_EQUIL){
      
      if(plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).lock()){
	plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).molecule.dissociate();
	plane.NEIGH_X(row,col,plane.cell(row,col).molecule.bon_nei).release();
      }
      else{
	return true;
      }
    }
    else{
      if(plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_row).lock()){
	plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_col).molecule.dissociate();
	plane.cell(plane.cell(row,col).molecule.bon_row,plane.cell(row,col).molecule.bon_row).release();
      }
      else{
	return true;
      }
    }
    plane.cell(row,col).molecule.dissociate();
    return false;
  }

  /***************
   * REPLICATION *
   ***************/
  cum_rate += plane.cell(row,col).molecule.rate_repli();
  if(p<cum_rate){
    if(Para::model_type!=Para::PARALLEL_WELL_MIXED_EQUIL){
      const unsigned par = plane.cell(row,col).molecule.bon_nei;

      /* Lock the partner molecule */
      if(plane.NEIGH_X(row,col,par).lock()){

	/* See scale_rates() in main.cpp for the explanation on why we use
	   neigh_7_select() here. */
	unsigned nei = wrapper.random.Integer(7u) + 1;
      
	nei = plane.neigh_7_select(nei,par);

	/* Lock the neighbor */
	if(plane.NEIGH_X(row,col,nei).lock()){
	  if(plane.NEIGH_X(row,col,nei).molecule.is_empty() &&
	     plane.cell(row,col).vesicle_index == plane.NEIGH_X(row,col,nei).vesicle_index){
	    plane.NEIGH_X(row,col,nei).molecule.replicate(plane.cell(row,col).molecule,plane.NEIGH_X(row,col,par).molecule);
	  }
	  plane.NEIGH_X(row,col,nei).release();
	}
	else{
	  /* Couldn't lock the nei. Abort. */
	  plane.NEIGH_X(row,col,par).release();
	  return true;
	}

	plane.NEIGH_X(row,col,par).release();
      }
      else{
	/* We couldn't lock the partner molecule. Abort. */
	return true;
      }
    }
    /* If the system is well-mixed */
    else{
      unsigned par_row = plane.cell(row,col).molecule.bon_row;
      unsigned par_col = plane.cell(row,col).molecule.bon_col;

      if(plane.cell(par_row,par_col).lock()){

	unsigned n_row,n_col;
	do{
	  n_row = wrapper.random.Integer(Para::sys_nrow) + 1;
	  n_col = wrapper.random.Integer(Para::sys_ncol) + 1;
	} while(n_row==row && n_col==col);
    
	if(plane.cell(n_row,n_col).lock()){
	  if(plane.cell(n_row,n_col).molecule.is_empty()){
	    plane.cell(n_row,n_col).molecule.replicate(plane.cell(row,col).molecule,plane.cell(par_row,par_col).molecule);
	  }
	  plane.cell(n_row,n_col).release();
	}
	else{
	  plane.cell(par_row,par_col).release();
	  return true;
	}
	plane.cell(par_row,par_col).release();
      }
      else{
	return true;
      }

    }
    /* Updating is done */
    return false;
  }
  /* Nothing happend up to here. HERE IS THE POINT WHERE MORE
     REACTIONS CAN FOLLOW. */
  return false;
}

/*
  I need not lock the square at (row,col) to update, because the
  algorithm never looks at the same square within one time
  step. Moreover, since we use atomic_volume, we actually need not lock
  anything.
*/
void MyCA::parallel_update_vesicle_1square(const unsigned row,const unsigned col)
{
  /* The state of vesicle of the current square may be modified if this
     square doesn't contain a molecule */
  if(!plane.cell(row,col).molecule.is_empty())
    return;

  /* Get the thread-specific random number generator. */
  RandomTBB::reference random_wrapper = random_tbb.local();

  /* Get a neighbor from Moore neighborhood */
  unsigned nei = random_wrapper.random.Integer(8u) +1;
  
  /* We need to modify the state only if the first square and its
     neighbor are from different vesicle. The first square can get the
     state of the neighbor square. */
   int index_this,index_nei;
  if((index_this = plane.cell(row,col).vesicle_index)
     != (index_nei = plane.NEIGH_X(row,col,nei).vesicle_index)){

    /* To introduce an empty space between vesicles once in a while */
    if(index_this!=-1 && index_nei!=-1)
      if(random_wrapper.random.Real()<0.001){
	index_nei=-1;
      }

    /* First, calculate the surface free energy differential */
     int energy_before=0, energy_after=0, index_tmp;

    /* Here, it's faster to use a for-loop, though this sounds crazy
       (may be, it's reasonable). Note that boundaries have
       vesicle_index=-1.  Thus, when boundaries are chosen as a
       neighbor, the neighbor will be considered as media (empty). Since
       copying process only modifies the vesicle state of the current
       cell, boundaries' state will be never modified. */
    for(unsigned i=1;i!=9;++i){
      index_tmp = plane.NEIGH_X(row,col,i).vesicle_index;
      energy_J(energy_before,index_this,index_tmp);
      energy_J(energy_after,index_nei,index_tmp);
    }

    /* Second, calculate the volume free energy differential */
    int energy_diff;
    if(index_this != -1 && index_nei != -1){
      energy_diff = energy_after - energy_before 
	+ 2*Para::lambda*(atomic_vesicles.get_target_volume(index_this)-atomic_vesicles.get_volume(index_this)+atomic_vesicles.get_volume(index_nei)-atomic_vesicles.get_target_volume(index_nei)+1);
    }else if(index_this != -1){
      energy_diff = energy_after - energy_before 
	+ 1 + 2*Para::lambda*(atomic_vesicles.get_target_volume(index_this)-atomic_vesicles.get_volume(index_this));
    }else{
      energy_diff = energy_after - energy_before 
	+ 1 + 2*Para::lambda*(atomic_vesicles.get_volume(index_nei)-atomic_vesicles.get_target_volume(index_nei));
    }

    /* Third, determine if we change the state */
    if(is_copy_state(energy_diff,random_wrapper)){

      /* copy the state */
      plane.cell(row,col).vesicle_index = index_nei;

      /* update volume */
      if(index_this != -1){
	/* Note that if the model implements parallelized CPM, volume is
	   atomic. See vesicle-infor.hpp */
	atomic_vesicles.add_volume(index_this,-1);
	/* In the serial version of the algorithm, if the volume becomes
	   zero here, we delete this vesicle by setting the volume to -1
	   (done by destroby_vesicle()). However, in the parallel
	   version, this can introduce a mistake. To see this, consider
	   a case in which one thread's index_this and another thread's
	   index_nei are identical. To avoid this mistake, we have to
	   postpone destroying vesicles till all CA squares are
	   updated. This will be done in
	   MyCA::parallel_update_vesicle_whole(). */
      }
      if(index_nei != -1)
	atomic_vesicles.add_volume(index_nei);
    }
  }
}

void MyCA::ParallelDivideCollectData::operator()(const tbb::blocked_range2d<size_t>& r) const
{
  MyCA& loc_myca = my_myca;
  tbb::atomic<bool>& loc_is_divide = my_is_divide;
  const int loc_volume_threshold = Para::volume_threshold;

  for(size_t row=r.rows().begin(); row!=r.rows().end(); ++row)
    for(size_t col=r.cols().begin(); col!=r.cols().end();++col){
      if(loc_myca.plane.cell(row,col).vesicle_index>-1){
	/* Set vesicles with volume>threshold to divide, and collect
	   information to calculate PCA. */
	if(loc_myca.atomic_vesicles.get_volume(loc_myca.plane.cell(row,col).vesicle_index)>=loc_volume_threshold){
	  /* Tell that this vesicle is going to divide. Here, I don't
	     have to use compare_and_swap() to solve the concurrency
	     problem because as long as we can set is_divide=true at any
	     moment duing collecting the data for calculating principle
	     components, the alrogithm after this
	     ---i.e. ParallelDivideCollectData::operator()--- will
	     work. */
	  if(!loc_myca.atomic_vesicles.is_divide(loc_myca.plane.cell(row,col).vesicle_index)){
	    loc_myca.atomic_vesicles.set_is_divide(loc_myca.plane.cell(row,col).vesicle_index);
	  }
	  /* This is just to know that there is a division event at
	     all. */
	  if(!loc_is_divide){
	    loc_is_divide = true;
	  }
	  loc_myca.atomic_vesicles.PCA_collect_data(loc_myca.plane.cell(row,col).vesicle_index,row,col);
	}
      }
    }
}

void MyCA::ParallelDivideDivideVes::operator()(const tbb::blocked_range2d<size_t>& r) const
{
  MyCA& loc_myca = my_myca;
  
  /* Go through a (subset of a) whole plane to divide vesicles. We
     assign some squares to new vesicles */
  int mother_ind;
  for(unsigned row=r.rows().begin();row!=r.rows().end();++row)
    for(unsigned col=r.cols().begin();col!=r.cols().end();++col){
      /* Here, I use 'mother_ind' because I cannot use
	 plane.cell(row,col).vesicle_index directly, for this might be
	 changed during the division. For the same reason, don't use
	 const reference either! */
      mother_ind = loc_myca.plane.cell(row,col).vesicle_index;
      if(mother_ind>-1){
	/* Check if this vesicle has to be divided (if so, is_divide()
	   return true). */
	if(loc_myca.atomic_vesicles.is_divide(mother_ind)){
	  /* is_daughter() tells whether this square belongs to a
	     daughter vesicle. If so, we change the state of this
	     square */
	  /* Note that it may be possible that the next 'is_daughter'
	     never returns true. This happens, e.g., when the volume
	     of the mother vesicle is very small. In that case, even
	     if the mother vesicle's is_divide flag is true, it will
	     not divide. */
	  if(loc_myca.atomic_vesicles.is_daughter(mother_ind,row,col)){
	    /* If a daughter hasn't been assigned to this mother
	       vesicle, get_daughter_ind() returns -1. In that case,
	       we assign one.    
	       
	       To make it thread-safe, it appears to be a good idea to
	       combine the initial checking of daughter_ind, which is
	       done below by get_daugther_ind(), and subsequent setting
	       of it, which is done below by set_daughter_ind(), through
	       using swap_and_compare(). However, I deliverately
	       separate these two operations here for speed-up.  I think
	       using swap_and_compare() every time is a waste of
	       time. In the following method, I use swap_and_compare()
	       only if a thread sees that get_daughter_ind() returns
	       -1. Once some thread successfully sets the daughter_ind
	       to a valid value (i.e. value>-1), from thereon all
	       threads don't have to do swap_and_compare() except for
	       those rare cases in which threads access to daughter_ind
	       during the process of seting daugther_ind. We just have
	       to careful (i.e. use swap_and_compare) only in those rare
	       cases. See also AtomicVesicles::set_daughter_ind() for
	       more on this point.
	    */

	    if(loc_myca.atomic_vesicles.get_daughter_ind(mother_ind)==-1){
	      /* Control flows comes here only in a rare
		 cases. Optimally, we like it if a flow comes only once
		 for one mother_ind. But, there is a chance that flows
		 come here multiple times. set_daughter_ind() can
		 handle those rare dangerous cases. */
	      loc_myca.atomic_vesicles.set_daughter_ind(mother_ind);
	    }
	    /* A daughter index has already been assined */
	    /* This square is now belonging to the daughter */
	    loc_myca.plane.cell(row,col).vesicle_index = loc_myca.atomic_vesicles.get_daughter_ind(mother_ind);
	      
	    /* update volume */
	    loc_myca.atomic_vesicles.add_volume(mother_ind,-1);
	    loc_myca.atomic_vesicles.add_volume(loc_myca.atomic_vesicles.get_daughter_ind(mother_ind));
	    
	    /* a complex molecule that spans between mother and
	       daughter vesicles must dissociate. */
	    if(loc_myca.plane.cell(row,col).molecule.is_complex()){
	      unsigned nei_row=0,nei_col=0;
	      loc_myca.plane.XY_NEIGH_X(row,col,loc_myca.plane.cell(row,col).molecule.bon_nei,nei_row,nei_col);
	      if(!loc_myca.atomic_vesicles.is_daughter(mother_ind,nei_row,nei_col)){
		loc_myca.plane.cell(row,col).molecule.dissociate();
		loc_myca.plane.cell(nei_row,nei_col).molecule.dissociate();
	      }
	    }
	  }
	}
      }
    }
}

void MyCA::parallel_divide_vesicle(const long Time)
{
  /* Collect information to calculate principle component of vesicles
     that must be divided. */
  atomic_vesicles.PCA_reset();
  tbb::atomic<bool> is_divide;
  is_divide=false;
  tbb::parallel_for(tbb::blocked_range2d<size_t>(1,nrow+1,1,ncol+1),ParallelDivideCollectData(*this,is_divide));

  /* Calculate principle component of vesicles that must be divided */
  atomic_vesicles.calculate_PCA_for_dividing_vesicles();

  /* If there are vesicles that have to divide, is_divide is true. */
  if(is_divide){
    /* set target_volume (the mother and daughter get the fraction of
       the original target_volume of the mother proportional to their
       volume */
    tbb::parallel_for(tbb::blocked_range2d<size_t>(1,nrow+1,1,ncol+1),ParallelDivideDivideVes(*this));
    atomic_vesicles.divide_target_volume();
  }
}

/* See also the header file myca.hpp  */
void MyCA::ParallelVesicleNeutralGrowthCount::operator()(const tbb::blocked_range2d<size_t>& r) const
{
  MyCA& loc_myca = my_myca;
  
  /* We count the number of molecules in vesicles. */
  for(unsigned row=r.rows().begin();row!=r.rows().end();++row)
    for(unsigned col=r.cols().begin();col!=r.cols().end();++col){
      if(loc_myca.plane.cell(row,col).vesicle_index!=-1 
	 && !loc_myca.plane.cell(row,col).molecule.is_empty()){
	loc_myca.atomic_vesicles.add_target_volume(loc_myca.plane.cell(row,col).vesicle_index);
      }
    }
}

/* See also the header file myca.hpp  */
void MyCA::ParallelVesicleFixedGrowth::operator()(const tbb::blocked_range2d<size_t>& r) const
{
  MyCA& loc_myca = my_myca;
  const int loc_target_volume = Para::fixed_target_volume;

  /* We count the number of molecules in vesicles. */
  for(unsigned row=r.rows().begin();row!=r.rows().end();++row)
    for(unsigned col=r.cols().begin();col!=r.cols().end();++col){
      if(loc_myca.plane.cell(row,col).vesicle_index!=-1 
	 && !loc_myca.plane.cell(row,col).molecule.is_empty()){
	
	/* MyCA::parallel_vesicle_fixed_growth() first sets the target
	   volume of every vesicle to zero. So, if a vesicle contains
	   one molecule and if its target volume is still zero, set it
	   to Para::fixed_target_volume. */
	if(loc_myca.atomic_vesicles.get_target_volume(loc_myca.plane.cell(row,col).vesicle_index) == 0)
	  loc_myca.atomic_vesicles.set_target_volume(loc_myca.plane.cell(row,col).vesicle_index,loc_target_volume);
      }
    }
  
}
