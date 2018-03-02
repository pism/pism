// THIS IS A TEST FILE TO SEPARATE THE GEOMETRY FROM THE REST OF PICO


namespace pism {
namespace ocean {


void Pico::test() {
	// TEST FUNCTION  
	m_log->message(2, "TEST...\n");

}


// used in IdentifyMask
const int Pico::imask_inner        = 2;
const int Pico::imask_outer        = 0;
const int Pico::imask_exclude      = 1;
const int Pico::imask_unidentified = -1;




// To be used solely in round_basins()
double Pico::most_frequent_element(const std::vector<double> &v)
  {   // Precondition: v is not empty
      std::map<double, double> frequencyMap;
      int maxFrequency = 0;
      double mostFrequentElement = 0;
      for (double x : v)
      {
          double f = ++frequencyMap[x];
          if (f > maxFrequency)
          {
              maxFrequency = f;
              mostFrequentElement = x;
          }
      }

      return mostFrequentElement;
  }

//! Round non-integer basin mask values to integers.

//! Basin mask can have non-integer values from PISM regridding for points that lie at
//! basin boundaries.
//! Find such point here and set them to the integer value that is most frequent next to it.
void Pico::round_basins() {

  // FIXME: THIS routine should be applied once in init, and roundbasins should
  // be stored as field (assumed the basins do not change with time).

  double id_fractional;
  std::vector<double> neighbours = {0,0,0,0};

  IceModelVec::AccessList list;
  list.add(cbasins);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // do not consider domain boundaries (they should be far from the shelves.)
    if ((i==0) | (j==0) | (i>(Mx-2)) | (j>(My-2))){
      id_fractional = 0.;
    } else {
      id_fractional = (cbasins)(i,j);
      neighbours[0] = (cbasins)(i+1,j+1);
      neighbours[1] = (cbasins)(i-1,j+1);
      neighbours[2] = (cbasins)(i-1,j-1);
      neighbours[3] = (cbasins)(i+1,j-1);

      // check if this is an interpolated number:
      // first condition: not an integer
      // second condition: has no neighbour with same value
      if ((id_fractional != round(id_fractional)) ||
          ((id_fractional != neighbours[0]) &&
          (id_fractional != neighbours[1]) &&
          (id_fractional != neighbours[2]) &&
          (id_fractional != neighbours[3]))){

        double most_frequent_neighbour = most_frequent_element(neighbours);
        (cbasins)(i,j) = most_frequent_neighbour;
        // m_log->message(2, "most frequent: %f at %d,%d\n",most_frequent_neighbour,i,j);
      }
    }

  }
}




//! Create masks that indicate ocean on continental shelf, ice rises as well as open ocean.

//! ocean_continental_shelf: ocean on the continental shelf without detached submarine islands
//! icerises: grounded ice not connected to the main ice body
//! ocean: ocean without holes in ice shelves, extends beyond continental shelf
//! lakes: subglacial lakes without access to the ocean
//! We here use a search algorithm, starting at the center or the boundary of the domain.
//! We iteratively look for regions which satisfy one of the three types named above.
void Pico::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {

  m_log->message(5, "starting identifyMASK routine\n");

  // Assume that the center of the domain belongs to main ice body.
  int seed_x = (Mx - 1)/2,
      seed_y = (My - 1)/2;

  double linner_identified = 0.0,
         all_inner_identified = 1.0,
         previous_step_identified = 0.0;

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S &m_topg = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list;
  list.add(inputmask);
  list.add(m_mask);
  list.add(m_topg);

  inputmask.set(imask_unidentified);

  // Find starting points for iteration.
  if ((masktype=="ocean_continental_shelf" || masktype=="icerises") && (seed_x >= m_grid->xs()) && (seed_x < m_grid->xs()+m_grid->xm()) && (seed_y >= m_grid->ys())&& (seed_y < m_grid->ys()+m_grid->ym())){
    inputmask(seed_x,seed_y)=imask_inner;
  }
  else if (masktype=="ocean" || masktype=="lakes"){
    //assume that some point on the domain boundary belongs to the open ocean
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if ((i==0) | (j==0) | (i>(Mx-2)) | (j>(My-2))){
        inputmask(i,j)=imask_inner;
      }
    }
  }

  // Iteratively find region which satisfies condition for coninental shelf ocean,
  // ice rise or open ocean.
  int iteration_round = 0;
  while(all_inner_identified > previous_step_identified){

    iteration_round+=1;
    previous_step_identified = all_inner_identified;

    for (Points p(*m_grid); p; p.next()) {

      const int i = p.i(), j = p.j();
      bool masktype_condition = false;

      if (masktype=="ocean_continental_shelf"){
        masktype_condition = (m_mask(i,j)!=maskocean || m_topg(i,j) >= continental_shelf_depth);}
      else if (masktype=="icerises"){
        masktype_condition = (m_mask(i,j)==maskgrounded);
      }
      else if (masktype=="ocean"){
        masktype_condition = (m_mask(i,j)==maskocean);
      }
      else if (masktype=="lakes"){
        masktype_condition = (m_mask(i,j)==maskocean || m_mask(i,j)==maskfloating);
      }

      if (masktype_condition && inputmask(i,j)==imask_unidentified &&
        (inputmask(i,j+1)==imask_inner || inputmask(i,j-1)==imask_inner ||
         inputmask(i+1,j)==imask_inner || inputmask(i-1,j)==imask_inner)){
         inputmask(i,j)=imask_inner;
         linner_identified+=1;
      }
      else if (masktype_condition == false){
        inputmask(i,j)=imask_outer;
      }
    }

    inputmask.update_ghosts();

    all_inner_identified = GlobalSum(m_grid->com, linner_identified);

  }

  // Set all unidentified grid cells to value for excluded areas (ice rises
  // or submarine islands)
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputmask(i,j)==imask_unidentified){
      inputmask(i,j)=imask_exclude;
    }

    if (masktype=="ocean_continental_shelf"){ //exclude ice covered parts
      if (m_mask(i,j)!=maskocean && inputmask(i,j) == imask_inner){
        inputmask(i,j) = imask_outer;
      }
    }
  }

}



//! Create mask that indicates indicidual ice shelves

// FIXME, this is ugly code, would be nicer and faster to use a breadth/depth first search here 
// some attempt made below, but I do not know how to access the neighboring Points...
void Pico::identify_shelf_mask() {

  m_log->message(5, "starting identify_shelf_mask routine \n");

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(shelf_mask);
  list.add(m_mask);
  list.add(lake_mask);
  if (exicerises_set) { list.add(icerise_mask); list.add(ocean_mask); }

  shelf_mask.set(0);

  std::vector<double> labels_counter (Mx*My,0);// labels_couter[i] = number of ice shelf cells with the number i, at maxim
  std::vector<double> labels_counter_global (Mx*My,0);// labels_couter[i] = number of ice shelf cells with the number i, a

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  // label all shelf cells
  while( global_continue_loop !=0 ) { 
    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition;
      if (exicerises_set) { // either floating or an ice rise..
        condition = ( (m_mask(i,j)==maskfloating  && lake_mask(i,j)!=1) || icerise_mask(i,j)==imask_exclude );
      }
      else {
        condition = (m_mask(i,j)==maskfloating && lake_mask(i,j)!=1);
      }

      if (condition) {
        
        if (shelf_mask(i,j)==0) { // if not labeled yet, set to the minimum label of any neighbor
          
          if ( shelf_mask(i-1,j)>0 && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>=shelf_mask(i-1,j) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>=shelf_mask(i-1,j) || shelf_mask(i,j-1)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i-1,j) || shelf_mask(i,j+1)==0)) {
            shelf_mask(i,j) = shelf_mask(i-1,j);
            labels_counter[static_cast<int>(shelf_mask(i-1,j))]++;
            local_continue_loop = 1;
          }
          if ( shelf_mask(i+1,j)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)> shelf_mask(i+1,j) || shelf_mask(i-1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>=shelf_mask(i+1,j) || shelf_mask(i,j-1)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i+1,j) || shelf_mask(i,j+1)==0)){

            shelf_mask(i,j) = shelf_mask(i+1,j);
            labels_counter[static_cast<int>(shelf_mask(i+1,j))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i,j-1)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)>shelf_mask(i,j-1) || shelf_mask(i-1,j)==0) && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>shelf_mask(i,j-1) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)>=shelf_mask(i,j-1) || shelf_mask(i,j+1)==0)){
            shelf_mask(i,j) = shelf_mask(i,j-1);
            labels_counter[static_cast<int>(shelf_mask(i,j-1))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i,j+1)>0 && (shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)>shelf_mask(i,j+1) || shelf_mask(i-1,j)==0) && (shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)>shelf_mask(i,j+1) || shelf_mask(i+1,j)==0) && (shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)>shelf_mask(i,j+1) || shelf_mask(i,j-1)==0)){

            shelf_mask(i,j) = shelf_mask(i,j+1);
            labels_counter[static_cast<int>(shelf_mask(i,j+1))]++;
            local_continue_loop = 1;
          }
          if (shelf_mask(i-1,j)==0 && shelf_mask(i+1,j)==0 && shelf_mask(i,j-1)==0 && shelf_mask(i,j+1)==0 ) {
            shelf_mask(i,j) = Mx*i + j; //just make sure that each shelf cell is initialized to a different value
            labels_counter[static_cast<int>(shelf_mask(i,j))]++;
            local_continue_loop = 1;
          }

        } else { // if there is a neighbor with label smaller than (i,j)
          
          if ( shelf_mask(i-1,j)>0 && shelf_mask(i-1,j)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j);
            labels_counter[static_cast<int>(shelf_mask(i-1,j))]++;
            local_continue_loop = 1;
          }          
          if ( shelf_mask(i+1,j)>0 && shelf_mask(i+1,j)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j);
            labels_counter[static_cast<int>(shelf_mask(i+1,j))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i,j-1)>0 && shelf_mask(i,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i,j-1);
            labels_counter[static_cast<int>(shelf_mask(i,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i,j+1)>0 && shelf_mask(i,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i,j+1);
            labels_counter[static_cast<int>(shelf_mask(i,j+1))]++;
            local_continue_loop = 1;
          } 
          // full eight neigbors, could also be done above but should not make a difference
          if ( shelf_mask(i-1,j-1)>0 && shelf_mask(i-1,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j-1);
            labels_counter[static_cast<int>(shelf_mask(i-1,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i-1,j+1)>0 && shelf_mask(i-1,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i-1,j+1);
            labels_counter[static_cast<int>(shelf_mask(i-1,j+1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i+1,j-1)>0 && shelf_mask(i+1,j-1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j-1);
            labels_counter[static_cast<int>(shelf_mask(i+1,j-1))]++;
            local_continue_loop = 1;
          } 
          if ( shelf_mask(i+1,j+1)>0 && shelf_mask(i+1,j+1)<shelf_mask(i,j) ){
            labels_counter[static_cast<int>(shelf_mask(i,j))]--;
            shelf_mask(i,j) = shelf_mask(i+1,j+1);
            labels_counter[static_cast<int>(shelf_mask(i+1,j+1))]++;
            local_continue_loop = 1;
          } 
        } // check whether labeled or neighboring a cell with smaller label
      } // if condition  
    } // grid
    
    shelf_mask.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while 
  
  
  for (int k=0;k<Mx*My;k++){ labels_counter_global[k] = GlobalSum(m_grid->com, labels_counter[k]);}
  
  // remove non-existing labels
  double new_label_current = 1;
  std::vector<double> new_labels (Mx*My,0); 
  
  for (int k=0;k<Mx*My;k++){
    if (labels_counter_global[k] != 0){
      new_labels[k] = new_label_current;
      new_label_current++;
    } // no else case, skip 
  }

  // set the new shelf mask labels
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int label = static_cast<int>(shelf_mask(i,j));   
    shelf_mask(i,j) = new_labels[label];
  }

  numberOfShelves = new_label_current; // total numer of shelves +1
  m_log->message(5, "Number of shelves = %d\n", numberOfShelves-1); // internally calculated with +1
}


/*  // NOTE: This was the try to use a depth-search to identify the connected components, does not work this way... 
  // NOTE: This does not work since accessing the physical neighbors of a grid point is not possible...
  // using own points instead (i,j) will probably be problematic with parallel computing...
   for (Points p(*m_grid); p; p.next()) { //FIXME correct? or do we need PointsWithGhosts?
    const int i = p.i(), j = p.j();
    // if current cell is floating and not labeled...
    if (m_mask(i,j) == maskfloating && shelf_mask(i,j)==0){
      m_log->message(5, "starting a depth-first search... \n");
      std::vector<Points> stack; // create a stack for the queue
      stack.push_back(p); // add current grid point to the stack
      while (!stack.empty()) { // as long as the stack is non-empty
        Points q = stack.back(); // get the last element and
        stack.pop_back(); // remove the last element on the stack
        const int k = q.i(), l = q.j(); // get the coordinates of q 
        if (m_mask(k,l) == maskfloating && shelf_mask(k,l)==0){
          shelf_mask(k,l) = label_current; // label as current
          // all all neigbors that are floating and not labeled to the stack
          //if (k>=1 && shelf_mask(k-1,l)==0 && m_mask(k-1,l)==maskfloating){
          Points q_new;
          q_new.i = k-1; 
          q_new.j = l;  
          stack.push_back(q_new); 
        }
      }
      label_current++;
    }
  }
*/



//! Compute for each ice shelf cell distance to grounding line and ice front

//! DistGL: distance to grounding line
//! DistIF: distance to calving front
//! Ice holes within the shelf are treated like ice shelf cells,
//! if exicerises_set, also ice rises are treated like ice shelf cells.

void Pico::compute_distances() {

  m_log->message(5, "starting compute_distances routine\n");

  double currentLabelGL = 1; // to find DistGL, 1 if floating and directly adjacent to a grounded cell
  double currentLabelIF = 1; // to find DistIF, 1 if floating and directly adjacent to an ocean cell

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(m_mask);
  list.add(DistIF);
  list.add(cbasins);
  list.add(DistGL);
  list.add(ocean_mask);

  if (exicerises_set) { list.add(icerise_mask); }

  DistGL.set(0);
  DistIF.set(0);

  // Find the grounding line and the ice front and
  // set DistGL to 1 if ice shelf cell is next to the grounding line,
  // set DistIF to 1 if ice shelf cell is next to the calving front.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (exicerises_set) {
      condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
    }
    else {
      condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
    }

    if (condition) { //if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

      // label the shelf cells adjacent to the grounding line with DistGL = 1
      bool neighbor_to_land;
      if (exicerises_set) {
        neighbor_to_land = (  icerise_mask(i,j+1)==imask_inner || icerise_mask(i,j-1)==imask_inner ||
          icerise_mask(i+1,j)==imask_inner || icerise_mask(i-1,j)==imask_inner ||
           icerise_mask(i+1,j+1)==imask_inner || icerise_mask(i+1,j-1)==imask_inner ||
           icerise_mask(i-1,j+1)==imask_inner || icerise_mask(i-1,j-1)==imask_inner );
      } else {
        neighbor_to_land = (  m_mask(i,j+1)<maskfloating || m_mask(i,j-1)<maskfloating ||
           m_mask(i+1,j)<maskfloating || m_mask(i-1,j)<maskfloating ||
          m_mask(i+1,j+1)<maskfloating || m_mask(i+1,j-1)<maskfloating ||
          m_mask(i-1,j+1)<maskfloating || m_mask(i-1,j-1)<maskfloating );
      }

      if (neighbor_to_land ){
        // i.e. there is a grounded neighboring cell (which is not ice rise!)
        DistGL(i,j) = currentLabelGL;
      } // no else

      // label the shelf cells adjacent to the calving front with DistIF = 1,
      // we do not need to exclude ice rises in this case.
      bool neighbor_to_ocean;
      neighbor_to_ocean = (ocean_mask(i,j+1)==imask_inner || ocean_mask(i,j-1)==imask_inner || ocean_mask(i+1,j)==imask_inner || ocean_mask(i-1,j)==imask_inner);

      if (neighbor_to_ocean) {
        DistIF(i,j) = currentLabelIF;
      }

    }
  }

  DistGL.update_ghosts();
  DistIF.update_ghosts();

  // DistGL calculation: Derive the distance from the grounding line for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1;
  while( global_continue_loop !=0 ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
      }

      if ( condition && DistGL(i,j)==0 &&
        (DistGL(i,j+1)==currentLabelGL || DistGL(i,j-1)==currentLabelGL ||
        DistGL(i+1,j)==currentLabelGL || DistGL(i-1,j)==currentLabelGL) ) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistGL(i,j) = currentLabelGL+1;
          local_continue_loop = 1;
      } //if

    } // for

    currentLabelGL++;
    DistGL.update_ghosts();

    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistGL

  // DistIF calculation: Derive the distance from the calving front for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1; // start loop
  while( global_continue_loop !=0  ) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (exicerises_set) {
        condition = (m_mask(i,j)==maskfloating || icerise_mask(i,j)==imask_exclude || ocean_mask(i,j)==imask_exclude);
      }
      else {
        condition = (m_mask(i,j)==maskfloating || ocean_mask(i,j)==imask_exclude);
      }

      if ( condition && DistIF(i,j)==0 &&
        (DistIF(i,j+1)==currentLabelIF || DistIF(i,j-1)==currentLabelIF ||
        DistIF(i+1,j)==currentLabelIF || DistIF(i-1,j)==currentLabelIF) ) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
          DistIF(i,j)=currentLabelIF+1;
          local_continue_loop = 1;
      } //if

    } // for

    currentLabelIF++;
    DistIF.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistIF

}





//! Compute the ocean_box_mask

//! Determine number of boxes for each basin based on max(DistGL).
//! Use a relative distance to the grounding line determine the ocean_box_mask
//! Finally, compute the extent of each ocean box in each basin.

void Pico::identify_ocean_box_mask(const Constants &cc) {

  m_log->message(5, "starting identify_ocean_box_mask routine\n");

  // Find the maximal DistGL and DistIF for each basin
  std::vector<double> max_distGL(numberOfShelves);
  std::vector<double> max_distIF(numberOfShelves);
  std::vector<double> lmax_distGL(numberOfShelves);
  std::vector<double> lmax_distIF(numberOfShelves);

  double lmax_distGL_ref = 0.0; 
  double max_distGL_ref = 0.0; 

  const IceModelVec2CellType &m_mask = *m_grid->variables().get_2d_cell_type("mask");

  for(int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){ max_distGL[shelf_id]=0.0; max_distIF[shelf_id]=0.0;lmax_distGL[shelf_id]=0.0; lmax_distIF[shelf_id]=0.0;}

  IceModelVec::AccessList list;
  list.add(shelf_mask);
  list.add(DistGL);
  list.add(DistIF);
  list.add(ocean_box_mask);
  list.add(lake_mask);
  list.add(m_mask);

  for (Points p(*m_grid); p; p.next()) {
  const int i = p.i(), j = p.j();
    int shelf_id = shelf_mask(i,j);

    if ( DistGL(i,j)> lmax_distGL[shelf_id] ) {
      lmax_distGL[shelf_id] = DistGL(i,j);
    }
    if ( DistIF(i,j)> lmax_distIF[shelf_id] ) {
      lmax_distIF[shelf_id] = DistIF(i,j);
    }
    if (DistGL(i,j)>lmax_distGL_ref){
      lmax_distGL_ref = DistGL(i,j);
    }
  }


  for (int l=0;l<numberOfShelves;l++){
    max_distGL[l] = GlobalMax(m_grid->com, lmax_distGL[l]);
    max_distIF[l] = GlobalMax(m_grid->com, lmax_distIF[l]);
  }
  max_distGL_ref = GlobalMax(m_grid->com, lmax_distGL_ref);  

  // Compute the number of boxes for each basin
  // based on maximum distance between calving front and grounding line (in DistGL)
  // this is done by interpolating between nmin=1 and nmax=numberOfBoxes
  // this will be equal to numberOfBoxes for a 'large' ice shelf

  std::vector<int> lnumberOfBoxes_perShelf(numberOfShelves);

  int n_min = 1; //
  double zeta = 0.5; // hard coded for now

  for (int l=0;l<numberOfShelves;l++){
    lnumberOfBoxes_perShelf[l] = 0;
    lnumberOfBoxes_perShelf[l] = n_min + static_cast<int>( 
    		round(pow((max_distGL[l]/max_distGL_ref), zeta) *(numberOfBoxes-n_min))); 
    lnumberOfBoxes_perShelf[l] = PetscMin(lnumberOfBoxes_perShelf[l],cc.default_numberOfBoxes);
    m_log->message(5, "lnumberOfShelves[%d]=%d \n", l, lnumberOfBoxes_perShelf[l]);
  }


  // Define the ocean boxes in ocean_box_mask
  // this is based on the relative distance to the grounding line (computed from DistGL and DistIF)
  // and the number of boxes for the basin
  ocean_box_mask.set(0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_mask(i,j)==maskfloating && DistGL(i,j)>0 && DistIF(i,j)>0 && ocean_box_mask(i,j)==0){
      int shelf_id = shelf_mask(i,j);
      int n = lnumberOfBoxes_perShelf[shelf_id];
      // relative distance between grounding line and ice front
      double r = DistGL(i,j)*1.0/(DistGL(i,j)*1.0+DistIF(i,j)*1.0);

      for(int k=0;k<n;++k){

        // define the ocean_box_mask using rule (n-k)/n< (1-r)**2 <(n-k+1)/n
        // FIXME: is there a more elegant way to ensure float?
        if ( ((n*1.0-k*1.0-1.0)/(n*1.0) <= pow((1.0-r),2)) && (pow((1.0-r), 2) <= (n*1.0-k*1.0)/n*1.0) ){

          // ensure that boxnumber of a cell cannot be bigger than the distance to the grounding line
          if (DistGL(i,j) < k+1) {
            ocean_box_mask(i,j) = DistGL(i,j);
          // if smaller or equal, use default case: set to current box number
          } else{
            ocean_box_mask(i,j) = k+1;
          }
        }//if

      } //for
    }
  } // for

  // set all floating cells which have no ocean_box_mask value to numberOfBoxes+1.
  // do not apply this to cells that are subglacial lakes, 
  // since they are not accessible for ocean waters and hence should not be treated in this model
  // For these, beckmann-goose melting will be applied, see calculate_basal_melt_missing_cells
  // those are the cells which are not reachable from grounding line or calving front.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (m_mask(i,j)==maskfloating && ocean_box_mask(i,j)==0 && lake_mask(i,j)!=1){ // floating, no sub-glacial lake
      ocean_box_mask(i,j) = numberOfBoxes + 1;
    }

  }

  // Compute the number of cells per box and shelf and save to counter_boxes.
  const int nBoxes = numberOfBoxes+2;

  // Compute the number of cells per box and shelf and save to counter_boxes.
  counter_boxes.resize(numberOfShelves, std::vector<double>(2,0));  
  std::vector<std::vector<int> > lcounter_boxes(numberOfShelves, std::vector<int>(nBoxes));

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    for (int l=0;l<nBoxes;l++){
      lcounter_boxes[shelf_id][l]=0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int box_id = static_cast<int>(round(ocean_box_mask(i,j)));
    if (box_id > 0){ // floating
      int shelf_id = shelf_mask(i,j);
      lcounter_boxes[shelf_id][box_id]++;
    }
  }

  for (int shelf_id=0;shelf_id<numberOfShelves;shelf_id++){
    counter_boxes[shelf_id].resize(nBoxes);
    for (int l=0;l<nBoxes;l++){
      counter_boxes[shelf_id][l] = GlobalSum(m_grid->com, lcounter_boxes[shelf_id][l]);
    }
  }

}



} // end of namespace ocean
} // end of namespace pism