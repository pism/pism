// THIS IS A TEST FILE TO SEPARATE THE GEOMETRY FROM THE REST OF PICO

#include "Pico.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Vars.hh"

namespace pism {
namespace ocean {

void Pico::test() {
  // TEST FUNCTION
  m_log->message(2, "TEST...\n");
}

// To be used solely in round_basins()
static double most_frequent_element(const std::vector<double> &v) { // Precondition: v is not empty
  std::map<double, double> frequencyMap;
  int maxFrequency           = 0;
  double mostFrequentElement = 0;
  for (double x : v) {
    double f = ++frequencyMap[x];
    if (f > maxFrequency) {
      maxFrequency        = f;
      mostFrequentElement = x;
    }
  }

  return mostFrequentElement;
}

//! Round non-integer basin mask values to integers.

//! Basin mask can have non-integer values from PISM regridding for points that lie at
//! basin boundaries.
//! Find such point here and set them to the integer value that is most frequent next to it.
void round_basins(IceModelVec2S &basin_mask) {

  // FIXME: THIS routine should be applied once in init, and roundbasins should
  // be stored as field (assumed the basins do not change with time).

  IceGrid::ConstPtr grid = basin_mask.grid();

  int
    Mx = grid->Mx(),
    My = grid->My();

  double id_fractional;
  std::vector<double> neighbours = { 0, 0, 0, 0 };

  IceModelVec::AccessList list(basin_mask);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // do not consider domain boundaries (they should be far from the shelves.)
    if ((i == 0) | (j == 0) | (i > (Mx - 2)) | (j > (My - 2))) {
      id_fractional = 0.0;
    } else {
      id_fractional = basin_mask(i, j);
      neighbours[0] = basin_mask(i + 1, j + 1);
      neighbours[1] = basin_mask(i - 1, j + 1);
      neighbours[2] = basin_mask(i - 1, j - 1);
      neighbours[3] = basin_mask(i + 1, j - 1);

      // check if this is an interpolated number:
      // first condition: not an integer
      // second condition: has no neighbour with same value
      if ((id_fractional != round(id_fractional)) ||
          ((id_fractional != neighbours[0]) && (id_fractional != neighbours[1]) && (id_fractional != neighbours[2]) &&
           (id_fractional != neighbours[3]))) {

        basin_mask(i, j) = most_frequent_element(neighbours);
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

//! look for connected region to a certain seed grid point or region that is set at
//! the beginning. Seed can be the center fo the domain or the domain boundary.
//! this connected region will be identified with the INNER mask value.
//! all regions that are not reached by connectedness of an INNER point are set to
//! EXCLUDE. These regions are the ones of interest at the end (i think), like lakes
//! that are not connected to the ocean, or ice rises that are not connceted to the main ice sheet.
//!
//!

void Pico::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {

  m_log->message(5, "starting identifyMASK routine\n");

  // Assume that the center of the domain belongs to main ice body.
  int seed_x = (m_Mx - 1) / 2, seed_y = (m_My - 1) / 2;

  double linner_identified = 0.0, all_inner_identified = 1.0, previous_step_identified = 0.0;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S &bed         = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list{ &inputmask, &mask, &bed };

  // the algorithm will work through all unidentified grid points, until all are set with
  // some other mask value
  inputmask.set(UNIDENTIFIED);

  // Find starting points for iteration.
  // for the two mask types of ocean continental shelf or ice rises, set the point at the center
  // of the domain to the value INNER
  // the extra conditions with m_grid are not fully clear to me.
  if ((masktype == "ocean_continental_shelf" || masktype == "icerises") && (seed_x >= m_grid->xs()) &&
      (seed_x < m_grid->xs() + m_grid->xm()) && (seed_y >= m_grid->ys()) && (seed_y < m_grid->ys() + m_grid->ym())) {
    inputmask(seed_x, seed_y) = INNER;
  // for the two mask types of ocean or lakes, set all domain boundary values to INNER
  } else if (masktype == "ocean" || masktype == "lakes") {
    //assume that some point on the domain boundary belongs to the open ocean
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      if ((i == 0) | (j == 0) | (i > (m_Mx - 2)) | (j > (m_My - 2))) {
        inputmask(i, j) = INNER;
      }
    }
  }

  // Iteratively find region which satisfies condition for coninental shelf ocean,
  // ice rise or open ocean.
  int iteration_round = 0;

  double continental_shelf_depth = m_config->get_double("ocean.pico.continental_shelf_depth");

  // continue as long as new grid cells are found for which values are set
  while (all_inner_identified > previous_step_identified) {

    iteration_round += 1;
    previous_step_identified = all_inner_identified;

    for (Points p(*m_grid); p; p.next()) {

      const int i = p.i(), j = p.j();
      bool masktype_condition = false;

      // masktype==true seems to indicate a grid point that could
      // potentially belong to the inner part of a certain region.

      if (masktype == "ocean_continental_shelf") {
        // true if not being open ocean or being less deep than continential shelf depth
        // so all the ice sheet plus shelves plus shallower ocean
        masktype_condition = (mask(i, j) != MASK_ICE_FREE_OCEAN || bed(i, j) >= continental_shelf_depth);
      } else if (masktype == "icerises") {
        // true for grounded ice
        masktype_condition = (mask(i, j) == MASK_GROUNDED);
      } else if (masktype == "ocean") {
        // true for open ocean
        masktype_condition = (mask(i, j) == MASK_ICE_FREE_OCEAN);
      } else if (masktype == "lakes") {
        // true for open ocean or for floating ice
        // combinded with the domain boundary being set to INNER before,
        // lake inputmask seem to be INNER for all open ocean and all ice shelf,
        // and exclude for floating ice that could not be reached through connectedness
        // from open ocean or ice shelves.
        masktype_condition = (mask(i, j) == MASK_ICE_FREE_OCEAN || mask(i, j) == MASK_FLOATING);
      }

      // in unidentified regions with true masktype, look for neighbours that belong
      // to the inner part of a region. set to inner if at least one neighbour is inner.
      // this is '4 dot connectedness'.
      // as only here inputmask is set to inner, inner points always require
      // masktype condition true
      if (masktype_condition && inputmask(i, j) == UNIDENTIFIED &&
          (inputmask(i, j + 1) == INNER || inputmask(i, j - 1) == INNER ||
           inputmask(i + 1, j) == INNER || inputmask(i - 1, j) == INNER)) {
        inputmask(i, j) = INNER;
        linner_identified += 1;

      // next condition needs inputmask(i, j) == UNIDENTIFIED
      // so set all UNIDENTIFIED point with masktype_condition false to OUTER
      } else if (masktype_condition == false) {
        inputmask(i, j) = OUTER;
      }
    }

    inputmask.update_ghosts();

    // add if new grid points were hit that were UNIDENTIFIED before.
    all_inner_identified = GlobalSum(m_grid->com, linner_identified);
  }

  // Set all unidentified grid cells to value for excluded areas (ice rises
  // or submarine islands)
  // I think this means excluded areas are basically that what we look for.
  // (we may change naming here). They were not reached by connectedness,
  // nor were they set to OUTER before.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputmask(i, j) == UNIDENTIFIED) {
      inputmask(i, j) = EXCLUDE;
    }

    // regions that are not ice free ocean
    if (masktype == "ocean_continental_shelf") { //exclude ice covered parts
      if (mask(i, j) != MASK_ICE_FREE_OCEAN && inputmask(i, j) == INNER) {
        inputmask(i, j) = OUTER;
      }
    }
  }
}


//! Create mask that indicates indicidual ice shelves

// FIXME, this is ugly code, would be nicer and faster to use a breadth/depth first search here
// some attempt made below, but I do not know how to access the neighboring Points...
void Pico::identify_shelf_mask() {

  m_log->message(5, "starting identify_shelf_mask routine \n");

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{ &m_shelf_mask, &mask, &m_lake_mask };

  if (m_exicerises_set) {
    list.add(m_icerise_mask);
    list.add(m_ocean_mask);
  }

  m_shelf_mask.set(0);

  std::vector<double> labels_counter(m_Mx * m_My,
                                     0); // labels_couter[i] = number of ice shelf cells with the number i, at maxim
  std::vector<double> labels_counter_global(m_Mx * m_My,
                                            0); // labels_couter[i] = number of ice shelf cells with the number i, a

  double global_continue_loop = 1;
  double local_continue_loop  = 0;

  // label all shelf cells
  while (global_continue_loop != 0) {
    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition;
      if (m_exicerises_set) { // either floating or an ice rise..
        condition = ((mask(i, j) == MASK_FLOATING && m_lake_mask(i, j) != 1) || m_icerise_mask(i, j) == EXCLUDE);
      } else {
        condition = (mask(i, j) == MASK_FLOATING && m_lake_mask(i, j) != 1);
      }

      if (condition) {

        if (m_shelf_mask(i, j) == 0) { // if not labeled yet, set to the minimum label of any neighbor

          if (m_shelf_mask(i - 1, j) > 0 &&
              (m_shelf_mask(i + 1, j) > 0 && m_shelf_mask(i + 1, j) >= m_shelf_mask(i - 1, j) ||
               m_shelf_mask(i + 1, j) == 0) &&
              (m_shelf_mask(i, j - 1) > 0 && m_shelf_mask(i, j - 1) >= m_shelf_mask(i - 1, j) ||
               m_shelf_mask(i, j - 1) == 0) &&
              (m_shelf_mask(i, j + 1) > 0 && m_shelf_mask(i, j + 1) >= m_shelf_mask(i - 1, j) ||
               m_shelf_mask(i, j + 1) == 0)) {
            m_shelf_mask(i, j) = m_shelf_mask(i - 1, j);
            labels_counter[static_cast<int>(m_shelf_mask(i - 1, j))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i + 1, j) > 0 &&
              (m_shelf_mask(i - 1, j) > 0 && m_shelf_mask(i - 1, j) > m_shelf_mask(i + 1, j) ||
               m_shelf_mask(i - 1, j) == 0) &&
              (m_shelf_mask(i, j - 1) > 0 && m_shelf_mask(i, j - 1) >= m_shelf_mask(i + 1, j) ||
               m_shelf_mask(i, j - 1) == 0) &&
              (m_shelf_mask(i, j + 1) > 0 && m_shelf_mask(i, j + 1) >= m_shelf_mask(i + 1, j) ||
               m_shelf_mask(i, j + 1) == 0)) {

            m_shelf_mask(i, j) = m_shelf_mask(i + 1, j);
            labels_counter[static_cast<int>(m_shelf_mask(i + 1, j))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i, j - 1) > 0 &&
              (m_shelf_mask(i - 1, j) > 0 && m_shelf_mask(i - 1, j) > m_shelf_mask(i, j - 1) ||
               m_shelf_mask(i - 1, j) == 0) &&
              (m_shelf_mask(i + 1, j) > 0 && m_shelf_mask(i + 1, j) > m_shelf_mask(i, j - 1) ||
               m_shelf_mask(i + 1, j) == 0) &&
              (m_shelf_mask(i, j + 1) > 0 && m_shelf_mask(i, j + 1) >= m_shelf_mask(i, j - 1) ||
               m_shelf_mask(i, j + 1) == 0)) {
            m_shelf_mask(i, j) = m_shelf_mask(i, j - 1);
            labels_counter[static_cast<int>(m_shelf_mask(i, j - 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i, j + 1) > 0 &&
              (m_shelf_mask(i - 1, j) > 0 && m_shelf_mask(i - 1, j) > m_shelf_mask(i, j + 1) ||
               m_shelf_mask(i - 1, j) == 0) &&
              (m_shelf_mask(i + 1, j) > 0 && m_shelf_mask(i + 1, j) > m_shelf_mask(i, j + 1) ||
               m_shelf_mask(i + 1, j) == 0) &&
              (m_shelf_mask(i, j - 1) > 0 && m_shelf_mask(i, j - 1) > m_shelf_mask(i, j + 1) ||
               m_shelf_mask(i, j - 1) == 0)) {

            m_shelf_mask(i, j) = m_shelf_mask(i, j + 1);
            labels_counter[static_cast<int>(m_shelf_mask(i, j + 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i - 1, j) == 0 && m_shelf_mask(i + 1, j) == 0 && m_shelf_mask(i, j - 1) == 0 &&
              m_shelf_mask(i, j + 1) == 0) {
            m_shelf_mask(i, j) = m_Mx * i + j; //just make sure that each shelf cell is initialized to a different value
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]++;
            local_continue_loop = 1;
          }

        } else { // if there is a neighbor with label smaller than (i,j)

          if (m_shelf_mask(i - 1, j) > 0 && m_shelf_mask(i - 1, j) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i - 1, j);
            labels_counter[static_cast<int>(m_shelf_mask(i - 1, j))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i + 1, j) > 0 && m_shelf_mask(i + 1, j) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i + 1, j);
            labels_counter[static_cast<int>(m_shelf_mask(i + 1, j))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i, j - 1) > 0 && m_shelf_mask(i, j - 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i, j - 1);
            labels_counter[static_cast<int>(m_shelf_mask(i, j - 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i, j + 1) > 0 && m_shelf_mask(i, j + 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i, j + 1);
            labels_counter[static_cast<int>(m_shelf_mask(i, j + 1))]++;
            local_continue_loop = 1;
          }
          // full eight neigbors, could also be done above but should not make a difference
          if (m_shelf_mask(i - 1, j - 1) > 0 && m_shelf_mask(i - 1, j - 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i - 1, j - 1);
            labels_counter[static_cast<int>(m_shelf_mask(i - 1, j - 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i - 1, j + 1) > 0 && m_shelf_mask(i - 1, j + 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i - 1, j + 1);
            labels_counter[static_cast<int>(m_shelf_mask(i - 1, j + 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i + 1, j - 1) > 0 && m_shelf_mask(i + 1, j - 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i + 1, j - 1);
            labels_counter[static_cast<int>(m_shelf_mask(i + 1, j - 1))]++;
            local_continue_loop = 1;
          }
          if (m_shelf_mask(i + 1, j + 1) > 0 && m_shelf_mask(i + 1, j + 1) < m_shelf_mask(i, j)) {
            labels_counter[static_cast<int>(m_shelf_mask(i, j))]--;
            m_shelf_mask(i, j) = m_shelf_mask(i + 1, j + 1);
            labels_counter[static_cast<int>(m_shelf_mask(i + 1, j + 1))]++;
            local_continue_loop = 1;
          }
        } // check whether labeled or neighboring a cell with smaller label
      }   // if condition
    }     // grid

    m_shelf_mask.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while


  for (int k = 0; k < m_Mx * m_My; k++) {
    labels_counter_global[k] = GlobalSum(m_grid->com, labels_counter[k]);
  }

  // remove non-existing labels
  double new_label_current = 1;
  std::vector<double> new_labels(m_Mx * m_My, 0);

  for (int k = 0; k < m_Mx * m_My; k++) {
    if (labels_counter_global[k] != 0) {
      new_labels[k] = new_label_current;
      new_label_current++;
    } // no else case, skip
  }

  // set the new shelf mask labels
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int label = static_cast<int>(m_shelf_mask(i, j));
    m_shelf_mask(i, j) = new_labels[label];
  }

  m_n_shelves = new_label_current;                                // total numer of shelves +1
  m_log->message(5, "Number of shelves = %d\n", m_n_shelves - 1); // internally calculated with +1
}


/*  // NOTE: This was the try to use a depth-search to identify the connected components, does not work this way...
  // NOTE: This does not work since accessing the physical neighbors of a grid point is not possible...
  // using own points instead (i,j) will probably be problematic with parallel computing...
   for (Points p(*m_grid); p; p.next()) { //FIXME correct? or do we need PointsWithGhosts?
    const int i = p.i(), j = p.j();
    // if current cell is floating and not labeled...
    if (mask(i,j) == MASK_FLOATING && shelf_mask(i,j)==0){
      m_log->message(5, "starting a depth-first search... \n");
      std::vector<Points> stack; // create a stack for the queue
      stack.push_back(p); // add current grid point to the stack
      while (!stack.empty()) { // as long as the stack is non-empty
        Points q = stack.back(); // get the last element and
        stack.pop_back(); // remove the last element on the stack
        const int k = q.i(), l = q.j(); // get the coordinates of q
        if (mask(k,l) == MASK_FLOATING && shelf_mask(k,l)==0){
          shelf_mask(k,l) = label_current; // label as current
          // all all neigbors that are floating and not labeled to the stack
          //if (k>=1 && shelf_mask(k-1,l)==0 && mask(k-1,l)==MASK_FLOATING){
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

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{ &mask, &m_DistIF, &m_cbasins, &m_DistGL, &m_ocean_mask };

  if (m_exicerises_set) {
    list.add(m_icerise_mask);
  }

  m_DistGL.set(0);
  m_DistIF.set(0);

  // Find the grounding line and the ice front and
  // set DistGL to 1 if ice shelf cell is next to the grounding line,
  // set DistIF to 1 if ice shelf cell is next to the calving front.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    bool condition;
    if (m_exicerises_set) {
      condition =
          (mask(i, j) == MASK_FLOATING || m_icerise_mask(i, j) == EXCLUDE || m_ocean_mask(i, j) == EXCLUDE);
    } else {
      condition = (mask(i, j) == MASK_FLOATING || m_ocean_mask(i, j) == EXCLUDE);
    }

    if (condition) { //if this is an ice shelf cell (or an ice rise) or a hole in an ice shelf

      // label the shelf cells adjacent to the grounding line with DistGL = 1
      bool neighbor_to_land;
      if (m_exicerises_set) {
        neighbor_to_land =
            (m_icerise_mask(i, j + 1) == INNER || m_icerise_mask(i, j - 1) == INNER ||
             m_icerise_mask(i + 1, j) == INNER || m_icerise_mask(i - 1, j) == INNER ||
             m_icerise_mask(i + 1, j + 1) == INNER || m_icerise_mask(i + 1, j - 1) == INNER ||
             m_icerise_mask(i - 1, j + 1) == INNER || m_icerise_mask(i - 1, j - 1) == INNER);
      } else {
        neighbor_to_land =
            (mask(i, j + 1) < MASK_FLOATING || mask(i, j - 1) < MASK_FLOATING || mask(i + 1, j) < MASK_FLOATING ||
             mask(i - 1, j) < MASK_FLOATING || mask(i + 1, j + 1) < MASK_FLOATING || mask(i + 1, j - 1) < MASK_FLOATING ||
             mask(i - 1, j + 1) < MASK_FLOATING || mask(i - 1, j - 1) < MASK_FLOATING);
      }

      if (neighbor_to_land) {
        // i.e. there is a grounded neighboring cell (which is not ice rise!)
        m_DistGL(i, j) = currentLabelGL;
      } // no else

      // label the shelf cells adjacent to the calving front with DistIF = 1,
      // we do not need to exclude ice rises in this case.
      bool neighbor_to_ocean;
      neighbor_to_ocean = (m_ocean_mask(i, j + 1) == INNER || m_ocean_mask(i, j - 1) == INNER ||
                           m_ocean_mask(i + 1, j) == INNER || m_ocean_mask(i - 1, j) == INNER);

      if (neighbor_to_ocean) {
        m_DistIF(i, j) = currentLabelIF;
      }
    }
  }

  m_DistGL.update_ghosts();
  m_DistIF.update_ghosts();

  // DistGL calculation: Derive the distance from the grounding line for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1;
  while (global_continue_loop != 0) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (m_exicerises_set) {
        condition = (mask(i, j) == MASK_FLOATING || m_icerise_mask(i, j) == EXCLUDE ||
                     m_ocean_mask(i, j) == EXCLUDE);
      } else {
        condition = (mask(i, j) == MASK_FLOATING || m_ocean_mask(i, j) == EXCLUDE);
      }

      if (condition && m_DistGL(i, j) == 0 &&
          (m_DistGL(i, j + 1) == currentLabelGL || m_DistGL(i, j - 1) == currentLabelGL ||
           m_DistGL(i + 1, j) == currentLabelGL || m_DistGL(i - 1, j) == currentLabelGL)) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
        m_DistGL(i, j) = currentLabelGL + 1;
        local_continue_loop = 1;
      } //if

    } // for

    currentLabelGL++;
    m_DistGL.update_ghosts();

    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistGL

  // DistIF calculation: Derive the distance from the calving front for
  // all ice shelf cells iteratively.
  // Ice holes within the shelf are treated like ice shelf cells,
  // if exicerises_set, also ice rises are treated like ice shelf cells.
  global_continue_loop = 1; // start loop
  while (global_continue_loop != 0) {

    local_continue_loop = 0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool condition; // this cell is floating or an hole in the ice shelf (or an ice rise)
      if (m_exicerises_set) {
        condition = (mask(i, j) == MASK_FLOATING || m_icerise_mask(i, j) == EXCLUDE ||
                     m_ocean_mask(i, j) == EXCLUDE);
      } else {
        condition = (mask(i, j) == MASK_FLOATING || m_ocean_mask(i, j) == EXCLUDE);
      }

      if (condition && m_DistIF(i, j) == 0 &&
          (m_DistIF(i, j + 1) == currentLabelIF || m_DistIF(i, j - 1) == currentLabelIF ||
           m_DistIF(i + 1, j) == currentLabelIF || m_DistIF(i - 1, j) == currentLabelIF)) {
        // i.e. this is an shelf cell with no distance assigned yet and with a neighbor that has a distance assigned
        m_DistIF(i, j) = currentLabelIF + 1;
        local_continue_loop = 1;
      } //if

    } // for

    currentLabelIF++;
    m_DistIF.update_ghosts();
    global_continue_loop = GlobalMax(m_grid->com, local_continue_loop);

  } // while: find DistIF
}


//! Compute the ocean_box_mask

//! Determine number of boxes for each basin based on max(DistGL).
//! Use a relative distance to the grounding line determine the ocean_box_mask
//! Finally, compute the extent of each ocean box in each basin.

void Pico::identify_ocean_box_mask(const BoxModel &cc) {

  m_log->message(5, "starting identify_ocean_box_mask routine\n");

  // Find the maximal DistGL and DistIF for each basin
  std::vector<double> max_distGL(m_n_shelves, 0.0);
  std::vector<double> max_distIF(m_n_shelves, 0.0);
  std::vector<double> lmax_distGL(m_n_shelves, 0.0);
  std::vector<double> lmax_distIF(m_n_shelves, 0.0);

  // distGL describes distance to the ice front
  // distIF describes distance to the calving front
  double lmax_distGL_ref = 0.0;
  double max_distGL_ref  = 0.0;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{ &m_shelf_mask, &m_DistGL, &m_DistIF, &m_ocean_box_mask, &m_lake_mask, &mask };

  // find the maximum distance to the grounding line and the maximum
  // distance to the calving for each ice shelf (identified by shelf_id).
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int shelf_id = m_shelf_mask(i, j);

    if (m_DistGL(i, j) > lmax_distGL[shelf_id]) {
      lmax_distGL[shelf_id] = m_DistGL(i, j);
    }
    if (m_DistIF(i, j) > lmax_distIF[shelf_id]) {
      lmax_distIF[shelf_id] = m_DistIF(i, j);
    }
    // find lmax_distGL_ref as the maximum distance to the grounding line
    // per processor domain.
    if (m_DistGL(i, j) > lmax_distGL_ref) {
      lmax_distGL_ref = m_DistGL(i, j);
    }
  }


  for (int l = 0; l < m_n_shelves; l++) {
    max_distGL[l] = GlobalMax(m_grid->com, lmax_distGL[l]);
    max_distIF[l] = GlobalMax(m_grid->com, lmax_distIF[l]);
  }
  // find the maximum grounding line distance the whole domain.
  max_distGL_ref = GlobalMax(m_grid->com, lmax_distGL_ref);

  // Compute the number of boxes for each basin
  // based on maximum distance between calving front and grounding line (in DistGL)
  // this is done by interpolating between nmin=1 and nmax=numberOfBoxes
  // this will be equal to numberOfBoxes for a 'large' ice shelf

  std::vector<int> lnumberOfBoxes_perShelf(m_n_shelves);

  int n_min   = 1;   //
  double zeta = 0.5; // hard coded for now

  int number_of_boxes = m_config->get_double("ocean.pico.number_of_boxes");

  for (int l = 0; l < m_n_shelves; l++) {
    lnumberOfBoxes_perShelf[l] = 0;

    // equation (9) of PICO paper https://doi.org/10.5194/tc-2017-70
    lnumberOfBoxes_perShelf[l] =
        n_min + static_cast<int>(round(pow((max_distGL[l] / max_distGL_ref), zeta) * (m_numberOfBoxes - n_min)));

    // never have more boxes than default number of boxes
    lnumberOfBoxes_perShelf[l] = PetscMin(lnumberOfBoxes_perShelf[l], number_of_boxes);
    m_log->message(5, "lnumberOfShelves[%d]=%d \n", l, lnumberOfBoxes_perShelf[l]);
  }


  // Define the ocean boxes in ocean_box_mask
  // this is based on the relative distance to the grounding line (computed from DistGL and DistIF)
  // and the number of boxes for the basin
  m_ocean_box_mask.set(0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    // only enter calculation if ice is floating, grounding line and ice front distance are larger zero,
    // and the ocean box mask has not been changed from zero before.
    if (mask(i, j) == MASK_FLOATING && m_DistGL(i, j) > 0 && m_DistIF(i, j) > 0 && m_ocean_box_mask(i, j) == 0) {
      int shelf_id = m_shelf_mask(i, j);
      int n        = lnumberOfBoxes_perShelf[shelf_id];
      // relative distance between grounding line and ice front
      // this distance is non-dimensional
      // equation (10) of PICO paper https://doi.org/10.5194/tc-2017-70
      double r = m_DistGL(i, j) * 1.0 / (m_DistGL(i, j) * 1.0 + m_DistIF(i, j) * 1.0);

      // loop over number of boxes for this shelf.
      for (int k = 0; k < n; ++k) {

        // define the ocean_box_mask using rule (n-k)/n< (1-r)**2 <(n-k+1)/n
        // FIXME: is there a more elegant way to ensure float?
        // equation (11) of PICO paper https://doi.org/10.5194/tc-2017-70
        // (minus 1 an squared of the version written in the paper.)
        // see also Appendix B of the paper, which motivates equation 11.
        if (((n * 1.0 - k * 1.0 - 1.0) / (n * 1.0) <= pow((1.0 - r), 2)) &&
            (pow((1.0 - r), 2) <= (n * 1.0 - k * 1.0) / n * 1.0)) {

          // ensure that boxnumber of a cell cannot be bigger than the distance to the grounding line
          if (m_DistGL(i, j) < k + 1) {
            m_ocean_box_mask(i, j) = m_DistGL(i, j);
            // if smaller or equal, use default case: set to current box number
          } else {
            m_ocean_box_mask(i, j) = k + 1;
          }
        } //if

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
    if (mask(i, j) == MASK_FLOATING && m_ocean_box_mask(i, j) == 0 &&
        m_lake_mask(i, j) != 1) { // floating, no sub-glacial lake
      m_ocean_box_mask(i, j) = m_numberOfBoxes + 1;
    }
  }

  // Compute the number of cells per box and shelf and save to counter_boxes.
  const int nBoxes = m_numberOfBoxes + 2;

  // Compute the number of cells per box and shelf and save to counter_boxes.
  // counter_boxes is used in Pico.cc to determine the area covered by a certain box by
  // counber_boxes * dx * dy
  counter_boxes.resize(m_n_shelves, std::vector<double>(2, 0));
  std::vector<std::vector<int> > lcounter_boxes(m_n_shelves, std::vector<int>(nBoxes));

  for (int shelf_id = 0; shelf_id < m_n_shelves; shelf_id++) {
    for (int l = 0; l < nBoxes; l++) {
      lcounter_boxes[shelf_id][l] = 0;
    }
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    int box_id = static_cast<int>(round(m_ocean_box_mask(i, j)));
    if (box_id > 0) { // floating
      int shelf_id = m_shelf_mask(i, j);
      lcounter_boxes[shelf_id][box_id]++;
    }
  }

  for (int shelf_id = 0; shelf_id < m_n_shelves; shelf_id++) {
    counter_boxes[shelf_id].resize(nBoxes);
    for (int l = 0; l < nBoxes; l++) {
      counter_boxes[shelf_id][l] = GlobalSum(m_grid->com, lcounter_boxes[shelf_id][l]);
    }
  }
}


} // end of namespace ocean
} // end of namespace pism
