#ifndef __PISMIO
#define __PISMIO

#include "nc_util.hh"
#include "LocalInterpCtx.hh"

//! A class containing IO functions used to read and write spatial variables.
class PISMIO : public NCTool {
public:
  PISMIO(IceGrid *my_grid);
  virtual ~PISMIO() {}

  using NCTool::open_for_writing;
  PetscErrorCode open_for_writing(string filename, bool append,
				  bool check_dims = false);
  PetscErrorCode get_grid(const char filename[]);
  PetscErrorCode create_dimensions() const;

  PetscErrorCode get_grid_info(grid_info &g) const;
  PetscErrorCode get_grid_info_2d(grid_info &g) const;

  PetscErrorCode get_var(int varid, Vec g, GridType dims, int t) const;
  PetscErrorCode put_var(int varid, Vec g, GridType dims) const;

  PetscErrorCode regrid_var(int varid, GridType dim_flag, LocalInterpCtx &lic, Vec g,
			    bool useMaskInterp) const;
  PetscErrorCode set_MaskInterp(MaskInterp *mi_in);

private:
  int compute_block_size(GridType dims, int* count) const;
  PetscErrorCode compute_start_and_count(int varid, int *pism_start, int *pism_count, GridType dims,
					 size_t* &nc_start, size_t* &nc_count, ptrdiff_t* &imap) const;
  bool check_dimensions() const;
  IceGrid* grid;
  MaskInterp  *myMaskInterp;
};

#endif	// __PISMIO
