// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petsc.h>
#include "iceModelVec.hh"
#include "columnSystem.hh"


columnSystemCtx::columnSystemCtx(int my_nmax) : nmax(my_nmax) {
  if (nmax < 1) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system too small\n");
    PetscEnd();
  }
  if (nmax > 1000000) {
    PetscPrintf(PETSC_COMM_WORLD,
      "columnSystem ERROR: nmax of system unreasonable (> 10^6)\n");
    PetscEnd();
  }
  Lp   = new PetscScalar[nmax-1];
  L    = Lp-1; // ptr arith.; note L[0]=Lp[-1] not allocated
  D    = new PetscScalar[nmax];
  U    = new PetscScalar[nmax-1];
  rhs  = new PetscScalar[nmax];
  work = new PetscScalar[nmax];
}


columnSystemCtx::~columnSystemCtx() {
  delete [] Lp;
  delete [] D;
  delete [] U;
  delete [] rhs;
  delete [] work;
}


PetscErrorCode columnSystemCtx::solveTridiagonalSystem(
                  PetscInt n, PetscScalar **x) {
  // modified slightly from Numerical Recipes version
  if (x == NULL) { SETERRQ(-999,"x is NULL in columnSystemCtx"); }
  if (*x == NULL) { SETERRQ(-998,"*x is NULL in columnSystemCtx"); }
  if (n < 1) { SETERRQ(-997,"instance size n < 1 in columnSystemCtx"); }
  if (n > nmax) { SETERRQ(-996,"instance size n too large in columnSystemCtx"); }
  PetscScalar b;
  b = D[0];
  if (b == 0.0) { return 1; }
  (*x)[0] = rhs[0]/b;
  for (int i=1; i<n; ++i) {
    work[i] = U[i-1]/b;
    b = D[i] - L[i] * work[i];
    if (b == 0.0) { return i+1; }
    (*x)[i] = (rhs[i] - L[i] * (*x)[i-1]) / b;
  }
  for (int i=n-2; i>=0; --i) {
    (*x)[i] -= work[i+1] * (*x)[i+1];
  }
  return 0;
}



ageSystemCtx::ageSystemCtx(int my_Mz)
      : columnSystemCtx(my_Mz), Mz(my_Mz) {
  callcount = 0;
  nuEQ = -1.0;
  dx = -1.0;
  dy = -1.0;
  dtTempAge = -1.0;
  zlevEQ = NULL;
  tau3 = NULL;
}


PetscErrorCode ageSystemCtx::ageSetConstants(
    PetscScalar my_dx, PetscScalar my_dy, PetscScalar my_dtTempAge,
    PetscScalar my_nuEQ,
    PetscScalar *my_zlevEQ) {

  if (callcount > 1) {
    SETERRQ(1,"ageSetConstants() should only be called once");
  }
  callcount++;

  if (my_dx <= 0.0) { SETERRQ(2,"invalid dx in ageSetConstants()"); }
  dx = my_dx;
  if (my_dy <= 0.0) { SETERRQ(3,"invalid dy in ageSetConstants()"); }
  dy = my_dy;
  if (my_dtTempAge <= 0.0) { SETERRQ(4,"invalid dtTempAge in ageSetConstants()"); }
  dtTempAge = my_dtTempAge;
  if (my_nuEQ <= 0.0) { SETERRQ(5,"invalid nuEQ in ageSetConstants()"); }
  nuEQ = my_nuEQ;

  if (my_zlevEQ == NULL) { SETERRQ(6,"zlevEQ == NULL in ageSetConstants()"); }
  zlevEQ = my_zlevEQ;
  
  return 0;
}


PetscErrorCode ageSystemCtx::ageColumnWiseSetUpAndSolve(
    PetscInt i, PetscInt j, PetscInt ks,
    PetscScalar *u, PetscScalar *v, PetscScalar *w,
    IceModelVec3 &tau3,
    PetscScalar **x) {

  PetscErrorCode ierr;
  if (callcount == 0) {
    SETERRQ(1,"ageColumnWiseSetUpAndSolve() says ageSetConstants() has not been called");
  }
  if (dx <= 0.0) { SETERRQ(2,"invalid dx in ageColumnWiseSetUpAndSolve()"); }
  if (dy <= 0.0) { SETERRQ(3,"invalid dy in ageColumnWiseSetUpAndSolve()"); }
  if (dtTempAge <= 0.0) { SETERRQ(4,"invalid dtTempAge in ageColumnWiseSetUpAndSolve()"); }
  if (nuEQ <= 0.0) { SETERRQ(5,"invalid nuEQ in ageColumnWiseSetUpAndSolve()"); }
  if (zlevEQ == NULL) { SETERRQ(6,"zlevEQ == NULL in ageColumnWiseSetUpAndSolve()"); }

  // set up system: 0 <= k < ks
  for (PetscInt k=0; k<ks; k++) {
    planeStar ss;  // note ss.ij = tau[k]
    ierr = tau3.getPlaneStarZ(i,j,zlevEQ[k],&ss); CHKERRQ(ierr);
    // do lowest-order upwinding, explicitly for horizontal
    rhs[k] =  (u[k] < 0) ? u[k] * (ss.ip1 -  ss.ij) / dx
                         : u[k] * (ss.ij  - ss.im1) / dx;
    rhs[k] += (v[k] < 0) ? v[k] * (ss.jp1 -  ss.ij) / dy
                         : v[k] * (ss.ij  - ss.jm1) / dy;
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    rhs[k] = ss.ij + dtTempAge * (1.0 - rhs[k]);

    // do lowest-order upwinding, *implicitly* for vertical
    PetscScalar AA = nuEQ * w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        L[k] = - AA;
        D[k] = 1.0 + AA;
        U[k] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        L[k] = 0.0;
        D[k] = 1.0 - AA;
        U[k] = + AA;
      }
    } else { // k == 0 case
      // note L[0] not an allocated location
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        D[0] = 1.0;
        U[0] = 0.0;
        rhs[0] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        D[0] = 1.0 - AA;
        U[0] = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < ks"
      
  // surface b.c. at ks
  if (ks>0) {
    L[ks] = 0;
    D[ks] = 1.0;   // ignore U[ks]
    rhs[ks] = 0.0;  // age zero at surface
  }
  // done setting up system

  // solve it
  ierr = solveTridiagonalSystem(ks+1,x); CHKERRQ(ierr);

  return 0;
}

