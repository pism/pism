
PetscErrorCode IceModel::assembleSSAMatrix(Vec vNu[2], Mat A) {
  const PetscInt  Mx=grid.Mx, My=grid.My, M=2*My;
  const PetscScalar   dx=grid.dx, dy=grid.dy;
  const PetscScalar   one = 1.0;
  PetscErrorCode  ierr;
  PetscScalar     **mask, **nu[2], **u, **v, **tauc;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vubarSSA, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vvbarSSA, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscInt J = 2*j;
      const PetscInt rowU = i*M + J;
      const PetscInt rowV = i*M + J+1;
      if (intMask(mask[i][j]) == MASK_SHEET) {
        ierr = MatSetValues(A, 1, &rowU, 1, &rowU, &one, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, 1, &rowV, &one, INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscInt im = (i + Mx - 1) % Mx, ip = (i + 1) % Mx;
        const PetscInt Jm = 2 * ((j + My - 1) % My), Jp = 2 * ((j + 1) % My);
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */

// these DO use thickness H because nu[][][] is actually viscosity times thickness
        const PetscScalar c00 = nu[0][i-1][j];
        const PetscScalar c01 = nu[0][i][j];
        const PetscScalar c10 = nu[1][i][j-1];
        const PetscScalar c11 = nu[1][i][j];

#if 0
        /*
        * These are stencils for
        *  - \nabla \cdot ( \nu \nabla u ) = 0
        *      and
        *  - \nabla \cdot ( \nu \nabla v ) = 0
        * They leave u and v uncoupled.
        */
        const PetscInt stencilSize = 5;
        const PetscInt colU[stencilSize] = {
          i*M+Jp,
          im*M+J,     i*M+J,      ip*M+J,
          i*M+Jm};
        const PetscScalar valU[stencilSize] = {
          -c11/dy2,
          -c00/dx2,   (c00+c01)/dx2+(c10+c11)/dy2,    -c01/dx2,
          -c10/dy2};
        const PetscInt colV[stencilSize] = {
          i*M+Jp+1,
          im*M+J+1,   i*M+J+1,    ip*M+J+1,
          i*M+Jm+1};
        const PetscScalar valV[stencilSize] = {
          -c11/dy2,
          -c00/dx2,   (c00+c01)/dx2+(c10+c11)/dy2,    -c01/dx2,
          -c10/dy2};
#else
#if 0
        /*
        * These are the stencils for constant thickness and constant viscosity.
        * They are not scaled for grid size so use 0 as RHS.
        * Only use for basic testing.
        */
        const PetscInt stencilSize = 9;
        const PetscInt colU[] = {
          im*M+Jp+1,  i*M+Jp,     ip*M+Jp+1,
          im*M+J,     i*M+J,      ip*M+J,
          im*M+Jm+1,  i*M+Jm,     ip*M+Jm+1 };
        const PetscScalar valU[] = {0.75, -1, -0.75,        -4, 10, -4,     -0.75, -1, 0.75};
        const PetscInt colV[] = {
          im*M+Jp,    i*M+Jp+1,   ip*M+Jp,
          im*M+J+1,   i*M+J+1,    ip*M+J+1,
          im*M+Jm,    i*M+Jm+1,   ip*M+Jm };
        const PetscScalar valV[] = {0.75, -4, -0.75,        -1, 10, -1,     -0.75, -4, 0.75};
#else
        const PetscInt stencilSize = 13;
        /* The locations of the stencil points for the U equation */
        const PetscInt colU[] = {
          /*       */ i*M+Jp,
          im*M+Jp+1,  i*M+Jp+1,   ip*M+Jp+1,
          im*M+J,     i*M+J,      ip*M+J,
          im*M+J+1,               ip*M+J+1,
          /*       */ i*M+Jm,
          im*M+Jm+1,  i*M+Jm+1,   ip*M+Jm+1};
        /* The values at those points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };

        /* The locations of the stencil points for the V equation */
        const PetscInt colV[] = {
          im*M+Jp,        i*M+Jp,     ip*M+Jp,
          /*           */ i*M+Jp+1,
          im*M+J,                     ip*M+J,
          im*M+J+1,       i*M+J+1,    ip*M+J+1,
          im*M+Jm,        i*M+Jm,     ip*M+Jm,
          /*           */ i*M+Jm+1 };
        /* The values at those points */
        PetscScalar valV[] = {
          (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
          /*               */ -4*c11/dy2,
          2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
          -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
          -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
          /*               */ -4*c10/dy2 };
#endif
#endif

        /* Dragging ice experiences friction at the bed determined by the
         *    basalDrag[x|y]() methods.  These may be a plastic, pseudo-plastic,
         *    or linear friction law according to basal->drag(), ultimately. */
        if (intMask(mask[i][j]) == MASK_DRAGGING) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basalDragx(tauc, u, v, i, j);
          valV[7] += basalDragy(tauc, u, v, i, j);
        }

        ierr = MatSetValues(A, 1, &rowU, stencilSize, colU, valU, INSERT_VALUES); CHKERRQ(ierr);
        ierr = MatSetValues(A, 1, &rowV, stencilSize, colV, valV, INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vubarSSA, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vvbarSSA, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vtauc, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[0], &nu[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vNu[1], &nu[1]); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  return 0;
}


