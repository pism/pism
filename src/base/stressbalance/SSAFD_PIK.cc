// Copyright (C) 2011 Torsten Albrecht and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "SSAFD.hh"

PetscErrorCode SSAFD_PIK::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = SSAFD::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
                    "  [... including PIK CFBC implementation]\n"); CHKERRQ(ierr);  
}

PetscErrorCode SSAFD_PIK::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode ierr;

  ierr = verbPrintf(3,grid.com, "SSAFD_PIK:assemble_matrix is called\n"); CHKERRQ(ierr);

  // put the new matrix assembly here

  const PetscScalar   dx=grid.dx, dy=grid.dy;
  // next constant not too sensitive, but must match value in assembleSSARhs():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream
  IceModelVec2V vel = velocity;         // a shortcut

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  PetscReal beta_shelves_drag_too = config.get("beta_shelves_drag_too");
  bool shelvesDragToo = config.get_flag("shelves_drag_too");

  /* matrix assembly loop */

  ierr = nuH.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);

  //IceModelVec2S &thk = *thickness; 
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  //ierr = bed->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PismMask mask_value = mask->value(i,j);

	  PetscScalar   Ho = (*thickness)(i,j),
					He = (*thickness)(i+1,j),
					Hw = (*thickness)(i-1,j),
					Hn = (*thickness)(i,j+1),
					Hs = (*thickness)(i,j-1);				
		
	  PetscTruth onIcefreeOcean=PETSC_FALSE;
	  if (Ho<=1.0) {  
		onIcefreeOcean=PETSC_TRUE;
	  }		
	
	  PetscTruth atBoundary=PETSC_FALSE;
	  if (Ho>100.0 && (He<=1.0 || Hw<=1.0 || Hs<=1.0 || Hn<=1.0)) {  //defined so far via ice thickness
		atBoundary=PETSC_TRUE;
	  }




	  if (onIcefreeOcean) { // vanish ice velocities on the ice free ocean
		MatStencil  row, col;
        row.j = i; row.i = j; row.c = 0;
        col.j = i; col.i = j; col.c = 0;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
        row.c = 1;
        col.c = 1;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);

      } else if (atBoundary) { 
	
		const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
		
		const PetscScalar c00 = nuH(i-1,j,0);
        const PetscScalar c01 = nuH(i,j,0);
        const PetscScalar c10 = nuH(i,j-1,1);
        const PetscScalar c11 = nuH(i,j,1);

		//in i-direction 
		PetscInt 	aMn=1, aPn=1,
	  				aMM=1, aPP=1,
	  				aMs=1, aPs=1;
	
		//in -j-direction
		PetscInt 	bPw=1, bPP=1, bPe=1,
	  				bMw=1, bMM=1, bMe=1;


		//direct adjacent neighbors
	   	if (He<=1.0) aPP=0;
	    if (Hw<=1.0) aMM=0;
		if (Hn<=1.0) bPP=0;
	    if (Hs<=1.0) bMM=0;
	
		//neighbors in the corners
	    PetscScalar Hne = (*thickness)(i+1,j+1),
					Hse = (*thickness)(i+1,j-1),
					Hnw = (*thickness)(i-1,j+1),
					Hsw = (*thickness)(i-1,j-1);

		//for the each single boundary to decide, which derivative to drop
		if (Hn<=1.0 || Hne<=1.0) aPn=0;
		if (He<=1.0 || Hne<=1.0) bPe=0;
		if (He<=1.0 || Hse<=1.0) bMe=0;
		if (Hs<=1.0 || Hse<=1.0) aPs=0;
		if (Hs<=1.0 || Hsw<=1.0) aMs=0;
	    if (Hw<=1.0 || Hsw<=1.0) bMw=0;
		if (Hw<=1.0 || Hnw<=1.0) bPw=0;
		if (Hn<=1.0 || Hnw<=1.0) aMn=0;

        const PetscInt sten = 14;//one more than default
        MatStencil  row, col[sten];

		// This complicated stencil is the result of adding factors (0 if boundary or 1 if not) for every single derivatice of the SSA equations across one of the 4 direct or 8 neighboring grid cell boundaries.
		// If this ice grid cell was surrounded only by other ice grid cells, we would get the standard stencil of the SSA.
		// Derivatives across the 4 direct neighbors are replaced by hydrostatic pressure on the rhs.

	    PetscScalar valU[] = {

	    /*                                           */ -bPP*c11/dy2,
	    (2*bPw*aMM*c00+aMn*bPP*c11)/d4,                 (2*bPP*(c00*aMM-c01*aPP)+c11*bPP*(aPn-aMn))/d4,                -(2*bPe*aPP*c01+aPn*bPP*c11)/d4,
	    -4*aMM*c00/dx2,                                  4*(aPP*c01+aMM*c00)/dx2+(bPP*c11+bMM*c10)/dy2,                 -4*aPP*c01/dx2,
	    (aMM*(bPP*c11-bMM*c10)+2*c00*aMM*(bMw-bPw))/d4, (2*(c01*aPP-c00*aMM)*(bPP-bMM)+(c11*bPP-c10*bMM)*(aPP-aMM))/d4, (aPP*(c10*bMM-c11*bPP)+2*c01*aPP*(bPe-bMe))/d4,
	    /*                                           */ -bMM*c10/dy2,
	    -(2*bMw*aMM*c00+aMs*bMM*c10)/d4,                (2*bMM*(aPP*c01-c00*aMM)+c10*bMM*(aMs-aPs))/d4,                 (2*bMe*aPP*c01+aPs*bMM*c10)/d4};


	    PetscScalar valV[] = {

	    (2*aMn*bPP*c11+bPw*aMM*c00)/d4,                 (bPP*(c00*aMM-aPP*c01)+2*c11*bPP*(aPn-aMn))/d4,                -(2*aPn*bPP*c11+bPe*aPP*c01)/d4,
	    /*                                           */ -4*bPP*c11/dy2,
	    (2*aMM*(bPP*c11-c10*bMM)+c00*aMM*(bMw-bPw))/d4, (2*(bPP*c11-c10*bMM)*(aPP-aMM)+(aPP*c01-c00*aMM)*(bPP-bMM))/d4, (2*aPP*(c10*bMM-c11*bPP)+c01*aPP*(bPe-bMe))/d4,
	    -aMM*c00/dx2,                                    4*(bPP*c11+bMM*c10)/dy2+(aPP*c01+aMM*c00)/dx2,                 -aPP*c01/dx2,
	    -(2*aMs*bMM*c10+bMw*aMM*c00)/d4,                (bMM*(aPP*c01-c00*aMM)+2*c10*bMM*(aMs-aPs))/d4,                 (2*aPs*bMM*c10+bMe*aPP*c01)/d4,
	    /*                                           */ -4*bMM*c10/dy2 };



    	if (include_basal_shear && (mask_value == MASK_DRAGGING_SHEET)) {
      		// Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
      		valU[5] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
      		valV[8] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);// without bc it would be the 7th point of the stencil
    	}

        // make shelf drag a little bit if desired
        if (shelvesDragToo && (mask_value == MASK_FLOATING)) {
          //ierr = verbPrintf(1,grid.com,"... SHELF IS DRAGGING ..."); CHKERRQ(ierr);
          valU[5] += beta_shelves_drag_too;
          valV[8] += beta_shelves_drag_too;// without bc it would be the 7th point of the stencil
        }

	 	// build "u" equation: NOTE TRANSPOSE
     	row.j = i; row.i = j; row.c = 0;
     	const PetscInt UI[] = {
       	/*       */ i,
       	i-1,        i,          i+1,
       	i-1,        i,          i+1,
       	i-1,        i,          i+1,
       	/*       */ i,
       	i-1,        i,          i+1};
     	const PetscInt UJ[] = {
       	/*       */ j+1,
       	j+1,        j+1,        j+1,
       	j,          j,          j,
       	j,          j,          j,
       	/*       */ j-1,
       	j-1,        j-1,        j-1};
     	const PetscInt UC[] = {
       	/*       */ 0,
       	1,          1,          1,
       	0,          0,          0,
       	1,          1,          1,
       	/*       */ 0,
       	1,          1,          1};
     	for (PetscInt m=0; m<sten; m++) {
       		col[m].j = UI[m]; col[m].i = UJ[m], col[m].c = UC[m];
     	}
     	ierr = MatSetValuesStencil(A,1,&row,sten,col,valU,INSERT_VALUES); CHKERRQ(ierr);

     // build "v" equation: NOTE TRANSPOSE
     row.j = i; row.i = j; row.c = 1;
     const PetscInt VI[] = {
       i-1,        i,          i+1,
       /*       */ i,
       i-1,        i,          i+1,
       i-1,        i,          i+1,
       i-1,        i,          i+1,
       /*       */ i};
     const PetscInt VJ[] = {
       j+1,        j+1,        j+1,
       /*       */ j+1,
       j,          j,          j,
       j,          j,          j,
       j-1,        j-1,        j-1,
       /*       */ j-1};
     const PetscInt VC[] = {
       0,          0,          0,
       /*       */ 1,
       0,          0,          0,
       1,          1,          1,
       0,          0,          0,
       /*       */ 1};
     for (PetscInt m=0; m<sten; m++) {
       col[m].j = VI[m]; col[m].i = VJ[m], col[m].c = VC[m];
     }
     ierr = MatSetValuesStencil(A,1,&row,sten,col,valV,INSERT_VALUES); CHKERRQ(ierr);
	
	//the rest is just a copy of SSAFD
      } else if (mask_value == MASK_SHEET) { // FIXME: replace with MASK_BC
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where boundary value to SSA is set
        MatStencil  row, col;
        row.j = i; row.i = j; row.c = 0;
        col.j = i; col.i = j; col.c = 0;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
        row.c = 1;
        col.c = 1;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients  nu H:
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */
        const PetscScalar c00 = nuH(i-1,j,0);
        const PetscScalar c01 = nuH(i,j,0);
        const PetscScalar c10 = nuH(i,j-1,1);
        const PetscScalar c11 = nuH(i,j,1);

        const PetscInt sten = 13;
        MatStencil  row, col[sten];

        /* start with the values at the points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };
        PetscScalar valV[] = {
          (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
          /*               */ -4*c11/dy2,
          2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
          -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
          -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
          /*               */ -4*c10/dy2 };

        /* Dragging ice experiences friction at the bed determined by the
         *    basalDrag[x|y]() methods.  These may be a plastic, pseudo-plastic,
         *    or linear friction law according to basal->drag(), which gets called
         *    by basalDragx(),basalDragy().  */
        if (include_basal_shear && (mask_value == MASK_DRAGGING_SHEET)) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
          valV[7] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
        }

        // make shelf drag a little bit if desired
        if (shelvesDragToo && (mask_value == MASK_FLOATING)) {
          //ierr = verbPrintf(1,grid.com,"... SHELF IS DRAGGING ..."); CHKERRQ(ierr);
          valU[5] += beta_shelves_drag_too;
          valV[7] += beta_shelves_drag_too;
        }


        // build "u" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 0;
        const PetscInt UI[] = {
          /*       */ i,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          i-1,                    i+1,
          /*       */ i,
          i-1,        i,          i+1};
        const PetscInt UJ[] = {
          /*       */ j+1,
          j+1,        j+1,        j+1,
          j,          j,          j,
          j,                      j,
          /*       */ j-1,
          j-1,        j-1,        j-1};
        const PetscInt UC[] = {
          /*       */ 0,
          1,          1,          1,
          0,          0,          0,
          1,                      1,
          /*       */ 0,
          1,          1,          1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = UI[m]; col[m].i = UJ[m], col[m].c = UC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valU,INSERT_VALUES); CHKERRQ(ierr);

        // build "v" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 1;
        const PetscInt VI[] = {
          i-1,        i,          i+1,
          /*       */ i,
          i-1,                    i+1,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          /*       */ i};
        const PetscInt VJ[] = {
          j+1,        j+1,        j+1,
          /*       */ j+1,
          j,                      j,
          j,          j,          j,
          j-1,        j-1,        j-1,
          /*       */ j-1};
        const PetscInt VC[] = {
          0,          0,          0,
          /*       */ 1,
          0,                      0,
          1,          1,          1,
          0,          0,          0,
          /*       */ 1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = VI[m]; col[m].i = VJ[m], col[m].c = VC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valV,INSERT_VALUES); CHKERRQ(ierr);


      }
    }
  }

  ierr = vel.end_access(); CHKERRQ(ierr);  
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);  
  ierr = nuH.end_access(); CHKERRQ(ierr);

  ierr = thickness->end_access(); CHKERRQ(ierr);
  //ierr = bed->end_access(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0;
}





PetscErrorCode SSAFD_PIK::assemble_rhs(Vec rhs) {
  	PetscErrorCode ierr;


   ierr = verbPrintf(3,grid.com, "SSAFD_PIK:assemble_rhs is called\n"); CHKERRQ(ierr);

	const PetscScalar   dx=grid.dx, dy=grid.dy;
	PISMVector2     **rhs_uv;

	// next constant not too sensitive, but must match value in assembleSSAMatrix():
	const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream;

	ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

	// get driving stress components
	ierr = compute_driving_stress(taud); CHKERRQ(ierr);

	ierr = taud.begin_access(); CHKERRQ(ierr);
	ierr = DAVecGetArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);
  	
	if (vel_bc && bc_locations) {
	  ierr = vel_bc->begin_access(); CHKERRQ(ierr);
	  ierr = bc_locations->begin_access(); CHKERRQ(ierr);
	}
	
	ierr = thickness->begin_access(); CHKERRQ(ierr);
    ierr = bed->begin_access(); CHKERRQ(ierr);


	for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
	  for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	 	PetscScalar   	Ho = (*thickness)(i,j),
						He = (*thickness)(i+1,j),
						Hw = (*thickness)(i-1,j),
						Hn = (*thickness)(i,j+1),
						Hs = (*thickness)(i,j-1);				
		
	  PetscTruth onIcefreeOcean=PETSC_FALSE;
	  if (Ho<=1.0) {  
		onIcefreeOcean=PETSC_TRUE;
	  }			
	  PetscTruth atBoundary=PETSC_FALSE;
	  if (Ho>100.0 && (He<=1.0 || Hw<=1.0 || Hs<=1.0 || Hn<=1.0)) {  
		atBoundary=PETSC_TRUE;
	  }
	
	
	  if (onIcefreeOcean) {
		  rhs_uv[i][j].u = 0.0;
	      rhs_uv[i][j].v = 0.0;

      } else if (atBoundary) {
	 	PetscInt aMM=1, aPP=1, bMM=1,bPP=1;
		//direct adjacent neighbors
	   	if (He<=1.0) aPP=0;
	    if (Hw<=1.0) aMM=0;
		if (Hn<=1.0) bPP=0;
	    if (Hs<=1.0) bMM=0;
	
	
	  	const double standard_gravity = config.get("standard_gravity");
		double ocean_rho = config.get("sea_water_density");
	
	  	const 	PetscScalar icepressure = ice.rho * standard_gravity * Ho;	  
	  			PetscScalar oceanpressure;
	
	    //if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
		PetscReal currentSeaLevel=0.0; //FIXME
	  	//ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);
	
	    const PetscScalar 	hgrounded = (*bed)(i,j) + Ho,
							hfloating = currentSeaLevel + (1.0 - ice.rho/ocean_rho) * Ho,
							Ho2 = Ho*Ho;
							
			 PetscScalar	ho=0.0,
							tdx=0.0, tdy=0.0;
	
	
	  	if (Ho>0.0 && (*bed)(i,j)<(currentSeaLevel-(ice.rho/ocean_rho) * Ho)) { //calving front boundary condition for floating shelf
	    	oceanpressure = 0.5 * ice.rho * standard_gravity * (1-(ice.rho/ocean_rho))*Ho2; // this is not really the oceanpressure, but the difference between oceanpressure and isotrop.normal stresses (=pressure) from within the ice
			ho=hfloating;
			//ierr = verbPrintf(3,grid.com, "SSAFD_PIK_INFO: oceanpressure at i=%d,j=%d equals=%e\n",i,j,oceanpressure); CHKERRQ(ierr);
		} else { 
			ho=hgrounded;
            if( (*bed)(i,j) >= currentSeaLevel){//boundary condition for cliff --> zero stress = oceanpressure
	      		oceanpressure = 0.5 * ice.rho * standard_gravity * Ho2; // this is not 'zero' because the isotrop.normal stresses (=pressure) from within the ice figures on RHS           	
            }else{//boundary condition for marine terminating glacier
	      		oceanpressure = 0.5 * ice.rho * standard_gravity * (Ho2-(ocean_rho/ice.rho)*(currentSeaLevel-(*bed)(i,j))*(currentSeaLevel-(*bed)(i,j)));
	        }
	    }
	
	
	
	  if (aPP==0 && aMM==1) tdx= icepressure*ho/dx; //here we take the direct gradient at the boundary (not centered)
	  else if (aMM==0 && aPP==1) tdx= -icepressure*ho/dx;
	  else if (aPP==0 && aMM==0) tdx= 0; //in case of some kind of ice nose, or ice bridge
	
	  if (bPP==0 && bMM==1) tdy= icepressure*ho/dy;
	  else if (bMM==0 && bPP==1) tdy= -icepressure*ho/dy;
	  else if (bPP==0 && bMM==0) tdy= 0;
	
	  rhs_uv[i][j].u = tdx-(aMM-aPP)*oceanpressure/dx;
      rhs_uv[i][j].v = tdy-(bMM-bPP)*oceanpressure/dy;
	
	
	
	//the rest is just a copy of SSAFD
	  } else if (vel_bc && (bc_locations->value(i,j) == MASK_SHEET)) { // FIXME: replace with MASK_BC
		  rhs_uv[i][j].u = scaling * (*vel_bc)(i,j).u;
	      rhs_uv[i][j].v = scaling * (*vel_bc)(i,j).v;
      } else {
	// usual case: use already computed driving stress
      rhs_uv[i][j].u = taud(i,j).u;
      rhs_uv[i][j].v = taud(i,j).v;
	  }
	}
	}
	
	
	if (vel_bc) {
	  ierr = bc_locations->end_access(); CHKERRQ(ierr);
	  ierr = vel_bc->end_access(); CHKERRQ(ierr);
	}
	
	ierr = thickness->end_access(); CHKERRQ(ierr);
    ierr = bed->end_access(); CHKERRQ(ierr);

	ierr = taud.end_access(); CHKERRQ(ierr);
	ierr = DAVecRestoreArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);

	ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);
	
  return 0;
}



