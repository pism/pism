// Copyright (C) 2011-2013 Torsten Albrecht
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


#include <cmath>
#include <petscdmda.h>
#include "iceModel.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"

// FIXME: remove unused commented-out code.

//! \file iMfractures.cc implementing calculation of fracture density with PIK options -fractures.

PetscErrorCode IceModel::calculateFractureDensity() {
  const PetscScalar dx = grid.dx, dy = grid.dy, Mx = grid.Mx, My = grid.My;
  PetscErrorCode ierr;
  
  IceModelVec2S vFDnew = vWork2d[0], vFAnew = vWork2d[1];
  
  // get SSA velocities and related strain rates and stresses
  IceModelVec2V *ssa_velocity;
  ierr = stress_balance->get_2D_advective_velocity(ssa_velocity); CHKERRQ(ierr);  
  ierr = stress_balance->compute_2D_principal_strain_rates(*ssa_velocity, vMask, strain_rates);
  ierr = stress_balance->compute_2D_stresses(*ssa_velocity, vMask, deviatoric_stresses);
  ierr = ssa_velocity->begin_access(); CHKERRQ(ierr);
  ierr = strain_rates.begin_access(); CHKERRQ(ierr);
  ierr = deviatoric_stresses.begin_access(); CHKERRQ(ierr);
    
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vFD.begin_access(); CHKERRQ(ierr);
  ierr = vFD.copy_to(vFDnew); CHKERRQ(ierr);
  ierr = vFDnew.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access();  CHKERRQ(ierr);
  
  const bool dirichlet_bc = config.get_flag("ssa_dirichlet_bc");
  if (dirichlet_bc) {
    ierr = vBCMask.begin_access();  CHKERRQ(ierr);
    ierr = vBCvel.begin_access();  CHKERRQ(ierr);
  }
  
  const bool write_fd = config.get_flag("write_fd_fields");
  if (write_fd) {
    ierr = vFG.begin_access();  CHKERRQ(ierr);
    ierr = vFH.begin_access();  CHKERRQ(ierr);
    ierr = vFE.begin_access();  CHKERRQ(ierr);
    ierr = vFT.begin_access();  CHKERRQ(ierr);
    ierr = vFA.begin_access();  CHKERRQ(ierr);
    ierr = vFA.copy_to(vFAnew); CHKERRQ(ierr);
    ierr = vFAnew.begin_access(); CHKERRQ(ierr);
  } 
  
  MaskQuery M(vMask);
  PetscScalar tempFD;
  
  //options
  /////////////////////////////////////////////////////////
  PetscScalar soft_residual = 1.0;
  PetscBool fracture_soft; 
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-fracture_softening", &soft_residual, &fracture_soft);
  //assume linear response function: E_fr = (1-(1-soft_residual)*phi) -> 1-phi
  //more: T. Albrecht, A. Levermann; Fracture-induced softening for large-scale ice dynamics; (2013), 
  //The Cryosphere Discussions 7; 4501-4544; DOI:10.5194/tcd-7-4501-2013
  
  //get four options for calculation of fracture density.
  //1st: fracture growth constant gamma
  //2nd: fracture initiation stress threshold sigma_cr 
  //3rd: healing rate constant gamma_h
  //4th: healing strain rate threshold
  //more: T. Albrecht, A. Levermann; Fracture field for large-scale ice dynamics; (2012), 
  //Journal of Glaciology, Vol. 58, No. 207, 165-176, DOI: 10.3189/2012JoG11J191.
  
  PetscInt  Nparamf=4;
  PetscReal inarrayf[4] = {1.0, 7.0e4, 0.0, 2.0e-10};
  PetscBool  fractures_set;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-fractures", inarrayf, &Nparamf, &fractures_set);
  CHKERRQ(ierr);
  if ((Nparamf > 4) || (Nparamf < 4)) {
     ierr = verbPrintf(1, grid.com, "PISM ERROR: option -fractures provided with more or fewer than 4 arguments ... ENDING ...\n");CHKERRQ(ierr);
    PISMEnd();
  }
  PetscReal gamma = inarrayf[0], 
    initThreshold = inarrayf[1], 
    gammaheal     = inarrayf[2], 
    healThreshold = inarrayf[3];
  
  ierr = verbPrintf(3, grid.com,"PISM-PIK INFO: fracture density is found with parameters:\n gamma=%.2f, sigma_cr=%.2f, gammah=%.2f, healing_cr=%.1e and soft_res=%f \n", gamma,initThreshold,gammaheal,healThreshold,soft_residual); 
  CHKERRQ(ierr);
    
  //const bool do_fracground = config.get_flag("do_frac_on_grounded"); 
  PetscBool do_fracground;   
  ierr = PetscOptionsHasName(NULL,"-do_frac_on_grounded",&do_fracground); CHKERRQ(ierr);
   
  PetscScalar fdBoundaryValue = 0.0;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-phi0", &fdBoundaryValue, PETSC_NULL);
     
  PetscBool constant_healing;   
  ierr = PetscOptionsHasName(NULL,"-constant_healing",&constant_healing); CHKERRQ(ierr);
   
  PetscBool fracture_weighted_healing;   
  ierr = PetscOptionsHasName(NULL,"-fracture_weighted_healing",&fracture_weighted_healing); CHKERRQ(ierr);
   
  PetscBool max_shear_stress;   
  ierr = PetscOptionsHasName(NULL,"-max_shear",&max_shear_stress); CHKERRQ(ierr);
   
  PetscBool lefm;   
  ierr = PetscOptionsHasName(NULL,"-lefm",&lefm); CHKERRQ(ierr);
   
  PetscBool constant_fd;   
  ierr = PetscOptionsHasName(NULL,"-constant_fd",&constant_fd); CHKERRQ(ierr);
   
  PetscBool fd2d_scheme;   
  ierr = PetscOptionsHasName(NULL,"-scheme_fd2d",&fd2d_scheme); CHKERRQ(ierr);
  
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
  
      tempFD=0;
      //SSA: v . grad memField
          
      PetscScalar uvel=(*ssa_velocity)(i,j).u;
      PetscScalar vvel=(*ssa_velocity)(i,j).v;
          
      if (fd2d_scheme){
              
        if (uvel>=dx*vvel/dy && vvel>=0.0) //1
          tempFD = uvel*(vFD(i,j)-vFD(i-1,j))/dx + vvel*(vFD(i-1,j)-vFD(i-1,j-1))/dy;
        else if (uvel<=dx*vvel/dy && uvel>=0.0) //2
          tempFD = uvel*(vFD(i,j-1)-vFD(i-1,j-1))/dx + vvel*(vFD(i,j)-vFD(i,j-1))/dy;
        else if (uvel>=-dx*vvel/dy && uvel<=0.0) //3
          tempFD = -uvel*(vFD(i,j-1)-vFD(i+1,j-1))/dx + vvel*(vFD(i,j)-vFD(i,j-1))/dy;
        else if (uvel<=-dx*vvel/dy && vvel>=0.0) //4
          tempFD = -uvel*(vFD(i,j)-vFD(i+1,j))/dx + vvel*(vFD(i+1,j)-vFD(i+1,j-1))/dy;
        else if (uvel<=dx*vvel/dy && vvel<=0.0) //5
          tempFD = -uvel*(vFD(i,j)-vFD(i+1,j))/dx - vvel*(vFD(i+1,j)-vFD(i+1,j+1))/dy;
        else if (uvel>=dx*vvel/dy && uvel<=0.0) //6
          tempFD = -uvel*(vFD(i,j+1)-vFD(i+1,j+1))/dx - vvel*(vFD(i,j)-vFD(i,j+1))/dy;
        else if (uvel<=-dx*vvel/dy && uvel>=0.0) //7
          tempFD = uvel*(vFD(i,j+1)-vFD(i-1,j+1))/dx - vvel*(vFD(i,j)-vFD(i,j+1))/dy;
        else if (uvel>=-dx*vvel/dy && vvel<=0.0) //8
          tempFD = uvel*(vFD(i,j)-vFD(i-1,j))/dx - vvel*(vFD(i-1,j)-vFD(i-1,j+1))/dy;
        else{
          ierr = verbPrintf(3,grid.com,"######### missing case of angle %f of %f and %f at %d, %d \n",atan(vvel/uvel)/M_PI*180.,uvel*3e7,vvel*3e7,i,j);
        }
      }
      else{
        tempFD += uvel * (uvel<0 ? vFD(i+1,j)-vFD(i,j):vFD(i,j)-vFD(i-1,j))/dx; 
        tempFD += vvel * (vvel<0 ? vFD(i,j+1)-vFD(i,j):vFD(i,j)-vFD(i,j-1))/dy;
      }
  
      vFDnew(i,j)-= tempFD * dt;
          
      //sources /////////////////////////////////////////////////////////////////
      ///von mises criterion
      PetscScalar txx = deviatoric_stresses(i , j, 0),
                  tyy = deviatoric_stresses(i , j, 1),
                  txy = deviatoric_stresses(i , j, 2),
        
      T1 = 0.5*(txx+tyy)+sqrt(0.25*PetscSqr(txx-tyy)+PetscSqr(txy)),//Pa
      T2 = 0.5*(txx+tyy)-sqrt(0.25*PetscSqr(txx-tyy)+PetscSqr(txy)),//Pa     
      sigmat = sqrt(PetscSqr(T1)+PetscSqr(T2)-T1*T2);
  
  
      ///max shear stress criterion (more stringent than von mises)
      if (max_shear_stress){
        PetscScalar maxshear=PetscAbs(T1);
        maxshear=PetscMax(maxshear,PetscAbs(T2));
        maxshear=PetscMax(maxshear,PetscAbs(T1-T2));
  
        sigmat=maxshear;
      }
      
      ///lefm mixed-mode criterion
      if (lefm){
        // PetscScalar sigmabeta = 45.0*M_PI/180.0; //crack orientation with respect to principal stress axes
        PetscScalar sigmamu = 0.1; //friction coefficient between crack faces
        
        //PetscScalar sigmac = 0.08; //initial crack depth 8cm
        //PetscScalar sigmac = 1/M_PI; //initial crack depth 32cm
        PetscScalar sigmac = 0.64/M_PI; //initial crack depth 20cm
        
        // PetscScalar sigmateta = 0.0;
        // if ((txx(i,j)-tyy(i,j))!=0.0)
        //   sigmateta = 0.5*atan(2.0*txy(i,j)/(txx(i,j)-tyy(i,j)));
        // else
        //   sigmateta = 45.0*M_PI/180.0;
        //sigmateta = sigmateta +45.*M_PI/180;
          
        PetscScalar sigmabetatest, sigmanor, sigmatau, Kone, Ktwo,
          KSI, KSImax, sigmatetanull;
        KSImax=0.0;
        //sigmabetanull=0.0;
        for (PetscInt l = 46; l <= 90; ++l) { //optimize for various precursor angles beta 
          sigmabetatest = l*M_PI/180.0;
          
          //rist_sammonds99
          sigmanor = 0.5*(T1+T2)-(T1-T2)*cos(2*sigmabetatest);
          sigmatau = 0.5*(T1-T2)*sin(2*sigmabetatest);
          //shayam_wu90
          if (sigmamu*sigmanor<0.0) {//compressive case
            if (abs(sigmatau) <= abs(sigmamu*sigmanor))
              sigmatau=0.0;
            else {
              if (sigmatau>0) //coulomb friction opposing sliding
                sigmatau+=(sigmamu*sigmanor);
              else
                sigmatau-=(sigmamu*sigmanor);
            }
          }
            
          //stress intensity factors
          Kone = sigmanor*sqrt(M_PI*sigmac);//normal
          Ktwo = sigmatau*sqrt(M_PI*sigmac);//shear
            
          if (Ktwo==0.0)
            sigmatetanull=0.0;
          else //eq15 in hulbe_ledoux10 or eq15 shayam_wu90
            sigmatetanull=-2.0*atan((sqrt( PetscSqr(Kone) + 8.0 * PetscSqr(Ktwo) ) - Kone )/(4.0*Ktwo));
          //sigmatetanull=sigmateta;
          
          //0 = cos(0.5*sigmateta)*(Kone*sin(sigmateta) + Ktwo*(3.0*cos(sigmateta)-1.0));
          KSI = cos(0.5*sigmatetanull)*(Kone*cos(0.5*sigmatetanull)*cos(0.5*sigmatetanull) - 0.5*3.0*Ktwo*sin(sigmatetanull));
          // mode I stress intensity
   
          //if (KSI>KSImax)
          //  sigmabetanull=sigmabetatest;
          KSImax=std::max(KSI,KSImax);
          }
          sigmat=KSImax;      
          //sigmat=PetscAbs(0.5*(T1+T2)-(T1-T2)*cos(2*sigmabetanull));        
      }  
  
      //////////////////////////////////////////////////////////////////////////////
  
      //fracture density
      PetscScalar fdnew = gamma*(strain_rates(i,j,0)-0.0)*(1-vFDnew(i,j));
      if (sigmat > initThreshold) {  
        vFDnew(i,j)+= fdnew*dt; 
      }
  
      //healing
      PetscScalar fdheal = gammaheal*(strain_rates(i,j,0)-healThreshold);
      if (vH(i,j)>0.0){
        if (constant_healing){
          fdheal = gammaheal*(-healThreshold);
          if (fracture_weighted_healing)
            vFDnew(i,j)+= fdheal*dt*(1-vFD(i,j));
          else
            vFDnew(i,j)+= fdheal*dt;
        }
        else if (strain_rates(i,j,0) < healThreshold) {
          //vFDnew(i,j)+= gammaheal*(vPrinStrain1(i,j)-healThreshold)*dt;//*(1-vFD(i,j)); 
          if (fracture_weighted_healing)
            vFDnew(i,j)+= fdheal*dt*(1-vFD(i,j));
          else
            vFDnew(i,j)+= fdheal*dt; 
        }          
      }
  
      //bounding        
      if (vFDnew(i,j)<0.0)
        vFDnew(i,j)=0.0;
      if (vFDnew(i,j)>1.0)
         vFDnew(i,j)=1.0;
  
      //################################################################################
      // write related fracture quantities to nc-file
      // if option -write_fd_fields is set
      if (write_fd && vH(i,j)>0.0) {
        //fracture toughness
        vFT(i,j)=sigmat; 
  
        // fracture growth rate
        if (sigmat > initThreshold) {
          vFG(i,j)=fdnew; 
          //vFG(i,j)=gamma*(vPrinStrain1(i,j)-0.0)*(1-vFDnew(i,j)); 
        } else {
          vFG(i,j)=0.0;
        }
  
        // fracture healing rate    
        if (vH(i,j)>0.0){
          if (constant_healing || (strain_rates(i,j,0) < healThreshold)){
            if (fracture_weighted_healing)
              vFH(i,j)=fdheal*(1-vFD(i,j)); 
            else
              vFH(i,j)=fdheal;  
          }
          else
            vFH(i,j)=0.0;
        }
  
        //fracture age since fracturing occured
        vFAnew(i,j) -= dt * uvel * (uvel<0 ? vFA(i+1,j)-vFA(i, j):vFA(i, j)-vFA(i-1, j))/dx; 
        vFAnew(i,j) -= dt * vvel * (vvel<0 ? vFA(i,j+1)-vFA(i, j):vFA(i, j)-vFA(i, j-1))/dy;
        vFAnew(i,j)+= dt/grid.convert(1.0, "year", "seconds");         
        if (sigmat > initThreshold) 
          vFAnew(i,j) = 0.0;
  
        // additional flow enhancement due to fracture softening
        PetscScalar phi_exp=3.0;//flow_law->exponent();
        PetscScalar softening = pow((1.0-(1.0-soft_residual)*vFDnew(i,j)),-phi_exp);              
        if (vH(i,j)>0.0) {
          vFE(i,j)=1.0/pow(softening,1/3.0);
        }      
        else
          vFE(i,j)=1.0; 
      }  
       
      //boundary condition  
      if (dirichlet_bc && !do_fracground) {
        if (vBCMask.as_int(i,j) == 1){
          if (vBCvel(i,j).u != 0.0 || vBCvel(i,j).v != 0.0)
            vFDnew(i,j)=fdBoundaryValue;
          if (write_fd) {
            vFAnew(i,j)=0.0;
            vFG(i,j)=0.0;
            vFH(i,j)=0.0;
            vFE(i,j)=1.0;
            vFT(i,j)=0.0;
          }
        }
      }  
      // ice free regions and boundary of computational domain    
      if (vH(i,j)==0.0 || i==0 || j==0 || i==Mx-1 || j==My-1){
        vFDnew(i,j)=0.0;
        if (write_fd) {
          vFAnew(i,j)=0.0;
          vFG(i,j)=0.0;
          vFH(i,j)=0.0;
          vFE(i,j)=1.0;
          vFT(i,j)=0.0;
        }
      } 
        
      if (constant_fd) {//no fd evolution
        vFDnew(i,j)=vFD(i,j);
      }
    } 
  }
  
  ierr = ssa_velocity->end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vFD.end_access(); CHKERRQ(ierr);
  ierr = vFDnew.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access();  CHKERRQ(ierr);
  
  if (dirichlet_bc) {
    ierr = vBCMask.end_access();  CHKERRQ(ierr);
    ierr = vBCvel.end_access();  CHKERRQ(ierr);
  }
    
  if (write_fd) {
    ierr = vFG.end_access();  CHKERRQ(ierr);
    ierr = vFH.end_access();  CHKERRQ(ierr);
    ierr = vFE.end_access();  CHKERRQ(ierr);
    ierr = vFT.end_access();  CHKERRQ(ierr);
    ierr = vFA.end_access();  CHKERRQ(ierr);
    ierr = vFAnew.end_access(); CHKERRQ(ierr);
    ierr = vFAnew.update_ghosts(vFA); CHKERRQ(ierr);
  }
  
  ierr = strain_rates.end_access(); CHKERRQ(ierr);
  ierr = deviatoric_stresses.end_access(); CHKERRQ(ierr);
  ierr = vFDnew.update_ghosts(vFD); CHKERRQ(ierr);
 
  return 0;
}

