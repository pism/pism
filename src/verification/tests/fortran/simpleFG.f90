program simpleFG
!   Copyright (C) 2005-2007 Ed Bueler
!
!   This file is part of PISM.
!
!   PISM is free software; you can redistribute it and/or modify it under the
!   terms of the GNU General Public License as published by the Free Software
!   Foundation; either version 3 of the License, or (at your option) any later
!   version.
!
!   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
!   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
!   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
!   details.
!
!   You should have received a copy of the GNU General Public License
!   along with PISM; if not, write to the Free Software
!   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Runs testsFG; see reference there.
!
! Compilation and linking example:
!  $ gfortran -c testsFG.f90
!  $ gfortran -c simpleFG.f90
!  $ gfortran testsFG.o simpleFG.o -o simpleFG
!
! Example dialog:
!	$ ./simpleFG
!	 Enter  t, r, z  (in a, km, m, resp.):
!	500 500 0
!   Result from Test F (functions of r only):
!	       H =   1925.29529003980       (m)
!	       M = -1.050983313290947E-002  (m/a)
!	 (functions of r,z):
!	       T =   265.122614956470       (deg K)
!	       U = -6.123215695366028E-017  (m/a)
!	       w =  1.153402815543494E-018  (m/a)
!	   Sigma =  0.264346036328259       (10^(-3) deg K/a)
!	 Sigma_c = -0.373725522145257       (10^(-3) deg K/a)
!
!	 Result from Test G (functions of t,r only):
!	       H =   2101.89975093468       (m)
!	       M =  4.073801873505128E-002  (m/a)
!	 (functions of t,r,z):
!	       T =   267.835031282999       (deg K)
!	       U = -2.607041335277112E-016  (m/a)
!	       w =  3.064052270958646E-018  (m/a)
!	   Sigma =   1.21539202938520       (10^(-3) deg K/a)
!	 Sigma_c =  -1.32366417749281       (10^(-3) deg K/a)
!
! ELB 7/27/07
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use testsFG

   real(kind) :: t, r, z, H, TT, U, w, Sig, M, Sigc

   print *, 'Enter  t, r, z  (in a, km, m, resp.):'
   read *, t, r, z

   call testF(r*1000.0_kind,z,H,TT,U,w,Sig,M,Sigc)

   print *, 'Result from Test F (functions of r only):'
   print *, '      H =', H, ' (m)'
   print *, '      M =', M*SperA, ' (m/a)'
   print *, '(functions of r,z):'
   print *, '      T =', TT, ' (deg K)'
   print *, '      U =', U*SperA, ' (m/a)'
   print *, '      w =', w*SperA, ' (m/a)'
   print *, '  Sigma =', Sig*SperA*1.0e3_kind, ' (10^(-3) deg K/a)'
   print *, 'Sigma_c =', Sigc*SperA*1.0e3_kind, ' (10^(-3) deg K/a)'
   print *, ' '

   call testG(t*SperA,r*1000.0_kind,z,H,TT,U,w,Sig,M,Sigc)

   print *, 'Result from Test G (functions of t,r only):'
   print *, '      H =', H, ' (m)'
   print *, '      M =', M*SperA, ' (m/a)'
   print *, '(functions of t,r,z):'
   print *, '      T =', TT, ' (deg K)'
   print *, '      U =', U*SperA, ' (m/a)'
   print *, '      w =', w*SperA, ' (m/a)'
   print *, '  Sigma =', Sig*SperA*1.0e3_kind, ' (10^(-3) deg K/a)'
   print *, 'Sigma_c =', Sigc*SperA*1.0e3_kind, ' (10^(-3) deg K/a)'

end program simpleFG
