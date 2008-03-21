

/*! \fn PetscErrorCode IceModel::temperatureStep()
    \brief Takes a semi-implicit time-step for the temperature equation.

In summary, the conservation of energy equation is
    \f[ \rho c_p \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma,\f] 
where \f$T(t,x,y,z)\f$ is the temperature of the ice.  This equation is the shallow approximation
of the full three-dimensional conservation of energy.  Note \f$dT/dt\f$ stands for the material
derivative, so advection is included.  Here \f$\rho\f$ is the density of ice, 
\f$c_p\f$ is its specific heat, and \f$k\f$ is its conductivity.  Also \f$\Sigma\f$ is the volume
strain heating.

In summary, the numerical method is first-order upwind for advection and centered-differences with
semi-implicitness for the vertical conduction term.  We work from the bottom 
of the column upward in building the system to solve (in the semi-implicit time-stepping scheme).
The excess energy above pressure melting is converted to melt-water, and that a fraction 
of this melt water is transported to the base according to the scheme in excessToFromBasalMeltLayer().

The method uses equally-spaced calculation but the methods getValColumn(), setValColumn() interpolate 
back and forth from this equally-spaced calculational grid to the (usually) non-equally space storage 
grid.

In this procedure four scalar fields are modified: vHmelt, vbasalMeltRate, Tb3, and Tnew3.
But vHmelt, vbasalMeltRate and Tb3 will never need to communicate ghosted values (i.e. horizontal 
stencil neighbors.  The ghosted values for T3 are updated from the values in Tnew3 in the
communication done by temperatureAgeStep().
 
@cond REFMAN
Here is a more complete discussion and derivation.

Consider a column of a slowly flowing and heat conducting material as shown in the left side 
of the next figure.  (The left side shows a general column of flowing and heat conduction 
material showing a small segment \f$V\f$.  The right side shows a more specific column of ice 
flowing and sliding over bedrock.)  This is an \em Eulerian view so the material flows through 
a column which remains fixed (and is notional).  The column is vertical.  We will 
assume when needed that it is rectangular in cross-section with cross-sectional area 
\f$\Delta x\Delta y\f$.

\image latex earlycols.png "Left: a general column of material.  Right: ice over bedrock." width=3in

[FIXME: CONTINUE TO MINE eqns3D.tex FOR MORE]

The application of the geothermal flux at the base of a column is a special case for which 
we give a finite difference argument.  This scheme follows the equation (2.114) in 
Morton and Mayers (2005).  We have the boundary condition
	\f[  -k \frac{\partial T}{\partial z} = G(t,x,y) \f]
where \f$G(t,x,y)\f$ is the applied geothermal flux, and it is applied at level \f$z=-B_0\f$ 
in the bedrock (which is the only case considered here).  We <em> add a virtual lower grid 
point </em> \f$z_{-1} = z_0 - \Delta  z\f$ and we approximate the above boundary condition 
at $z_0$ by the centered-difference
	\f[  -k \frac{T_{1} - T_{-1}}{2 \Delta z} = G. \f]
Here \f$T_k = T_{ijk}^{l+1}\f$.  We also apply the discretized conduction equation at $z_0$.  These two combined equations 
yield a simplified form
	\f[(1 + 2 KR) T_0 - 2 K T_1 = T_0 + \frac{2\Delta t}{\rho c_p \Delta z} G \f]
where \f$K = k \Delta t (\rho c \Delta z^2)^{-1}\f$.

@endcond
 */
 
 

/*! \page bombproof About BOMBPROOF, a proposed numerical scheme for temperature and age
\latexonly

The goal is a completely bombproof finite difference numerical scheme for the 
advection-conduction-reaction problem for the temperature within the ice.  Note crucially 
that the region in which the temperature equation needs to be solved is changing in time.
``Bombproof'' means, roughly, as stable as possible in as many ways as possible.


\endlatexonly
 */ 


