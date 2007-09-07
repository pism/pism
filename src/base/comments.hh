
/*! \fn PetscErrorCode IceModel::temperatureStep()
    \brief Takes a semi-implicit time-step for the temperature equation.
    
@cond CONTINUUM
In summary, the conservation of energy equation is
    \f[ \rho c_p \frac{dT}{dt} = k \frac{\partial^2 T}{\partial z^2} + \Sigma,\f] 
where \f$T(t,x,y,z)\f$ is the temperature of the ice.  This equation is the shallow approximation
of the full three-dimensional conservation of energy.  Note \f$dT/dt\f$ stands for the material
derivative, so advection is included.  Here \f$\rho\f$ is the density of ice, 
\f$c_p\f$ is its specific heat, and \f$k\f$ is its conductivity.  Also \f$\Sigma\f$ is the volume
strain heating.
@endcond

@cond NUMERIC
In summary, the numerical method is first-order upwind for advection and centered-differences with
semi-implicitness for the vertical conduction term.  Note that we work from the bottom 
of the column upward in building the system to solve (in the semi-implicit time-stepping scheme).
Note that the excess energy above pressure melting is converted to melt-water, and that a fraction 
of this melt water is transported to the base according to the scheme in excessToFromBasalMeltLayer().
@endcond

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
@endcond
 */



