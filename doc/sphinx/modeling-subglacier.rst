.. _sec-modeling-subglacier:

Modeling choices:  The subglacier
=================================

.. contents::

.. _sec-basestrength:

Controlling basal strength
--------------------------

When using option :opt:`-stress_balance ssa+sia`, the SIA+SSA hybrid stress balance, a model for basal resistance is required.  This model for basal resistance is based, at least conceptually, on the hypothesis that the ice sheet is underlain by a layer of till [Clarke05]_.  The user can control the parts of this model:

- the so-called sliding law, typically a power law, which relates the ice base (sliding) velocity to the basal shear stress, and which has a coefficient which is or has the units of a yield stress,
- the model relating the effective pressure on the till layer to the yield stress of that layer, and
- the model for relating the amount of water stored in the till to the effective pressure on the till.

This subsection explains the relevant options.

The primary example of ``-stress_balance ssa+sia`` usage is in section :ref:`sec-start` of this Manual, but the option is also used in sections :ref:`subsect-MISMIP`, :ref:`subsect-MISMIP3d`, and :ref:`sec-jako`.

In PISM the key coefficient in the sliding is always denoted as yield stress :math:`\tau_c`, which is :var:`tauc` in PISM output files.  This parameter represents the strength of the aggregate material at the base of an ice sheet, a poorly-observed mixture of pressurized liquid water, ice, granular till, and bedrock bumps.  The yield stress concept also extends to the power law form, and thus most standard sliding laws can be chosen by user options (below).  One reason that the yield stress is a useful parameter is that it can be compared, when looking at PISM output files, to the driving stress (:var:`taud_mag` in PISM output files).  Specifically, where :var:`tauc` :math:`<` :var:`taud_mag` you are likely to see sliding if option ``-stress_balance ssa+sia`` is used.

A historical note on modeling basal sliding is in order.  Sliding can be added directly to a SIA stress balance model by making the sliding velocity a local function of the basal value of the driving stress.  Such an SIA sliding mechanism appears in ISMIP-HEINO [Calovetal2009HEINOfinal]_ and in EISMINT II experiment H [EISMINT00]_, among other places.  This kind of sliding is *not* recommended, as it does not make sense to regard the driving stress as the local generator of flow if the bed is not holding all of that stress [BBssasliding]_, [Fowler01]_.  Within PISM, for historical reasons, there is an implementation of SIA-based sliding only for verification test E; see section :ref:`sec-verif`.  PISM does *not* support this SIA-based sliding mode in other contexts.

Choosing the sliding law
^^^^^^^^^^^^^^^^^^^^^^^^

In PISM the sliding law can be chosen to be a purely-plastic (Coulomb) model, namely,

.. math::
   :name: eq-plastic

   |\boldsymbol{\tau}_b| \le \tau_c \quad \text{and} \quad \boldsymbol{\tau}_b = - \tau_c \frac{\mathbf{u}}{|\mathbf{u}|} \quad\text{if and only if}\quad |\mathbf{u}| > 0.

Equation :eq:`eq-plastic` says that the (vector) basal shear stress :math:`\boldsymbol{\tau}_b` is at most the yield stress :math:`\tau_c`, and that only when the shear stress reaches the yield value can there be sliding.  The sliding law can, however, also be chosen to be the power law

.. math::
   :name: eq-pseudoplastic

   \boldsymbol{\tau}_b =  - \tau_c \frac{\mathbf{u}}{u_{\text{threshold}}^q |\mathbf{u}|^{1-q}},

where :math:`u_{\text{threshold}}` is a parameter with units of velocity (see below).  Condition :eq:`eq-plastic` is studied in [SchoofStream]_ and [SchoofTill]_ in particular, while power laws for sliding are common across the glaciological literature (e.g.~see [CuffeyPaterson]_, [GreveBlatter2009]_).  Notice that the coefficient :math:`\tau_c` in :eq:`eq-pseudoplastic` has units of stress, regardless of the power :math:`q`.

In both of the above equations :eq:`eq-plastic` and :eq:`eq-pseudoplastic` we call :math:`\tau_c` the *yield stress*.  It corresponds to the variable ``tauc`` in PISM output files.  We call the power law :eq:`eq-pseudoplastic` a "pseudo-plastic" law with power :math:`q` and threshold velocity :math:`u_{\text{threshold}}`.  At the threshold velocity the basal shear stress :math:`\boldsymbol{\tau}_b` has exact magnitude :math:`\tau_c`.  In equation :eq:`eq-pseudoplastic`, :math:`q` is the power controlled by ``-pseudo_plastic_q``, and the threshold velocity :math:`u_{\text{threshold}}` is controlled by ``-pseudo_plastic_uthreshold``.  The plastic model :eq:`eq-plastic` is the :math:`q=0` case of :eq:`eq-pseudoplastic`.

See Table tab-sliding-power-law_ for options controlling the choice of sliding law. The purely plastic case is the default; just use ``-stress_balance ssa+sia`` to turn it on.   (Or use ``-stress_balance ssa`` if a model with no vertical shear is desired.)

.. warning::

   Options ``-pseudo_plastic_q`` and ``-pseudo_plastic_uthreshold`` have no effect if ``-pseudo_plastic`` is not set.

.. list-table:: Sliding law command-line options
   :name: tab-sliding-power-law
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-pseudo_plastic`
     - Enables the pseudo-plastic power law model. If this is not set the sliding law is purely-plastic, so ``pseudo_plastic_q`` and ``pseudo_plastic_uthreshold`` are inactive.
   * - :opt:`-plastic_reg` (m/a)
     - Set the value of :math:`\eps` regularization of the plastic law, in the formula :math:`\boldsymbol{\tau}_b = - \tau_c \mathbf{u}/\sqrt{|\mathbf{u}|^2 + \eps^2}`. The default is :math:`0.01` m/a. This parameter is inactive if ``-pseudo_plastic`` is set.
   * - :opt:`-pseudo_plastic_q`
     - Set the exponent :math:`q` in :eq:`eq-pseudoplastic`.  The default is :math:`0.25`.
   * - :opt:`-pseudo_plastic_uthreshold` (m/a)
     - Set :math:`u_{\text{threshold}}` in :eq:`eq-pseudoplastic`.  The default is :math:`100` m/a.

Equation :eq:`eq-pseudoplastic` is a very flexible power law form.  For example, the linear case is :math:`q=1`, in which case if :math:`\beta=\tau_c/u_{\text{threshold}}` then the law is of the form

.. math::

   \boldsymbol{\tau}_b = - \beta \mathbf{u}

(The ":math:`\beta`" coefficient is also called :math:`\beta^2` in some sources (see [MacAyeal]_, for example).)  If you want such a linear sliding law, and you have a value :math:`\beta=```beta`` in :math:`\text{Pa}\,\text{s}\,\text{m}^{-1}`, then use this option combination:

.. code-block:: none

   -pseudo_plastic \
   -pseudo_plastic_q 1.0 \
   -pseudo_plastic_uthreshold 3.1556926e7 \
   -yield_stress constant -tauc beta

This sets :math:`u_{\text{threshold}}` to 1 :math:`\text{m}\,\text{s}^{-1}` but using units :math:`\text{m}\,\text{a}^{-1}`.

More generally, it is common in the literature to see power-law sliding relations in the form

.. math::

   \boldsymbol{\tau}_b = - C |\mathbf{u}|^{m-1} \mathbf{u},

where :math:`C` is a constant, as for example in sections :ref:`subsect-MISMIP` and :ref:`subsect-MISMIP3d`.  In that case, use this option combination:

.. code-block:: none

   -pseudo_plastic \
   -pseudo_plastic_q m \
   -pseudo_plastic_uthreshold 3.1556926e7 \
   -yield_stress constant \
   -tauc C

Determining the yield stress
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Other than setting it to a constant, which only applies in some special cases, the above discussion does not determine the yield stress :math:`\tau_c`.  As shown in Table tab-yieldstress_, there are two schemes for determining :math:`\tau_c` in a spatially-variable manner:

- ``-yield_stress mohr_coulomb`` (the default) determines the yields stress by models of till material property (the till friction angle) and of the effective pressure on the saturated till, or
- ``-yield_stress constant`` allows the yield stress to be supplied as time-independent data, read from the input file.


In normal modelling cases, variations in yield stress are part of the explanation of the locations of ice streams [SchoofStream]_.  The default model ``-yield_stress mohr_coulomb`` determines these variations in time and space.  The value of :math:`\tau_c` is determined in part by a subglacial hydrology model, including the modeled till-pore water amount ``tillwat`` (subsection :ref:`subsect-subhydro`), which then determines the effective pressure :math:`N_{til}` (see below).  The value of :math:`\tau_c` is also determined in part by a material property field :math:`\phi=```tillphi``, the "till friction angle".  These quantities are related by the Mohr-Coulomb criterion [CuffeyPaterson]_:

.. math::
   :name: eq-mohrcoulomb

   \tau_c = c_{0} + (\tan\phi)\,N_{til}.

Here :math:`c_0` is called the "till cohesion", whose default value in PISM is zero (see [SchoofStream]_, formula (2.4)) but which can be set by option :opt:`-till_cohesion`.

Option combination ``-yield_stress constant -tauc X`` can be used to fix the yield stress to have value :math:`\tau_c=```X`` at all grounded locations and all times if desired.  This is unlikely to be a good modelling choice for real ice sheets.


.. list-table:: Command-line options controlling how yield stress is determined
   :name: tab-yieldstress
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-yield_stress mohr_coulomb`
     - The default.  Use equation :eq:`eq-mohrcoulomb` to determine :math:`\tau_c`.  Only effective if ``-stress_balance ssa`` or ``-stress_balance ssa+sia`` is also set.
   * - :opt:`-till_cohesion`
     - Set the value of the till cohesion (:math:`c_{0}`) in the plastic till model.  The value is a pressure, given in Pa.
   * - :opt:`-tauc_slippery_grounding_lines`
     - If set, reduces the basal yield stress at grounded-below-sea-level grid points one cell away from floating ice or ocean.  Specifically, it replaces the normally-computed :math:`\tau_c` from the Mohr-Coulomb relation, which uses the effective pressure from the modeled amount of water in the till, by the minimum value of :math:`\tau_c` from Mohr-Coulomb, i.e.~using the effective pressure corresponding to the maximum amount of till-stored water.  Does not alter the reported amount of till water, nor does this mechanism affect water conservation. 
   * - :opt:`-plastic_phi` (degrees)
     - Use a constant till friction angle. The default is :math:`30^{\circ}`.
   * - :opt:`-topg_to_phi` (*list of 4 numbers*)
     - Compute :math:`\phi` using equation :eq:`eq-phipiecewise`.
   * - :opt:`-yield_stress constant`
     - Keep the current values of the till yield stress :math:`\tau_c`.  That is, do not update them by the default model using the stored basal melt water.  Only effective if ``-stress_balance ssa`` or ``-stress_balance ssa+sia`` is also set. 
   * - :opt:`-tauc`
     - Directly set the till yield stress :math:`\tau_c`, in units Pa, at all grounded locations and all times.  Only effective if used with ``-yield_stress constant``, because otherwise :math:`\tau_c` is updated dynamically. 

We find that an effective, though heuristic, way to determine :math:`\phi=```tillphi`` in :eq:`eq-mohrcoulomb` is to make it a function of bed elevation [AschwandenAdalgeirsdottirKhroulev]_, [Martinetal2011]_, [Winkelmannetal2011]_.  This heuristic is motivated by hypothesis that basal material with a marine history should be weak [HuybrechtsdeWolde]_.  PISM has a mechanism setting :math:`\phi`=``tillphi`` to be a *piecewise-linear* function of bed elevation.  The option is

.. code-block:: none

   -topg_to_phi phimin,phimax,bmin,bmax

.. math::

   \newcommand{\phimin}{\phi_{\mathrm{min}}}
   \newcommand{\phimax}{\phi_{\mathrm{max}}}
   \newcommand{\bmin}{b_{\mathrm{min}}}
   \newcommand{\bmax}{b_{\mathrm{max}}}

Thus the user supplies 4 parameters: :math:`\phimin`, :math:`\phimax`, :math:`\bmin`, :math:`\bmax`, where :math:`b` stands for the bed elevation.  To explain these, we define :math:`M = (\phimax - \phimin) / (\bmax - \bmin)`.  Then

.. math::
   :name: eq-phipiecewise

   \phi(x,y) =
   \begin{cases}
     \phimin, & b(x,y) \le \bmin, \\
     \phimin + (b(x,y) - \bmin) \,M, & \bmin < b(x,y) < \bmax, \\
     \phimax, & \bmax \le b(x,y).
   \end{cases}

It is worth noting that an earth deformation model (see section :ref:`subsect-beddef`) changes :math:`b(x,y)=\mathrm{topg}` used in :eq:`eq-phipiecewise`, so that a sequence of runs such as

.. code-block:: none

   pismr -i foo.nc -bed_def lc -stress_balance ssa+sia -topg_to_phi 10,30,-50,0 ... -o bar.nc
   pismr -i bar.nc -bed_def lc -stress_balance ssa+sia -topg_to_phi 10,30,-50,0 ... -o baz.nc

will use *different* ``tillphi`` fields in the first and second runs.  PISM will print a warning during initialization of the second run:

.. code-block:: none

   * Initializing the default basal yield stress model...
     option -topg_to_phi seen; creating tillphi map from bed elev ...
   PISM WARNING: -topg_to_phi computation will override the 'tillphi' field
                 present in the input file 'bar.nc'!

Omitting the ``-topg_to_phi`` option in the second run would make PISM continue with the same ``tillphi`` field which was set in the first run.

Determining the effective pressure
----------------------------------

When using the default option ``-yield_stress mohr_coulomb``, the effective pressure on the till :math:`N_{til}` is determined by the modeled amount of water in the till.  Lower effective pressure means that more of the weight of the ice is carried by the pressurized water in the till and thus the ice can slide more easily.  That is, equation :eq:`eq-mohrcoulomb` sets the value of :math:`\tau_c` proportionately to :math:`N_{til}`.  The amount of water in the till is, however, a nontrivial output of the hydrology (subsection :ref:`subsect-subhydro`) and conservation-of-energy (section :ref:`subsect-energy`) submodels in PISM.

Following [Tulaczyketal2000]_, based on laboratory experiments with till extracted from an ice stream in Antarctica, [BuelervanPelt2015]_ propose the following parameterization which is used in PISM.  It is based on the ratio :math:`s=W_{til}/W_{til}^{max}` where :math:`W_{til}=` ``tillwat`` is the effective thickness of water in the till and :math:`W_{til}^{max}=` ``hydrology_tillwat_max`` is the maximum amount of water in the till (see subsection :ref:`subsect-subhydro`):

.. math::
   :name: eq-computeNtil

   N_{til} = \min\left\{P_o, N_0 \left(\frac{\delta P_o}{N_0}\right)^s \, 10^{(e_0/C_c) \left(1 - s\right).}\right\}

Here :math:`P_o` is the ice overburden pressure, which is determined entirely by the ice thickness and density, and the remaining parameters are set by options in Table tab-effective-pressure_.  While there is experimental support for the default values of :math:`C_c`, :math:`e_0`, and :math:`N_0`, the value of :math:`\delta=```till_effective_fraction_overburden`` should be regarded as uncertain, important, and subject to parameter studies to assess its effect.

FIXME: EVOLVING CODE:  If the ``tauc_add_transportable_water`` configuration flag is set (either in the configuration file or using the :opt:`-tauc_add_transportable_water` option), then the above formula becomes FIXME

.. list-table:: Command-line options controlling how till effective pressure :math:`N_{til}` in equation :eq:`eq-mohrcoulomb` is determined
   :name: tab-effective-pressure
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-till_reference_void_ratio`
     - :math:`= e_0` in :eq:`eq-computeNtil`, dimensionless, with default value 0.69 [Tulaczyketal2000]_
   * - :opt:`-till_compressibility_coefficient`
     - :math:`= C_c` in :eq:`eq-computeNtil`, dimensionless, with default value 0.12 [Tulaczyketal2000]_
   * - :opt:`-till_effective_fraction_overburden`
     - :math:`= \delta` in :eq:`eq-computeNtil`, dimensionless, with default value 0.02 [BuelervanPelt2015]_
   * - :opt:`-till_reference_effective_pressure`
     - :math:`= N_0` in :eq:`eq-computeNtil`, in Pa, with default value 1000.0 [Tulaczyketal2000]_

.. _sec-subhydro:

Subglacial hydrology
--------------------

At the present time, two simple subglacial hydrology models are implemented *and documented* in PISM, namely ``-hydrology null`` and ``-hydrology routing``; see Table tab-hydrologychoice_ and [BuelervanPelt2015]_.  In both models, some of the water in the subglacial layer is stored locally in a layer of subglacial till by the hydrology model.  In the  ``routing`` model water is conserved by horizontally-transporting the excess water (namely ``bwat``) according to the gradient of the modeled hydraulic potential.  In both hydrology models a state variable ``tillwat`` is the effective thickness of the layer of liquid water in the till; it is used to compute the effective pressure on the till (see the previous subsection).  The pressure of the transportable water ``bwat`` in the ``routing`` model does not relate directly to the effective pressure on the till.

.. list-table:: Command-line options to choose the hydrology model
   :name: tab-hydrologychoice
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`-hydrology null`
     - The default model with only a layer of water stored in till.  Not mass conserving in the map-plane but much faster than ``-hydrology routing``.  Based on "undrained plastic bed" model of [Tulaczyketal2000b]_.  The only state variable is ``tillwat``.
   * - :opt:`-hydrology routing`
     - A mass-conserving horizontal transport model in which the pressure of transportable water is equal to overburden pressure.  The till layer remains in the model, so this is a "drained and conserved plastic bed" model.  The state variables are ``bwat`` and ``tillwat``.

See Table tab-hydrology_ for options which apply to all hydrology models.  Note that the primary water source for these models is the energy conservation model which computes the basal melt rate ``basal_melt_rate_grounded``.  There is, however, also option :opt:`-hydrology_input_to_bed_file` which allows the user to *add* water directly into the subglacial layer, in addition to the computed ``basal_melt_rate_grounded`` values.  Thus ``-hydrology_input_to_bed_file`` allows the user to model drainage directly to the bed from surface runoff, for example.  Also option :opt:`-hydrology_bmelt_file` allows the user to replace the computed ``basal_melt_rate_grounded`` rate by values read from a file, thereby effectively decoupling the hydrology model from the ice dynamics (esp.~conservation of energy).

.. list-table:: Subglacial hydrology command-line options which apply to all hydrology models
   :name: tab-hydrology
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`hydrology_bmelt_file`
     - Specifies a NetCDF file which contains a time-independent field ``basal_melt_rate_grounded`` which has units of water thickness per time.  This rate *replaces* the conservation-of-energy computed rate ``basal_melt_rate_grounded``.
   * - :opt:`hydrology_const_bmelt` (m/s)
     - If ``-hydrology_use_const_bmelt`` is set then use this to set the constant rate (water thickness per time).
   * - :opt:`hydrology_input_to_bed_file`
     - Specifies a NetCDF file which contains a time-dependent field ``inputtobed`` which has units of water thickness per time.  This rate is *added to* the ``basal_melt_rate_grounded`` rate.
   * - :opt:`hydrology_input_to_bed_period` (a)
     - The period, in years, of ``-hydrology_input_to_bed_file`` data.
   * - :opt:`hydrology_input_to_bed_reference_year` (a)
     - The reference year for periodizing the ``-hydrology_input_to_bed_file`` data.
   * - :opt:`hydrology_tillwat_max` (m)
     - Maximum effective thickness for water stored in till.
   * - :opt:`hydrology_tillwat_decay_rate` (m/a)
     - Water accumulates in the till at the basal melt rate ``basal_melt_rate_grounded``, minus this rate.
   * - :opt:`-hydrology_use_const_bmelt`
     - Replace the conservation-of-energy basal melt rate ``basal_melt_rate_grounded`` with a constant.

The default model: ``-hydrology null``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this model the water is *not* conserved but it is stored locally in the till up to a specified amount; option :opt:`-hydrology_tillwat_max` sets that amount.  The water is not conserved in the sense that water above the ``hydrology_tillwat_max`` level is lost permanently.  This model is based on the "undrained plastic bed" concept of [Tulaczyketal2000b]_; see also [BBssasliding]_.

In particular, denoting ``tillwat`` by :math:`W_{til}`, the till-stored water layer effective thickness evolves by the simple equation

.. math::
   :name: eq-tillwatevolve

   \frac{\partial W_{til}}{\partial t} = \frac{m}{\rho_w} - C

where :math:`m=` :var:`basal_melt_rate_grounded}` (kg :math:`\text{m}^{-2}\,\text{s}^{-1}`), :math:`\rho_w` is the density of fresh water, and :math:`C` :var:`hydrology_tillwat_decay_rate`.  At all times bounds :math:`0 \le W_{til} \le W_{til}^{max}` are satisfied.

This ``-hydrology null`` model has been extensively tested in combination with the Mohr-Coulomb till (subsection :ref:`subsect-basestrength` above) for modelling ice streaming (see [AschwandenAdalgeirsdottirKhroulev]_ and [BBssasliding]_, among others).

The mass-conserving model: ``-hydrology routing``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this model the water *is* conserved in the map-plane.  Water does get put into the till, with the same maximum value ``hydrology_tillwat_max``, but excess water is horizontally-transported.  An additional state variable ``bwat``, the effective thickness of the layer of transportable water, is used by ``routing``.  This transportable water will flow in the direction of the negative of the gradient of the modeled hydraulic potential.  In the ``routing`` model this potential is calculated by assuming that the transportable subglacial water is at the overburden pressure [Siegertetal2009]_.  Ultimately the transportable water will reach the ice sheet grounding line or ice-free-land margin, at which point it will be lost.  The amount that is lost this way is reported to the user.

In this model ``tillwat`` also evolves by equation :eq:`eq-tillwatevolve`, but several additional parameters are used in determining how the transportable water ``bwat`` flows in the model; see Table tab-hydrologyrouting_.  Specifically, the horizontal subglacial water flux is determined by a generalized Darcy flux relation [Clarke05]_, [Schoofetal2012]_

.. math::
   :name: eq-flux

   \bq = - k\, W^\alpha\, |\grad \psi|^{\beta-2} \grad \psi

where :math:`\bq` is the lateral water flux, :math:`W=` ``bwat`` is the effective thickness of the layer of transportable water, :math:`\psi` is the hydraulic potential, and :math:`k`, :math:`\alpha`, :math:`\beta` are controllable parameters (Table tab-hydrologyrouting_).

In the ``routing`` model the hydraulic potential :math:`\psi` is determined by

.. math::
   :name: eq-hydraulicpotential

   \psi = P_o + \rho_w g (b + W)

where :math:`P_o=\rho_i g H` is the ice overburden pressure, :math:`g` is gravity, :math:`\rho_i` is ice density, :math:`\rho_w` is fresh water density, :math:`H` is ice thickness, and :math:`b` is the bedrock elevation.

For most choices of the relevant parameters and most grid spacings, the ``routing`` model is at least two orders of magnitude more expensive computationally than the ``null`` model.  This follows directly from the CFL-type time-step restriction on lateral flow of a fluid with velocity on the order of centimeters to meters per second (i.e.~the subglacial liquid water ``bwat``).  (By comparison, much of PISM ice dynamics time-stepping is controlled by the much slower velocity of the flowing ice.)  Therefore the user should start with short runs of order a few model years.  The option :opt:`-report_mass_accounting` is also recommended, so as to see the time-stepping behavior at ``stdout``.  Finally, ``daily`` or even ``hourly`` reporting for scalar and spatially-distributed time-series to see hydrology model behavior, especially on fine grids (e.g.~:math:`< 1` km).

.. list-table:: Command-line options specific to hydrology model ``routing``
   :name: tab-hydrologyrouting
   :header-rows: 1

   * - Option
     - Description
   * - :opt:`hydrology_hydraulic_conductivity` :math:`k`
     - :math:`=k` in formula :eq:`eq-flux`.
   * - :opt:`hydrology_null_strip` (km)
     - In the boundary strip water is removed and this is reported.  This option specifies the width of this strip, which should typically be one or two grid cells.
   * - :opt:`hydrology_gradient_power_in_flux` :math:`\beta`
     - :math:`=\beta` in formula :eq:`eq-flux`.
   * - :opt:`hydrology_thickness_power_in_flux` :math:`\alpha`
     - :math:`=\alpha` in formula :eq:`eq-flux`.
   * - :opt:`-report_mass_accounting`
     - At each major (ice dynamics) time-step, the duration of hydrology time steps is reported, along with the amount of subglacial water lost to ice-free land, to the ocean, and into the "null strip".

.. FIXME -hydrology distributed is not documented except by [BuelervanPelt2015]_

.. _sec-beddef:

Earth deformation models
------------------------

The option :opt:`bed_def` ``[iso, lc]`` turns one of the two available bed deformation models.

The first model ``-bed_def iso``, is instantaneous pointwise isostasy.  This model assumes that the bed at the starting time is in equilibrium with the load.  Then, as the ice geometry evolves, the bed elevation is equal to the starting bed elevation minus a multiple of the increase in ice thickness from the starting time: :math:`b(t,x,y) = b(0,x,y) - f [H(t,x,y) - H(0,x,y)]`.  Here :math:`f` is the density of ice divided by the density of the mantle, so its value is determined by setting the values of ``bed_deformation.mantle_density`` and ``constants.ice.density`` in the configuration file; see subsection :ref:`sec-pism-defaults`.  For an example and verification, see Test H in Verification section.

The second model ``-bed_def lc`` is much more physical.  It is based on papers by Lingle and Clark [LingleClark]_ and Bueler and others [BLKfastearth]_.  It generalizes and improves the most widely-used earth deformation model in ice sheet modeling, the flat earth Elastic Lithosphere Relaxing Asthenosphere (ELRA) model [Greve2001]_.  It imposes  essentially no computational burden because the Fast Fourier Transform is used to solve the linear differential equation [BLKfastearth]_.  When using this model in PISM, the rate of bed movement (uplift) and the viscous plate displacement are stored in the PISM output file and then used to initialize the next part of the run.  In fact, if gridded "observed" uplift data is available, for instance from a combination of actual point observations and/or paleo ice load modeling, and if that uplift field is put in a NetCDF variable with standard name ``tendency_of_bedrock_altitude`` in the input file, then this model will initialize so that it starts with the given uplift rate.

Here are minimal example runs to compare these models:

.. code-block:: none

   mpiexec -n 4 pisms -eisII A -y 8000 -o eisIIA_nobd.nc
   mpiexec -n 4 pisms -eisII A -bed_def iso -y 8000 -o eisIIA_bdiso.nc
   mpiexec -n 4 pisms -eisII A -bed_def lc -y 8000 -o eisIIA_bdlc.nc

Compare the :var:`topg`, :var:`usurf`, and :var:`dbdt` variables in the resulting output files. See also the comparison done in [BLKfastearth]_.

To include "measured" uplift rates during initialization, use the option :opt:`-uplift_file` to specify the name of the file containing the field :var:`dbdt` (CF standard name: ``tendency_of_bedrock_altitude``).

Use the :opt:`-topg_delta_file` option to apply a correction to the bed topography field read in from an input file. This sets the bed topography :math:`b` at the beginning of a run as follows:

.. math::
   :name: eq-bedcorrection

   b = b_{0} + \Delta b.

Here :math:`b_{0}` is the bed topography (:var:`topg`) read in from an input file and :math:`\Delta b` is the :var:`topg_delta` field read in from the file specified using this option.

A correction like this can be used to get a bed topography field at the end of a paleo-climate run that is closer to observed present day topography. The correction is computed by performing a "preliminary" run and subtracting modeled bed topography from present day observations. A subsequent run with this correction should produce a bed elevations that are closer to observed values.

.. _sec-bedsmooth:

Parameterization of bed roughness in the SIA
--------------------------------------------

Schoof [Schoofbasaltopg2003]_ describes how to alter the SIA stress balance to model ice flow over bumpy bedrock topgraphy. One computes the amount by which bumpy topography lowers the SIA diffusivity. An internal quantity used in this method is a smoothed version of the bedrock topography. As a practical matter for PISM, this theory improves the SIA's ability to handle bed roughness because it parameterizes the effects of "higher-order" stresses which act on the ice as it flows over bumps. For additional technical description of PISM's implementation, see the `Browser <pism-browser_>`_ page *Using Schoof's (2003) parameterized bed roughness technique in PISM*.

This parameterization is "on" by default when using ``pismr``. There is only one associated option: :opt:`-bed_smoother_range` gives the half-width of the square smoothing domain in meters. If zero is given, ``-bed_smoother_range 0`` then the mechanism is turned off. The mechanism is on by default using executable ``pismr``, with the half-width set to 5 km (``-bed_smoother_range 5.0e3``), giving Schoof's recommended smoothing size of 10 km [Schoofbasaltopg2003]_.

This mechanism is turned off by default in executables ``pisms`` and ``pismv``.

Under the default setting ``-o_size medium``, PISM writes fields :var:`topgsmooth` and :var:`schoofs_theta` from this mechanism. The thickness relative to the smoothed bedrock elevation, namely :var:`topgsmooth`, is the difference between the unsmoothed surface elevation and the smoothed bedrock elevation. It is *only used internally by this mechanism*, to compute a modified value of the diffusivity; the rest of PISM does not use this or any other smoothed bed. The field :var:`schoofs_theta` is a number :math:`\theta` between :math:`0` and :math:`1`, with values significantly below zero indicating a reduction in diffusivity, essentially a drag coefficient, from bumpy bed.
