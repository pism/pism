.. DO NOT EDIT. Update CITATION.cff to add more authors and see `authorship.py`.

.. include:: global.txt

.. include:: <isonum.txt>

.. _sec-pism-authors:

Authorship
==========

   *Copyright* |copy| *2004 -- 2023 the PISM authors.*

   *This file is part of PISM. PISM is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 3 of the License, or (at your option) any later
   version. PISM is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
   PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
   have received a copy of the GNU General Public License along with PISM. If not, write
   to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301 USA*

PISM is a joint project between developers in the ice sheet modeling group at the
University of Alaska (UAF), developers at the Potsdam Institute for Climate Impact
Research (PIK), and several additional developers listed here.

.. list-table::
   :header-rows: 1
   :widths: 1,1

   * - Name and affiliation
     - Areas of contribution

   * - `Constantine Khrulev <https://orcid.org/0000-0003-2170-7454>`_ (University of Alaska Fairbanks)
     - primary source code author, software design, numerics, input/output, higher-order stress balance, mass conservation, verification, documentation, testing, user support, maintenance, ...

   * - `Andy Aschwanden <https://orcid.org/0000-0001-8149-2315>`_ (University of Alaska Fairbanks)
     - project leader, thermodynamics, testing, user support, visualization, etc

   * - `Ed Bueler <https://orcid.org/0000-0001-8232-4145>`_ (University of Alaska Fairbanks)
     - former project leader, verification, earth deformation, numerics, thermodynamics, documentation

   * - `Jed Brown <https://orcid.org/0000-0002-9945-0639>`_ (University of Colorado Boulder)
     - source code original author, SSA numerics, parallelism using PETSc, etc

   * - David Maxwell (University of Alaska Fairbanks)
     - inverse methods, SSA finite element solver, Python bindings

   * - `Torsten Albrecht <https://orcid.org/0000-0001-7459-2860>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - ice shelf physics and numerics

   * - `Ronja Reese <https://orcid.org/0000-0001-7625-040X>`_ (University of Northumbria at Newcastle)
     - sub-shelf mass balance (PICO), ice shelf numerics, etc

   * - `Matthias Mengel <https://orcid.org/0000-0001-6724-9685>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - ice shelf physics (PICO)

   * - `Maria Martin <https://orcid.org/0000-0002-1443-0891>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - SeaRISE-Antarctica, Antarctic mean annual temperature parameterization, sub-shelf processes

   * - `Ricarda Winkelmann <https://orcid.org/0000-0003-1248-3217>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - sub-shelf mass balance (PICO), Antarctica processes, coupling, modeling, etc

   * - `Maria Zeitz <https://orcid.org/0000-0001-8129-1329>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - surface processes (dEBM-simple)

   * - `Anders Levermann <https://orcid.org/0000-0003-4432-4704>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - ice shelf dynamics theory, calving, etc

   * - `Johannes Feldmann <https://orcid.org/0000-0003-4210-0221>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - marine ice sheet processes, grounding line dynamics

   * - `Julius Garbe <https://orcid.org/0000-0003-3140-3307>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - surface mass and energy balance (dEBM-simple)

   * - `Marianne Haseloff <https://orcid.org/0000-0002-6576-5384>`_ (University of Wisconsin-Madison)
     - PISM-PIK development, SSA and SIA theory and numerics

   * - `Julien Seguinot <https://orcid.org/0000-0002-5315-0761>`_ (University of Bergen)
     - improvements to the temperature index (PDD) model

   * - `Sebastian Hinck <https://orcid.org/0000-0002-6428-305X>`_ (Deutsches Klimarechenzentrum GmbH)
     - improvements to the Lingle-Clark bed deformation model

   * - `Thomas Kleiner <https://orcid.org/0000-0001-7825-5765>`_ (Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research)
     - bug fixes in the enthalpy model

   * - `Elizabeth Fischer <https://orcid.org/0000-0002-3987-2696>`_ (University of Alaska Fairbanks)
     - bug fixes, coupling to a GCM

   * - `Anders Damsgaard <https://orcid.org/0000-0002-9284-1709>`_ (Aarhus University)
     - automatic testing

   * - Craig Lingle (University of Alaska Fairbanks)
     - original SIA model, earth deformation

   * - `Ward van Pelt <https://orcid.org/0000-0003-4839-7900>`_ (Uppsala Universitet)
     - subglacial hydrology analysis and design

   * - `Florian Ziemen <https://orcid.org/0000-0001-7095-5740>`_ (Deutsches Klimarechenzentrum GmbH)
     - bug fix in the interpolation code

   * - Nathan Shemonski (Zoox, Inc)
     - documentation, scripts and code related to the pre-v0.1 EISMINT Greenland setup

   * - `Kenneth Mankoff <https://orcid.org/0000-0001-5453-2019>`_ (Goddard Institute for Space Studies)
     - documentation

   * - `Joseph Kennedy <https://orcid.org/0000-0002-9348-693X>`_ (University of Alaska Fairbanks)
     - documentation

   * - Kyle Blum (University of Alaska Fairbanks)
     - bug fix in the SeaRISE-Antarctica setup

   * - Marijke Habermann (University of Alaska Fairbanks)
     - inverse modeling

   * - `Daniella DellaGiustina <https://orcid.org/0000-0002-5643-1956>`_ (University of Arizona)
     - regional modeling

   * - `Regine Hock <https://orcid.org/0000-0001-8336-9441>`_ (University of Oslo)
     - surface mass and energy balance

   * - `Moritz Kreuzer <https://orcid.org/0000-0002-8622-6638>`_ (Potsdam Institute for Climate Impact Research (PIK))
     - Python 3.8+ compatibility

   * - Enrico Degregori (Deutsches Klimarechenzentrum GmbH)
     - PICO optimization

   * - Simon Schoell (Potsdam Institute for Climate Impact Research (PIK))
     - calving bug fix
