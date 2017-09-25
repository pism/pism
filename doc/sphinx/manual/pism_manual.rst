.. include:: ../prologue.rst

.. default-role:: math

PISM User's Manual
==================

.. image:: /_static/pism-logo.png

.. figure:: figures/gris-flow-600m.png

   Magnitude of horizontal surface velocity from PISM run on a horizontal grid resolution
   of 600 m. Visualization by QGIS.

.. rubric:: Authorship

PISM is a joint project between developers in the ice sheet modeling group at the
University of Alaska (UAF), developers at the Potsdam Institute for Climate Impact
Research (PIK), and several additional developers listed here. Current (or recent)
affiliation and area of contributions shown:

.. list-table::
   :header-rows: 1

   * - Torsten Albrecht (PIK)
     - ice shelf physics and numerics 
   * - **Andy Aschwanden** (UAF)
     - scripts, visualization, thermodynamics, SeaRISE-Greenland
   * - Jed Brown (ANL)
     - source code original author, SSA numerics, PETSc underpinnings 
   * - **Ed Bueler** (UAF)
     - principal investigator, verification, earth deformation, SIA numerics,
       thermodynamics, documentation
   * - Dani DellaGiustina (UAF)
     - regional tools and modeling 
   * - Johannes Feldman (PIK)
     - marine ice sheet processes 
   * - Bob Fischer (GFDL)
     - coupling design 
   * - Marijke Habermann (UAF)
     - inversion
   * - Marianne Haseloff (UBC)
     - ice streams: physics and numerics
   * - Regine Hock (UAF)
     - surface mass and energy balance 
   * - **Constantine Khroulev** (UAF)
     - source code primary author, input/output, software design, climate couplers,
       parallelization, testing, user support, most documentation, most bug fixes,
       regional tools, ...
   * - Anders Levermann (PIK)
     - calving, ice shelf processes 
   * - Craig Lingle (UAF)
     - original SIA model, earth deformation 
   * - Maria Martin (PIK)
     - SeaRISE-Antarctica, Antarctica processes 
   * - Mattias Mengel (PIK)
     - marine ice sheet processes 
   * - David Maxwell (UAF)
     - inversion, SSA finite elements, python bindings 
   * - Ward van Pelt (IMAU)
     - hydrology analysis and design 
   * - Julien Seguinot (INK)
     - bug fixes, temperature index model 
   * - Ricarda Winkelmann (PIK)
     - Antarctica processes, coupling, and modeling  
   * - Florian Ziemen (UAF)
     - bug fixes, sliding 

Email the **highlighted** UAF developers at |pism-email|.

\textsl{Copyright (C) 2004--2014 The PISM Authors}
\medskip

\noindent \textsl{This file is part of PISM.  PISM is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.  PISM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  You should have received a copy of the GNU General Public License\index{GPL (\emph{GNU Public License})} along with PISM; see \emph{\texttt{COPYING}}.  If not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA}
\end{quote}
\vspace{0.5in}

\centerline{\textsc{Acknowledgements}}
\bigskip

\small
NASA Modeling, Analysis, and Prediction (MAP) program\index{Organizations!NASA!Modeling, Analysis, and Prediction program} grant \# NNX13AM16G and NASA Cryospheric Sciences program grant \# NNX13AK27G support the development of PISM from 2013 to 2017.  NASA MAP grant \# NNX09AJ38G supported the development of PISM from 2009 to 2013.  Development from 2002 to 2008 was supported by the NASA Cryospheric Sciences program\index{Organizations!NASA!Cryospheric Sciences program}.

The Snow, Ice, and Permafrost group\index{Organizations!Geophysical Institute!Snow, Ice, and Permafrost group} at the Geophysical Institute is the home for the University of Alaska PISM developers; find us in Elvey 410D.  The Arctic Region Supercomputing Center\index{Organizations!Arctic Region Supercomputing Center (ARSC)} has provided significant computational resources and technical help in the development of PISM.

Thanks for comments/questions from many PISM users around the world, including these not already listed as PISM authors:


    Tolly A{\dh}algeirsd{\'o}ttir, Antje Fitzner, Nick Golledge, Tore Hattermann, Moritz
    H\"utten, Thomas Kleiner, Leo van Kampenhout, Marianne Madsen, Malou Maris, Tim Morey,
    Mirena Olaizola, Christian Rodehacke, Nathan Shemonski, Sebastian Simonsen, Anne
    Solgaard, Ben Sperisen, Synne H\o{}yer Svendsen, Martin Truffer, Shuting Yang, Ryan
    Woodard

for helpful comments and questions on PISM and this *Manual*. Dave Covey, Don Bahls, and
Greg Newby have supported our hardware, software, and computations. Bob Bindschadler,
Sophie Nowicki, Jesse Johnson, and others in the SeaRISE group have motivated and assisted
PISM development in many ways.

Introduction
------------

Welcome!  All information about PISM is online at the home page

    |pism-docs|

Please see the extensive lists of PISM publications and applications at that page.

This User's Manual gives examples of how to run PISM using publicly-available data for:
the whole Greenland ice sheet, the Jakobshavn outlet glacier in Greenland, the Ross ice
shelf in Antarctica, and a number of simplified geometry tests. It documents all the PISM
options. It summarizes the continuum models used by PISM, and it illustrates how PISM's
numerical approximations are verified.

See the PISM Installation Manual\footnote{PDF for latest stable release at
\url{https://github.com/pism/pism/blob/dev/INSTALL.md}.} for how to
download\index{PISM!download source code} the PISM source code and install
it\index{PISM!install}, along with needed libraries. The PISM Climate Forcing
Manual\footnote{PDF for latest stable release at
\url{http://www.pism-docs.org/wiki/lib/exe/fetch.php?media=pism_forcing.pdf}.} extends the
User's Manual to cover additional couplings to atmosphere and ocean models and data.

Users who want to understand more deeply how PISM is designed, or who want to extend it,
will need to go beyond what is described here. See the \emph{Source Code
Browser}\index{PISM!\emph{Source Code Browser (HTML)}}, which is online for the latest
stable version.\footnote{At \url{http://www.pism-docs.org/wiki/doku.php?id=browser}.} It
can be generated from source code as described in the PISM Installation Manual. It gives a
complete view of the class/object structure of the PISM source code.


\begin{center}
\framebox{\parbox{5.0in}{ \emph{WARNING}:\index{PISM!warning}  PISM is an ongoing research project.  Ice sheet modeling requires many choices.  Please don't trust the results of PISM or any other ice sheet model without a fair amount of exploration.  Also, please don't expect all your questions to be answered here.  Write to us with questions at \PISMEMAIL.} }
\end{center}


\clearpage\newpage
\input{getting-started.tex}

\clearpage\newpage
\input{highlevelview.tex}

\clearpage\newpage
\input{init-boot.tex}

\clearpage\newpage
\input{modeling-computational.tex}

\clearpage\newpage
\input{modeling-dynamics.tex}

\clearpage\newpage
\input{modeling-subglacier.tex}

\clearpage\newpage
\input{modeling-marine.tex}

\clearpage\newpage
\input{practical-usage.tex}

\clearpage\newpage
\input{verification.tex}

\clearpage\newpage
\input{simplified-geometry.tex}

\clearpage\newpage
\input{validation.tex}

%\clearpage\newpage
%\input{storglaciaren.tex}

\clearpage\newpage
\input{jako.tex}

