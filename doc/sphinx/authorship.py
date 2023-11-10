#!/usr/bin/env python3

import yaml
import sys

print(""".. DO NOT EDIT. Update CITATION.cff to add more authors and see `authorship.py`.

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
     - Areas of contribution""")

for author in yaml.safe_load(sys.stdin)["authors"]:
    name = f"{author['given-names']} {author['family-names']}"
    try:
        name = f"`{name} <{author['orcid']}>`_"
    except KeyError:
        pass
    print(f"""
   * - {name} ({author['affiliation']})
     - {author['contribution-areas']}""")
