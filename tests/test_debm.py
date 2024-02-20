# Copyright (C) 2024 Andy Aschwanden, Constantine Khroulev
#
# This file is part of pism.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

import PISM
context = PISM.Context().ctx

def test_CalovGreveIntegrand():
    """
    Test the CalovGreveIntegrand
    """

    sigma = np.array([2.0, 0.0, 1.0])
    temperature = np.array([0.0, 2.0, -1.0])

    
    debm = PISM.SurfaceDEBMSimplePointwise(context)

    cgi = np.vectorize(debm.CalovGreveIntegrand)(sigma, temperature)

    assert_array_almost_equal(np.array([0.7979, 2.0000, 0.0833]), cgi, decimal=4)
