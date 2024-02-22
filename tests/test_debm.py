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

def test_hour_angle():
    """
    Test the calculation of the hour angle
    """

    phi = np.array([0.0, np.pi / 4.0, np.pi / 2.0])
    latitude = np.array([-np.pi / 2.0, 0.0, np.pi / 4.0])
    declination = np.array([np.pi / 8.0, 0.0, 0.0])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    hour_angle = np.vectorize(debm.hour_angle)(phi, latitude, declination)
    assert_array_almost_equal(np.array([0.0, 0.78539816, 0.0]), hour_angle, decimal=4)


def test_solar_longitude():
    """
    Test solar longitude
    """
    year_fraction = np.array([0.0, 1.0 / 12.0, 1.0])
    eccentricity = np.array([0.9, 1.0, 1.0])
    perhelion_longitude = np.array([np.pi / 8.0, -np.pi, 0.0])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    solar_longitude = np.vectorize(debm.solar_longitude)(
        year_fraction, eccentricity, perhelion_longitude
    )
    assert_array_almost_equal(
        np.array([-2.4174, -0.1785, 3.6222]), solar_longitude, decimal=4
    )


def test_solar_declination_present_day():
    """
    Test solar declination present day
    """
    year_fraction = np.array([0.0, 1.0 / 12.0, 1.0])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    solar_declination = np.vectorize(debm.solar_declination_present_day)(year_fraction)
    assert_array_almost_equal(
        np.array([-0.402449  , -0.30673297, -0.402449  ]), solar_declination, decimal=4
    )
    
def test_solar_declination_paleo():
    """
    Test solar declination paleo
    """

    obliquity = np.array([np.pi / 4])
    solar_longitude = np.array([np.pi * 3 / 4])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    solar_declination = np.vectorize(debm.solar_declination_paleo)(obliquity, solar_longitude)
    assert_array_almost_equal(
        np.array([0.55536037]), solar_declination, decimal=4
    )

def test_distance_factor_present_day():
    """
    Test distance factor present day
    """
    year_fraction = np.array([0.0, 1.0 / 12.0, 1.0])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    d = np.vectorize(debm.distance_factor_present_day)(year_fraction)
    assert_array_almost_equal(np.array([1.03505, 1.03081244, 1.03505   ]), d, decimal=4)


def test_distance_factor_paleo():
    """
    Test distance factor paleo
    """

    eccentricity = np.array([0.0167, 0.03])
    obliquity = np.deg2rad(np.array([23.14, 22.10]))
    perihelion_longitude = np.deg2rad(np.array([102.94719, -44.3]))

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    d = np.vectorize(debm.distance_factor_paleo)(eccentricity, perihelion_longitude, obliquity)

def test_insolation():
    """
    Test insolation
    """
    
    solar_constant = np.array([1361.0])
    distance_factor = np.array([1.1])
    hour_angle = np.array([0.8])
    latitude = np.array([np.pi / 4])
    declination = np.array([np.pi / 8])

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    insolation = np.vectorize(debm.insolation)(solar_constant, distance_factor, hour_angle, latitude, declination)

    assert_array_almost_equal(np.array([1282.10500694]), insolation, decimal=4)


def test_orbital_parameters():
    """
    Test orbital parameters
    """
    time = 2022.25
    
    debm = PISM.SurfaceDEBMSimplePointwise(context)
    orbital_parameters = debm.orbital_parameters(time)

    assert_almost_equal(-0.4024197463537702, orbital_parameters.declination, decimal=4)
    assert_almost_equal(1.0350505747618992, orbital_parameters.distance_factor, decimal=4)




def test_albedo():
    """
    Test albedo
    """

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    melt_rate = 1. / 3.15569259747e7
    albedo = debm.albedo(melt_rate, 2)
    
    assert_almost_equal(0.7972189470997155, albedo, decimal=4)

    
def test_atmospshere_transmissivity():
    """
    Test atmosphere transmissivity
    """

    debm = PISM.SurfaceDEBMSimplePointwise(context)
    elevation = np.array([0.0, 1000.0, 2000.0])
    transmissivity = np.vectorize(debm.atmosphere_transmissivity)(elevation)
    assert_array_almost_equal(np.array([0.65 , 0.682, 0.714]), transmissivity, decimal=4)

    
def test_melt():
    """
    Test melt
    """

    year_fraction = 0.
    dt = 1.0 / 12 
    temp = 323.0
    temp_sd = 12.0
    surface_elevation = 1000
    latitude = np.pi/4 * 3
    albedo = 0.47
    
    debm = PISM.SurfaceDEBMSimplePointwise(context)
    orbital_parameters = debm.orbital_parameters(year_fraction)
    declination = orbital_parameters.declination
    distance_factor = orbital_parameters.distance_factor
    melt_info = debm.melt(declination, distance_factor, dt, temp_sd, temp, surface_elevation, latitude, albedo)

    assert_almost_equal(3.41058271637613e-08, melt_info.insolation_melt, decimal=4)
    assert_almost_equal(1.3995243305283176e-07, melt_info.temperature_melt, decimal=4)
    assert_almost_equal(-9.003261572844215e-09, melt_info.offset_melt, decimal=4)
    assert_almost_equal(1.6505499864374886e-07, melt_info.total_melt, decimal=4)
