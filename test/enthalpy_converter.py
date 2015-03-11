import PISM
import numpy as np

config = PISM.Context().config

# list of converters
converters = {"Default": PISM.EnthalpyConverter(config),
              "Glint2 (Kirchhoff)": PISM.KirchhoffEnthalpyConverter(config),
              "verification (cold)": PISM.ColdEnthalpyConverter(config),
              "linear-in-temp C(T)": PISM.varcEnthalpyConverter(config)}


def try_all_converters(test):
    print ""
    for name, converter in converters.items():
        print "Testing '%s' converter..." % name
        test(name, converter)


def reversibility_test():
    "Converting from (E, P) to (T, omega, P)"

    def run(name, EC):
        # for a fixed pressure...
        H = 1000.0
        P = EC.pressure(H)

        # cold ice
        # enthalpy form
        E = EC.enthalpy(250.0, 0.0, P)
        # temperature form
        T = EC.temperature(E, P)
        omega = EC.water_fraction(E, P)
        # we should get the same E back
        assert E == EC.enthalpy(T, omega, P)

        # temperate ice
        # enthalpy form
        T_m = EC.melting_temperature(P)
        E = EC.enthalpy(T_m, 0.1, P)
        # temperature form
        T = EC.temperature(E, P)
        omega = EC.water_fraction(E, P)
        # we should get the same E back
        assert E == EC.enthalpy(T, omega, P)

    try_all_converters(run)


def temperate_temperature_test():
    "For temperate ice, an increase of E should not change T."

    def run(name, EC):
        H = 1000.0
        P = EC.pressure(H)
        T_m = EC.melting_temperature(P)
        E = EC.enthalpy(T_m, 0.005, P)

        if not EC.is_temperate(E, P):
            # skip the test if our converter still treats this ice as
            # cold
            return

        for delta_E in [0, 100, 1000]:
            assert EC.temperature(E + delta_E, P) == T_m

    try_all_converters(run)


def cts_computation_test():
    "E_cts should be the same no matter how you compute it."

    def run(name, EC):
        H = 1000.0
        P = EC.pressure(H)
        T_m = EC.melting_temperature(P)

        assert EC.enthalpy(T_m, 0.0, P) == EC.enthalpy_cts(P)

    try_all_converters(run)


def water_fraction_at_cts_test():
    "Water fraction at CTS is zero."

    def run(name, EC):
        H = 1000.0
        P = EC.pressure(H)
        E = EC.enthalpy_cts(P)

        assert EC.water_fraction(E, P) == 0

    try_all_converters(run)


def enthalpy_of_water_test():
    """Test the dependence of the enthalpy of water at T_m(p) on p."""

    config = PISM.Context().config
    c_w = config.get("water_specific_heat_capacity")

    EC = converters["Glint2 (Kirchhoff)"]

    depth0 = 0.0
    p0 = EC.pressure(depth0)
    T0 = EC.melting_temperature(p0)
    omega0 = 1.0
    E0 = EC.enthalpy(T0, omega0, p0)

    depth1 = 1000.0
    p1 = EC.pressure(depth1)
    T1 = EC.melting_temperature(p1)
    omega1 = 1.0
    E1 = EC.enthalpy(T1, omega1, p1)

    # if we change the pressure of water from p0 to p1 while keeping
    # it at T_m(p), its enthalpy should change and this change should
    # be equal to c_w * (T_m(p1) - T_m(p0))
    assert np.fabs((E1 - E0) - c_w * (T1 - T0)) < 1e-9


def plot_converter(name, EC):
    """Test an enthalpy converter passed as the argument."""

    H = 5000.0                  # ice thickness

    Z = np.linspace(0, H, int(H / 10))  # vertical levels

    p = np.zeros_like(Z)
    T_melting = np.zeros_like(Z)
    E_cts = np.zeros_like(Z)
    E = np.zeros_like(Z)
    T = np.zeros_like(Z)
    E_wet = np.zeros_like(Z)
    omega = np.zeros_like(Z)

    E[:] = 97000.0
    delta_E = 3000.0

    for i, z in enumerate(Z):
        depth = H - z
        p[i] = EC.pressure(depth)
        T_melting[i] = EC.melting_temperature(p[i])
        E_cts[i] = EC.enthalpy_cts(p[i])
        T[i] = EC.temperature(E[i], p[i])
        # dependence on pressure for high omega and T_melting
        E_wet[i] = EC.enthalpy(T_melting[i], 0.95, p[i])
        omega[i] = EC.water_fraction(E_cts[i] + delta_E, p[i]) * 100.0

    plt.figure(figsize=(15, 10))
    plt.subplot(2, 2, 1)
    plt.title("%s enthalpy converter" % name)
    plt.plot(Z, E, label="constant enthalpy", lw=2)
    plt.plot(Z, E_cts, label="enthalpy corresponding to CTS", lw=2)
    plt.legend(loc='best')
    plt.ylabel("J/kg")
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.title("%s enthalpy converter" % name)
    plt.plot(Z, omega, label="water fraction for E = E_cts + C", lw=2)
    plt.legend(loc='best')
    plt.ylabel("percent")
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(Z, T_melting, label="melting temperature", lw=2)
    plt.plot(Z, T, label="temperature corresponding to constant E", lw=2)
    plt.legend(loc='best')
    plt.ylabel("Kelvin")
    plt.grid()
    plt.xlabel("height above the base of the ice, m")

    plt.subplot(2, 2, 4)
    plt.plot(Z, E_wet, label="temperate ice enthalpy with high omega", lw=2)
    plt.legend(loc='best')
    plt.xlabel("height above the base of the ice, m")
    plt.ylabel("J/kg")
    plt.grid()


def compare():
    for name, converter in converters.items():
        plot_converter(name, converter)

    plt.show()

if __name__ == "__main__":
    import pylab as plt
    compare()
