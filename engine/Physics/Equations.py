import numpy as np

# General constants and values
e = 1.6021765e-19  # C
k = 1.3806503e-16  # erg/K
T0 = 300  # K

EV_TO_J = e
J_TO_ERG = 1e7
EV_TO_ERG = EV_TO_J * J_TO_ERG

# Constants for mobility equations
# Si
Ae_Si = 6.43e6
Be_Si = 7.13e-12
Ap_Si = 1.8e6
Bp_Si = 1.04e-12
# Ge
Ae_Ge = 18.7e6
Be_Ge = 30.6e-12
Ap_Ge = 8.02e6
Bp_Ge = 21.3e-12
# GaAs
Ae_GaAs = 53.5e6
Be_GaAs = 14.8e-12
Ap_GaAs = 1.8e6
Bp_GaAs = 2.38e-12

STATE_DENSITY_CONSTANT = 2.51e19


def state_density(m, T):
    """
    Calculates density of states for particles in semiconductor crystal cell.

    Args:
        m(float): effective mass in units of m0 (electron masses).
        T(float): current temperature in Kelvins.

    Returns:
        float: value of state density in concentration control units.
    """
    return STATE_DENSITY_CONSTANT * m ** 1.5 * (T / T0) ** 1.5


def electron_concentration(Nc, Ef, Eg, T):
    """
    Calculates electron concentration in semiconductor crystal cell.

    Args:
        Nc(float): density of states for electrons in concentration control units.
        Ef(float): Fermi level in eV.
        Eg(float): energy gap (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of electron concentration in concentration control units.
    """
    return Nc * np.exp(EV_TO_ERG * (Ef - Eg) / k / T)


def hole_concentration(Nv, Ef, T):
    """
    Calculates hole concentration in semiconductor crystal cell.

    Args:
        Nv(float): density of states for holes in concentration control units.
        Ef(float): Fermi level in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of hole concentration in concentration control units.
    """
    return Nv * np.exp(-EV_TO_ERG * Ef / k / T)


def pos_donor_concentration(Nd0, Ef, Ed, T):
    """
    Calculates positive donor ions concentration in semiconductor crystal cell.

    Args:
        Nd0(float): density of states for electrons in concentration control units.
        Ef(float): Fermi level in eV.
        Ed(float): donor energy level (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of positive donor ions concentration in concentration control units.
    """
    return Nd0 / (1 + np.exp(EV_TO_ERG * (-Ed + Ef) / k / T))


def neg_acceptor_concentration(Na0, Ef, Ea, T):
    """
    Calculates negative acceptor ions concentration in semiconductor crystal cell.

    Args:
        Na0(float): density of states for holes in concentration control units.
        Ef(float): Fermi level in eV.
        Ea(float): acceptor energy level (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of negative acceptor ions concentration in concentration control units.
    """
    return Na0 / (1 + np.exp(EV_TO_ERG * (Ea - Ef) / k / T))


def electron_concentration_d(Nc, Ef, Eg, T):
    """
    Calculates negative acceptor ions concentration derivative in semiconductor crystal cell.

    Args:
        Nc(float): density of states for electrons in concentration control units.
        Ef(float): Fermi level in eV.
        Eg(float): energy gap (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of electron concentration derivative in concentration control units.
    """
    return Nc * np.exp(EV_TO_ERG * (Ef - Eg) / k / T) * EV_TO_ERG / k / T


def hole_concentration_d(Nv, Ef, T):
    """
    Calculates hole concentration derivative in semiconductor crystal cell.

    Args:
        Nv(float): density of states for holes in concentration control units.
        Ef(float): Fermi level in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of hole concentration derivative in concentration control units.
    """
    return -Nv * np.exp(-EV_TO_ERG * Ef / k / T) * EV_TO_ERG / k / T


def pos_donor_concentration_d(Nd0, Ef, Ed, T):
    """
    Calculates positive donor ions concentration derivative in semiconductor crystal cell.

    Args:
        Nd0(float): density of states for electrons in concentration control units.
        Ef(float): Fermi level in eV.
        Ed(float): donor energy level (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of positive donor ions concentration derivative in concentration control units.
    """
    return Nd0 * EV_TO_ERG / k / T * np.exp(EV_TO_ERG * (-Ed + Ef) / k / T) / (1 + np.exp(EV_TO_ERG * (-Ed + Ef) / k / T)) ** 2


def neg_acceptor_concentration_d(Na0, Ef, Ea, T):
    """
    Calculates negative acceptor ions concentration derivative in semiconductor crystal cell.

    Args:
        Na0(float): density of states for holes in concentration control units.
        Ef(float): Fermi level in eV.
        Ea(float): acceptor energy level (from the top of valent zone) in eV.
        T(float): current temperature in Kelvins.

    Returns:
        float: value of negative acceptor ions concentration derivative in concentration control units.
    """
    return Na0 * EV_TO_ERG / k / T * np.exp(EV_TO_ERG * (Ea - Ef) / k / T) / (1 + np.exp(EV_TO_ERG * (Ea - Ef) / k / T)) ** 2


def get_mobility(Ndp, Nan, T, A, B):
    """
    Calculates particle mobility in semiconductor crystal cell.
    By formulae: mobility = Ae/(T**1.5 + Be*(Ndp + Nan)*/T**1.5)

    Args:
        Ndp(float): positive donor ions concentration in concentration control units.
        Nan(float): negative acceptor ions concentration in concentration control units.
        T(float): current temperature in Kelvins.
        A(float): coefficient, describing particle in equation.
        B(float): coefficient, describing particle in equation.

    Returns:
        float: value of electron mobility in concentration control units.
    """
    return A / (T ** 1.5 + B * (Ndp + Nan) * T ** -1.5)


def get_conductivity(ne, mn, np, mp):
    """
    Calculates conductivity of semiconductor crystal cell.

    Args:
        ne(float): electron concentration in concentration control units.
        mn(float): electron mobility in concentration control units.
        np(float): hole concentration in concentration control units.
        mp(float): hole mobility in concentration control units.

    Returns:
        float: value of conductivity in concentration control units.
    """
    return e * (ne * mn + np * mp)
