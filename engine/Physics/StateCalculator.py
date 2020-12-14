import numpy as np

import engine.Physics.Equations as Equations
import engine.Methematics.NumericalMethods as Math


# General parameters of semiconductor materials
# Eg - energy gap, eV
# mn - electron effective mass, amounts of m0
# mp - hole effective mass, amounts of m0
elements = {
    "Si": {"Eg": 1.120, "mn": 0.36, "mp": 0.81,
           "Ae": Equations.Ae_Si, "Ap": Equations.Ap_Si, "Be": Equations.Be_Si, "Bp": Equations.Bp_Si},
    "GaAs": {"Eg": 1.424, "mn": 0.85, "mp": 0.53,
             "Ae": Equations.Ae_GaAs, "Ap": Equations.Ap_GaAs, "Be": Equations.Be_GaAs, "Bp": Equations.Bp_GaAs},
    "Ge": {"Eg": 0.661, "mn": 0.22, "mp": 0.34,
           "Ae": Equations.Ae_Ge, "Ap": Equations.Ap_Ge, "Be": Equations.Be_Ge, "Bp": Equations.Bp_Ge}
}

class StateCalculator:
    # Can it be just a single value?
    # Then we don't need to initialize variables as an arrays, but we need to call methods multiple times..
    """Allows to set, get, and keep parameters during computation"""

    def __fill_parameters(self):
        """Resets saved parameters and initialize them to be zeroes"""
        self._Ef = np.zeros(len(self._Na_range))
        self._ne = np.zeros(len(self._Na_range))
        self._np = np.zeros(len(self._Na_range))
        self._Ndp = np.zeros(len(self._Na_range))
        self._Nan = np.zeros(len(self._Na_range))
        self._mn = np.zeros(len(self._Na_range))
        self._mp = np.zeros(len(self._Na_range))
        self._s = np.zeros(len(self._Na_range))

    def __init__(self, material=None, T=300):
        """Initialize calculator, default material is silicon"""
        if material is None:
            self._material = elements["Si"]
        else:
            self._material = elements[material]
        self._T = T
        self._Na_range = np.zeros(1)
        self._Nd_range = np.zeros(1)
        self._Ea = 0
        self._Ed = 0
        self.__fill_parameters()

    def set_material(self, material=None):
        """Changes using material in calculator to another"""
        self._material = elements[material]

    def set_temperature(self, T):
        """Changes using temperature in calculator to another"""
        self._T = T

    def set_additions(self, Na, Ea, Nd, Ed):
        """Changes using concentrations and energies of additions for material in calculator to other"""
        self._Na_range = np.array(Na)
        self._Nd_range = np.array(Nd)
        self._Ea = Ea
        self._Ed = Ed
        self.__fill_parameters()

    def update(self):
        """Processes set parameters and finds characteristics of interest"""
        _Nc = Equations.state_density(self._material["mn"], self._T)
        _Nv = Equations.state_density(self._material["mp"], self._T)
        for i in range(len(self._Na_range)):
            def balance_equation(x):
                return (
                        Equations.hole_concentration(_Nv, x, self._T) -
                        Equations.electron_concentration(_Nc, x, self._material["Eg"], self._T) -
                        Equations.neg_acceptor_concentration(self._Na_range[i], x, self._Ea, self._T) +
                        Equations.pos_donor_concentration(self._Nd_range[i], x, self._material["Eg"], self._Ed,
                                                          self._T))

            def balance_equation_d(x):
                return (
                        Equations.hole_concentration_d(_Nv, x, self._T) -
                        Equations.electron_concentration_d(_Nc, x, self._material["Eg"], self._T) -
                        Equations.neg_acceptor_concentration_d(self._Na_range[i], x, self._Ea, self._T) +
                        Equations.pos_donor_concentration_d(self._Nd_range[i], x, self._material["Eg"], self._Ed,
                                                            self._T))

            # _Ef = Math.newton_method(self._material["Eg"]/2, balance_equation, balance_equation_d, 60)
            _Ef = Math.dichotomy_method(-self._material["Eg"], 2*self._material["Eg"], balance_equation, 100)

            _Ef = _Ef[0]
            self._Ef[i] = _Ef
            self._ne[i] = Equations.electron_concentration(_Nc, _Ef, self._material["Eg"], self._T)
            self._np[i] = Equations.hole_concentration(_Nv, _Ef, self._T)
            self._Ndp[i] = Equations.pos_donor_concentration(self._Nd_range[i], _Ef, self._material["Eg"], self._Ed,
                                                             self._T)
            self._Nan[i] = Equations.neg_acceptor_concentration(self._Na_range[i], _Ef, self._Ea, self._T)
            self._mn[i] = Equations.get_mobility(
                self._Ndp[i], self._Nan[i], self._T, self._material["Ae"], self._material["Be"])
            self._mp[i] = Equations.get_mobility(
                self._Ndp[i], self._Nan[i], self._T, self._material["Ap"], self._material["Bp"])
            self._s[i] = Equations.get_conductivity(self._ne[i], self._mn[i], self._np[i], self._mp[i])

    def get_fermi_levels(self):
        """Returns value of Fermi level for set parameters"""
        return self._Ef

    def get_electron_concentration(self):
        """Returns value of electron concentration for set parameters"""
        return self._ne

    def get_hole_concentration(self):
        """Returns value of hole concentration for set parameters"""
        return self._np

    def get_pos_donor_concentration(self):
        """Returns value of positive charged donor ions concentration for set parameters"""
        return self._Ndp

    def get_neg_acceptor_concentration(self):
        """Returns value of negative charged acceptor ions concentration for set parameters"""
        return self._Nan

    def get_electron_mobility(self):
        """Returns value of electron mobility for set parameters"""
        return self._mn

    def get_hole_mobility(self):
        """Returns value of hole mobility for set parameters"""
        return self._mp

    def get_conductivity(self):
        """Returns value of conductivity for set parameters"""
        print((1 / self._s)[-1])
        return 1 / self._s

    def get_resistivity(self):
        """Returns value of resistivity for set parameters"""
        return 1 / self.get_conductivity()
