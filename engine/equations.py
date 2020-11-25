import numpy as np
import matplotlib.pyplot as plt

#Eg - energy gap, eV
#mn -
elements = {
    "Si": {"Eg": 1.120, "mn": 0.36, "mp": 0.81},
    "GaAs": {"Eg": 1.424, "mn": 0.85, "mp": 0.53},
    "Ge": {"Eg": 0.661, "mn": 0.22, "mp": 0.34}
}
e = 1.6021765e-19 #C
k = 1.3806503e-16 #erg/K

EV_TO_J = e
J_TO_ERG = 1e7
EV_TO_ERG = EV_TO_J * J_TO_ERG

T0 = 300 #K

T0np = 1
T0n = 1
T0pp = 1
T0p = 1
Ae = 6.43e6
Be = 7.13e-12
Ap = 1.8e6
Bp = 1.04e-12

alfa_nd = 2.51e19
beta_nd = EV_TO_ERG / k

def dichotomy_method(x_left, x_right, f, error=1e-5):
    a = x_left
    b = x_right
    c = (a + b) / 2

    f_a = f(a)
    f_b = f(b)
    f_c = f(c)

    REPEAT_AMOUNT = 1
    while abs(f_c) > error:
        if f_c * f_b < 0:
            a = c
            f_a = f_c
        if f_c * f_a < 0:
            b = c
            f_b = f_c
        c = (a + b) / 2
        f_c = f(c)
        REPEAT_AMOUNT += 1

    print("repeats: " + str(REPEAT_AMOUNT))
    return c


class Equations:
    """Physical equations for semiconductor parameters"""
    @staticmethod
    def get_state_density(m, T):
        return m**1.5 * (T/T0)**1.5
    @staticmethod
    def get_electron_concentration(Nc, Ef, Eg, T):
        return Nc * np.exp(beta_nd * (Ef - Eg) / T)
    @staticmethod
    def get_hole_concentration(Nv, Ef, T):
        return Nv * np.exp(-beta_nd * Ef / T)
    @staticmethod
    def get_pos_donor_concentration(Nd0, Ef, Ed, T):
        return Nd0 / (1 + np.exp(beta_nd * (Ed - Ef) / T))
    @staticmethod
    def get_neg_acceptor_concentration(Na0, Ef, Ea, T):
        return Na0 / (1 + np.exp(beta_nd * (Ea - Ef) / T))
    @staticmethod
    def get_electron_mobility(Ndp, Nan, T):
        return Ae / (T**1.5 + Be*(Ndp + Nan)*T**-1.5)#1 / ((T/T0np)**1.5 + (Ndp + Nan)*(T0n/T)**1.5)
    @staticmethod
    def get_hole_mobility(Ndp, Nan, T):
        return Ap / (T**1.5 + Bp*(Ndp + Nan)*T**-1.5)#1 / ((T/T0pp)**1.5 + (Ndp + Nan)*(T0p/T)**1.5)
    @staticmethod
    def get_conductivity(ne, mn, np, mp):
        return e * (ne * mn + np * mp)
    @staticmethod
    def get_fermi_level(_ndp_, _np_, _nan_, _ne_, right_lim, delta=1e-5):
        """Calculates value of fermi level in semiconductor
        \n_ndp_ function of positive donor concentration of fermi level
        \n_np_  function of hole concentration of fermi level
        \n_nan_ function of negative acceptor concentration of fermi level
        \n_ne_ function of electron concentration of fermi level
        \nright_lim value of right bound for fermi level searching
        """
        f = lambda x: _np_(x) - _ne_(x) + _ndp_(x) - _nan_(x)
        return dichotomy_method(0, right_lim, f, delta)


class ConcentrationCalculator:
    """"""

    def __fill_parameters(self):
        self._Ef = np.zeros(len(self._Na_range))
        self._ne = np.zeros(len(self._Na_range))
        self._np = np.zeros(len(self._Na_range))
        self._Ndp = np.zeros(len(self._Na_range))
        self._Nan = np.zeros(len(self._Na_range))
        self._mn = np.zeros(len(self._Na_range))
        self._mp = np.zeros(len(self._Na_range))
        self._s = np.zeros(len(self._Na_range))

    def __init__(self, material=None):
        if material is None:
            self._material = elements["Si"]
        else:
            self._material = material
        self._T = T0
        self._Na_range = np.zeros(1)
        self._Nd_range = np.zeros(1)
        self._Ea = 0
        self._Ed = 0
        self.__fill_parameters()

    def set_material(self, material=None):
        self._material = material

    def set_temperature(self, T):
        self._T = T

    def set_additions(self, Na, Ea, Nd, Ed):
        self._Na_range = np.array(Na) / alfa_nd
        self._Nd_range = np.array(Nd) / alfa_nd
        self._Ea = Ea
        self._Ed = Ed
        self.__fill_parameters()

    def update(self):
        _Nc = Equations.get_state_density(self._material["mn"], self._T)
        _Nv = Equations.get_state_density(self._material["mp"], self._T)
        print("Nc:", _Nc * alfa_nd)
        print("Nv:", _Nv * alfa_nd)
        for i in range(len(self._Na_range)):
            _ndp_ = lambda x: Equations.get_pos_donor_concentration(self._Nd_range[i], x, self._material["Eg"] - self._Ed, self._T)
            _np_ = lambda x: Equations.get_hole_concentration(_Nc, x, self._T)
            _nan_ = lambda x: Equations.get_neg_acceptor_concentration(self._Na_range[i], x, self._Ea, self._T)
            _ne_ = lambda x: Equations.get_electron_concentration(_Nv, x, self._material["Eg"], self._T)

            _Ef = Equations.get_fermi_level(_ndp_, _np_, _nan_, _ne_, self._material["Eg"], 1e-12)
            self._Ef[i] = _Ef
            self._ne[i] = _ne_(_Ef)
            self._np[i] = _np_(_Ef)
            self._Ndp[i] = _ndp_(_Ef)
            self._Nan[i] = _nan_(_Ef)
            self._mn[i] = Equations.get_electron_mobility(self._Ndp[i], self._Nan[i], self._T)
            self._mp[i] = Equations.get_hole_mobility(self._Ndp[i], self._Nan[i], self._T)
            self._s[i] = Equations.get_conductivity(self._ne[i], self._mn[i], self._np[i], self._mp[i])

    def get_fermi_levels(self):
        return self._Ef

    def get_electron_concentration(self):
        return self._ne * alfa_nd

    def get_hole_concentration(self):
        return self._np * alfa_nd

    def get_pos_donor_concentration(self):
        return self._Ndp * alfa_nd

    def get_neg_acceptor_concentration(self):
        return self._Nan * alfa_nd

    def get_electron_mobility(self):
        return self._mn

    def get_hole_mobility(self):
        return self._mp

    def get_conductivity(self):
        return self._s

    def get_resistivity(self):
        return 1 / self._s

# class Material:
#     """Imitates doped semiconductor"""
#
#     def get_electron_concentration(self):
#         return self._ne_(self._Ef)
#     def get_hole_concentration(self):
#         return self._np_(self._Ef)
#     def get_pos_donor_concentration(self):
#         return self._ndp_(self._Ef)
#     def get_neg_donor_concentration(self):
#         return self._nan_(self._Ef)
#
#     def __init__(self, element=None):
#         if element is None:
#             self._Element = elements["Si"]
#         else:
#             self._Element = element
#         self._Ed = 0
#         self._Nd0 = 0
#         self._Ea = 0
#         self._Na0 = 0
#         self._T = T0
#
#     def set_addition(self, Ed=0, Nd0=0, Ea=0, Na0=0):
#         self._Ed = Ed
#         self._Nd0 = Nd0
#         self._Ea = Ea
#         self._Na0 = Na0
#
#     def set_temperature(self, T=T0):
#         self._T = T
#
#     def get_pos_donor_n(self):
#         return self._Nd0 / (1 + np.exp((self._Element["Eg"] - self._Ed - self.Ef)/k/self._T))


if __name__ == '__main__':
    # app = QtWidgets.QApplication(sys.argv)
    # window = Ui()
    # sys.exit(app.exec_())

    calc = ConcentrationCalculator()
    calc.set_temperature(300)
    COUNTS = 100
    Nas = np.linspace(1e14, 1e19, COUNTS, endpoint=True)
    Ea = 0.045
    Nds = np.zeros(COUNTS)
    Ed = 0.045
    calc.set_additions(Nas, Ea, Nds, Ed)
    calc.update()

    #np.set_printoptions(precision=5, suppress=True)
    print("Fermi level", calc.get_fermi_levels())
    print("Electron concentration", calc.get_electron_concentration())
    print("Hole concentration", calc.get_hole_concentration())
    print("Positive acceptor concentration", calc.get_pos_donor_concentration())
    print("Negative acceptor concentration", calc.get_neg_acceptor_concentration())
    print("Electron mobility", calc.get_electron_mobility())
    print("Hole mobility", calc.get_hole_mobility())
    print("Conductivity", calc.get_conductivity())
    print("Resistivity", calc.get_resistivity())

    x = Nas
    y = calc.get_resistivity()

    fig, ax1 = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(16, 12)
    )

    ax1.plot(x, y, '-k')

    ax1.set(xlabel='x', ylabel='P(x)', title='Integral solution')
    ax1.grid(which='major', axis='x')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    plt.show()
