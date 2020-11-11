import numpy as np
import matplotlib.pyplot as plt

# General parameters of semiconductor materials
# Eg - energy gap, eV
# mn - electron effective mass, amounts of m0
# mp - hole effective mass, amounts of m0
elements = {
    "Si": {"Eg": 1.120, "mn": 0.36, "mp": 0.81},
    "GaAs": {"Eg": 1.424, "mn": 0.85, "mp": 0.53},
    "Ge": {"Eg": 0.661, "mn": 0.22, "mp": 0.34}
}

# General constants and values
e = 1.6021765e-19  # C
k = 1.3806503e-16  # erg/K

EV_TO_J = e
J_TO_ERG = 1e7
EV_TO_ERG = EV_TO_J * J_TO_ERG

T0 = 300  # K

# Constants for mobility equations
# T0np = 1
# T0n = 1
# T0pp = 1
# T0p = 1
Ae = 6.43e6
Be = 7.13e-12
Ap = 1.8e6
Bp = 1.04e-12

# Constants for reducing magnitude of values during computations
alfa_nd = 2.51e19
beta_nd = EV_TO_ERG / k


def newton_method(x0, f, f_d, error=1e-5):
    x = x0

    REPEAT_AMOUNT = 1
    while abs(f(x)) > error:
        x = x - f(x) / f_d(x)
        REPEAT_AMOUNT += 1

    # print("repeats: " + str(REPEAT_AMOUNT))
    return (x, REPEAT_AMOUNT)


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

    return (c, REPEAT_AMOUNT)


class Equations:  # TODO 1: shall it be replaced with a module with functions?
    """Physical equations for semiconductor parameters"""

    @staticmethod
    def get_state_density(m, T):
        return m ** 1.5 * (T / T0) ** 1.5

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
    def get_electron_concentration_d(Nc, Ef, Eg, T):
        return Nc * np.exp(beta_nd * (Ef - Eg) / T) * beta_nd / T

    @staticmethod
    def get_hole_concentration_d(Nv, Ef, T):
        return -Nv * np.exp(-beta_nd * Ef / T) * beta_nd / T

    @staticmethod
    def get_pos_donor_concentration_d(Nd0, Ef, Ed, T):
        return Nd0 * beta_nd / T * np.exp(beta_nd * (Ed - Ef) / T) / (1 + np.exp(beta_nd * (Ed - Ef) / T)) ** 2

    @staticmethod
    def get_neg_acceptor_concentration_d(Na0, Ef, Ea, T):
        return Na0 * beta_nd / T * np.exp(beta_nd * (Ea - Ef) / T) / (1 + np.exp(beta_nd * (Ea - Ef) / T)) ** 2

    @staticmethod
    def get_electron_mobility(Ndp, Nan, T):
        return Ae / (T ** 1.5 + Be * (
                    Ndp + Nan) * T ** -1.5)  # 1/((T/T0np)**1.5 + (Ndp + Nan)*(T0n/T)**1.5) #TODO 2: Shall we use our own equations for mobility estimation?

    @staticmethod
    def get_hole_mobility(Ndp, Nan, T):
        return Ap / (T ** 1.5 + Bp * (
                    Ndp + Nan) * T ** -1.5)  # 1/((T/T0pp)**1.5 + (Ndp + Nan)*(T0p/T)**1.5) #TODO 2: Shall we use our own equations for mobility estimation?

    @staticmethod
    def get_conductivity(ne, mn, np, mp):
        return alfa_nd * e * (ne * mn + np * mp)


#     @staticmethod
#     def get_fermi_level(_ndp_, _np_, _nan_, _ne_, right_lim, delta=1e-5):#TODO 3: We probably should use newton method instead of dichotomy, but im lazy and dont want to calculate derivative)
#         """Calculates value of fermi level in semiconductor using DICHOTOMY METHOD of equation solving
#         \n_ndp_ function of positive donor concentration of fermi level
#         \n_np_  function of hole concentration of fermi level
#         \n_nan_ function of negative acceptor concentration of fermi level
#         \n_ne_ function of electron concentration of fermi level
#         \nright_lim value of right bound for fermi level searching
#         """
#         f = lambda x: _np_(x) - _ne_(x) + _ndp_(x) - _nan_(x)
#         return dichotomy_method(0, right_lim, f, delta)


class ConcentrationCalculator:  # TODO 4: Shall we always set concentration parameter as an array?
    # Can it be just a single value?
    # Then we dont need to initialize variables as an arrays, but we need to call methods multiple times..
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
            self._material = material
        self._T = T
        self._Na_range = np.zeros(1)
        self._Nd_range = np.zeros(1)
        self._Ea = 0
        self._Ed = 0
        self.__fill_parameters()

    def set_material(self, material=None):
        """Changes using material in calculator to another"""
        self._material = material

    def set_temperature(self, T):
        """Changes using temperature in calculator to another"""
        self._T = T

    def set_additions(self, Na, Ea, Nd, Ed):
        """Changes using concentrations and energies of additions for material in calculator to other"""
        self._Na_range = np.array(Na) / alfa_nd
        self._Nd_range = np.array(Nd) / alfa_nd
        self._Ea = Ea
        self._Ed = Ed
        self.__fill_parameters()

    def update(self):
        """Processes set parameters and finds characteristics of interest"""
        _Nc = Equations.get_state_density(self._material["mn"], self._T)
        _Nv = Equations.get_state_density(self._material["mp"], self._T)
        for i in range(len(self._Na_range)):
            _ndp_ = lambda x: Equations.get_pos_donor_concentration(self._Nd_range[i], x,
                                                                    self._material["Eg"] - self._Ed, self._T)
            _np_ = lambda x: Equations.get_hole_concentration(_Nv, x, self._T)
            _nan_ = lambda x: Equations.get_neg_acceptor_concentration(self._Na_range[i], x, self._Ea, self._T)
            _ne_ = lambda x: Equations.get_electron_concentration(_Nc, x, self._material["Eg"], self._T)

            _ndp_d_ = lambda x: Equations.get_pos_donor_concentration_d(self._Nd_range[i], x,
                                                                        self._material["Eg"] - self._Ed, self._T)
            _np_d_ = lambda x: Equations.get_hole_concentration_d(_Nv, x, self._T)
            _nan_d_ = lambda x: Equations.get_neg_acceptor_concentration_d(self._Na_range[i], x, self._Ea, self._T)
            _ne_d_ = lambda x: Equations.get_electron_concentration_d(_Nc, x, self._material["Eg"], self._T)

            balance_equation = lambda x: _np_(x) - _ne_(x) + _ndp_(x) - _nan_(x)
            balance_equation_d = lambda x: _np_d_(x) - _ne_d_(x) + _ndp_d_(x) - _nan_d_(x)

            _Ef = newton_method(0, balance_equation, balance_equation_d, 1e-10)
            # _Ef = dichotomy_method(0, self._material["Eg"], balance_equation, 1e-10)
            print(balance_equation(_Ef[0]))
            print(balance_equation_d(_Ef[0]))
            print(_Ef[1])
            _Ef = _Ef[0]
            self._Ef[i] = _Ef
            self._ne[i] = _ne_(_Ef)
            self._np[i] = _np_(_Ef)
            self._Ndp[i] = _ndp_(_Ef)
            self._Nan[i] = _nan_(_Ef)
            self._mn[i] = Equations.get_electron_mobility(self._Ndp[i], self._Nan[i], self._T)
            self._mp[i] = Equations.get_hole_mobility(self._Ndp[i], self._Nan[i], self._T)
            self._s[i] = Equations.get_conductivity(self._ne[i], self._mn[i], self._np[i], self._mp[i])

    def get_fermi_levels(self):
        """Returns value of Fermi level for set parameters"""
        return self._Ef

    def get_electron_concentration(self):
        """Returns value of electron concentration for set parameters"""
        return self._ne * alfa_nd

    def get_hole_concentration(self):
        """Returns value of hole concentration for set parameters"""
        return self._np * alfa_nd

    def get_pos_donor_concentration(self):
        """Returns value of positive charged donor ions concentration for set parameters"""
        return self._Ndp * alfa_nd

    def get_neg_acceptor_concentration(self):
        """Returns value of negative charged acceptor ions concentration for set parameters"""
        return self._Nan * alfa_nd

    def get_electron_mobility(self):
        """Returns value of electron mobility for set parameters"""
        return self._mn

    def get_hole_mobility(self):
        """Returns value of hole mobility for set parameters"""
        return self._mp

    def get_conductivity(self):
        """Returns value of conductivity for set parameters"""
        return self._s

    def get_resistivity(self):
        """Returns value of resistivity for set parameters"""
        return 1 / self.get_conductivity()


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
