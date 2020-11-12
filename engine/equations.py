import numpy as np
import matplotlib.pyplot as plt

from engine.Physics.StateCalculator import StateCalculator

if __name__ == '__main__':
    # app = QtWidgets.QApplication(sys.argv)
    # window = Ui()
    # sys.exit(app.exec_())

    calc = StateCalculator()
    calc.set_temperature(300)
    COUNTS = 100
    Nas = np.linspace(1e14, 1e19, COUNTS, endpoint=True)
    Ea = 0.045
    Nds = np.zeros(COUNTS)
    Ed = 0.045
    calc.set_additions(Nas, Ea, Nds, Ed)
    calc.update()

    # np.set_printoptions(precision=5, suppress=True)
    # print("Fermi level", calc.get_fermi_levels())
    # print("Electron concentration", calc.get_electron_concentration())
    # print("Hole concentration", calc.get_hole_concentration())
    # print("Positive acceptor concentration", calc.get_pos_donor_concentration())
    # print("Negative acceptor concentration", calc.get_neg_acceptor_concentration())
    # print("Electron mobility", calc.get_electron_mobility())
    # print("Hole mobility", calc.get_hole_mobility())
    # print("Conductivity", calc.get_conductivity())
    # print("Resistivity", calc.get_resistivity())

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
