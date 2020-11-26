import sys
import engine.Physics.StateCalculator as StateCalculator
import gui.last_var as last_var
import numpy as np
import pyqtgraph as pg

from PyQt5 import QtCore, QtGui, QtWidgets , Qt


class MyWin(QtWidgets.QMainWindow):
    temperature = 100
    material = "Si"
    count = 100
    Namin = 0.1
    Namax = 10
    Nadelta = count
    Ndmin = 0.1
    Ndmax = 10
    Nddelta = count
    Ea = 0.045
    Ed = 0.045
    flag = False

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.calc = StateCalculator.StateCalculator(self.material)
        self.ui = last_var.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lineEdit_3.textChanged.connect(self.onChanged_Ndmax)
        self.ui.lineEdit_4.textChanged.connect(self.onChanged_Ed)
        self.ui.lineEdit_5.textChanged.connect(self.onChanged_Ndmin)
        self.ui.lineEdit_6.textChanged.connect(self.onChanged_Ea)
        self.ui.lineEdit_7.textChanged.connect(self.onChanged_Namin)
        self.ui.lineEdit_8.textChanged.connect(self.onChanged_Namax)
        self.ui.lineEdit_9.textChanged.connect(self.onChanged_Temp)
        self.ui.lineEdit_10.textChanged.connect(self.onChanged_count)
        self.ui.pushButton.clicked.connect(self.work_with_button)
        self.ui.checkBox.stateChanged.connect(self.set_log_axis)

    def set_log_axis(self):
        if self.ui.checkBox.isChecked():
            self.flag = True
        else: self.flag = False

    def work_with_button(self):

        self.material = self.ui.ComboBox.currentText()
        self.calc.set_material(str(self.material))
        self.calc.set_temperature(float(self.temperature))
        self.Na = np.linspace(float(self.Namin)*1e18, float(self.Namax)*1e18, int(self.count))
        self.Nd = np.linspace(float(self.Ndmin)*1e18, float(self.Ndmax)*1e18, int(self.count))
        self.calc.set_additions(self.Na, float(self.Ea), self.Nd, float(self.Ed))
        self.calc.update()
        self.full_all_field()
        self.draw_graphics()

    def draw_graph(self, elem, name_graph, y_values):
        pen = pg.mkPen(color=(0, 0, 0), width=1)
        styles = {"color": "#000", "font-size": "20px"}
        elem.setLabel("left", name_graph, **styles)
        elem.setLabel("bottom", "Acceptor concentration Na (1e18  cmˆ-3)", **styles)
        elem.addLegend()
        elem.setLogMode(False, self.flag)
        elem.plot(self.Na / 1e18, y_values, pen=pen, clear=True)

    def draw_graphics(self):

        self.draw_graph(self.ui.tab_4, "Electron concentration (cmˆ-3)", self.Electron_concentration)
        self.draw_graph(self.ui.tab, "Electron mobility (cmˆ2 * Vˆ1 * sˆ-1)", self.electron_mobility)
        self.draw_graph(self.ui.tab_2, "Hole concentration (1e18 cmˆ-3)", self.Hole_concentration / 1e18)
        self.draw_graph(self.ui.tab_3, "Negative acceptor concentration (1e18 cmˆ-3)",
                        self.Negative_acceptor_concentration / 1e18)
        self.draw_graph(self.ui.tab_5, "Positive acceptor concentration (cmˆ-3)", self.Positive_acceptor_concentration)
        self.draw_graph(self.ui.tab_6, "Conductivity (Ωˆ-1 * cmˆ-1)", self.Conductivity)
        self.draw_graph(self.ui.tab_7, "Resistivity (Ω * cm)", self.Resistivity)

    def full_all_field(self):
        self.Electron_concentration = self.calc.get_electron_concentration()
        self.Hole_concentration = self.calc.get_hole_concentration()
        self.Positive_acceptor_concentration = self.calc.get_pos_donor_concentration()
        self.Negative_acceptor_concentration = self.calc.get_neg_acceptor_concentration()
        self.electron_mobility = self.calc.get_electron_mobility()
        self.Hole_mobility = self.calc.get_hole_mobility()
        self.Conductivity = self.calc.get_conductivity()
        self.Resistivity = self.calc.get_resistivity()

    def onChanged_Temp(self, text):
        self.temperature = text

    def onChanged_Ed(self, text):
        self.Ed = text

    def onChanged_Ea(self, text):
        self.Ea = text

    def onChanged_Namin(self, text):
        self.Namin = text

    def onChanged_Namax(self, text):
        self.Namax = text

    def onChanged_count(self, text):
        self.count = text

    def onChanged_Ndmin(self, text):
        self.Ndmin = text

    def onChanged_Ndmax(self, text):
        self.Ndmax = text


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = MyWin()
    myapp.show()
    sys.exit(app.exec_())