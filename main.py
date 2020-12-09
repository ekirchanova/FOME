import sys
import engine.Physics.StateCalculator as StateCalculator
import gui.last_var as last_var
import numpy as np
import pyqtgraph as pg
import csv

from PyQt5 import QtCore, QtGui, QtWidgets , Qt


class MyWin(QtWidgets.QMainWindow):
    temperature = 300
    material = "Si"
    count = 100
    Namin = 0
    Namax = 0
    Ndmin = 0
    Ndmax = 0
    Ea = 0
    Ed = 0
    Nmin = 0.1
    Nmax = 1
    E = 0.045
    flag_logx = False
    flag_logy = False
    flag_n = True
    flag_p = False
    index = 3

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.calc = StateCalculator.StateCalculator(self.material)
        self.ui = last_var.Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lineEdit_3.textChanged.connect(self.onChanged_Nmax)
        self.ui.lineEdit_4.textChanged.connect(self.onChanged_E)
        self.ui.lineEdit_5.textChanged.connect(self.onChanged_Nmin)
        self.ui.lineEdit_9.textChanged.connect(self.onChanged_Temp)
        self.ui.lineEdit_10.textChanged.connect(self.onChanged_count)

        self.ui.pushButton.clicked.connect(self.work_with_button)
        self.ui.checkBox.stateChanged.connect(self.set_log_x_axis)
        self.ui.checkBox_2.stateChanged.connect(self.set_log_y_axis)
        self.ui.checkBox_3.stateChanged.connect(self.set_n)
        self.ui.checkBox_4.stateChanged.connect(self.set_p)
        self.ui.pushButton_2.clicked.connect(self.save_graphs)

    def save_graph(self,name_file, y_values):
        arxiv = []
        if (self.flag_p ):
            arxiv = zip(self.Na, y_values)
            name_file = name_file + "_Na.csv"
        elif (self.flag_n ):
            arxiv = zip(self.Nd, y_values)
            name_file = name_file + "_Nd.csv"
        with open(name_file, 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerows(arxiv)

    def save_graphs(self):
        self.save_graph("Electron_concentration", self.Electron_concentration)
        self.save_graph("Electron_mobility", self.Electron_mobility)
        self.save_graph("Hole_mobility", self.Hole_mobility)
        self.save_graph("Hole_concentration", self.Hole_concentration)
        self.save_graph("Negative_acceptor_concentration",
                        self.Negative_acceptor_concentration)
        self.save_graph("Positive_acceptor_concentration", self.Positive_acceptor_concentration)
        self.save_graph("Conductivity", self.Conductivity)
        self.save_graph("Resistivity", self.Resistivity)

    def set_n(self):
        if self.ui.checkBox_3.isChecked():
            self.flag_n = True
            self.flag_p = False
            self.ui.checkBox_4.setChecked(0)

    def set_p(self):
        if self.ui.checkBox_4.isChecked():
            self.flag_p = True
            self.flag_n = False
            self.ui.checkBox_3.setChecked(0)

    def set_log_x_axis(self):
        if self.ui.checkBox.isChecked():
            self.flag_logx = True
        else: self.flag_logx = False

    def set_log_y_axis(self):
        if self.ui.checkBox.isChecked():
            self.flag_logy = True
        else: self.flag_logy = False

    def work_with_flag(self):
        if (self.flag_p):
            self.Namin = self.Nmin
            self.Namax = self.Nmax
            self.Ea = self.E
            self.Ndmin = self.Nmin
            self.Ndmax = self.Nmax
            self.Ed = 0
        elif(self.flag_n):
            self.Ndmin = self.Nmin
            self.Ndmax = self.Nmax
            self.Ed = self.E
            self.Namin = self.Nmin
            self.Namax = self.Nmax
            self.Ea = 0

    def work_with_button(self):
        self.work_with_flag()
        self.index = self.ui.tabWidget.currentIndex()
        self.material = self.ui.ComboBox.currentText()
        self.calc.set_material(str(self.material))
        self.calc.set_temperature(float(self.temperature))
        self.Na = np.linspace(float(self.Namin)*1e18, float(self.Namax)*1e18, int(self.count))
        self.Nd = np.linspace(float(self.Ndmin)*1e18, float(self.Ndmax)*1e18, int(self.count))
        self.calc.set_additions(self.Na, float(self.Ea), self.Nd, float(self.Ed))
        self.calc.update()
        self.full_all_field()
        print(self.Hole_concentration)
        self.draw_graphics()

    def draw_graph(self, elem, name_graph, y_values):
        pen = pg.mkPen(color=(0, 0, 0), width=1)
        styles = {"color": "#000", "font-size": "20px"}
        elem.setLabel("left", name_graph, **styles)
        elem.addLegend()
        elem.setLogMode(self.flag_logx, self.flag_logy)
        if(self.flag_p):
            elem.setLabel("bottom", "Acceptor concentration Na (1e18  cmˆ-3)", **styles)
            elem.plot(self.Na/1e18 , y_values, pen=pen, clear=True)
        elif(self.flag_n):
            elem.setLabel("bottom", "Donor concentration Nd (1e18  cmˆ-3)", **styles)
            elem.plot(self.Nd/1e18, y_values, pen=pen, clear=True)


    def draw_graphics(self):
        self.draw_graph(self.ui.tab_4, "Electron concentration (cmˆ-3)", self.Electron_concentration/ 1e18)
        self.draw_graph(self.ui.tab_8, "Electron mobility (1e18 cmˆ2 * Vˆ1 * sˆ-1)", self.Electron_mobility)
        self.draw_graph(self.ui.tab_9, "Hole mobility (cmˆ2 * Vˆ1 * sˆ-1)", self.Hole_mobility)
        self.draw_graph(self.ui.tab_2, "Hole concentration (1e18 cmˆ-3)", self.Hole_concentration / 1e18)
        if(self.flag_n):
            self.draw_graph(self.ui.tab_3, "Positive donor concentration (1e18 cmˆ-3)", self.Positive_acceptor_concentration/1e18)
        elif(self.flag_p):
            self.draw_graph(self.ui.tab_3, "Negative acceptor concentration (1e18 cmˆ-3)",
                            self.Negative_acceptor_concentration / 1e18)

        self.draw_graph(self.ui.tab_6, "Conductivity (Ωˆ-1 * cmˆ-1)", self.Conductivity)
        self.draw_graph(self.ui.tab_7, "Resistivity (Ω * cm)", self.Resistivity)

    def full_all_field(self):
        self.Electron_concentration = self.calc.get_electron_concentration()
        self.Hole_concentration = self.calc.get_hole_concentration()
        self.Positive_acceptor_concentration = self.calc.get_pos_donor_concentration()
        self.Negative_acceptor_concentration = self.calc.get_neg_acceptor_concentration()
        self.Electron_mobility = self.calc.get_electron_mobility()
        self.Hole_mobility = self.calc.get_hole_mobility()
        self.Conductivity = self.calc.get_conductivity()
        self.Resistivity = self.calc.get_resistivity()

    def onChanged_Temp(self, text):
        self.temperature = text

    def onChanged_E(self, text):
        self.E = text

    def onChanged_Nmin(self, text):
        self.Nmin = text

    def onChanged_Nmax(self, text):
        self.Nmax = text


    def onChanged_count(self, text):
        self.count = text




if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    myapp = MyWin()
    myapp.show()
    sys.exit(app.exec_())
