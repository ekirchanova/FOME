#!/usr/bin/env python3

import time, random
import math
from collections import deque

start = time.time()


class RealtimePlot:
    def __init__(self, axes, max_entries=100):
        self.axis_x = deque(maxlen=max_entries)
        self.axis_y = deque(maxlen=max_entries)
        self.axes = axes
        self.max_entries = max_entries

        self.lineplot, = axes.plot([], [], "ro-")
        self.axes.set_autoscaley_on(True)

    def add(self, x, y):
        self.axis_x.append(x)
        self.axis_y.append(y)
        self.lineplot.set_data(self.axis_x, self.axis_y)
        self.axes.set_xlim(self.axis_x[0], self.axis_x[-1] + 1e-15)
        self.axes.relim();
        self.axes.autoscale_view()  # rescale the y-axis

    def animate(self, figure, callback, interval=50):
        import matplotlib.animation as animation
        def wrapper(frame_index):
            self.add(*callback(frame_index))
            self.axes.relim();
            self.axes.autoscale_view()  # rescale the y-axis
            return self.lineplot

        animation.FuncAnimation(figure, wrapper, interval=interval)


def main():
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots()
    display = RealtimePlot(axes)
    display.animate(fig, lambda frame_index: (time.time() - start, random.random() * 100))
    plt.show()

    fig, axes = plt.subplots()
    display = RealtimePlot(axes)
    while True:
        display.add(time.time() - start, random.random() * 100)
        plt.pause(0.001)


if __name__ == "__main__": main()



# # embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases
# #
# # Copyright (C) 2005 Florent Rougon
# #               2006 Darren Dale
# #               2015 Jens H Nielsen
# #
# # This file is an example program for matplotlib. It may be used and
# # modified with no restriction; raw copies as well as modified versions
# # may be distributed without limitation.
#
# from __future__ import unicode_literals
# import sys
# import os
# import random
# import matplotlib
# # Make sure that we are using QT5
# matplotlib.use('Qt5Agg')
# from PyQt5 import QtCore, QtWidgets
#
# from numpy import arange, sin, pi
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.figure import Figure
#
# progname = os.path.basename(sys.argv[0])
# progversion = "0.1"
#
#
# class MyMplCanvas(FigureCanvas):
#     """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
#
#     def __init__(self, parent=None, width=5, height=4, dpi=100):
#         fig = Figure(figsize=(width, height), dpi=dpi)
#         self.axes = fig.add_subplot(111)
#
#         self.compute_initial_figure()
#
#         FigureCanvas.__init__(self, fig)
#         self.setParent(parent)
#
#         FigureCanvas.setSizePolicy(self,
#                                    QtWidgets.QSizePolicy.Expanding,
#                                    QtWidgets.QSizePolicy.Expanding)
#         FigureCanvas.updateGeometry(self)
#
#     def compute_initial_figure(self):
#         pass
#
#
# class MyStaticMplCanvas(MyMplCanvas):
#     """Simple canvas with a sine plot."""
#
#     def compute_initial_figure(self):
#         t = arange(0.0, 3.0, 0.01)
#         s = sin(2*pi*t)
#         self.axes.plot(t, s)
#
#
# class MyDynamicMplCanvas(MyMplCanvas):
#     """A canvas that updates itself every second with a new plot."""
#
#     def __init__(self, *args, **kwargs):
#         MyMplCanvas.__init__(self, *args, **kwargs)
#         timer = QtCore.QTimer(self)
#         timer.timeout.connect(self.update_figure)
#         timer.start(1000)
#
#     def compute_initial_figure(self):
#         self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')
#
#     def update_figure(self):
#         # Build a list of 4 random integers between 0 and 10 (both inclusive)
#         l = [random.randint(0, 10) for i in range(4)]
#         self.axes.cla()
#         self.axes.plot([0, 1, 2, 3], l, 'r')
#         self.draw()
#
#
# class ApplicationWindow(QtWidgets.QMainWindow):
#     def __init__(self):
#         QtWidgets.QMainWindow.__init__(self)
#         self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
#         self.setWindowTitle("application main window")
#
#         self.file_menu = QtWidgets.QMenu('&File', self)
#         self.file_menu.addAction('&Quit', self.fileQuit,
#                                  QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
#         self.menuBar().addMenu(self.file_menu)
#
#         self.help_menu = QtWidgets.QMenu('&Help', self)
#         self.menuBar().addSeparator()
#         self.menuBar().addMenu(self.help_menu)
#
#         self.help_menu.addAction('&About', self.about)
#
#         self.main_widget = QtWidgets.QWidget(self)
#
#         l = QtWidgets.QVBoxLayout(self.main_widget)
#         sc = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
#         dc = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
#         l.addWidget(sc)
#         l.addWidget(dc)
#
#         self.main_widget.setFocus()
#         self.setCentralWidget(self.main_widget)
#
#         self.statusBar().showMessage("All hail matplotlib!", 2000)
#
#     def fileQuit(self):
#         self.close()
#
#     def closeEvent(self, ce):
#         self.fileQuit()
#
#     def about(self):
#         QtWidgets.QMessageBox.about(self, "About",
#                                     """embedding_in_qt5.py example
# Copyright 2005 Florent Rougon, 2006 Darren Dale, 2015 Jens H Nielsen
#
# This program is a simple example of a Qt5 application embedding matplotlib
# canvases.
#
# It may be used and modified with no restriction; raw copies as well as
# modified versions may be distributed without limitation.
#
# This is modified from the embedding in qt4 example to show the difference
# between qt4 and qt5"""
#                                 )
#
#
# qApp = QtWidgets.QApplication(sys.argv)
#
# aw = ApplicationWindow()
# aw.setWindowTitle("%s" % progname)
# aw.show()
# sys.exit(qApp.exec_())
# #qApp.exec_()




import numpy as np
import matplotlib.pyplot as plt

from engine.Physics.StateCalculator import StateCalculator

# if __name__ == '__main__':
#     # app = QtWidgets.QApplication(sys.argv)
#     # window = Ui()
#     # sys.exit(app.exec_())
#
#     calc = StateCalculator()
#     calc.set_temperature(300)
#     COUNTS = 100
#     Nas = np.linspace(1e14, 1e19, COUNTS, endpoint=True)
#     Ea = 0.045
#     Nds = np.zeros(COUNTS)
#     Ed = 0.045
#     calc.set_additions(Nas, Ea, Nds, Ed)
#     calc.update()
#
#     # np.set_printoptions(precision=5, suppress=True)
#     # print("Fermi level", calc.get_fermi_levels())
#     # print("Electron concentration", calc.get_electron_concentration())
#     # print("Hole concentration", calc.get_hole_concentration())
#     # print("Positive acceptor concentration", calc.get_pos_donor_concentration())
#     # print("Negative acceptor concentration", calc.get_neg_acceptor_concentration())
#     # print("Electron mobility", calc.get_electron_mobility())
#     # print("Hole mobility", calc.get_hole_mobility())
#     # print("Conductivity", calc.get_conductivity())
#     # print("Resistivity", calc.get_resistivity())
#
#     x = Nas
#     y = calc.get_resistivity()
#
#     fig, ax1 = plt.subplots(
#         nrows=1,
#         ncols=1,
#         figsize=(16, 12)
#     )
#
#     ax1.plot(x, y, '-k')
#
#     ax1.set(xlabel='x', ylabel='P(x)', title='Integral solution')
#     ax1.grid(which='major', axis='x')
#     ax1.set_xscale('log')
#     ax1.set_yscale('log')
#
#     plt.show()
