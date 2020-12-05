from PyQt5 import QtCore, QtGui, QtWidgets
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1300, 650)
        MainWindow.setUnifiedTitleAndToolBarOnMac(False)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(20, 20, 60, 16))
        self.label.setObjectName("label")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(390, 10, 900, 550))
        self.tabWidget.setStyleSheet("background-color: rgb(254, 181, 80);\n"
"background-color: rgb(177, 198, 222);")
        self.tabWidget.setObjectName("tabWidget")
        self.tab = pg.PlotWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        self.tab.setBackground("f0f0f0")
        self.tab_4 = pg.PlotWidget()
        self.tab_4.setObjectName("tab_4")
        self.tabWidget.addTab(self.tab_4, "")
        self.tab_4.setBackground("f0f0f0")
        self.tab_2 = pg.PlotWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.tab_2.setBackground("f0f0f0")
        self.tab_3 = pg.PlotWidget()
        self.tab_3.setObjectName("tab_3")
        self.tabWidget.addTab(self.tab_3, "")
        self.tab_3.setBackground("f0f0f0")
        self.tab_5 = pg.PlotWidget()
        self.tab_5.setObjectName("tab_5")
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_5.setBackground("f0f0f0")
        self.tab_6 = pg.PlotWidget()
        self.tab_6.setObjectName("tab_6")
        self.tabWidget.addTab(self.tab_6, "")
        self.tab_6.setBackground("f0f0f0")
        self.tab_7 = pg.PlotWidget()
        self.tab_7.setObjectName("tab_6")
        self.tabWidget.addTab(self.tab_7, "")
        self.tab_7.setBackground("f0f0f0")
        self.ComboBox = QtWidgets.QComboBox(self.centralwidget)
        self.ComboBox.setGeometry(QtCore.QRect(145, 15, 131, 31))
        self.ComboBox.setObjectName("ComboBox")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(20, 70, 101, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(20, 120, 121, 16))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(20, 170, 71, 16))
        self.label_4.setObjectName("label_4")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(280, 70, 21, 16))
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(280, 120, 91, 16))
        self.label_10.setObjectName("label_10")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(280, 170, 91, 16))
        self.label_11.setObjectName("label_11")

        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(280, 370, 91, 16))
        self.label_12.setObjectName("label_12")
        self.lineEdit_4 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_4.setGeometry(QtCore.QRect(150, 70, 121, 21))
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.lineEdit_5 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_5.setGeometry(QtCore.QRect(150, 120, 121, 21))
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.lineEdit_6 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_6.setGeometry(QtCore.QRect(150, 220, 121, 21))
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.lineEdit_7 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_7.setGeometry(QtCore.QRect(150, 270, 121, 21))
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.lineEdit_8 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_8.setGeometry(QtCore.QRect(150, 320, 121, 21))
        self.lineEdit_8.setObjectName("lineEdit_8")
        self.lineEdit_3 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_3.setGeometry(QtCore.QRect(150, 170, 121, 21))
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.lineEdit_9 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_9.setGeometry(QtCore.QRect(150, 420, 121, 21))
        self.lineEdit_9.setObjectName("lineEdit_9")
        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(20, 270, 71, 16))
        self.label_21.setObjectName("label_21")
        self.lineEdit_10 = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit_10.setGeometry(QtCore.QRect(150, 370, 121, 21))
        self.lineEdit_10.setObjectName("lineEdit_10")
        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        self.label_24.setGeometry(QtCore.QRect(20, 220, 121, 21))
        self.label_24.setObjectName("label_24")
        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(280, 220, 91, 16))
        self.label_14.setObjectName("label_14")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(280, 270, 91, 16))
        self.label_15.setObjectName("label_15")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(20, 320, 71, 16))
        self.label_5.setObjectName("label_5")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(20, 370, 71, 16))
        self.label_7.setObjectName("label_7")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(280, 320, 91, 16))
        self.label_16.setObjectName("label_16")
        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        self.label_25.setGeometry(QtCore.QRect(20, 420, 91, 16))
        self.label_25.setObjectName("label_25")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(280, 420, 41, 16))
        self.label_13.setObjectName("label_13")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(150, 470, 121, 61))
        self.pushButton.setObjectName("pushButton")
        self.checkBox = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox.setGeometry(QtCore.QRect(280, 470, 81, 31))
        self.checkBox.setObjectName("checkBox")
        self.checkBox_2 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_2.setGeometry(QtCore.QRect(280, 490, 81, 31))
        self.checkBox_2.setObjectName("checkBox")
        self.checkBox_3 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_3.setGeometry(QtCore.QRect(280, 510, 81, 31))
        self.checkBox_3.setObjectName("checkBox")
        self.checkBox_4 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_4.setGeometry(QtCore.QRect(280, 530, 81, 31))
        self.checkBox_4.setObjectName("checkBox")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1183, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(20, 470, 121, 61))
        self.pushButton_2.setObjectName("pushButton")

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "Material"))
       # self.tabWidget.setStyleSheet(background-color: white')
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Mobility"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), _translate("MainWindow", " e Concentration"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "hole Concentration"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Negative acceptor concentration"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), _translate("MainWindow", "Positive acceptor concentration"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), _translate("MainWindow", "Conductivity"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_7), _translate("MainWindow", "Resistivity"))
        self.ComboBox.setCurrentText(_translate("MainWindow", "choose material"))
        materials = ["Si","GaAs","Ge"]
        self.ComboBox.addItems(materials)
        self.label_2.setText(_translate("MainWindow", "Donor\'s energy"))
        self.lineEdit_4.setText("0.045")
        self.lineEdit_5.setText("0.1")
        self.lineEdit_3.setText("10")
        self.lineEdit_6.setText("0.045")
        self.lineEdit_7.setText("0.1")
        self.lineEdit_8.setText("10")
        self.lineEdit_10.setText("100")
        self.lineEdit_9.setText("100")

        self.label_3.setText(_translate("MainWindow", "min Nd"))
        self.label_4.setText(_translate("MainWindow", "max Nd"))
        self.label_9.setText(_translate("MainWindow", "eV"))
        self.label_10.setText(_translate("MainWindow", "*10ˆ18 cmˆ-3"))
        self.label_11.setText(_translate("MainWindow", "*10ˆ18 cmˆ-3"))
        self.label_12.setText(_translate("MainWindow", "pts"))
        self.label_21.setText(_translate("MainWindow", "min Nd"))
        self.label_24.setText(_translate("MainWindow", "Acceptor's energy"))
        self.label_14.setText(_translate("MainWindow", "eV"))
        self.label_15.setText(_translate("MainWindow", "*10ˆ18 cmˆ-3"))
        self.label_5.setText(_translate("MainWindow", "max Na"))
        self.label_7.setText(_translate("MainWindow", "Count steps"))
        self.label_16.setText(_translate("MainWindow", "*10ˆ18 cmˆ-3"))
        self.label_25.setText(_translate("MainWindow", "Temperature"))
        self.label_13.setText(_translate("MainWindow", "K"))
        self.pushButton.setText(_translate("MainWindow", "DRAW PLOT"))
        self.checkBox.setText(_translate("MainWindow", "LogAxis X"))
        self.checkBox_2.setText(_translate("MainWindow", "LogAxis Y"))
        self.checkBox_3.setText(_translate("MainWindow", "n"))
        self.checkBox_3.setChecked(1)
        self.checkBox_4.setText(_translate("MainWindow", "p"))
        self.pushButton_2.setText(_translate("MainWindow", "SAVE"))