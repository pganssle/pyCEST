"""
Main GUI handler
"""

from __future__ import unicode_literals
import sys
import os

from matplotlib.backends import qt_compat
use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE
if use_pyside:
    from PySide import QtGui, QtCore
else:
    from PyQt4 import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import pyCEST
progname = 'pyCEST ({})'.format(pyCEST.__version__)

class MainCanvas(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        self.axes.hold(False)

        super(MainCanvas, self).__init__(fig)

        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(progname)

        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction('&Load', self.fileLoad,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_L)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtGui.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QVBoxLayout(self.main_widget)
        sc = MainCanvas(self.main_widget, width=10, height=8, dpi=100)
        l.addWidget(sc)

        self.canvas = sc

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)


    def fileLoad(self):
        dataPath = QtGui.QFileDialog.getExistingDirectory(self,
            "Select a Bruker folder")

        if dataPath:
            load_bruker(self.canvas, dataPath)

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtGui.QMessageBox.about(self, 'About',
            progname)

from pyCEST.old import nylib
from matplotlib import cm

from . import roi

def load_bruker(canvas, dataPath):
    field = float(nylib.BrukerPar(dataPath, 'acqp', 'BF1='))

    data = nylib.Paravision2dseqNew(dataPath)
    freq = nylib.BrukerPar(dataPath, 'method', 'MT_Offsets_NoM0=') / field

    ax = canvas.axes

    ax.imshow(data[0, :, :], cmap=cm.gray)
    ROIH = roi.ROIHandler(ax)
    ROIH.connect()

    canvas.draw()

def make_gui(argv):
    qApp = QtGui.QApplication(argv)
    aw = ApplicationWindow()
    aw.setWindowTitle(progname)

    aw.show()
    sys.exit(qApp.exec_())