import sys

import numpy as np

from pyCEST.gui import gui
from pyCEST.gui import roi
from pyCEST.old import nylib

from matplotlib import pyplot as plt
from matplotlib import cm
import skimage as ski

if __name__ == "__main__":
    '''
    # Load the data
    dataPath = 'Data/62'
    field = float(nylib.BrukerPar(dataPath, 'acqp', 'BF1='))

    data = nylib.Paravision2dseqNew(dataPath)
    freq = nylib.BrukerPar(dataPath, 'method', 'MT_Offsets_NoM0=') / field

    plt.figure()
    plt.imshow(data[0, :, :], cmap=cm.gray)
    ax = plt.gca()

    ROIH = roi.ROIHandler(ax)
    ROIH.connect()

    plt.show()
    '''

    gui.make_gui(sys.argv)