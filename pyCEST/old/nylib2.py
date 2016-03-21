import cjlib
from PIL import Image, ImageDraw
import numpy as np
from matplotlib.pyplot import *

def roipolyny(data, XY):
    img = Image.new('L', (data.shape[-1], data.shape[-2]), 0)
    ImageDraw.Draw(img).polygon(XY, outline=1, fill=1)
    mask = np.array(img)
    return mask

def applyMask(data, dataMask):
    dataLength = data.shape[0]
    Ints = np.zeros((dataLength))
    for ii in range(dataLength):
        intsTemp = np.ma.array(data[ii], mask=(dataMask==False))
        Ints[ii] = intsTemp.mean()
    return Ints

def getRoi(data, inputTime = -1):
    figure(1)
    clf()
    figure(1).subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0)
    cjlib.mimage( data[0] )
    
    XY = ginput(inputTime)
    mask = roipolyny(data, XY)
    print(XY)
    
    # plot the roi
    XY=np.array(XY)
    XY = np.concatenate ((XY, [XY[0,:]]), axis=0 )
    plot(XY.transpose()[0], XY.transpose()[1], 'r-*', linewidth = 4)
    figure(1).canvas.draw()
    
    return XY, mask

def sector_mask(shape,centre,radius,angle_range):
    """
    taken from http://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array
    Return a boolean mask for a circular sector. The start/stop angles in  
    `angle_range` should be given in clockwise order.
    """

    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
    tmin,tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
            tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
    theta = np.arctan2(x-cx,y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)

    # circular mask
    circmask = r2 <= radius*radius

    # angular mask
    anglemask = theta <= (tmax-tmin)

    return circmask*anglemask
