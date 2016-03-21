from __future__ import division
from numpy import *
import numpy as np
from scipy.optimize import fmin_l_bfgs_b, leastsq, fmin_tnc, fminbound, fmin_bfgs, fmin_powell, fmin, fminbound, fmin_slsqp, fmin_cg, nnls , curve_fit
from scipy import interpolate, ogrid, mgrid, ndimage, rand
#from scipy.stsci.convolve import convolve2d
from numpy import array, fromfile, frombuffer, concatenate, zeros, arange, mean, prod, tile, mgrid, exp, polyfit, log, arccos, squeeze, uint16, vectorize, double, memmap, linspace, nonzero, unique, conj, pi, arctan, imag, real, diff, sum, dot, isinf, sqrt, inf, logical_and, roll, polyfit, polyval, median, interp, isnan, cos, sin
from numpy.linalg import pinv, svd
#import numpy.numarray.nd_image
import scipy.ndimage.measurements
import scipy.ndimage.morphology
from matplotlib.pyplot import imshow, clim, gray, colorbar, ginput, axis, gca, axes, figure, clf, plot, legend, xticks, yticks, gcf, title, xlabel, ylabel, grid
from matplotlib.mlab import prctile
import matplotlib.colors 
from progressbar import ProgressBar, Percentage, Bar, ETA
import math 
import struct
import re
import traceback
from matplotlib.widgets import Slider, Button, RadioButtons
import time
import glob
import gzip
import gc
import sys, traceback
try:
    import cPickle as pickle
except ImportError:
    import pickle
#import nifti

##==============================================================================
##
##  General Things
##
##==============================================================================

def exist(var):
    print("exist BROKEN")
    return globals().has_key(var) | vars().has_key(var) | locals().has_key(var)

##==============================================================================
##
##  Masking and Regions
##
##==============================================================================
def applymask(d, mask):

    if d.shape[-1] != mask.shape[-1] or d.shape[-2] != mask.shape[-2]:
        return None
    
    # Determine the size of the mask
    nrows = d.shape[-1]
    ncols = d.shape[-2]
    nslices = d.shape[-3]

    # Compute the mean
    rows = arange(nrows)
    cols = arange(ncols)

    #print(rows)
    #print(cols)

    s = 0
    #print([ d[s,c,r] for r in rows for c in cols if mask[c,r] ])

    mm = array([ mean( [ d[s,c,r] for r in rows for c in cols if mask[c,r] ] ) for s in arange(nslices) ])
    ss = array([ std( [ d[s,c,r] for r in rows for c in cols if mask[c,r] ] ) for s in arange(nslices) ])

    return r_[mm], r_[ss]

def roipoly():
    XY = ginput(-1)

    blah1, blah2, ncols, nrows = gca().dataLim.bounds

    mask = array([ pnpoly(i,j,XY) for i in range(nrows) for j in range(ncols)  ]).reshape(ncols,nrows)

    return mask

def roipolyXY():
    XY = ginput(-1)

    blah1, blah2, ncols, nrows = gca().dataLim.bounds
    
    nrows, ncols = int(nrows), int(ncols)
    
    mask = array([ pnpoly(j,i,XY) for i in range(nrows) for j in range(ncols)  ]).reshape(ncols,nrows)

    return XY,mask

def createmask(d, percentileThreshold = 75 ):
    mask = d > prctile(d, percentileThreshold )
    mask2 = scipy.ndimage.morphology.binary_fill_holes( scipy.ndimage.morphology.binary_erosion( mask, iterations=3 ) )
    L = scipy.ndimage.measurements.label( mask2 )[0]
    hist = np.histogram( L, np.max(L) )[0]
    second_largest = nonzero( hist == np.sort(hist)[-2] )[0]
    mask = L == second_largest

    return mask


##==============================================================================
##
##  Reading Files 
##
##==============================================================================

#  Get the diffusion gradients from the PAR file
def diffgrads(filename, outputFilename=None):

    doublev = np.vectorize(np.double)

    if filename.find('.rec') > -1:
        filename = filename.replace('.rec', '.par')
    elif filename.find('.REC') > -1:
        filename = filename.replace('.REC', '.PAR')
    
    # Open and read in the lines in the PAR file
    try:
        fp = open(filename)
        lines = fp.readlines()
        fp.close()
    except IOError:
        print("diffgrads - Error: Could not read in " + filename)
        return None

    #  Get all the lines at the bottom of the file
    blah = [ aa for aa in lines if re.search("([-.\d]+\s+){4,}", aa)]

    # Convert everything to numbers
    nums = array( [doublev( re.findall('[-.e\d]+', l) ) for l in blah ] )

    # Find the appropriate columns where the first number is a 1
    grads = array( [ A[-4:-1] for A in nums if A[0] == 1 ]  )

    # If there is an output filename, then output the gradients to a file.
    if not outputFilename == None:
        fp = open(outputFilename, 'wt')
        for ii,gg in enumerate( grads ):
            fp.write( "%d: %f, %f, %f\n" % (ii, gg[0], gg[1], gg[2]) )    
        fp.close()

    return grads

#def rec2nii(filename, filename_out, scale=1):
#    d = readrec(filename, outfmt='float32')
#
#    filename_par = filename.replace('REC', 'PAR')
#    filename_par = filename_par.replace('rec', 'par')
#    fp = open(filename_par)
#    par = fp.readlines()
#    fp.close()
#        
#    re_number = re.compile("[\d]+")
#    re_spaces = re.compile("\s+")
#    for ii,l in enumerate( par ):
#        numbers = re_spaces.split(l.lstrip().rstrip())        
#        if re_number.match(numbers[0]) and re_number.match(numbers[1]) and len(numbers) > 10:
#            slthk = float(numbers[22])+float(numbers[23]) 
#            xsize = float(numbers[28]) 
#            ysize = float(numbers[29])
#            break;
#
#    tnii = nifti.NiftiImage(d*scale )
#    tnii.setPixDims((xsize, ysize, slthk))
#    tnii.setFilename(filename_out)
#    tnii.save()
#
#def rec2niiSP(filename, filename_out, a, b, scale=1):
#    d = readrec(filename, outfmt='float32')
#
#    filename_par = filename.replace('REC', 'PAR')
#    filename_par = filename_par.replace('rec', 'par')
#    fp = open(filename_par)
#    par = fp.readlines()
#    fp.close()
#        
#    re_number = re.compile("[\d]+")
#    re_spaces = re.compile("\s+")
#    for ii,l in enumerate( par ):
#        numbers = re_spaces.split(l.lstrip().rstrip())        
#        if re_number.match(numbers[0]) and re_number.match(numbers[1]) and len(numbers) > 10:
#            slthk = float(numbers[22])+float(numbers[23]) 
#            xsize = float(numbers[28]) 
#            ysize = float(numbers[29])
#            break;
#
#    tnii = nifti.NiftiImage(d[a,b]*scale )
#    tnii.setPixDims((xsize, ysize, slthk))
#    tnii.setFilename(filename_out)
#    tnii.save()

def getTE(filename):
    
    # Set the correct filenames
    if filename.find('.par') > -1:
        filename_par = filename
        filename_rec = filename.replace('.par', '.rec')
    elif filename.find('.PAR') > -1:
        filename_par = filename
        filename_rec = filename.replace('.PAR', '.REC')
    elif filename.find('.rec') > -1:
        filename_rec = filename
        filename_par = filename.replace('.rec', '.par')
    elif filename.find('.REC') > -1:
        filename_rec = filename
        filename_par = filename.replace('.REC', '.PAR')
    else:
        filename_par = filename + '.PAR'
        filename_rec = filename + '.REC'

    # Read in the PAR
    if filename_par.find('.gz') > -1:
        fp = gzip.GzipFile(filename_par, 'rt')
    else:
        fp = open(filename_par, 'rt')
    par = fp.readlines()
    fp.close()

    re_number = re.compile("[\d]+")
    #re_float_number = re.compile("[-]\d[-.e\d]*")
    re_float_number = re.compile("[-+]\d[\d*][.][\d*][eE][-+][d\*]")
    re_spaces = re.compile("\s+")

    # Get the reconstruction resolution
    numbers = re_number.findall( par[-3] )
    xsize,ysize = int(numbers[9]), int(numbers[10])

    # Get the RI, RS and SS to rescale the images into the 
    # proper floating point values
    # DV = PV * RS + RI   FP = DV / (RS * SS)
    # Convert the pixel value to floating point value
    te = []
    for ii,l in enumerate( par ):
        numbers = re_spaces.split(l.lstrip().rstrip())
        if re_number.match(numbers[0]) and re_number.match(numbers[1]) and len(numbers) > 10:
            te.append( float( numbers[30] ) )

    te = array(te)

    # Get the number of phases
    expr = re.compile("Max. number of echoes")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    nechoes = int(numbers[0])

    return unique(te)

def readrec(filename, byteswap=False, sq=True, outfmt='float32', magnitudeOnly=False):
    
    # Set the correct filenames
    if filename.find('.par') > -1:
        filename_par = filename
        filename_rec = filename.replace('.par', '.rec')
    elif filename.find('.PAR') > -1:
        filename_par = filename
        filename_rec = filename.replace('.PAR', '.REC')
    elif filename.find('.rec') > -1:
        filename_rec = filename
        filename_par = filename.replace('.rec', '.par')
    elif filename.find('.REC') > -1:
        filename_rec = filename
        filename_par = filename.replace('.REC', '.PAR')
    else:
        filename_par = filename + '.PAR'
        filename_rec = filename + '.REC'

    # Read in the PAR
    if filename_par.find('.gz') > -1:
        fp = gzip.GzipFile(filename_par, 'rt')
    else:
        fp = open(filename_par, 'rt')
    par = fp.readlines()
    fp.close()

    re_number = re.compile("[\d]+")
    #re_float_number = re.compile("[-]\d[-.e\d]*")
    re_float_number = re.compile("[-+]\d[\d*][.][\d*][eE][-+][d\*]")
    re_spaces = re.compile("\s+")

    # Get the reconstruction resolution
    numbers = re_number.findall( par[-3] )
    xsize,ysize = int(numbers[9]), int(numbers[10])

    # Get the RI, RS and SS to rescale the images into the 
    # proper floating point values
    # DV = PV * RS + RI   FP = DV / (RS * SS)
    # Convert the pixel value to floating point value
    ri,rs,ss = [], [], []
    slices = []
    echoes = []
    dynamics = []
    phases = []
    types = []
    gorients=[]
    for ii,l in enumerate( par ):
        numbers = re_spaces.split(l.lstrip().rstrip())
        if re_number.match(numbers[0]) and re_number.match(numbers[1]) and len(numbers) > 10:
            ri.append( float( numbers[11] ) )
            rs.append( float( numbers[12] ) )
            ss.append( float( numbers[13] ) )
            slices.append( int( numbers[0] ) )
            echoes.append( int( numbers[1] ) )
            dynamics.append( int( numbers[2] ) )
            phases.append( int( numbers[3] ) )
            types.append( int( numbers[4] ) )
            gorients.append( [ numbers[-4:-1] ] )

    ri, rs, ss = array(ri), array(rs), array(ss)
    types = array(types)
    echoes = array(echoes)
    dynamics = array(dynamics)
    slices = array(slices)
    phases = array(phases)

    types_unique = np.unique( types )
    ntypes = len(types_unique)

    # Get the number of slices
    expr = re.compile("Max. number of slices")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    zsize = int(numbers[0])

    # Get the number of phases
    expr = re.compile("Max. number of cardiac")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    nphases = int(numbers[0])

    # Get the number of phases
    expr = re.compile("Max. number of dynamics")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    ndynamics = int(numbers[0])

    # Get the number of phases
    expr = re.compile("Max. number of echoes")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    nechoes = int(numbers[0])

    # Get the number of diffusion gradients
    expr = re.compile("Max. number of gradient orients")
    numbers = re_number.findall( filter(expr.search, par)[0] )
    norients = int(numbers[0])
    if ( norients > 1):
        norients += 1

    #  Read in the file
    if filename_rec.endswith('gz'):
        fp = gzip.GzipFile(filename_rec, 'rb')
        print('cjlib::readrec: GZIP NOT WORKING YET')
    else:
        fp = open(filename_rec, 'rb')

    #d = zeros( (norients, ntypes, ndynamics,nechoes,zsize, xsize*ysize), dtype=outfmt )
    d = zeros( (ntypes, norients, ndynamics,nechoes, zsize,xsize*ysize), dtype=outfmt )

    gorient_prev = [-1, -1, -1]
    if outfmt == 'float32':
        gorienti= -1
        for ii in range(ndynamics*norients*nechoes*ntypes*zsize):
            if not gorient_prev == gorients[ii]:
                gorienti = gorienti + 1
                gorient_prev = gorients[ii]
            ti = nonzero( types_unique == types[ii] )[0]
            d[ti,gorienti, dynamics[ii]-1,echoes[ii]-1,slices[ii]-1] = array( (np.fromfile(fp, dtype=np.dtype('uint16'), count=xsize*ysize) * rs[ii] + ri[ii] ) / (rs[ii] * ss[ii] ), outfmt)
    else:
        for ii in range(ndynamics*norients*nechoes*ntypes*zsize):
            if not gorient_prev == gorients[ii]:
                gorienti = gorienti + 1
                gorient_prev = gorients[ii]
            ti = nonzero( types_unique == types[ii] )[0]
            d[ti,gorienti,dynamics[ii]-1, echoes[ii]-1, slices[ii]-1] = np.fromfile(fp, dtype='uint16', count=xsize*ysize)

    fp.close()

    print(d.shape)

    d = d.reshape( (ntypes, norients, ndynamics, nechoes, zsize, ysize, xsize ) )

    if magnitudeOnly:
        d = d[0]

    return squeeze( d )

def readfa(filename, shape):

    fa = readraw(filename, shape, intype='float32')

    return fa

def readvec(filename, shape):

    if len(shape) == 3:
        shape.append(3)

    vec = readraw(filename, shape, intype='float32')

    return vec

def readdec(filename, shape):
    return readcmap(filename, shape )

def readcmap(filename, shape):

    if len(shape) == 3:
        shape.append(3)

    vec = readraw(filename, shape, intype='uint8')

    return vec


def readrawMEM(filename, shape, intype='int16', byteSwap=False):
    """ readraw - To read in a raw file and reformat it to the right shape """

    #  Read in the file
    return memmap(filename, dtype=intype, mode='r', shape=shape)

def readraw(filename, shape, intype='int16', byteSwap=False):
    """ readraw - To read in a raw file and reformat it to the right shape """

    #  Read in the file
    if filename.endswith('gz'):
        fp = gzip.open(filename, 'rb')
    else:
        fp = open(filename, 'rb')

    d = fromfile(file=fp, dtype=intype, count=shape.prod()).reshape(shape)

    d.byteswap(byteSwap)

    return d


def readanalyze(filename):

    # Set the corimgt filenames
    if filename.endswith('hdr'):
        filename_hdr = filename
        filename_img = filename.replace('hdr', 'img')
    elif filename.endswith('HDR'):
        filename_hdr = filename
        filename_img = filename.replace('HDR', 'IMG')
    elif filename.endswith('img'):
        filename_img = filename
        filename_hdr = filename.replace('img', 'hdr')
    elif filename.endswith('IMG'):
        filename_img = filename
        filename_hdr = filename.replace('IMG', 'HDR')
    else:
        filename_hdr = filename + '.HDR'
        filename_img = filename + '.IMG'

    # Read in the header
    fp = open(filename_hdr, 'rb')
    line = fp.read(40)

    endian = ">"

    sizeof_hdr,data_type,db_name,extents,session_error,regular,hkey_un0 = struct.unpack(">L10s18sLHss", line)

    if sizeof_hdr != 348:
        endian = "<"
        sizeof_hdr,data_type,db_name,extents,session_error,regular,hkey_un0 = struct.unpack("%sL10s18sLHss" % endian, line)

    line = fp.read(68)
    dim = zeros(8)
    pixdim = zeros(8)
    dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], unused8, unused9, unused10, unused11, unused12, unused13, unused14, datatype, bitpix, dim_un0, pixdim[0], pixdim[1], pixdim[2], pixdim[3], pixdim[4], pixdim[5], pixdim[6], pixdim[7] = struct.unpack("%s8h7hhhh8f"%endian, line)

    line = fp.read(36)
    vox_offset, funused1, funused2, funused3, cal_max, cal_min, compressed, verified, glmax, glmin = struct.unpack("%s8f2h"%endian, line)

    line = fp.read(200)
    descrip, aux_file, orient, originator, generated, scannum, patient_id, exp_date, exp_time, hist_un0, views, vols_added, start_field, field_skip, omax, omin, smax, smin  = struct.unpack("%s80s24ss10s10s10s10s10s10s3s8L"%endian, line)

    fp.close()

    # read in the data
    byteswap = False
    if endian == ">":
        byteswap = True

    if datatype == 2:
        type = 'uint8'
    elif datatype == 4:
        type = 'int16'
    elif datatype == 8:
        type = 'int32'
    elif datatype == 16:
        type = 'float32'
    elif datatype == 32:
        type = 'complex64'
    elif datatype == 64:
        type = 'float64'
    else:
        print("Unknown type " + datatype)
        return

    d = readraw(filename_img, (dim[1:dim[0]+1][::-1]), byteSwap=byteswap, intype=type )

    return d.squeeze()

def parse_procpar(lines, val):
    tomatch = re.compile(val)
    for ii, l in enumerate(lines):
        if tomatch.match(l):
            return lines[ii+1].split(" ")[1]

def readvimg(filename):

    fdfs = glob.glob(filename+"/*.fdf")

    for ii,f in enumerate(fdfs):
        d = readfdf( f )
        if ii == 0:
            data = zeros( ( len(fdfs), d.shape[1], d.shape[0] ) )
        data[ii] = d

    return data    

def readfdf(filename):
    
    fp = open( filename, 'rb' )

    xsize = -1
    ysize = -1
    zsize = 1
    bigendian = -1
    done = False

    while not done :

        line = fp.readline()    

        if( len( line ) >= 1 and line[0] == chr(12) ):
            break

        if( len( line ) >= 1 and line[0] != chr(12) ):

            if( line.find('bigendian') > 0 ):
                endian = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('echos') > 0 ):
                nechoes = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('echo_no') > 0 ):
                echo_no = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('nslices') > 0 ):
                nslices = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('slice_no') > 0 ):
                sl = line.split('=')[-1].rstrip('\n; ').strip(' ')

            if( line.find('matrix') > 0 ):
                m = re.findall('(\d+)', line.rstrip())
                
                if len(m) == 2:
                    xsize, ysize = int(m[0]), int(m[1])
                elif len(m) == 3:
                    xsize, ysize, zsize = int(m[0]), int(m[1]), int(m[2])

    fp.seek(-xsize*ysize*zsize*4,2)

    if bigendian == 1:
        fmt = ">%df" % (xsize*ysize*zsize)
    else:
        fmt = "<%df" % (xsize*ysize*zsize)

    data = struct.unpack(fmt, fp.read(xsize*ysize*zsize*4))
    data = array( data ).reshape( [xsize, ysize, zsize ] ).squeeze()

    fp.close()

    return data

##==============================================================================
##
##  Show Images
##
##==============================================================================

def vimage(d):
    """ Very slow at this point """

    mimage(d[:,:,1])
    main_ax = gca()
    
    axslice  = axes([0.1, 0.1, 0.8, 0.05])
    axis('off')
    slice_slider = Slider(axslice, 'Slice', 1, d.shape[2], valinit=1)

    def update(val):
        axes(main_ax)
        mimage(d[:,:,int(slice_slider.val)])

    slice_slider.on_changed(update)


def mimage(d, cmap=gray):
    imshow( d )
    axis('image') # needed so that ginput doesn't resize the image
    clim([ prctile(d,1) , prctile(d, 99) ])
    xticks([])
    yticks([])
#    gca().get_axes().set_position([0,0,1,1]) #commented out by ny temporarily

def mimagecb(d, cmap=None):
    mimage(d, cmap)
    ax = gcf().add_axes([0.85, 0.05, 0.04, 0.9])
    colorbar(cax=ax, orientation='vertical')

def mmontage( d, cmap=gray ):

    slices, rows, cols = d.shape

    N = round( math.sqrt( slices ) )

    im_cols = N
    im_rows = N
    
    if im_rows*im_cols < slices:
        im_rows = im_rows + 1

    extra = im_cols * im_rows - slices

    ii = 0
    d2 = zeros( (im_rows*rows, im_cols*cols), dtype=d.dtype)

    s1 = time.time()
    for ri in arange(0, im_rows):
        for ci in arange(0, im_cols):
            if ii == slices: break
            d2[ri*rows:(ri+1)*rows,ci*cols:(ci+1)*cols] = d[ii,:,:]
            ii = ii + 1

    mimage(d2, cmap)

def mmontageNEW( d, rc=None, cmap=None ):

    slices, rows, cols = d.shape

    N = round( math.sqrt( slices ) )

    if rc == None:
        im_cols = N
        im_rows = N
        
        if im_rows*im_cols < slices:
            im_rows = im_rows + 1
    
        extra = im_cols * im_rows - slices
    else:
        im_rows, im_cols = rc

    ii = 0
    d2 = zeros( (im_rows*rows, im_cols*cols), dtype=d.dtype)

    s1 = time.time()
    for ri in arange(0, im_rows):
        for ci in arange(0, im_cols):
            if ii == slices: break
            d2[ri*rows:(ri+1)*rows,ci*cols:(ci+1)*cols] = d[ii,:,:]
            ii = ii + 1
    mimage(d2, cmap, )

def mmontagecb( d, cmap=None ):
    mmontage(d, cmap)
    ax = gcf().add_axes([0.85, 0.05, 0.04, 0.9])
    colorbar(cax=ax, orientation='vertical')

def showdec( vec, fa):
    """ Display a diffusion encoded colormap based on a slice of vec and fa """
    imshow(abs(vec) * tile(fa, (3,1,1)).transpose(1,2,0) )

def showdecmore(vec, fa):
    """ Supposedly for a 3d volume, but doesn't really work"""
    print('this doesn\'t work properly')
    mmontage( abs(vec) * tile(fa, (3,1,1,1)).transpose(1,2,3,0))

def showcmap( d ):
    slices, rows, cols, rgb = d.shape

    N = round( math.sqrt( slices ) )

    im_cols = N
    im_rows = N
    
    if im_rows*im_cols < slices:
        im_rows = im_rows + 1

    extra = im_cols * im_rows - slices

    ii = 0
    d2 = zeros( (im_rows*rows, im_cols*cols, 3), dtype=d.dtype)

    s1 = time.time()
    for ri in arange(0, im_rows):
        for ci in arange(0, im_cols):
            if ii == slices: break
            d2[ri*rows:(ri+1)*rows,ci*cols:(ci+1)*cols,0] = d[ii,:,:,2]
            d2[ri*rows:(ri+1)*rows,ci*cols:(ci+1)*cols,1] = d[ii,:,:,0]
            d2[ri*rows:(ri+1)*rows,ci*cols:(ci+1)*cols,2] = d[ii,:,:,1]
            ii = ii + 1

    imshow( d2 )
    xticks([])
    yticks([])

def plotcest( freq, data ):

    plot( freq, data )
    xlabel('Frequency Offset from Water (ppm)')
    ylabel('Relative Signal (%)')
    grid('on')
    xl = gca().get_xlim()
    if xl[0] < xl[-1]:
        gca().set_xlim( xl[::-1] )

##==============================================================================
##
##  Fitting
##
##==============================================================================

def t1abs( a1, d1, a2, d2, TR):

    if a1 > np.pi:
        a1 = a1 / 180 * pi

    if a2 > np.pi:
        a2 = a2 / 180 * pi

    r = d2 / d1

    denom = (sin(a2)-r*sin(a1))/(cos(a1)*sin(a2)-r*cos(a2)*sin(a1))

    T1 = -TR / log( denom )

    return T1

def t2fit(te, mm):
    """ Do a mono-exponential decay curve fit to given data """

    # Calculate some mins and max
    noise_max = 2*min(mm)
    pd_max = 2*max(mm)

    coeffs = polyfit( te, log(mm), 1 )    
    t2_guess = -1 / coeffs[0]

    # Lambda to calculate the residuals between the simulated and the measured data
    residuals = lambda x, te, mm: sum(   (  (x[0]*array([math.exp(-tt/x[1]) for tt in te])+x[2] ) - mm  )**2   )

    # Set the initial parameters: PD T2 offset
    p0 = array([mm[0], t2_guess, mm[-1]/2])

    # Call the optimization program
    plsq = fmin_l_bfgs_b(residuals, p0, args=(te, mm), bounds=[(0, pd_max), (0.1, 1200), (0, noise_max) ], approx_grad=True)
    #plsq = fmin_tnc(residuals, p0, args=(te, mm), bounds=[(0, pd_max), (t2_guess/2, t2_guess*2), (0, noise_max) ], approx_grad=True, messages=0)
    
    # Return the appropriate values
    return plsq[0]

def t2fit_leastsq(te, mm):
    """ Do a mono-exponential decay curve fit to given data """

    # Calculate some mins and max
    noise_max = 2*min(mm)
    pd_max = 2*max(mm)

    coeffs = polyfit( te, log(mm), 1 )    
    t2_guess = -1 / coeffs[0]

    # Lambda to calculate the residuals between the simulated and the measured data

    #residuals = lambda x, te, mm: (x[0]*array([math.exp(-tt/x[1]) for tt in te])+x[2] ) - mm  
    def residuals( x, te, mm ): 

        if any(x<0):
            x[x<0] = 0.001 
        res = (x[0]*array([math.exp(-tt/x[1]) for tt in te])+x[2] ) - mm  
        return res

    p0 = array([exp(coeffs[1]), t2_guess, mm[-1]/2])

    # Call the optimization program
    x,ier = leastsq(residuals, p0, args=(te, mm))
    
    # Return the appropriate values
    return x,ier

##==============================================================================
##
##  Filtering
##
##==============================================================================

def gauss_kern(size, sizey=None):
     """ Returns a normalized 2D gauss kernel array for convolutions """
     size = int(size)
     if not sizey:
         sizey = size
     else:
         sizey = int(sizey)
     x, y = mgrid[-size:size+1, -sizey:sizey+1]
     g = exp(-(x**2/float(size)+y**2/float(sizey)))
     return g / g.sum()

def ave_kern(size, sizey=None, sizez=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """

    if( sizey == None ):
        sizey = size

    if( sizez == None ):
        A = np.ones( (size, sizey) )
    else:
        A = np.ones( (size, sizey, sizez) )

    return A / sum(A) 


def image_smooth(im, kern, N=1):
    #aa = convolve2d( im, kern)
    aa = im
    for ii in range(N):
        aa = np.numarray.nd_image.convolve( aa, kern )
    return aa

rad2deg = lambda f: f/np.pi*180.0

def aniso(v, kappa=-1, N=1):

    if kappa == -1:
        kappa = prctile(v, 40)

        vf = v.copy()

        for ii in range(N):
                dE = -vf + roll(vf,-1,0)
                dW = vf - roll(vf,1,0)

                dN = -vf + roll(vf,-1,1)
                dS = vf - roll(vf,1,1)

        if len(v.shape) > 2:
            dU = -vf + roll(vf,-1,2)
            dD = vf - roll(vf,1,2)

            vf = (vf + 
                  3./28. * (((exp(- (abs(dE) / kappa)**2 ) * dE) - (exp(- (abs(dW) / kappa)**2 ) * dW)) +
                            ((exp(- (abs(dN) / kappa)**2 ) * dN) - (exp(- (abs(dS) / kappa)**2 ) * dS))))

        if len(v.shape) > 2:
                        vf += 1./28. * ((exp(- (abs(dU) / kappa)**2 ) * dU) - (exp(- (abs(dD) / kappa)**2 ) * dD))

        return vf


##==============================================================================
##
##  Calculating
##
##==============================================================================

def dam( f1, f2, radians=True):
    """ Double angle B1 quantification method.
    1. R Stollberger and P Wach, "Imaging of the active B1 field in vivo," J Magnetic Resonance in Medicine: Official Journal of the Society of Magnetic Resonance in Medicine / Society of Magnetic Resonance in Medicine 35, no. 2 (February 1996): 246-251.

    """
    flip = arccos( f2[1] / (2 * f1[1] ) )

    if radians == False:
        flip = rad2deg( flip )

    return flip

def mtr( off, on, prthresh=70):
    mask = off > prctile( off, prthresh)

    mtr = (off - on) / off * mask * 100

    return mtr

# Rephase a real/imaginary pair of images 
def rephaseV( real_in, imag_in, N=3, kernel_size=4 ):
    gk = ave_kern(kernel_size, kernel_size, kernel_size)
    real_s, imag_s = real_in, imag_in

    for dummy in range(N):
        real_s = image_smooth( real_s, gk) 
        imag_s = image_smooth( imag_s, gk) 

    co = real_in + 1j * imag_in
    cs = real_s + 1j * imag_s
    
    phasor = np.conj( cs / abs(cs) )

    rephased = co * phasor

    return rephased

def rephase( real_in, imag_in, N=3, kernel_size=4 ):
    gk = ave_kern(kernel_size)
    real_s, imag_s = real_in, imag_in

    for dummy in range(N):
        if len( real_s.shape ) > 2: 
            for sl in range( real_s.shape[0] ):
                real_s[sl] = image_smooth( real_s[sl], gk) 
                imag_s[sl] = image_smooth( imag_s[sl], gk) 
        else:
            real_s = image_smooth( real_s, gk) 
            imag_s = image_smooth( imag_s, gk) 

    co = real_in + 1j * imag_in
    cs = real_s + 1j * imag_s
    
    phasor = np.conj( cs / abs(cs) )

    rephased = co * phasor

    return rephased

def unwrap_phase( real_in, imag_in, N=10, kernel_size=4 ):
    tt = rephase( real_in, imag_in, N=N, kernel_size=kernel_size )
    return arctan( imag(tt) / real(tt) )

#  THis is used to calculate a frequency offset map
#  from a WASSR acquisition (APT/CEST).
#  Ref: Mina Kim et al., "Water saturation shift referencing (WASSR) for chemical exchange saturation transfer (CEST) experiments" MRM 61:6 2009 pp 1441-1450.
def calculateOffsetMap( freq, data, percentileThreshold=60, extrapolate=1.0, freqFactor=10 ):
    ##  If extrapolate is > 1.0 then the curve will be extrapolated.

    sort_inds = freq.argsort()
    _freq = freq[sort_inds]
    _data = data[sort_inds]

    ##  Check to make sure the frequencies are all increasing
    if len( _data.shape ) != 4:
        print("calculateOffsetMap: The data must be a volume.")
        raise

    ##  Calculate the threshold on which to work
    thresh = prctile( _data[-1], percentileThreshold )    

    ##  Set the array of frequncy offsets
    minp = zeros(_data[0].shape)

    ##  Find all the coordinates in the data that are above the threshold
    coords = array( np.nonzero( _data[-1] > thresh ) ).transpose()

    pbar = ProgressBar(widgets=['Calc Offset Map ', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    x, A, w = -1, -1, -1

    #cjinds = logical_and( _freq<inf, _freq!=0 )
    cjinds = nonzero( _freq<inf )[0]
    for ii,coord in enumerate(coords):
        
        ##  Get the data.
        s,r,c = coord
        mm = _data[:,s,r,c]

        x = wassrLL( _freq[cjinds], mm[cjinds] )
        minp[s,r,c] = x

        pbar.update(ii)

    pbar.finish()

    return minp

##  Fix the asymmetry map based on the water offset calculated above
def correctB0( target_freqs, freq, data, minp, percentileThreshold=60 ):

    ##  Calculate the threshold
    thresh = prctile( data[0], percentileThreshold )    

    ##  Create the pos and neg offset information
    output = zeros( concatenate( (array([len(target_freqs)]), data[0].shape), axis=0 ) )

    coords = array( np.nonzero( data[0] > thresh ) ).transpose()

    pbar = ProgressBar(widgets=['Calc Asym', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()
    for ii,coord in enumerate(coords):
        ##  Get the data.
        s,r,c = coord

        for jj, tf in enumerate( target_freqs ):

            inds = nonzero( abs( freq - tf ) < 300 )[0]

            # Do the positive side
            mm = data[inds,s,r,c]
            tt = interpolate.splrep( freq[inds], mm, k=1 )

            output[jj,s,r,c] = interpolate.splev( tf + minp[s,r,c], tt )

        pbar.update(ii)

    pbar.finish()

    return output

##  Fix the asymmetry map based on the water offset calculated above
def calculateAsym( freq, data, minp, percentileThreshold=60 ):

    ##  Calculate the threshold
    thresh = prctile( data[0], percentileThreshold )    

    pos448 = zeros(data[0].shape)
    neg448 = zeros(data[0].shape)

    coords = array( np.nonzero( data[0] > thresh ) ).transpose()

    pbar = ProgressBar(widgets=['Calc Asym', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    inds_pos = [ii for ii in range(len(freq)) if freq[ii] > 0 and freq[ii] < inf]
    inds_pos = [inds_pos[a] for a in freq[inds_pos].argsort()]
    
    inds_neg = [ii for ii in range(len(freq)) if freq[ii] < 0]
    inds_neg = [inds_neg[a] for a in freq[inds_neg].argsort()]

    inds_neg = inds_neg[::-1]

    for ii,coord in enumerate(coords):
        ##  Get the data.
        s,r,c = coord

        # Do the positive side
        mm = data[inds_pos,s,r,c]
        tt = interpolate.splrep(freq[inds_pos], mm, k=5)
        ns = interpolate.splev(448+minp[s,r,c],tt)
        pos448[s,r,c] = ns

        # Do the negative side
        mm = data[inds_neg,s,r,c]
        tt = interpolate.splrep(-freq[inds_neg], mm, k=3)
        ns = interpolate.splev(448-minp[s,r,c],tt)
        neg448[s,r,c] = ns

        pbar.update(ii)

    pbar.finish()

    return (neg448, pos448)

def wassrNEW(freq, data, freqFactor=10):
    yy = data - min(data)
    w0 = mean( freq[ nonzero( yy < median(yy)/2)[0] ] ) + (rand(1)-0.5)/1e3

    def symscore( guess, freq, data, tt ):
        if len(guess) == 1:
            #tt = interpolate.splrep(freq, data, k=3)
            yc = interpolate.splev( guess-freq, tt )
            return mean( abs( (data-yc)**2 ) )

    tt = interpolate.splrep(freq, data, k=3)
    x = fmin_bfgs( symscore, 2*w0[0], args=(freq, data, tt), full_output=False, disp=False, retall=False )[0]/2

    return x

def wassrLL(freq, data):

    def llminp(x, freq, data ):
        A = x[0]
        x0 = x[1]
        w = x[2]    

        Ap = 1 - w**2/(w**2 + 4*(x0 - freq)**2)
        x0p = -A*w**2*(-8*x0 + 8*freq)/(w**2 + 4*(x0 - freq)**2)**2
        wp = A*(-2*w/(w**2 + 4*(x0 - freq)**2) + 2*w**3/(w**2 + 4*(x0 - freq)**2)**2)
        return [Ap, x0p, wp]

    def f(freq, A, x0, w):

        return A * ( 1 - ( w**2 / ( w**2 + 4 * (x0-freq)**2 ) ) )

    def llmin(x, freq, data):
        A = x[0]
        x0 = x[1]
        w = x[2]    

        ll = A * ( 1 - ( w**2 / ( w**2 + 4 * (x0-freq)**2 ) ) )

        return sum( (data-ll)**2 )

    A = max(data)
    x0 = freq[ nonzero( data == min(data) )[0][0] ] # Hz
    w = 30 # Hz

    try:
        x = curve_fit(f, freq, data, p0=(A,x0,w), sigma=None)[0]
        
    except RuntimeError:
        x = fmin( llmin, (A,x0,w), args=(freq, data), full_output=False, disp=False, retall=False)

    return x[1]

def wassr(freq, data):
    yy = data - min(data)
    w0 = mean( freq[ nonzero( yy < median(yy)/2)[0] ] ) + (rand(1)-0.5)/1e3

    def symscore( guess, freq, data, freq_rev, data_rev):
        try:
            nf = (guess-freq_rev)
            tt = interpolate.splrep(guess-freq_rev, data_rev, k=3)
            yc = interpolate.splev( freq, tt )
            score = mean([ abs((data[ii]-yc[ii])**2) for ii,j in enumerate( yc ) if nf[ii] > freq[0] and nf[ii] < freq[-1] ] )
        except:
            score = 999999999;

        return score

    x = fmin( symscore, 2*w0[0], args=(freq, data, freq[::-1], data[::-1]), full_output=False,disp=False,retall=False )[0]/2
    #x = fminbound( symscore, -50, 50, args=(freq, data, freq[::-1], data[::-1]), full_output=False,disp=False)[0]/2

    return x


def wassrblah(freq, data, freqFactor=10):
    yy = data - min(data)
    w0 = mean( freq[ nonzero( yy < median(yy)/2)[0] ] ) + (rand(1)-0.5)/1e3

    def symscore( guess, freq, data, tt ):
        if len(guess) == 1:
            #tt = interpolate.splrep(freq, data, k=3)
            yc = interpolate.splev( guess-freq, tt )
            return mean( abs( (data-yc)**2 ) )

    superFreq = linspace( freq.min(), freq.max(), freq.shape[0]*freqFactor )
    tt = interpolate.splrep( freq, data, k=3 )
    superData = interpolate.splev( superFreq, tt )

    tt = interpolate.splrep(superFreq, superData, k=3)
    x = fmin_bfgs( symscore, 2*w0[0], args=(superFreq, superData, tt), epsilon=1.0,full_output=False, disp=False, retall=False )[0]/2

    #fig100 = figure(100)
    #clf()
    #plot( freq, data )
    #plot( superFreq, superData )
    #plot( [x, x], [0, data.max()] )
    #fig100.canvas.draw()

    return x

def wassrCorrect(freq, data, offset):

    freqSuper = linspace( freq[0], freq[-1], freq.shape[0]*50 )
    tt = interpolate.splrep( freq, data, k=3)
    dataSuper = interpolate.splev( freqSuper, tt )

    tt = interpolate.splrep(freqSuper, dataSuper, k=3)
    mmout = interpolate.splev( offset - freqSuper, tt )

    tt = interpolate.splrep(freqSuper, mmout, k=3)
    mmout = interpolate.splev( freq, tt )

    return mmout

def cestFit( freq, mm, fitinds, newfreq=None ):
    def lorentzian(freq, A, x0, w, b, k): 
        tt = A * ( 1 - ( w**2 / ( k* ( w**2 + (x0-freq)**2 ) ) ) ) + b 
        tt[ nonzero(tt < 0 ) ] = 0
        return tt

    def lorentzian_min(x, freq, data):
        return sum( (lorentzian(freq, x[0], x[1], x[2], x[3], x[4]) - (data) )**2 )  

    def lorentzianFit(freq, data):

        A = max(data)
        x0 = freq[ nonzero( data == min(data) )[0][0] ] # Hz
        w = 0.1 # Hz
        b = 0 
        k = 1

        #x = curve_fit(lorentzian, freq, data, p0=(A,x0,w), sigma=None)[0]
        x = fmin_l_bfgs_b(lorentzian_min, fprime=[], x0=(A,x0,w, b, k), args=(freq,data), approx_grad=True, bounds=((None, None),(None, None),(None, None),(None, None),(1,1)), factr=10.0)[0]

        #  This will return A, x0, w, b
        return x
    
        ## I got this off of stackoverflow.com
        ## http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-a-an-extrapolated-result-beyond-the-input-rang
    def extrap1d(interpolator):
        xs = interpolator.x
        ys = interpolator.y

        def pointwise(x):
            if x < xs[0]:
                return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            elif x > xs[-1]:
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
                return interpolator(x)

        def ufunclike(xs):
            return array(map(pointwise, array(xs)))

        return ufunclike

    A,x0,w,b,k = lorentzianFit( freq[fitinds], mm[fitinds] )

    #print('Fitted parameters:')
    #print('\tA:', A)
    #print('\tx0:', x0)
    #print('\tw:', w)
    #print('\tb:', b)
    #print('\tk:', k)

    if newfreq == None:
        newfreq = linspace(freq[1], freq[-1], 2*len(freq))        

    #tt = interpolate.splrep( freq[3:-2], mm[3:-2], k=1)
    #mm_fixed = interpolate.splev( newfreq + x0, tt)
    ff = interpolate.interp1d( freq-x0, mm, kind='linear')
    ff = extrap1d( ff )
    mm_fixed = ff( newfreq )

#    tt = interpolate.splrep( freq[1:], lorentzian(freq[1:], A, x0, w, b,k), k=1)
#    lorentzian_fixed = interpolate.splev( newfreq + x0, tt)
    gg = extrap1d( interpolate.interp1d( freq-x0, lorentzian(freq, A, x0, w, b,k), kind='slinear' ) )
    lorentzian_fixed = gg( newfreq )
#    lorentzian_fixed = lorentzian(newfreq, A, x0, w, b,k)

#    figure(100)
#    clf()
#    plot( newfreq, mm_fixed, 'bo')
#    figure(100).canvas.draw()

    return newfreq, mm_fixed, lorentzian_fixed, A, x0, w, b, k



def swi( filename ):

    #  Read in the data
    d = readrec( filename )

    # Determine the file types
    fp = open( filename.replace('REC', 'PAR'), 'rt' )
    par = fp.readlines()
    fp.close()

    re_spaces = re.compile("\s+")

    # Get the reconstruction resolution
    numbers = re_number.findall( par[-3] )
    xsize,ysize = int(numbers[9]), int(numbers[10])

    # Get the RI, RS and SS to rescale the images into the 
    # proper floating point values
    # DV = PV * RS + RI   FP = DV / (RS * SS)
    # Convert the pixel value to floating point value
    ri,rs,ss = [], [], []
    types = []
    for ii,l in enumerate( par ):
        numbers = re_spaces.split(l.lstrip().rstrip())
        if re_number.match(numbers[0]) and re_number.match(numbers[1]) and len(numbers) > 10:
            print(numbers)
            ri.append( float( numbers[11] ) )
            rs.append( float( numbers[12] ) )
            ss.append( float( numbers[13] ) )
            types.append( int( numbers[4] ) )

    types = np.unique( types )
 
    print(types)


def swiRI( rr, ii, N=5 ):
    for num in range(N):
        rr = scipy.ndimage.filters.uniform_filter( rr, [3,3,3] )
        ii = scipy.ndimage.filters.uniform_filter( ii, [3,3,3] )

    return abs(rr+1j*ii) * ( ( pi - angle(rr+1j*ii))/pi )**4

def polyCorrect(freq, data, offset):

    freqSuper = linspace( freq[0], freq[-1], len(freq)*50 )
    dataSuper = interpolate.spline(freq, data, freqSuper, order=3, kind='smoothest')

    mmout = interpolate.spline(freqSuper, dataSuper, offset-freqSuper, order=3, kind='smoothest')

    mmout = interpolate.spline(freqSuper, mmout, freq, order=3, kind='smoothest')

    return mmout

def nnls_fit( te, y, t2 ):
    A = exp(- outer( te,  r_[ [1/t2a for t2a in t2], 0]) )

    if False:
        H = 0.0*diag(1*ones((A.shape[1],)))

        #H = diag(1*ones((A.shape[1],)))
        #H = H + diag(-1*ones((A.shape[1],)), k=1)[:-1,:-1]
        yt = zeros(( A.shape[1] ))
        Att = concatenate( (A, H), axis=0 )
        ytt = concatenate( (y, yt), axis=0 )

        x = scipy.optimize.nnls(Att, ytt)[0]
    else:
        x = scipy.optimize.nnls(A, y)[0]

    ## Compute the fitted data
    y_fit = inner(A, x)

    ## Compute the chi2
    chi2 = sqrt( sum( ( y - y_fit)**2 ) )

    return x, y_fit, chi2
#    return x, 0,0

def t2map_nnls( te, data, t2, percentileThreshold=65 ):
    print('create x and chi2')
    x = zeros(( r_[len(t2)+1, data[0].shape] ))
    chi2 = zeros(( data[0].shape ))

    print('compute threshold')
    threshold = prctile( data[0], percentileThreshold )

    print('create coords')
    coords = array( np.nonzero( data[0] > threshold ) ).transpose()

    print('t2 is', t2)

    pbar = ProgressBar(widgets=['Calc T2star', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    for ii,coord in enumerate(coords):
        si,ci,ri = coords[ii]
        x[:,si,ci,ri],y_fit,chi2[si,ci,ri] = nnls_fit(te, squeeze(data[:,si,ci,ri]), t2)

#        mm = data[:,si,ci,ri]
#        mm[ nonzero(mm==0) ] = 1
#        a,b,c = t2fit(te, mm)
#
#        figure(1)
#        clf()
#        plot( te, data[:,si,ci,ri], 'x' )
#        plot( te, y_fit, label='nnls' )
#        plot( te, a*exp(-te/b)+c, label='monoexp' )
#
#        tt_t2_monoexp = b
#        tt_t2_nnls = dot(x[:,si,ci,ri], r_[t2,1] ) / sum(x[:,si,ci,ri])
#        tt_t2_nnls = exp( dot(x[:,si,ci,ri], log(r_[t2,1]) ) / sum(x[:,si,ci,ri]) )
#
#        title('Monoexp = %3.2f   NNLS = %3.2f ' % ( tt_t2_monoexp, tt_t2_nnls ) )
#        legend()
#        figure(1).canvas.draw()
#    
#        if( abs(tt_t2_monoexp - tt_t2_nnls ) > 5 ):
#            inds = nonzero( x[:-1,si,ci,ri] > 0 )[0]
#            print(x[inds,si,ci,ri], t2[inds])
#            blah = raw_input('waitng: ')

        pbar.update(ii)

    pbar.finish()

    return x, chi2

def t2map_logfit( te, data, percentileThreshold=65 ):
    print('create x and chi2')
    t2 = zeros(( data[0].shape ))

    print('compute threshold')
    threshold = prctile( data[0], percentileThreshold )

    print('create coords')
    coords = array( np.nonzero( data[0] > threshold ) ).transpose()

    pbar = ProgressBar(widgets=['Calc T2star', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    for ii,coord in enumerate(coords):
        si,ci,ri = coords[ii]
        
        mm = data[:,si,ci,ri]
        mm[ nonzero(mm==0) ] = 1
        blah = log( mm )
        tt = polyfit( te, blah, 1 )
        if tt[0] == 0:
            t2[si,ci,ri] = 0
        else:
            t2[si,ci,ri] = -1 / tt[0]

        pbar.update(ii)

    pbar.finish()

    return t2

def t2map_monoexp( te, data, percentileThreshold=65 ):
    print('create x and chi2')
    t2 = zeros(( data[0].shape ))

    print('compute threshold')
    threshold = prctile( data[0], percentileThreshold )

    print('create coords')
    coords = array( np.nonzero( data[0] > threshold ) ).transpose()

    pbar = ProgressBar(widgets=['Calc T2star', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    for ii,coord in enumerate(coords):
        si,ci,ri = coords[ii]
        
        mm = data[:,si,ci,ri]
        mm[ nonzero( mm == 0 ) ] = 1
        t2[si,ci,ri] = t2fit(te, mm)[1]

        pbar.update(ii)

    pbar.finish()

    return t2



def meanT2( x, t2v ):
    t2 = [ dot( x[:,s,r,c], r_[t2v,1] )/sum(x[:,s,r,c]) for s in range(x.shape[1]) for r in range(x.shape[2]) for c in range(x.shape[3]) ]
    t2 = array(t2).reshape( x.shape[1:] )
    t2[ nonzero( x[:-1].sum(axis=0) == 0 ) ] = 0
    return t2

def diff_calc( data, grads, b, pthresh=70 ):
    
    zeroi = nonzero(sum(grads, axis=1)==0)[0]
    nonzeroi = nonzero(sum(grads, axis=1)!=0)[0]

    ## Create the G matrix
    gradsnonzero = grads[nonzeroi]
    G = zeros((len(nonzeroi), 6))
    for ii in range( len(nonzeroi) ):
        G[ii,0] = b * gradsnonzero[ii,0] * gradsnonzero[ii,0]
        G[ii,1] = b * gradsnonzero[ii,1] * gradsnonzero[ii,1]
        G[ii,2] = b * gradsnonzero[ii,2] * gradsnonzero[ii,2]
        G[ii,3] = b * gradsnonzero[ii,0] * gradsnonzero[ii,1] * 2
        G[ii,4] = b * gradsnonzero[ii,0] * gradsnonzero[ii,2] * 2
        G[ii,5] = b * gradsnonzero[ii,1] * gradsnonzero[ii,2] * 2

    Ginv = pinv(G)

    threshold = prctile( data[0], pthresh )
    
    evalues = zeros( concatenate(([3], data.shape[1:])) )
    evectors = zeros( concatenate(([3], data.shape[1:])) )

    coords = array( np.nonzero( data[0] > threshold ) ).transpose()

    pbar = ProgressBar(widgets=['DTI Calculation: ', Percentage(), Bar(), ETA()], maxval=coords.shape[0]).start()

    for ii,coord in enumerate(coords):
        si,ci,ri = coords[ii]
    
        tt = data[:,si,ci,ri]
        tt = tt[nonzeroi] / mean(tt[zeroi])
        tt = -log(tt)
        tt[ isinf(tt) ] = 0
        tt[ isnan(tt) ] = 0
        
        tv = dot( Ginv , tt )

        TV = zeros((3,3))    
        TV[0,0] = tv[0]
        TV[1,1] = tv[1]
        TV[2,2] = tv[2]

        TV[0,1] = tv[3]
        TV[1,0] = tv[3]

        TV[0,2] = tv[4]
        TV[2,0] = tv[4]

        TV[1,2] = tv[5]
        TV[2,1] = tv[5]

        try:
            #w,v = scipy.linalg.eig(TV)
            #inds = argsort(w)[::-1]

            #evalues[:,si,ci,ri] = w[inds]
            #evectors[:,si,ci,ri] = v[:,inds[0]]

            u,s,v = scipy.linalg.svd(TV)

            evalues[:,si,ci,ri] = s
            evectors[:,si,ci,ri] = u[:,0]
        except:
            print('here')
            pass

        pbar.update(ii)

    pbar.finish()
    return evalues, evectors

def diff_fa( evalues ):
    sum_sq_values = sum(evalues**2, axis=0)
    sde = (evalues[0]-evalues[1])**2 + (evalues[1]-evalues[2])**2 + (evalues[0] - evalues[2])**2 
    fa = sqrt( sde / (2 * sum_sq_values) )
    fa[ nonzero( sum_sq_values == 0 ) ] = 0    

    return fa

def diff_cmap( evectors, fa):
    ev = zeros(evectors.shape)
    ev[0] = evectors[0] * fa
    ev[1] = evectors[1] * fa
    ev[2] = evectors[2] * fa
    cmap = abs( ev.swapaxes(3,0).swapaxes(0,1).swapaxes(1,2) )
    return cmap

def readgrad(fn):
    if not fn.endswith('grad'):
        print('readgrad: Meant to be used on a gradient file.')
        return None
    
    fp = open(fn, 'rt')
    lines = fp.readlines()
    fp.close()

    #re_float_number = re.compile("-{0,1}\d{0,}[.]\d{0,}")
    #numbers = re_float_number.findall(lines[0])
    a = lambda l: re.findall("-{0,1}\d{0,}[.]\d{0,}", l)
    numbers = map( vectorize(float), map( a, lines ) )

    return array( numbers )

##
##  Image Alignment
##
def align(fn_static, fn_dyn):
    static_ap_fov, static_fh_fov, static_rl_fov, static_nx, static_ny, static_nz, static_ap_offcenter, static_fh_offcenter, static_rl_offcenter = getinfo(fn_static)

    dyn_ap_fov, dyn_fh_fov, dyn_rl_fov, dyn_nx, dyn_ny, dyn_nz, dyn_ap_offcenter, dyn_fh_offcenter, dyn_rl_offcenter = getinfo(fn_dyn)
    d_dyn = readrec(fn_dyn.replace('PAR','REC'))

    dyn_z,dyn_x,dyn_y = ogrid[-dyn_fh_fov/2+dyn_fh_offcenter:dyn_fh_fov/2+dyn_fh_offcenter:dyn_fh_fov/dyn_nz,-dyn_rl_fov/2:dyn_rl_fov/2:dyn_rl_fov/dyn_nx,-dyn_ap_fov/2:dyn_ap_fov/2:dyn_ap_fov/dyn_ny]
    static_z,static_x,static_y = mgrid[-static_fh_fov/2+static_fh_offcenter:static_fh_fov/2+static_fh_offcenter:static_fh_fov/static_nz,-static_rl_fov/2:static_rl_fov/2:static_rl_fov/static_nx,-static_ap_fov/2:static_ap_fov/2:static_ap_fov/static_ny]
    x0 = dyn_x[0,0,0]
    y0 = dyn_y[0,0,0]
    z0 = dyn_z[0,0,0]
    dx = dyn_x[0,1,0] - x0
    dy = dyn_y[0,0,1] - y0
    dz = dyn_z[1,0,0] - z0
    ivals = (static_x - x0)/dx
    jvals = (static_y - y0)/dy
    kvals = (static_z - z0)/dz
    coords = array([kvals, ivals, jvals])
    return ndimage.map_coordinates(d_dyn, coords)
    

def getinfo(fn):

    fp = open(fn, 'rt')
    lines = fp.readlines()
    fp.close()
    
    re_number = re.compile("-{0,1}\d+\.{0,1}\d+")
    
    # Get the number of slices
    expr = re.compile("Max. number of slices")
    numbers = re_number.findall( filter(expr.search, lines)[0] )
    nz = int(numbers[0])
    
    # Get the number of voxels
    expr = re.compile("Scan resolution")
    numbers = re_number.findall( filter(expr.search, lines)[0] )
    nx,ny = int(numbers[0]), int(numbers[1])
    
    # Get the offcenter information
    expr = re.compile("Off Centre midslice")
    numbers = re_number.findall( filter(expr.search, lines)[0] )
    ap_offcenter,fh_offcenter,rl_offcenter = float(numbers[0]), float(numbers[1]), float(numbers[2])
    
    # Get the angulation information
    expr = re.compile("Angulation midslice")
    numbers = re_number.findall( filter(expr.search, lines)[0] )
    ap_angulation,fh_angulation,rl_angulation = float(numbers[0]), float(numbers[1]), float(numbers[2])
    
    # Get the FOV information
    expr = re.compile("FOV \(ap,fh")
    numbers = re_number.findall( filter(expr.search, lines)[0] )
    ap_fov,fh_fov,rl_fov = float(numbers[0]), float(numbers[1]), float(numbers[2])
    
    return ap_fov, fh_fov, rl_fov, nx, ny, nz, ap_offcenter, fh_offcenter, rl_offcenter
    
def interp3(z,y,x,data, Z,Y,X):
    """  This is to be used to interpolate a 3D dataset similar 
             to the interpn command in Matlab
    """
    x0 = x[0,0,0]
    y0 = y[0,0,0]
    z0 = z[0,0,0]
    dx = x[0,1,0] - x0
    dy = y[0,0,1] - y0
    dz = z[1,0,0] - z0
    ivals = (newx - x0)/dx
    jvals = (newy - y0)/dy
    kvals = (newz - z0)/dz
    coords = array([kvals, ivals, jvals])
    newf = ndimage.map_coordinates(d7, coords)

    return newf

cdict =   {'red':   ((0, 0.5, 0.5), 
                     (0.33, 1, 1), 
                     (0.65,1, 1), 
                     (1, 0.5, 0.5)),
               'green': ((0., 1, 1), 
                     (0.33,1, 1), 
                     (0.65,0,0), 
                     (1, 0, 0)), 
               'blue':  ((0., 0.5, 0.5), 
                     (0.33, 0, 0), 
                     (1, 0, 0))} 


cmap_gr = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

##==============================================================================
##
##  Pickler
##
##==============================================================================

# Do a pickle save but put the name 
# of the variable before each of the data
# values.

def savep(fn, **kwargs):
    fp = open(fn, 'wb')

    for key in kwargs.keys():
        pickle.dump(key, fp)
        pickle.dump(kwargs[key], fp)

    fp.close()

# Do a pickle load, but the data is paired
# so the first value is the name of the variable
# and the second is the actual value.
def loadp(fn):
    try:
        fp = open(fn, 'rb')
    except:
        print('cjlib.pload: Could not open file %s' % fn)
        return

    while True:
        #  Get the variable name
        try:    
            this_is_a_temporary_name = pickle.load(fp)
        except:
            break

        #  Get the data
        xxxxxx={}
        exec("""xxxxxx = { '%s': pickle.load(fp) }""" % this_is_a_temporary_name)
        import inspect
        inspect.currentframe(1).f_globals.update( xxxxxx )

    fp.close()

