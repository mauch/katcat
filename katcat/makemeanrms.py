#!/bin/python
import os

from astropy.io import fits
import numpy as np
import numba
import random
from scipy import interpolate


MAD_TO_SD = 1.4826

@numba.jit(nopython=True)
def mad(data):
    med = np.nanmedian(data)
    dev = np.abs(data - med)
    return med, MAD_TO_SD * np.nanmedian(dev)

@numba.jit(nopython=True)
def mode(data):
    med, sd = mad(data)
    cut = data - med < 0.5 * sd
    cut &= data - med > -1.0 * sd
    data = data[cut]
    # Find the x=0 point of a histogram derivative
    # - should be the mode for a Gaussian
    # Average 2000 elements per bin
    # no less than 10 and no more than 50 bins
    # also- just return median for < 1000 data points
    if len(data) < 1000:
        return med
    nbin = len(data) / 1000
    nbin = max(10, nbin)
    nbin = min(50, nbin)
    hist, bin_edges = np.histogram(data, bins=nbin)
    bins = bin_edges[1:-1]
    deriv = hist[1:] - hist[:-1]
    # Use np.lstsq to compute line of best fit to the derivative
    # and find its y=0 zero point
    A = np.vstack((bins, np.ones(len(deriv)),)).T
    m, c = np.linalg.lstsq(A, deriv.astype(np.float64))[0]
    y0 = -c / m
    return y0

@numba.jit(nopython=True)
def itersdmean(data):
    diff = 1
    i = 0
    data = data[data != 0.0]
    if len(data) == 0:
        return 0.0, 0.0
    med, sd = mad(data)
    while diff > 0:
        i += 1
        if i > 50:
            return sd, mode(data)
        old_sd = sd
        old_len = len(data)
        cut = np.abs(data - med) < 5.0 * sd
        if np.all(~cut):
            return sd, med
        data = data[cut] 
        med, sd = mad(data)
        diff = old_len - len(data)
        if sd == 0.0:
            mean = med
            return old_sd, mode(data)
    mean = mode(data)
    return sd, mean

def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fit_gauss(data, initx0, initsigma):
    try:
        if len(data) < 1000:
            return initx0, initsigma
        hist, bin_edges = np.histogram(data[data < initx0 + 0.2*initsigma], bins=50)
        bins = (bin_edges[1:] + bin_edges[:-1]) / 2.
        out = curve_fit(gaussian, bins, hist, p0=[max(hist), initx0, initsigma])
        return out[0][1], out[0][2]
    except RuntimeError:
        return initx0, initsigma

@numba.jit(nopython=True,parallel=True)
def medfilt(data, region, x_sub, y_sub, sdarr, meanarr):
    nxpix = len(x_sub)
    nypix = len(y_sub)
    halfbox = region // 2
    for x_lp in range(nxpix):
        x = x_sub[x_lp]
        x_start = max(x - halfbox, 0)
        x_stop = min(x + halfbox, data.shape[0])
    
        for y_lp in numba.prange(nypix):
            y = y_sub[y_lp]
            y_start = max(y - halfbox, 0)
            y_stop = min(y + halfbox, data.shape[1])
            this_data = data[x_start:x_stop,y_start:y_stop].flatten().copy()
            this_sd, this_mean = itersdmean(this_data)
            sdarr[x_lp, y_lp] = this_sd
            meanarr[x_lp, y_lp] = this_mean


def medfiltobit(data, region, x_sub, y_sub, sdarr, meanarr):

    desc = data.Desc.Dict
    xmax = desc['inaxes'][0]
    ymax = desc['inaxes'][1]
    nxpix = len(x_sub)
    nypix = len(y_sub)
    halfbox = region // 2
    for x_lp, x in enumerate(x_sub):
        x += 1
        x_start = max(x - halfbox, 1)
        x_stop = min(x + halfbox, xmax)
        for y_lp, y in enumerate(y_sub):
            y += 1
            y_start = max(y - halfbox, 1)
            y_stop = min(y + halfbox, ymax)
            blc = [x_start, y_start, 1, 1, 1]
            trc = [x_stop, y_stop, 1, 1, 1]
            p = Image.PReadPlane(data, err, blc=blc, trc=trc)
            sdarr[x_lp, y_lp] = p.RMS
            meanarr[x_lp, y_lp] = p.Mean

def makemeanmedmap(infile):
    """
    Generate a mean and a median map from infile.
    """
    inname = os.path.splitext(os.path.split(infile)[1])[0]
    meanfile = inname + '_MEAN.fits'
    rmsfile = inname + '_RMS.fits'
    if os.path.exists(meanfile) and os.path.exists(rmsfile):
        return meanfile, rmsfile
    fitsfile = fits.open(infile)
    data = fitsfile[0].data[0, 0].astype(np.float)
    mask = np.where(np.isnan(data))
    data = np.nan_to_num(data)
    sddata = np.zeros(data.shape, dtype=np.float)
    meandata = np.zeros(data.shape, dtype=np.float)
    x_sub = range(0, data.shape[0], 16)
    y_sub = range(0, data.shape[1], 16)
    x_coord, y_coord = np.meshgrid(x_sub, y_sub)
    x_coord, y_coord = x_coord.flatten(), y_coord.flatten()
    sddata_sub = np.zeros((len(x_sub), len(y_sub)), dtype=np.float)
    meandata_sub = np.zeros((len(x_sub), len(y_sub)), dtype=np.float)
    medfilt(data, 551, x_sub, y_sub, sddata_sub, meandata_sub)
    f = interpolate.RectBivariateSpline(x_sub, y_sub, sddata_sub)
    sddata = f(range(data.shape[0]),range(data.shape[1]))
    sddata[mask] = np.nan
    fitsfile[0].data[0, 0] = sddata
    fitsfile[0].data = fitsfile[0].data[0:1, 0:1, :, :]
    fitsfile.writeto(rmsfile, overwrite=True)
    f = interpolate.RectBivariateSpline(x_sub, y_sub, meandata_sub)
    meandata = f(range(data.shape[0]),range(data.shape[1]))
    meandata[mask] = np.nan
    fitsfile[0].data[0:1, 0:1] = meandata
    fitsfile.writeto(meanfile, overwrite=True)
    return meanfile, rmsfile



