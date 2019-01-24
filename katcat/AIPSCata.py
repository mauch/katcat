#from AIPS import AIPS

from AIPSLiteTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry import AIPSData as Wizardry
import OErr,Image,OTObit,re,History,ImageUtil
from math import *
import numpy
from scipy import interpolate
import pyfits
import os
import matplotlib
import string,random
import pylab


#Find a key anywhere in a header and return its 'value'
#Assume key/value pair is key = value
def get_key_thorough (header,key):
    for card in header.cards:
        data = str(card)
        if key in data:
            data=data.split()
            try:
                return data[data.index(key)+2]
            except:
                return data[data.index(key+'=')+1]
    return None

# 
# Script to extract planes from MFImage cube
def Extract (inIm, outName, err):
    """
    Extract planes of an MFImage cube to individual fits images

    * inIm     = input image filename
    * outName  = root of output fits image, disk=0
                 each plane will be outName.planenn.fits
    * err      = Obit error/message stack

    Returns the number of planes created
    """
    # Get freq info
    Pyfits_Im = pyfits.open(inIm)
    h = Pyfits_Im[0].header
    nterm = int(get_key_thorough(h,"NTERM"))
    nspec = int(get_key_thorough(h,"NSPEC"))
    freq = []
    for ifreq in range(0,nspec):
        t = "FREQ%4.4d"%(ifreq+1)
        freq.append(float(get_key_thorough(h,t)))

    #Loop over planes extracting
    #Use Obit to extract planes
    for iplane in range(nterm,nspec):
        # Create output
        outFile = outName+".plane%2.2d.fits"%(iplane-nterm+1)
        #Get the required data
        plane_data = Pyfits_Im[0].data[0,iplane]
        new_hdu = pyfits.PrimaryHDU(data=plane_data[numpy.newaxis,numpy.newaxis,:],header=h)
        new_hdu.header['CTYPE3'] = 'FREQ'
        new_hdu.header['NAXIS3'] = 1
        new_hdu.header['CRVAL3'] = freq[iplane-nterm]
        new_hdulist=pyfits.HDUList([new_hdu])
        new_hdulist.writeto(outFile,clobber=True)

    return nspec-nterm+1
    # end loop
    # end Extract

#
# Primary beam correction for single image at single frequency
#
def pbcorr_obit(image_filename,output_image,err,mingain=0.05,scratch_dir='/tmp'):

    #Load the original image into Obit
    im = Image.newPFImage("in",image_filename ,0 ,True , err)
    #Construct an output Obit Image
    im_pbcorr = Image.newPFImage("wt" , scratch_dir+'/crop_pb.fits' ,0 ,False , err)
    im.Clone(im_pbcorr, err)
    #PBCOrrect
    #ImageUtil.PPBCorr(im, im, im_pbcorr, err, antSize=12.0, PBmin=0.1)
    ImageUtil.PPBImage(im,im_pbcorr,err,minGain=mingain,antSize=13.5)
    inimage=pyfits.open(image_filename)
    beam=pyfits.open(scratch_dir+'/crop_pb.fits')
    pbcorrdata=inimage[0].data/beam[0].data
    new_hdu = pyfits.PrimaryHDU(data=pbcorrdata,header=inimage[0].header)
    pyfits.HDUList([new_hdu]).writeto(output_image,clobber=True)
    im_pbcorr = Image.newPFImage("wt" , output_image ,0 ,False , err)
    #cleanup
    os.remove(scratch_dir+'/crop_pb.fits')
    return im_pbcorr

#
# Primary beam correction for multi plane image (Work on fits image)
#
def pbcorr_multiplane(image_filename,output_image,err,scratch_dir='/tmp/'):

    nplanes=Extract(image_filename,scratch_dir+"/TEMP",err)
    basename = os.path.splitext(os.path.basename(image_filename))[0]
    #Got number of planes, now split them out and primary beam correct them.
    imArr = []
    for plane in range(nplanes-1):
        imArr.append(pbcorr_obit(scratch_dir+"/TEMP.plane%2.2d.fits"%(plane+1),scratch_dir+"/PBCORR.plane%2.2d.fits"%(plane+1),err,scratch_dir=scratch_dir))
        #Delete the image
        os.remove(scratch_dir+"/TEMP.plane%2.2d.fits"%(plane+1))
        #imArr.append(Image.newPFImage("P%2.2d"%(plane+1),"TEMP.plane%2.2d.fits"%(plane+1),0,True,err))
    output = Image.newPFImage("out" ,scratch_dir+'/pbcor_'+basename+'.fits' ,0 ,False , err)
    import SpectrumFit
    sf=SpectrumFit.PCreate("Fitter",2)
    #Form final image by fitting the spectrum in the pb corrected planes at each pixel.
    sf.ImArr(imArr,output,err)
    #Crop the image at the 15% level of the primary beam
    crop_pb(scratch_dir+'/pbcor_'+basename+'.fits',output_image,err,minlevel=0.05,scratch_dir=scratch_dir)
    #cleanup
    for plane in range(nplanes-1):
        os.remove(scratch_dir+"/PBCORR.plane%2.2d.fits"%(plane+1))
    os.remove(scratch_dir+'/pbcor_'+basename+'.fits')

#
# Crop an image at a given perventage of a primary beam pattern
#
def crop_pb(image,outimage,err,minlevel=0.15,scratch_dir='/tmp/'):
    #Crop the image at the 15% level of the primary beam
    wtim=Image.newPFImage("weight" ,scratch_dir+'/crop_pb.fits' ,0 ,False , err)
    im=Image.newPFImage("orig" ,image ,0 ,False , err)
    im.Clone(wtim,err)
    ImageUtil.PPBImage(im,wtim,err,minGain=minlevel)
    #Crop using pyfits
    output=pyfits.open(image)
    wtim=pyfits.open(scratch_dir+'/crop_pb.fits')
    output[0].data[0,0]=output[0].data[0,0]*wtim[0].data[0,0]/wtim[0].data[0,0]
    output.writeto(outimage,clobber=True)
    #cleanup
    os.remove(scratch_dir+'/crop_pb.fits')


# Get image plane
#
def getplane(imagedata,plane=1):

	imagedata2=imagedata[:]
	imagedata2[3]=imagedata2[3]+1
	
	image=makeaipsimage(imagedata)
	image2=makeaipsimage(imagedata2)

	subim = AIPSTask('subim')
	subim.indata=image
	subim.outdata=image2
	subim.blc=AIPSList([0,0,plane])
	subim.trc=AIPSList([0,0,plane])
	subim.go()
	image.zap()

	return imagedata2

#
# Get Antenna diameters
#
def getdiameter(tele='GMRT'):
    """
    return the tlescope diameter in meter
    """
    diameter = 0
    #
    if (tele =='GMRT'):
	# GMRT 30 antennas 45 m diameter
        diameter = 45
    if (tele =='WSRT'):
	# WSRT 14 antennas 25 m diameter
        diameter = 25
    if (tele =='VLA'):
	# VLA  27 antennas 25 m diameter
        diameter = 25
    if (tele =='VLBA'):
	# VLBA  10 antennas 25 m diameter
        diameter = 25
    if tele == 'ATCA':
	# ATCA  6 antennas 22 m diameter
	diameter = 22.
    if tele == 'KAT-7':
        # KAT-7 7 antennas 12 m diameter 
        diameter = 12.
    if tele == 'MeerKAT':
        diameter = 13.5
    #AUTO KAT-7
    diameter = 12.
    return(diameter)

#
# Get the field of view
#
def fov(freq,antenna='None'):
    """
    freq in Hz
    returns a vector of FOV (deg, arcmin,arcsec)
    """
    light_c=299792458.0        # m s-1
    wl=light_c/freq

    baseline = getdiameter(antenna)
    # half power beam width = 1.02 * 
    # first null              1.22 *
    
    fwhm=1.22*wl/baseline*360.0*60.*60.*1000./(2.0*pi)

    return(fwhm/(1000.*60.*60.),fwhm/(1000.*60.),fwhm/1000.)





#
#Return a random filename
#
def random_filename(chars=string.hexdigits, length=16, prefix='', suffix='', verify=True, attempts=10):
    for attempt in range(attempts):
        filename = ''.join([random.choice(chars) for i in range(length)])
        filename = prefix + filename + suffix
        if not verify or not os.path.exists(filename):
            return filename



#
# Get the range of RA and Dec spanned by an image.
#
def imrange(image):
	"""
	gives back the image range in degrees
	"""
	imhead=image.header()
	ra     = [imhead['crval'][0],imhead['crpix'][0],imhead['cdelt'][0],imhead['inaxes'][0]]
	dec    = [imhead['crval'][1],imhead['crpix'][1],imhead['cdelt'][1],imhead['inaxes'][1]]
	pramax  = ra[0] - (ra[3]-ra[1]) * ra[2]
	pramin  = ra[0] + (ra[3]-ra[1]) * ra[2]
	if pramin > pramax:
		ramax = pramin
		ramin = pramax
	else:
		ramax = pramax
		ramin = pramin

	pdecmin = dec[0] - (dec[3]-dec[1]) * dec[2]
	pdecmax = dec[0] + (dec[3]-dec[1]) * dec[2]
	if pdecmin > pdecmax:
		decmax = pdecmin
		decmin = pdecmax
	else:
		decmax = pdecmax
		decmin = pdecmin
	
	return(ramin,ramax,decmin,decmax)


#
# Given an aips data structure ['imname','imtype',disk,seq],
# setup an aipsimage or an aips wizardry image object.
#
def makeaipsimage(data):

        image = AIPSImage(data[0],data[1],data[2],data[3])
	return image

def makeobitimage(data,err, check=True):
	
	image = Image.newPAImage(data[0],data[0],data[1],data[2],data[3],check,err)
	return image
#
# Given an x,y pixel position in an aips image (aipsdata),
# bilinearlly interpolate the value of the image pixel at
# x,y.
#
def pixinterpold(pixx,pixy,imagedata):

    err=OErr.OErr()
    image = makeobitimage(imagedata,err)
    image.Open(Image.READONLY,err)
    image.Read(err)
    pixels=image.FArray
    
    values=numpy.zeros((25,))
    points=numpy.zeros((25,2,))
    nloop=-1
    for x in range(int(pixx-3),int(pixx+2)):
        for y in range(int(pixy-3),int(pixy+2)):
            nloop+=1
            values[nloop] = pixels.get(x,y)
            points[nloop] = [x,y]
    val=interpolate.griddata(points,values,[pixx-1,pixy-1],method='cubic')
    return val[0]

#
# Given an x,y pixel position in an aips image (aipsdata),
# bilinearlly interpolate the value of the image pixel at
# x,y.
#
def pixinterp(pixx,pixy,imagedata):
    
    values=numpy.zeros((25,))
    points=numpy.zeros((25,2,))
    nloop=-1
    for x in range(int(pixx-3),int(pixx+2)):
        for y in range(int(pixy-3),int(pixy+2)):
            nloop+=1
            values[nloop] = imagedata[y, x]
            points[nloop] = [y, x]
    val=interpolate.griddata(points,values,[pixy-1 ,pixx-1],method='cubic')
    return val[0]
#
# Fit  the peak flux of a source that has 
# known maj,min,pa and position by fitting
# ONLY the peak flux of a known elliptical gaussian
# in an image. Try JMFIT first and IMFIT second-
# return 0.0 if neither of them work.
#
def getknownflux(imagedata,maj,min,pa,x,y):
    
    image=makeaipsimage(imagedata)
    imhead=image.header()
    incrx=imhead['cdelt'][0] * 3600.0
    incry=imhead['cdelt'][1] * 3600.0
    
    if (abs(incrx)!=abs(incry)):
        print "WARNING: incrx not = incry -> undefined behaviour of getflux!!!!!"
    
    majpix=abs(maj/incrx)
    minpix=abs(min/incrx)
         
    boxhalfsize=abs((maj/incrx)/2.)+3.
    imfit=AIPSTask('jmfit')
    imfit.indata=image
    imfit.blc=AIPSList([int(numpy.floor(abs(x-(boxhalfsize+1)))),int(numpy.floor(abs(y-(boxhalfsize+1))))])
    imfit.trc=AIPSList([int(numpy.ceil(abs(x+(boxhalfsize+1)))),int(numpy.ceil(abs(y+(boxhalfsize+1))))])
    imfit.ngauss=1
    imfit.ctype=AIPSList([1])
    imfit.gmax=AIPSList([pixinterp(x,y,imagedata)])
    imfit.gpos=AIPSList([[x,y]])
    imfit.gwidth=AIPSList([[majpix,minpix,pa]])
    imfit.domax=AIPSList([1])
    imfit.dopos=AIPSList([[1,1]])
    # JMFIT always fails if the PA is held fixed- so allow it to vary.
    imfit.dowidth=AIPSList([[1,1,1]])
    imfit.bwsmear=0
    imfit.radius=boxhalfsize*100
    imfit.niter=2000
    #imfit.docrt=1
    imfit.fitout=''
    imfit.dooutput=-1
    imfit.offset=0
    imfit.domodel=-1
    imfit.outvers=-1
    imfit.pbparm=AIPSList([0])
    try:
        imfit.go()
        flux=imfit.fmax[1]
        sizehi = imfit.fwidth[1][1]
        sizelo = imfit.fwidth[1][2]
        sizepa = imfit.fwidth[1][3]
        return flux,-1.*sizehi*incrx,sizelo*incry,sizepa
    except:
        print "JMFIT failed- try IMFIT!"
        imfit2=AIPSTask('imfit')
        imfit2.indata=image
        imfit2.blc=AIPSList([int(numpy.floor(abs(x-(boxhalfsize+1)))),int(numpy.floor(abs(y-(boxhalfsize+1))))])
        imfit2.trc=AIPSList([int(numpy.ceil(abs(x+(boxhalfsize+1)))),int(numpy.ceil(abs(y+(boxhalfsize+1))))])
        imfit2.ngauss=1
        imfit2.ctype=AIPSList([1])
        imfit2.gmax=AIPSList([0])
        imfit2.gpos=AIPSList([[x,y]])
        imfit2.gwidth=AIPSList([[majpix,minpix,pa]])
        imfit2.domax=AIPSList([1])
        imfit2.dopos=AIPSList([[0,0]])
        imfit2.dowidth=AIPSList([[0,0,1]])
        imfit2.bwsmear=0
        imfit2.radius=boxhalfsize*100
        imfit2.niter=2000
        #imfit2.docrt=1
        imfit2.fitout=''
        imfit2.dooutput=-1
        imfit2.offset=0
        imfit2.domodel=-1
        imfit2.outvers=-1
        imfit2.pbparm=AIPSList([0])
        try:
            imfit2.go()
            flux=imfit2.fmax[1]
        except:
            print "IMFIT has failed as well!!!!! Returning Flux=0.0"
            flux=0.0
      
    return flux

#
# Add a gaussian to an image defined in imagedata
# at position xpix,ypix with amp and width=[maj,min,pa].
# Return the new imagedata (optionally destroy the input image).
# 
def addsource(xpix,ypix,amp,width,imagedata,zap=0):

    newimagedata = imagedata[:]
    newimagedata[3] += 1
    image=makeaipsimage(imagedata)
    newimage=makeaipsimage(newimagedata)
    immod = AIPSTask('immod')
    immod.indata = image
    immod.outdata = newimage
    immod.opcode = 'GAUS'
    immod.ngaus = 1
    immod.fmax = AIPSList([float(amp),0])
    immod.fpos = AIPSList([[float(xpix),float(ypix)]])
    immod.fwidth = AIPSList([width])
    immod.flux = 0
    immod.factor = 1
    immod.go()
    if zap==1: image.zap()
    return(newimagedata)

#
# Load a FITS image from PWD to the AIPS disk with porperties
# imagedata (use Obit).
#
def loadfitspwd(infilename,imagedata,user,err):

    #inimage=Image.newPFImage(imagedata[0],infilename,0,True, err)
    #outimage=Image.newPAImage(imagedata[0],imagedata[0],imagedata[1],imagedata[2],imagedata[3],False,err)
    #Image.PCopy(inimage,outimage,err)
    
    return OTObit.imlod(infilename,0,imagedata[0],imagedata[1],imagedata[2],imagedata[3],err)


#    image = makeaipsimage(imagedata)
#    fitld = AIPSTask('fitld')
#    fitld.douvcomp = -1
#    fitld.digicor  = -1
#    fitld.datain   = 'PWD:'+infilename
#    fitld.outdata  = image
#    fitld.msgkill  = -2
#    fitld.bchan    = 1
#    fitld.echan    = 1
    
#    fitld.go()

#
# Do a primary beam correction on an image using the
# pbparm array. Create outimage. inverse=1 means
# remove the primary beam.
#
def pbcorr(inimagedata,outimagedata,pbparm,inverse=0):

    pbaipsparm=AIPSList(pbparm)
    inimage=makeaipsimage(inimagedata)
    outimage=makeaipsimage(outimagedata)
    pbcor          = AIPSTask('pbcor')
    pbcor.outdata  = outimage
    pbcor.indata   = inimage
    pbcor.doinvers = inverse
    pbcor.pbparm   = pbaipsparm
    pbcor.go()    
    return(outimage)

#
# Return the rms in an image region using IMEAN in aips.
# Keep on iterating estimates until mean estimate changes by
# 0.1%(default) or 10(default) iterations have occured.
#
def getrms(inimagedata,blc,trc,tolerance=0.001,maxiter=100):

    inimage=makeaipsimage(inimagedata)
    imean          = AIPSTask('imean')
    imean.indata   = inimage
    imean.blc      = AIPSList(blc)
    imean.trc      = AIPSList(trc)
    imean.doinvers = -1
    imean.nboxes   = 0
    imean.pixrange = AIPSList([0,0])
    imean.functype = ''
    imean.pixavg   = 0
    imean.pixstd   = 0
    imean.docat    = 0
    imean.dotv     = -1
    imean.grchan   = 0
    imean.go()

    #Set the tolerance to be a fraction of current measured average.
    tolerance=fabs(imean.pixavg*tolerance)
    diff=tolerance+1.0
    count=0
    while diff>tolerance and count<maxiter:
        prevmean=imean.pixavg
        imean.go()
        diff=fabs(imean.pixavg - prevmean)
        count+=1
    if count==maxiter:
        print 'No convergence in IMEAN after '+str(maxiter)+' iterations!'

    return fabs(imean.pixstd)

#
# Return the mean in an image region using IMEAN in aips.
# Keep on iterating estimates until mean estimate changes by
# 0.1%(default) or 10(default) iterations have occured.
#
def getmean(inimagedata,blc,trc,tolerance=0.0001,maxiter=50):

    inimage=makeaipsimage(inimagedata)
    imean          = AIPSTask('imean')
    imean.indata   = inimage
    imean.blc      = AIPSList(blc)
    imean.trc      = AIPSList(trc)
    imean.doinvers = -1
    imean.nboxes   = 0
    imean.pixrange = AIPSList([0,0])
    imean.functype = ''
    imean.pixavg   = 0
    imean.pixstd   = 0
    imean.docat    = 0
    imean.dotv     = -1
    imean.grchan   = 0
    imean.go()

    #Set the tolerance to be a fraction of current measured average.
    tolerance=fabs(imean.pixavg*tolerance)
    diff=tolerance+1.0
    count=0
    while diff>tolerance and count<maxiter:
        prevmean=imean.pixavg
        imean.go()
        diff=fabs(imean.pixavg - prevmean)
        count+=1
    if count==maxiter:
        print 'No convergence in IMEAN after '+str(maxiter)+' iterations!'

    return fabs(imean.pixavg)

#
# Wrapper for the task math in AIPS.
#
def immath(inimagedata,outimagedata,opcode,cparm):

    inimage=makeaipsimage(inimagedata)
    outimage=makeaipsimage(outimagedata)

    aipsmaths         = AIPSTask('maths')
    aipsmaths.indata  = inimage
    aipsmaths.outdata = outimage
    aipsmaths.opcode  = opcode
    aipsmaths.cparm   = AIPSList(cparm)
    aipsmaths.go()

    return outimage 

# Stolen from Cotton!
def EVLAImFITS(inImage, filename, outDisk, err, fract=None, quant=None, \
          exclude=["AIPS HI","AIPS PL","AIPS SL"], include=["AIPS CC"],
          headHi=False, logfile=""):
    """
    Write AIPS image as FITS
    
    Write a Image data set as a FITAB format file
    History also copied

    * inImage    = Image data to copy
    * filename   = name of FITS file, any whitespace characters replaced with underscore
    * outDisk    = FITS directory number
    * err        = Python Obit Error/message stack
    * fract      = Fraction of RMS to quantize
    * quant      = quantization level in image units, has precedence over fract
      None or <= 0 => use fract.
    * exclude    = List of table types NOT to copy
      NB: "AIPS HI" isn't really a table and gets copied anyway
    * include    = List of table types to copy
    * headHi     = if True move history to header, else leave in History table
    """
    ################################################################
    #
    # Checks
    if not Image.PIsA(inImage):
        raise TypeError,"inImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    # Deblank filename
    fn = re.sub('\s','_',filename)
    # Set output
    outImage = Image.newPFImage("FITS Image DATA", fn, outDisk, False, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating FITS data")
    # Check for valid pixels
    if inImage.Desc.Dict["maxval"]<=inImage.Desc.Dict["minval"]:
        fract=None; quant=None
    # Copy
    if fract or quant:
        Image.PCopyQuantizeFITS (inImage, outImage, err, fract=fract, quant=quant)
    else:
        Image.PCopy (inImage, outImage, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying Image data to FITS")
    # Copy History
    inHistory  = History.History("inhistory",  inImage.List, err)
    outHistory = History.History("outhistory", outImage.List, err)
    History.PCopy(inHistory, outHistory, err)
    # Add this programs history
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit imtab",err)
    if fract:
        outHistory.WriteRec(-1,"imtab   / Quantized at "+str(fract)+" RMS",err)
    outHistory.WriteRec(-1,"imtab   / FITS file "+fn+", disk "+str(outDisk),err)
    outHistory.Close(err)
    # History in header?
    if headHi:
        # Copy back to header
        inHistory  = History.History("inhistory",  outImage.List, err)
        History.PCopy2Header (inHistory, outHistory, err)
        # zap table
        outHistory.Zap(err)
    OErr.printErrMsg(err, "Error with history")
    # Copy Tables
    Image.PCopyTables (inImage, outImage, exclude, include, err)
    del outImage
    # end EVLAImFITS

#
# Wrapper for the AIPS task imsad:
# inimage:   input image to catalogue
# cutoff:    cutoff level for gaussians
# outname:   name of sad output file (written to PWD)
# highcut:   list of upper cutoff levels above cutoff for gaussians
# blc:       blc of fit region
# trc:       trc of fit region
# gain:      AIPS gain
# icut:      retry level for AIPS
# dowidth:   fit widths (>0 => yes)
# maxwidth:  maximum size if fitting widths ( 0 => none)
# resid:     keep residual image (>0 => yes)
# residdata: data for residual image if resid>0
#
def imsad(indata,cutoff,outname='',uppercut=[],blc=[0,0],trc=[0,0],gain=0.3,icut=0.1,dowidth=1,maxwidth=5,resid=0,residdata=[]):

    task             = AIPSTask('sad')
    inimage = makeaipsimage(indata)
    task.indata      = inimage
    task.blc         = AIPSList(blc)
    task.trc         = AIPSList(trc)
    if len(uppercut)>0:
        task.cparm   = AIPSList(uppercut+[cutoff,0])
    else:
        task.cparm   = AIPSList([cutoff,0])
    task.doresid     = resid
    if resid>0:
        residimage   = makeaipsimage(residdata)
        task.outdata = residimage
    task.ngauss      = 80000
    task.icut        = icut
    task.sort        = 'RA'
    if outname=='':
        task.docrt   = 0
        task.fitout  = ''
    else:
        task.docrt   = 132
        task.fitout  = outname
    task.outvers     = 0
    task.doall       = 1
    if dowidth==1:
        width = [[1,1,1],[1,1,1],[1,1,1],[1,1,1]]
    else:
        width = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    task.dowidth     = AIPSList(width)
    task.gain        = gain
    task.dparm       = AIPSList([cutoff, 0, 0.0, maxwidth, -1, 0, cutoff, 0, 0])
    task.prtlev      = 0
    task.pbparm      = AIPSList([0,0,0,0,0,0,0])
    task.go()

#
# Wrapper for the AIPS task imsad with an rms image:
# inimage:   input image to catalogue
# inrmsimage: input rms map for imsad
# cutoff:    cutoff level for gaussians
# outname:   name of sad output file (written to PWD)
# highcut:   list of upper cutoff levels above cutoff for gaussians
# blc:       blc of fit region
# trc:       trc of fit region
# gain:      AIPS gain
# icut:      retry level for AIPS
# dowidth:   fit widths (>0 => yes)
# maxwidth:  maximum size if fitting widths ( 0 => none)
# resid:     keep residual image (>0 => yes)
# residdata: data for residual image if resid>0
#
def imsadrms(indata,inrmsimage,cutoff,outname='',uppercut=[],blc=[0,0],trc=[0,0],gain=0.3,icut=0.1,dowidth=1,maxwidth=5,resid=0,residdata=[]):

    task             = AIPSTask('sad')
    inimage = makeaipsimage(indata)
    in2image = makeaipsimage(inrmsimage)
    task.indata      = inimage
    task.in2data     = in2image
    task.blc         = AIPSList(blc)
    task.trc         = AIPSList(trc)
    if len(uppercut)>0:
        task.cparm   = AIPSList(uppercut+[cutoff,0])
    else:
        task.cparm   = AIPSList([cutoff,0])
    task.doresid     = resid
    if resid>0:
        residimage   = makeaipsimage(residdata)
        task.outdata = residimage
    task.ngauss      = 80000
    task.icut        = icut
    task.sort        = 'RA'
    if outname=='':
        task.docrt   = 0
        task.fitout  = ''
    else:
        task.docrt   = 132
        task.fitout  = outname
    task.outvers     = 0
    task.doall       = 1
    if dowidth==1:
        width = [[1,1,1],[1,1,1],[1,1,1],[1,1,1]]
    else:
        width = [[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
    task.dowidth     = AIPSList(width)
    task.gain        = gain
    task.dparm       = AIPSList([cutoff, 0, 0.0, maxwidth, -1, 0, cutoff, 0, 2, 0])
    task.prtlev      = 0
    task.pbparm      = AIPSList([0,0,0,0,0,0,0])
    task.go()

#
# Wrapper for the Obit task findsou:
# inimage:   input image to catalogue
# cutoff:    cutoff level for gaussians
# outname:   name of output file (written to PWD)
# highcut:   list of upper cutoff levels above cutoff for gaussians
# blc:       blc of fit region
# trc:       trc of fit region
# gain:      gain
# icut:      retry level
# dowidth:   fit widths (>0 => yes)
# maxwidth:  maximum size if fitting widths ( 0 => none)
# resid:     keep residual image (>0 => yes)
# residdata: data for residual image if resid>0
#
def findsou(indata,cutoff,err,outname='',uppercut=[],blc=[0,0],trc=[0,0],gain=0.3,icut=0.1,dowidth=1,maxwidth=5,resid=0,residdata=['RESID', 'IM', 1, 1]):

    task             = OTObit.ObitTask('FndSou')
    task._max_dict['NGauss'] = 100000.0
    inimage          = makeobitimage(indata, err)
    OTObit.setname(inimage, task)
    task.BLC         = blc
    task.TRC         = trc
    task.CutOff      = cutoff
    task.doResid     = bool(resid)
    if resid>0:
        residimage   = makeobitimage(residdata, err, check=False)
        OTObit.setoname(residimage, task)
    task.NGauss      = 100000
    task.NPass       = 1
    task.Retry       = icut
    task.Sort        = 'X'
    if outname=='':
        task.OutPrint  = ''
    else:
        task.OutPrint  = outname
    task.doMult       = True
    task.doWidth     = bool(dowidth)
    task.Gain        = gain
    task.Parms       = [cutoff, float(maxwidth), 2.0, 1.0, 2.0, -1, 0.1]
    task.prtLv      = 1
    task.nThreads   = 8
    task.go()
#
# Copy an aips cataloge image from arg1 to arg2
#
def copyimage(imagedata1,imagedata2):

    image1=makeaipsimage(imagedata1)
    image2=makeaipsimage(imagedata2)

    copy         = AIPSTask('subim')
    copy.indata  = image1
    copy.outdata = image2
    copy.blc     = AIPSList([0.0])
    copy.trc     = AIPSList([0,0])
    copy.xinc    = 1
    copy.yinc    = 1
    copy.opcode  = ''
    copy.go()

#
# Wrapper for the AIPS task rmsd
# takes inimage,outimage,
# boxsize- the size of the box to compute the rms in pixels.
# opcode- 'RMS' or 'MEAN' or 'BLAN'
# xinc,yinc- the amount of increment to use per step.
#
def rmsd(inimagedata,outimagedata,boxsize,opcode,xinc,yinc,scalr1=0,scalr2=0,scalr3=0,flux=0):

    inimage=makeaipsimage(inimagedata)
    outimage=makeaipsimage(outimagedata)

    rmsd         = AIPSTask('rmsd')
    rmsd.indata  = inimage
    rmsd.outdata = outimage
    rmsd.blc     = AIPSList([0,0])
    rmsd.trc     = AIPSList([0,0])
    rmsd.imsize  = AIPSList(boxsize)
    rmsd.opcode  = opcode
    rmsd.optype  = ''
    rmsd.scalr1  = scalr1
    rmsd.scalr2  = scalr2
    rmsd.scalr3  = scalr3
    rmsd.flux    = flux
    rmsd.xinc    = int(xinc)
    rmsd.yinc    = int(yinc)
    rmsd.go()

    return outimage

#
# Convienience task to write out a FITS file from an AIPS image
#
def fitsout(inimagedata,outfile):
    
    err = OErr.OErr()

    x = Image.newPAImage(inimagedata[0], inimagedata[0], inimagedata[1], inimagedata[2], inimagedata[3], True, err)
    
    #inimage=makeaipsimage(inimagedata)
    if os.path.exists(outfile):
        print 'Cannot create '+outfile
        raise IOError
    xf = EVLAImFITS (x, outfile, 0, err)
    #fittp         = AIPSTask('fittp')
    #fittp.indata  = inimage
    #fittp.dataout = outfile
    #fittp.go()

    # Correct an AIPS bug with NaN pixels.
    outimage=pyfits.open(outfile,mode='update')
    #print numpy.nanmin(outimage[0].data)
    outimage[0].header.set('datamin',numpy.nanmin(outimage[0].data))
    outimage[0].header.set('datamax',numpy.nanmax(outimage[0].data))
    outimage.flush()
    outimage.close()
                             

#
# Find the maximum pixel in an AIPS image. Do it using
# pyfits which is more reliable than AIPS at dealing
# with blanked pixels.
#
def getimmax(inimagedata,outputdir):

    filename=outputdir+random_filename(length=8,suffix='.FITS')
    fitsout(inimagedata,filename)
    fitsimage=pyfits.open(filename)
    maxpix=numpy.nanmax(fitsimage[0].data)
    fitsimage.close()
    os.remove(filename)
    return maxpix

#
# Combine two image using the aips task comb.
#
def imcombine(in1imagedata,in2imagedata,outimagedata,
              opcode,aparm,align=1,bparm=[0,0,0,0,0,0]):

    in1image=makeaipsimage(in1imagedata)
    in2image=makeaipsimage(in2imagedata)
    outimage=makeaipsimage(outimagedata)
    comb         = AIPSTask('comb')
    comb.indata  = in1image
    comb.in2data = in2image
    comb.outdata = outimage
    comb.doalign = align
    comb.blc     = AIPSList([0,0])
    comb.trc     = AIPSList([0,0])
    comb.opcode  = opcode
    comb.aparm   = AIPSList(aparm)
    comb.bparm   = AIPSList(bparm)
    comb.go()
    return(outimage)

#
# Get the maximum pixel in an aips defined image.
# pb=Remove a primary beam correction?
# returns the maximum pixel in image units.
#
def maxpix(imagedata,pb=0):

    image=makeaipsimage(imagedata)
    imhead=image.header()
    obsfreq=imhead['crpix'][2]
    telescope=imhead['teles']
    
    outimagedata = imagedata[:]
    if pb>0:
        pbeamparm = [0.001,1] + list(primarybeamcor(obsfreq,telescope))
        outname=random_filename(length=8,chars=string.ascii_uppercase)
        outimagedata[0] = outname
        pbcorr(imagedata,outimagedata,pbeamparm,1)

    outimage=makeaipsimage(outimagedata)
    outhead=outimage.header()
    return(outhead['maxval'])

# Make a 'radio_source' array out of raw sad output
# (which is a 'file' object).
def makeradioarray(sadoutput):

    out_array=[]
    for source in sadoutput.data:
        num=int(str(source).split()[2])
        ra=float(str(source).split()[0])
        dec=float(str(source).split()[1])
        thispos=position(ra,dec,'dd')
        pflux=float(str(source).split()[10])
        flux=float(str(source).split()[12])
        fmaj=float(str(source).split()[16])
        fmin=float(str(source).split()[17])
        fpa=float(str(source).split()[18])
        out_array.append(radio_source(num,thispos,pflux,flux,[fmaj,fmin,fpa]))

    return(out_array)

def func(x,a,b):
    return a*x + b
 

# Iteratively work out the image dynamic range- do this by going to the
# positions in 'brightsources' and checking first that the source is an isolated
# point source- if this source is, then look around the source (to some fraction
# 'boxsize' of the primary beam), in increasing rms steps until all the nearby
# peaks are removed. The function will return a dynamic range determination
# for each bright source in 'sources'
def getdr(brightcatalogue,imagedata,rms,err,blc=[0,0],trc=[0,0]):

    from CATparams import *
    from catalogue import *    
    from scipy.optimize import curve_fit
    from scipy import interpolate
        
    image=makeaipsimage(imagedata)
    imhead=image.header()
    beammajpix=imhead['beamMaj']*imhead['cdelt'][0]
    obsfreq=imhead['crval'][2] if imhead['ctype'][2] in ['FREQ     ','SPECLNMF '] else imhead['crval'][3]
    telescope=imhead['teles']
    
    # Get the size in pixels to search around each bright source.
    widtharcsec = fov(obsfreq,telescope)[2]*localwidth
    widthpixels = widtharcsec/(fabs(imhead['cdelt'][0])*3600.0)

    # Run IMSAD with cutoff=3 to get initial source estimates.
    sadfilename = random_filename(length=8,suffix='.sad')
    imsad(imagedata,4.0*rms,uppercut=[8.0*rms],
          outname=sadfilename,maxwidth=maxfitwidth*beammajpix,gain=gain,icut=highcut*ifact,blc=blc,trc=trc)
    sadfile = open(sadfilename, 'r')
    allsources = file(sadfile)
    sadfile.close()
    os.remove(sadfilename)

    allxpix = [float(str(source).split()[-4]) for source in allsources.data]
    allypix = [float(str(source).split()[-3]) for source in allsources.data]
    allposns = zip(allxpix,allypix)
    allflux = [float(str(source).split()[10]) for source in allsources.data]
    alldr=[]
    allradius=[]
    # Get the surface density of 3sigma sources in the image.
    # Number of x and y pixels in the image- take 95% of total image area.
    pixradius = (min(imhead['inaxes'][0:2])*0.95)/2.0
    # count the number of 3sigma sources in circle spanned by pixdiameter.
    # image centerpixel
    allbins = numpy.arange(0,pixradius,widthpixels)
    alloffsets = [pixoffset(*(list(posn[0:2]) + list(imhead['crpix'][0:2]))) for posn in allposns]
    alldata,allbins=numpy.histogram(numpy.array(alloffsets),numpy.array(allbins))
    allareas=(2.0*acos(0.0)*allbins[1:]**2.0)-(2.0*acos(0.0)*(allbins[:-1])**2.0)
    alldataperarea=alldata/allareas
    try:	
        fitp,fitc=curve_fit(func,(numpy.array(allbins[1:])-(widthpixels/2.0))[numpy.nonzero(alldata)],numpy.log10(alldataperarea[numpy.nonzero(alldata)]))    
    except:
        fitp=[0,0]
    for rawbrightsource,brightsource in zip(brightcatalogue.rawdata,brightcatalogue.sourcedata):
        # is it resolved and isolated
        if brightsource.isresolved(brightcatalogue.survey) and (brightsource.isisolated(brightcatalogue.sourcedata,widtharcsec)==0):            

            xpix = float(rawbrightsource[-4])
            ypix = float(rawbrightsource[-3])
            flux = float(rawbrightsource[10])
            maxpeak=0.0

            # Find the radius around this source with increase surface density
            alllocaloff=[]
            offsetcentre=pixoffset(xpix,ypix,*list(imhead['crpix'][0:2]))

            for num,thispos in enumerate(allposns):
                offset=pixoffset(xpix,ypix,thispos[0],thispos[1])
                if offset<10*widthpixels and offset > 0.01:
                    alllocaloff.append(offset)

            alllocalsorted=numpy.sort(alllocaloff)
            bins=numpy.array([0]+list(alllocalsorted[4::5]))
            areas=(2.0*acos(0.0)*bins[1:]**2.0)-(2.0*acos(0.0)*(bins[:-1])**2.0)
            dataperarea=5.0/areas

            backgrounddata=10**func(offsetcentre,*fitp)
            thiszero=numpy.nonzero((dataperarea-(3.0*backgrounddata))<0)
            
            if len(bins)==1:
                thisradius=widthpixels
            elif len(thiszero[0])==0:
                thisradius=dataperarea[-1]+((10*widthpixels)-dataperarea[-1])/2.0
                allradius=allradius+[thisradius]
            elif thiszero[0][0]>0:         #If this zero[0][0]==0 then there is no overdensity.
                firstzero=thiszero[0][0]
                xval=(dataperarea-(3.0*backgrounddata))[firstzero-1:firstzero+1]
                yval=bins[firstzero:firstzero+2]
                zeropoint=interpolate.interp1d(xval[::-1],yval[::-1])
                thisradius=zeropoint(0)
                allradius=allradius+[thisradius]
            else:
                thisradius=widthpixels


            
            # 'peaklim' sets the 5*sigma* limit above which sources can be considered "real"
            # *Sigma* is here defined in terms of a "peak" being fitted close to  a bright
            # source is a 5*sigma* noise peak.
            peaklim=(flux/mindr)*5.0
            for num,thispos in enumerate(allposns):
                offset=pixoffset(xpix,ypix,thispos[0],thispos[1])
                if offset<thisradius and offset > 0.01:
                    alllocaloff.append(offset)
                    if allflux[num]>maxpeak and allflux[num]<peaklim:
                        maxpeak=allflux[num]
            if maxpeak>0.0:
                thisdr=flux/(maxpeak/5.0)
                alldr.append(thisdr)
    allavradius=numpy.average(allradius)
    allmedradius=numpy.median(allradius)
    
    if len(alldr)>0:
        print '**** Median local dynamic range close to',str(len(alldr)),'bright sources is:',str(numpy.median(alldr))
        if len(allradius)>0:
            print '**** Median radius of artifacts around bright sources is:',str(numpy.median(allradius))
            return allmedradius,numpy.median(alldr),len(alldr)
        else:
            return beammajpix,numpy.median(alldr),len(alldr)
    
    else:
        print widthpixels*3.0
        return beammajpix,0.0,0

#Quick wrapper for the aips task zap.#
def zap(inimagedata):

    image=makeaipsimage(inimagedata)
    image.zap()


#Given an image- calculate the local rms and mean in 'rmsboxsize' steps
# and subtract the mean and divide by the local rms and return this image
def makermsimage(inimagedata,outimagedata,rmsboxsize,debug=1,stepsize=[0,0]):

    image=makeaipsimage(inimagedata)
    imhead=image.header()
    beammaj=imhead['beamMaj'] * 3600.0
    beammin=imhead['beamMin'] * 3600.0
    incr1 = imhead['cdelt'][0] * 3600.0
    incr2 = imhead['cdelt'][1] * 3600.0

    beammajpix =  abs(beammaj/incr1)
    beamminpix = abs(beammin/incr2)

    
    rmsimagename = random_filename(length=8)
    rmsimagedata = [rmsimagename, 'RMS', 1, 1]
    meanimagedata = [rmsimagename, 'MEAN', 1, 1]
    meansubimagedata = [rmsimagename, 'MSUB', 1, 1]

    if stepsize[0]==0:
        stepsize[0]=int(beammajpix*5)
    if stepsize[1]==0:
        stepsize[1]=int(beammajpix*5)
    
    rmsd(inimagedata,rmsimagedata,rmsboxsize,'MAD',*stepsize)
    rmsd(inimagedata,meanimagedata,rmsboxsize,'MEAN',*stepsize)

    imcombine(inimagedata,meanimagedata,meansubimagedata,'SUM',[1,-1,0,0,0,0,0,-1,0,0])
    imcombine(meansubimagedata,rmsimagedata,outimagedata,'DIV',[1,0,0,0,0,0,0,-1,0,0])

    if debug>0:
        try:
            outfile=random_filename()
            fitsout(meanimagedata,outfile)
            os.rename(outfile,'MEANIMAGE.FITS')
            outfiles.append('MEANIMAGE.FITS')
        except:
            print '**** Skipping MEANMAGE.FITS'
        try:
            outfile=random_filename()
            fitsout(rmsimagedata,outfile)
            os.rename(outfile,'RMSIMAGE.FITS')
            outfiles.append('RMSIMAGE.FITS')
        except:
            print '**** Skipping RMSIMAGE.FITS'
        


    zap(meansubimagedata)

    return(rmsimagedata,meanimagedata)
