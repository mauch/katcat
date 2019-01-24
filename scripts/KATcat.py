#! /usr/bin/env python
#
# Based on Tom Mauch's WSRT_CATALOGER
#
# changed:
# HRK made its AIPSLite compatible
# TXM Oct 6 2010- use pyfits and PWD.
# TXM Oct 13 2010- revamped script to use a bright source model to avoid artifacts.
# TXM Nov 14 2010- added diagnostic position accuracy (by comparison with FIRST)
#                  and flux density accuracy (by comparison with NVSS and WENSS)
#                  plots.
# TXM Nov 14 2010- Added vo output and moves all output files into OUTPUT directory.
# TXM Dec 9  2010- Added spectral index plots, added FIRST annotation file output
#                  added variation of offset with phase center distance plot,
#                  added FIRST,WENSS,NVSS catalogues to vo table output
#                  added facility to do position comparison with the NVSS catalogue
#                  added parameter input file CATparams.py
#                  fixed reversed RA,Dec axes in position offset plots
#                  added user supplied input catalogue for comarison.
#                  fixed numerous bugs related to potting
# TXM Jan 6  2011- Added option for NOT doing the primary beam correction
#                  when computing the bright source model.
# TXM Feb 1  2011- Added option to search the SUMSS catalogue for southern
#                  images.
# TXM May 9  2011- Added position offsets and flux offsets to output. Changed default FIRST
#                  position offset to 10arcsec. 
#                  Fixed porblem with SAD not deleting sources fitted with a very extended
#                  width.
#                  Added maxfitwidth to the parameters (default 5) to check what maximum
#                  fitted
#                  width to accept.
# TXM May 16 2011- Addead header information to GMRTfluxdata.txt output.
#                  Changes position offset comparison to only show point
#                  source objects with flux >3mJy. 
# TXM Jun 02 2011- Added extra statistics on flux offset for point sources in NVSS+GMRT only.
#                  Added a new plot 'SNUMS_fluxplot' to show the catalogue source numbers
#                  in the flux density comparison plot.
# TXM Jun 06 2011- Changed fitting to JMFIT in AIPS and now fit point sources with PSF.
# TXM Jun 24 2011- Fixed some JMFIT bug and tried ot make it compatible
#                  (in a reduced form) with non-GMRT data.
# TXM Aug 25 2011- COmpletely overhauled the code to make it
#                  more modular and easy to update.
#                  Many minor changes everywhere- but majer changes are the
#                  implementation of the 'catalogue' and 'external_cat'
#                  classes which manipulate the sad catalogue and the
#                  external comparson catalogues.
# TXM Sep 29 2011- Added a routine to determine the local dynamic range around bright sources by
#                  looking at the brightness of fitted peaks clost to them. Also added the offsets
#                  of the source positions in pixels.
# TXM Jan 10 2012- Fixed a bug in searchcats.py when reding the fitted data from NVSS.
#
# TXM Feb 13 2012- Added a determination of the median radius around bright sources of artifacts.
#
import os
import shutil
import subprocess
import sys
import time
import math
import warnings
import numpy as np
from time import strptime
from time import sleep
from math import *
from string import *

# Obit Setup
import OErr, OSystem, UV, AIPS, FITS, OTObit
import ObitTalkUtil

# Option parser
from optparse import OptionParser
#
# Use AIPSLite
import katcat.AIPSLite as AIPSLite
from katcat.AIPSLiteTask import AIPSTask, AIPSList
#from AIPSTask import *
#
# Ignore Warnings
#
warnings.simplefilter('ignore')
#
#VO Table code
#
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
#
#
#  FITS manipulation code
#
import pyfits
#
# Various local libraries for this software.
#
from katcat.searchcats import *
from katcat.makeplots import *
#from katcat.CATparams import *
from katcat.catalogue import *
from katcat.AIPSCata import *
import katcat.AIPSSetup as AIPSSetup
from katcat import makemeanrms
#
# Parse command line options
#
usage = "%prog [options] h5fileToImage"
description = "Make a catalogue from a KAT-7 FITS image"
parser = OptionParser( usage=usage, description=description)
parser.add_option("--parms", default=None, help="Overwrite the default cataloguing parameters using a parameter file.")
parser.add_option("--outputdir", default='./', help="Specify the output data directory.")
parser.add_option("--scratchdir", default='./', help="Specify the scratch directory.")
parser.add_option("--pbcor", action="store_true", default=False, help="Write out primary beam corrected flux densities.")
(options, args) = parser.parse_args()


#  Want to check that the listed file is present where specified.
#  And have we got the right number of arguments?
if len(args) < 1:
    parser.print_help()
    sys.exit()
infile = args[0]
#infilename = os.path.split(infile)[1]
infilename=infile
inname = os.path.splitext(os.path.split(infilename)[1])[0].replace('.IClean','')
if not os.path.exists(infile):
    print "File not found."
    sys.exit(-1)


#Get usr defined paramaters
if options.parms:
    exec(open(options.parms).read())
else:
    from katcat.CATparams import *

# Update header with beam info
fitsfile = pyfits.open(infile,mode='update')
if not 'CLEANBMJ' in fitsfile[0].header:
    fitsfile[0].header['CLEANBMJ']=fitsfile[0].header.get('BMAJ',-1)
    fitsfile[0].header['CLEANBMN']=fitsfile[0].header.get('BMIN',-1)
    fitsfile[0].header['CLEANBPA']=fitsfile[0].header.get('BPA',-1)
    fitsfile.flush()
if user_bmaj>0.:
    fitsfile[0].header['CLEANBMJ']=user_bmaj / 3600.
    fitsfile[0].header['CLEANBMN']=user_bmin / 3600.
    fitsfile[0].header['CLEANBPA']=user_bpa
    fitsfile.flush()


#check frequency header
if fitsfile[0].header['CTYPE3']=='FREQ':
    multiplane=False
else:
    multiplane=True

fitsfile.close()

# Make mean and rms fits images
print "### Generating MEAN and RMS maps ..."
meanfitsimage, rmsfitsimage = makemeanrms.makemeanmedmap(infile)
print "### MEAN and RMS maps in %s and %s."%(meanfitsimage, rmsfitsimage)

ObitSys = AIPSSetup.AIPSSetup(options.scratchdir)
# Use disk 1 always, and the AIPS number defined in CATparams
err        = OErr.OErr()
debug       = 0
userdisk    = 1
AIPS.userno = OSystem.PGetAIPSuser()
outfiles    = []           # Array to list the output files to be copied to the OUTPUT directory at the end.
scratchlist = [options.scratchdir+'/aipsdisk',options.scratchdir+'/da00']           # List of files to remove after the code has run.


#  Read in FITS file from the disk as an AIPS image.
imagedata=['CATFILE', 'FITS', 1, userdisk]
image=makeaipsimage(imagedata)
loadfitspwd(infilename,imagedata,AIPS.userno,err)

# Rean in Mean and RMS maps
meanimagedata = ['CATFILE', 'MEAN', 1, userdisk]
rmsimagedata = ['CATFILE', 'RMS', 1, userdisk]
loadfitspwd(meanfitsimage, meanimagedata, AIPS.userno, err)
loadfitspwd(rmsfitsimage, rmsimagedata, AIPS.userno, err)

rmsarray = pyfits.open(rmsfitsimage)[0].data[0,0]
meanarray = pyfits.open(meanfitsimage)[0].data[0,0]

#subim out the 1st plane to work on
imagedata=getplane(imagedata,1)
image=makeaipsimage(imagedata)
imhead=image.header()

try:
    telescope = imhead['teles']
except:
    telescope ='NONE'

#  Get image pixel information
naxis = imhead['naxis']
if naxis >= 2:
    pixx = imhead['crpix'][0]
    pixy = imhead['crpix'][1]
else:
    print 'There are not 2 or more axes in the .FITS image'
    sys.exit(-1)

# Get beam info from the aips image.
beammaj = imhead['beamMaj'] * 3600.0
beammin = imhead['beamMin'] * 3600.0
bpa = imhead['beamPA']
if beammaj<=0.0 or beammin<=0.0:
    print "**** NO BEAM INFROMATION IN HEADER!"
    print "**** QUIT!"
    sys.exit(1)

# Put some image data into well named variables.
# Image pixel increment in incr1,2; length of axis size1,2;
# center pixel cent1,2; axis info in axistype1,2.;
# observing frequency obsfreq
incr1 = imhead['cdelt'][0] * 3600.0
incr2 = imhead['cdelt'][1] * 3600.0
size1 = imhead['inaxes'][0]
size2 = imhead['inaxes'][1]
cent1 = int(size1/2.0)
cent2 = int(size2/2.0)
axistype1 = imhead['ctype'][0]
axistype2 = imhead['ctype'][1]
if 'FREQ     ' in imhead['ctype'][2] or 'SPECLNMF ' in imhead['ctype'][2]:
    freq = imhead['crval'][2]
else:
    freq = imhead['crval'][3]
obsfreq= freq
# Size of the beam in pixels.
beammajpix = abs(beammaj/incr1)
beamminpix = abs(beammin/incr2)

# Get the size of the rms box for local rms determination in pixels.
rmsboxsize = [int(rmsbox * beammajpix),-1]
brightrmsboxsize = [int(rmsbox * beammajpix)*3,-1]

# Make Primary beam image for fluxes and original image for rms values
if options.pbcor:
    #Construct a primary beam corrected image
    pbcorfitsname=options.scratchdir+'/'+random_filename(length=8,suffix='.FITS')
    if multiplane:
        pbcorr_multiplane(infile,pbcorfitsname,err,scratch_dir=options.scratchdir)
    else:
        pbcorr_obit(infile,pbcorfitsname,err,mingain=0.05,scratch_dir=options.scratchdir)
    pbcorimagedata = ['PBCOR', 'PBCOR', 1 , userdisk]
    loadfitspwd(pbcorfitsname,pbcorimagedata,AIPS.userno,err)
    #subim out the 1st plane to work on
    pbcorimagedata=getplane(pbcorimagedata,1)
    scratchlist.append(pbcorfitsname)
    #Crop the primary image to the pbcor image size
    crop_infile=options.scratchdir+'/'+random_filename(length=8,suffix='.FITS')
    crop_pb(infile,crop_infile,err,minlevel=0.1,scratch_dir=options.scratchdir)
    imagedata=['CATFILE', 'CROP', 1, userdisk]
    image=makeaipsimage(imagedata)
    loadfitspwd(crop_infile,imagedata,AIPS.userno,err)
    #subim out the 1st plane to work on
    imagedata=getplane(imagedata,1)
    scratchlist.append(crop_infile)
else:
    pbcorimagedata=imagedata

# Subtract the mean image from the image to use as the work image
meansubimagedata = ['CATFILE', 'MEANS', 1, 1]
imcombine(imagedata,meanimagedata,meansubimagedata,'SUM',[1,-1,0,0,0,0,0,-1,0,0])
imagedata = meansubimagedata

# Get the image rms in a 20% circular region in the center of the image.
rmsradius = int(min((size1,size2))*0.1)
rmsblc = [-1,rmsradius]
rmstlc = [cent1,cent2]
imagerms = getrms(imagedata,rmsblc,rmstlc)
imagemean = getmean(imagedata,rmsblc,rmstlc)

# Make global mean and rms images
globalmeandata = ['GLOBMEAN', 'GLOB', 1, 1]
globalrmsdata = ['GLOBRMS', 'GLOB', 1, 1]
immath(imagedata,globalmeandata,'POLY',[imagemean,0,0,0])
immath(imagedata,globalrmsdata,'POLY',[imagerms,0,0,0])

globalrmsarray = np.ones_like(rmsarray) * imagerms
globalmeanarray = np.ones_like(meanarray) * imagemean

# Get the maximum pixel in the pbcor image.
maxpix=getimmax(pbcorimagedata,options.outputdir+'/')
keepraw=1
modelrmsimagedata = ['RMSMODEL','RMS',1,1]


# Create an rms image and a mean image from the input image,
# using rmsd in AIPS. Box size is defined (in beams) in CATparams.py
#immath(imagedata,modelrmsimagedata,'POLY',[-imagemean,1,0,0])

# Run IMSAD with cutoff=highcut to get list of bright sources
# Bsource area
# Use 50% of the image
imarea=0.25
boffset1 = int(size1*imarea)
boffset2 = int(size2*imarea)
bblc = [cent1-boffset1,cent2-boffset2]
btrc = [cent1+boffset1,cent2+boffset2]
sadfilename = options.outputdir + '/' + random_filename(length=8,suffix='.sad')
imsad(imagedata,highcut*imagerms, outname=sadfilename,
      maxwidth=maxfitwidth*beammajpix, gain=gain, icut=highcut*ifact)
sadfile = open(sadfilename, 'r')
brightsources = file(sadfile)
sadfile.close()
scratchlist.append(sadfilename)

# Make a 'catalogue' object out of the bright sources and the global images.
brightcatalogue = sad_catalogue(infile,imagedata,globalrmsarray,
                                globalmeanarray,brightsources,maxfitwidth,imagerms)
 
# Get the image local dynamic range around the bright sources.
# Pad the image area
#pimarea=min(imarea+0.1,0.5)
#poffset1 = int(size1*imarea)
#poffset2 = int(size2*imarea)
#pblc = [cent1-boffset1,cent2-boffset2]
#ptrc = [cent1+boffset1,cent2+boffset2]
#localwidth,localdr,numsources=getdr(brightcatalogue,imagedata,imagerms,err,blc=pblc,trc=ptrc)
#numsources=0
#if numsources==0:
#localdr=maxpix/imagerms
numsources=10
# Create an rms image and a mean image from the input image,
# using rmsd in AIPS. Box size is defined (in beams) in CATparams.py
#rmsimagedata,meanimagedata = makermsimage(imagedata,[random_filename(length=8),'TMP', 1, 1],rmsboxsize,debug)
        
#rmspbcorimagedata,meanpbcorimagedata = makermsimage(pbcorimagedata,[random_filename(length=8),'TMP', 1, 1],rmsboxsize,debug)

localdr=2500.0
localwidth=100.0
# Add the bright source model to the rms image. We use a gaussian with
# width that is 3% of the primary beam of the telescope and amplitude
# given by the local dynamic range at the posoiton of
# each bright source.
widtharcsec = fov(obsfreq,telescope)[2]*localwidth
widthpixels = fabs(localwidth/2.0)
for source,catsource in zip(brightcatalogue.rawdata,brightcatalogue.sourcedata):
    #print widthpixels,catsource.pflux,localdr,localwidth,imagerms
    xpix = float(source[-4])
    ypix = float(source[-3])
    fittedmax = float(source[10])
    amplitude = ((catsource.pflux/localdr)*noisefact)/1000.0#-imagerms
    #print amplitude,catsource.pflux,localdr,noisefact,imagerms
    width = [widthpixels,widthpixels,0.0]
    rmsimagedata = addsource(xpix,ypix,amplitude,width,rmsimagedata)
#exit(-1)
outfile=options.outputdir + '/' + random_filename()
fitsout(rmsimagedata,outfile)
os.rename(outfile,options.outputdir + '/' + inname + '_RMSCORR.FITS')
outfiles.append(options.outputdir + '/' + inname + '_RMSCORR.FITS')

# Write out the rms model image.
#if debug>0 or dorms>0:
#    try:
#        outfile=options.outputdir + '/' + random_filename()
#        fitsout(rmsimagedata,outfile)
#        os.rename(outfile,options.outputdir + '/' + inname + '_RMS.FITS')
#        outfiles.append(options.outputdir + '/' + inname + '_RMS.FITS')
#        outfile=options.outputdir + '/' + random_filename()
#        fitsout(meanimagedata,outfile)
#        os.rename(outfile,options.outputdir + '/' + inname + '_MEAN.FITS')
#        outfiles.append(options.outputdir + '/' + inname + '_MEAN.FITS')
#    except:
#        print '**** Skipping MODELIMAGE.FITS'

# Now divide by the rms
# in the imput image to prepare an image that will be catalogued.
#sadimagedata = ['CATFILE','RMSD',1,1]
#imcombine(meansubimagedata,rmsimagedata,sadimagedata,'DIV',[1,0,0,0,0,0,0,-1,0,0])

# Write out the mean subtracted and rms divided images.
outfile=options.outputdir + '/' + random_filename()
#fitsout(sadimagedata, outfile)
#os.rename(outfile, options.outputdir + '/' + inname + '_SN.fits')
#outfiles.append(options.outputdir + '/' + inname + '_SN.fits')
#outfile=options.outputdir + '/' + random_filename()
fitsout(meansubimagedata, outfile)
os.rename(outfile,options.outputdir + '/' + inname + '_MSUB.fits')
outfiles.append(options.outputdir + '/' + inname + '_MSUB.fits')

# Now produce the catalogue using the AIPS task SAD. Which is run
# on the mean subtracted, rms divided image and gives sources in
# in S/N units. These positions are then refitted using JMFIT or IMFIT
# (or directly from the model image) later on to obtain fluxes in Janskys.
sadcatout = options.outputdir + '/' + random_filename(length=8,suffix='.sad')

residdata = ['IMAGE', 'RED', 1, 1]
imsadrms(meansubimagedata,rmsimagedata,cutoff,outname=sadcatout,uppercut=[1.0*cutoff],
      gain=gain,icut=ifact*cutoff,resid=1,maxwidth=maxfitwidth*beammajpix,residdata=residdata)
#residdata = ['IMAGE', 'RED', 1, 1]
#findsou(sadimagedata,cutoff,err, outname=sadcatout,uppercut=[],
#      gain=gain,icut=ifact*cutoff,maxwidth=maxfitwidth*beammajpix,
#      resid=1,residdata=residdata)
outfile=options.outputdir + '/' + random_filename()
fitsout(residdata, outfile)
os.rename(outfile,options.outputdir + '/' + inname + '_RESID.fits')
outfiles.append(options.outputdir + '/' + inname + '_RESID.fits')
#fitsout(sadimagedata, outfile)
#os.rename(outfile,options.outputdir + '/' + inname + '_CATIMAGE.fits')
#outfiles.append(options.outputdir + '/' + inname + '_CATIMAGE.fits')

# Place the sad output into 'sadoutput' object.
sadfile = open(sadcatout, 'r')
sadoutput = file(sadfile)
sadfile.close()

# Copy the sad file to its final name if the user wants it.
#if keepsad>0:
sadfilename = options.outputdir + '/' + inname +'.sad'
os.rename(sadcatout,sadfilename)
outfiles.append(sadfilename)
print '**** AIPS: SAD output stored in: '+sadfilename
#else:
#    scratchlist.append(sadcatout)
#    print '**** SAD output file removed'




# Construct a catalogue object from the sad output
# AT this point JMfit is run to get accurate flux densities. **Not any more ***
catalogue=sad_catalogue(infile,imagedata,rmsarray,meanarray,sadoutput,maxfitwidth,imagerms,fluxlimit=fluxcutoff)

# Write a raw catalogue if the user wants it.
if keepraw>0:
    rawfilename = options.outputdir + '/' + inname +'.raw'
    catalogue.writeraw(rawfilename)
    outfiles.append(rawfilename)
    print '**** Raw output data stored in: '+rawfilename

# Write a formatted catalogue if the user wants is.
if keepcat>0:
    catfilename = options.outputdir + '/' + inname + '.cat'
    catalogue.writecat(catfilename)
    outfiles.append(catfilename)
    print '**** Formatted output catalogue stored in: '+catfilename

# Write an annotation file for kvis if the user wants it.
if keepann>0:
    annfilename = options.outputdir + '/' + inname + '.ann'
    catalogue.makeann(annfilename,'BLUE')
    outfiles.append(annfilename)
    print '**** Kvis annotation file for catalogue stored in: '+annfilename


#Comparing with external catalogues.
externalsurveys=[]
firstsurvey,nvsssurvey,wensssurvey,sumsssurvey=survey('FIRST'),survey('NVSS'),survey('WENSS'),survey('SUMSS')
firstdata=external_cat(firstsurvey,catalogue,FIRSTcompare)
nvssdata=external_cat(nvsssurvey,catalogue,NVSScompare)
wenssdata=external_cat(wensssurvey,catalogue,WENSScompare)
sumssdata=external_cat(sumsssurvey,catalogue,SUMSScompare)
for i in firstdata,nvssdata,wenssdata,sumssdata:
    if len(i.rawcatdata)>0:
        externalsurveys.append(i)

#VO output..
if keepvo>0:
    votable = VOTableFile()
    resource = Resource()
    votable.resources.append(resource)
    resource.tables.append(catalogue.makevo(votable,description=telescope+' Data'))
    for surveydata in externalsurveys:
        resource.tables.append(surveydata.makevo(votable,description=surveydata.survey.name+' Data'))
    votable.to_xml(options.outputdir+'/'+inname+'.vo')
    outfiles.append(options.outputdir+'/'+inname+'.vo')
    print '**** vo catalogue containg all available source lists is in '+inname+'.vo'

#Make the required plots for each survey.
for surveydata in externalsurveys:
    outfiles=outfiles + surveydata.docomparison(options.outputdir+'/' + inname + '_')

#Print out the best flux comparison (and attempt 2 frequency flux compare)
outfiles=outfiles + fluxcompare(externalsurveys,options.outputdir+'/')

#Output image stats
outfiles=outfiles + [catalogue.getimstats(options.outputdir+'/' + inname + '_STATS.txt',externalsurveys,[localdr,numsources])]

#Output sad header file
#outfiles=outfiles + [catalogue.printheader(options.outputdir+'/CATALOGUE_HEADER.txt')]

#Output file information
#outfiles=outfiles + [catalogue.printfiles(options.outputdir+'/FILE_INFORMATION.txt',externalsurveys,sad=keepsad,cat=keepcat,raw=keepraw,ann=keepann,vo=keepvo)]

#Copy the output to the output directory
now=time.localtime()
for i in scratchlist:
    if os.path.isdir(i): shutil.rmtree(i)
    else: os.remove(i)

# deletes the AIPS environment
print '****\n**** CATALOGER done !\n****'
print '****\n**** Check in '+ options.outputdir +' for the results.\n****'

