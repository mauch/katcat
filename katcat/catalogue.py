#
# Code to deal with catalogue manipulations
# Main classes are:
# 'external_cat': deals with catalogues for comparison
# 'survey':       contains pimary statistics about various surveys
# 'catalogue':    deals with the output of the sad catalogue
#
# Written by txm 24/8/2011

from sys import *
import math
from math import *
from astrop import *

from AIPSTask import *
from AIPSData import *
from Wizardry import AIPSData as Wizardry
from astropy.io.votable.tree import *
from AIPSCata import *

from searchcats import *

import matplotlib
matplotlib.use('Agg')
import pylab
import scipy

import CATparams

import string
import re
import os


#
# A class to deal with the different All Sky surveys
#    
class survey:

    def __init__(self,name,in_bmaj=0.0,in_bmin=0.0,in_bpa=0.0,in_frequency=0.0,in_filename=None,in_sigma=0.0,in_telescope='',in_id=''):
        self.name=name
        self.in_bmaj=in_bmaj
        self.in_bmin=in_bmin
        self.in_bpa=in_bpa
        self.in_filename=in_filename
        # Define some vital statistics of each survey.
        self.frequency=in_frequency
        self.sigma=in_sigma
        self.telescope=in_telescope
        if self.name=='FIRST':
            self.id='FIRST'
            self.sigma=0.2
            self.frequency=1400.0
        elif self.name=='NVSS':
            self.id='NVSS'
            self.sigma=0.5
            self.frequency=1400.0
        elif self.name=='SUMSS':
            self.id='SUMSS'
            self.sigma=2.0
            self.frequency=843.0
        elif self.name=='WENSS':
            self.id='WENSS'
            self.sigma=3.6
            self.frequency=325.0
        else:
            self.id=in_id
            self.sigma=in_sigma
            self.frequency=in_frequency
        
    #
    # given a catimagedata object, return the survey data covering this field.
    #
    def getdata(self,catimagedata):

        if self.id=='SAD':
            self.sigma=self.in_sigma
            rawcatdata=catimagedata.rawdata
        elif self.id=='FIRST':
            rawcatdata=getfirstdata(catimagedata.pos.ra_hms,catimagedata.pos.dec_hms,catimagedata.imsize,'')
            if rawcatdata==[]:
                print "**** No FIRST data available at image position..."
        elif self.id=='NVSS':
            rawcatdata=getnvssdata(catimagedata.pos.ra_hms,catimagedata.pos.dec_hms,catimagedata.imsize)
            if rawcatdata==[]:
                print "**** No NVSS data available at image position..."
        elif self.id=='WENSS':
            rawcatdata=getincdata(catimagedata.pos.ra_hms,catimagedata.pos.dec_hms,catimagedata.imsize,1)
            if rawcatdata==[]:
                print "**** No WENSS data available at image position..."
        elif self.id=='SUMSS':
            rawcatdata=getincdata(catimagedata.pos.ra_hms,catimagedata.pos.dec_hms,catimagedata.imsize,2)
            if rawcatdata==[]:
                print "**** No SUMSS data available at image position..."
        else:
            rawcatdata=[]

        return rawcatdata

    #
    # Given a raw catalogue associated with 'survey.name', return an array of object radio_source
    # constructed from the raw data
    # 
    def makearray(self,rawdata):

        if self.id=='SAD':
            catdata=[]
            for source in rawdata.rawdata:
                # Fitted Major+Minor are equal to the beam if too small.
                fitmaj=max(float(source[16]),self.in_bmaj)
                fitmin=max(float(source[17]),self.in_bmin)
                fitpa=float(source[18])
                fitxpix=float(source[34])
                fitypix=float(source[35])
                # The raw SAD output contains source fluxes in SNR
                # First get the flux in the image units (Jy/beam).
                peakflux,fitmaj,fitmin,fitpa=getknownflux(rawdata.aipsdata,fitmaj,fitmin,fitpa,fitxpix,fitypix)
                
                #peakflux=pixinterp(fitxpix,fitypix,rawdata.aipsdata)
                # If the fitting has failed, use the SAD value (snr*rms)+mean
                if peakflux==0.0:
                    localrms = pixinterp(fitypix,fitxpix,rawdata.rmsdata)
                    localmean = pixinterp(fitypix,fitxpix,rawdata.meandata)
                    peakflux = (float(source[10])*localrms) + localmean

                # Jy -> mJy
                peakflux = peakflux*1000.0
                #Integrated flux = I*source_area/beam_area
                intflux=peakflux*(fitmaj*fitmin)/(self.in_bmaj*self.in_bmin)
                catdata.append(radio_source(int(source[2]),position(float(source[0]),float(source[1]),'dd'),
                                            peakflux,intflux,[fitmaj,fitmin,fitpa,0]))
        
        elif self.id=='FIRST':
            catdata=[]
            for sourcenum,source in enumerate(rawdata):
                catdata.append(radio_source(sourcenum,position([float(source[1]),float(source[2]),float(source[3])],
                                                     [source[4][0],float(source[4]),float(source[5]),float(source[6])],'hms'),
                                            float(source[8]),float(source[9]),
                                            [float(source[14]),float(source[15]),float(source[16]),0]))

        elif self.id=='NVSS':
            catdata=[]
            for sourcenum,source in enumerate(rawdata):
                thispos=position([float(source[1]),float(source[2]),float(source[3])],
                                 [source[4][0],float(source[4]),float(source[5]),float(source[6])],'hms')
                maj=float(source[18])
                min=float(source[19])
                pa=float(source[20])
                catdata.append(radio_source(sourcenum,thispos,float(source[7]),float(source[7]),[maj,min,pa,0]))
        
        elif self.id=='WENSS':
            catdata=[]
            for sourcenum,source in enumerate(rawdata):
                thispos=position([float(source[1]),float(source[2]),float(source[3])],
                                 [source[4][0],float(source[4]),float(source[5]),float(source[6])],'hms')
                if float(source[12])==0:
                    maj=self.resolution(thispos)[0]
                    min=self.resolution(thispos)[1]
                    pa=0.0
                else:
                    maj,min,pa=elconv(self.resolution(thispos),[float(source[12]),float(source[13]),float(source[14])])
                catdata.append(radio_source(sourcenum,thispos,float(source[10]),float(source[11]),[maj,min,pa,0]))


        elif self.id=='SUMSS':
            catdata=[radio_source(sourcenum,position([float(source[1]),float(source[2]),float(source[3])],
                                                     [source[4][0],float(source[4]),float(source[5]),float(source[6])],'hms'),
                                  float(source[10]),float(source[12]),[float(source[14]),float(source[15]),float(source[16]),0])
                     for sourcenum,source in enumerate(rawdata)]

        return catdata
            

    #
    # Return the resolution of the named survey.
    #
    def resolution(self,posn):
        if self.id=='SAD':
            return([self.in_bmaj,self.in_bmin,self.in_bpa])
        elif self.id=='FIRST':
            if (posn.dec_dd>4.5558333):
                beam    = [5.4,5.4,0.0]
            else:
                beam    = [6.4,5.4,0.0]
                if (posn.dec_dd<-2.506944):
                    if (posn.ra_dd>(21*15)) or (posn.ra_dd<(3*15)):
                        beam[0] = 6.8
            return(beam)

        elif self.id=='NVSS':
            return([45.0,45.0,0.0])

        elif self.id=='SUMSS':
            min=45.0
            pa=0.0
            maj=45.0/sind(fabs(posn.dec_dd))
            return([maj,min,pa])
            
        elif self.id=='WENSS':
            min=54.0
            pa=0.0
            maj=54.0/sind(fabs(posn.dec_dd))
            return([maj,min,pa])

    #
    # Return the VO table metadata associated with the raw surey data and the raw data formatted accordingly.
    #
    def vo_metadata(self,rawdata,votable):
        if self.id=='SAD':

            vofields=[Field(votable,ID="sourcenum", name="Source Number", datatype="int", ucd="meta.id"),
                      Field(votable,ID="radeg", name="Right Ascension",  datatype="double", ref="J2000", ucd="pos.eq.ra",unit="deg",width="12",precision="12"),
                      Field(votable,ID="decdeg", name="Declination", datatype="double", ref="J2000", ucd="pos.eq.dec",unit="deg",width="12",precision="12"),
                      Field(votable,ID="raerr", name="RA Error", datatype="float", ucd="stat.error;pos.eq.ra", unit="arcsec"),
                      Field(votable,ID="decerror", name="Dec Error", datatype="float", ucd="stat.error;pos.eq.dec", unit="arcsec"),
                      Field(votable,ID="pflux", name="Peak Flux", datatype="float", ucd="phot.flux.density.sb;em.radio.750-1500MHz", unit="mJy/beam"),
                      Field(votable,ID="pfluxerr", name="Peak Flux Error", datatype="float", ucd="stat.err;phot.flux.density.sb;em.radio.750-1500MHz", unit="mJy/beam"),
                      Field(votable,ID="totflux", name="Total Flux Density", datatype="float", ucd="phot.flux.density;em.radio.750-1500MHz", unit="mJy"),
                      Field(votable,ID="totfluxerr", name="Total Flux Density Error", datatype="float", ucd="stat.error;phot.flux.density;em.radio.750-1500MHz", unit="mJy"),
                      Field(votable,ID="fitmaj", name="Fitted Major Axis FWHM", datatype="float", ucd="phys.angSize;stat.fit", unit="arcsec"),
                      Field(votable,ID="fitmin", name="Fitted Minor Axis FWHM", datatype="float", ucd="phys.angSize;stat.fit", unit="arcsec"),
                      Field(votable,ID="fitpa", name="Fitted Position Angle", datatype="float", ucd="pos.posAng;stat.fit", unit="deg"),
                      Field(votable,ID="fitmajerr", name="Major Axis FWHM Error", datatype="float", ucd="stat.err;phys.angSize;stat.fit", unit="arcsec"),
                      Field(votable,ID="fitminerr", name="Minor Axis FWHM Error", datatype="float", ucd="stat.err;phys.angSize;stat.fit", unit="arcsec"),
                      Field(votable,ID="fitpaerr", name="Position Angle Error", datatype="float", ucd="stat.err;pos.posAng;stat.fit", unit="deg"),
                      Field(votable,ID="xpix", name="X Pixel", datatype="float", ucd="pos.cartesian.x"),
                      Field(votable,ID="ypix", name="Y Pixel", datatype="float", ucd="pos.cartesian.y"),
                      Field(votable,ID="imname", name="Image Name", datatype="char", arraysize='*'),
                      Field(votable,ID="locrms", name="Local image rms", datatype="float", ucd="stat.stdev", unit="mjy/beam"),
                      Field(votable,ID="locmean", name="Local image mean", datatype="float", ucd="stat.sdev", unit="mjy/beam")]

            formatdata=[]
            for rawsource,source in zip(rawdata.rawdata,rawdata.sourcedata):
                localrms = pixinterp(float(rawsource[34]),float(rawsource[35]),rawdata.rmsdata)
                localmean = pixinterp(float(rawsource[34]),float(rawsource[35]),rawdata.meandata)
                outarray=[]
                outarray.append(int(rawsource[2]))
                outarray.append(source.position.ra_dd)
                outarray.append(source.position.dec_dd)
                raerr,decerr=source.posnerrorcondon(localrms,self)
                outarray.append(raerr)
                outarray.append(decerr)
                perr,interr=source.fluxerrorcondon(localrms,self)
                outarray.append(source.pflux)
                outarray.append(perr)
                outarray.append(source.flux)
                outarray.append(interr)
                outarray.append(source.maj)
                outarray.append(source.min)
                outarray.append(source.pa)
                majerr,minerr,paerr=source.sizeerrorcondon(localrms,self)
                outarray.append(majerr)
                outarray.append(minerr)
                outarray.append(paerr)
                outarray.append(float(rawsource[34]))
                outarray.append(float(rawsource[35]))
                outarray.append(rawdata.fitsfile)
                outarray.append(localrms*1000.0)
                outarray.append(localmean*1000.0)
                formatdata=formatdata+[outarray]

        elif self.id=='FIRST':
            vofields=[Field(votable,ID="firstoffset", name="Offset from field center", datatype="float", ucd="pos.angDistance", unit="arcsec"),
                      Field(votable,ID="firstradeg", name="Right Ascension",  datatype="double", ref="J2000", ucd="pos.eq.ra",unit="deg",width="12",precision="12"),
                      Field(votable,ID="firstdecdeg", name="Declination", datatype="double", ref="J2000", ucd="pos.eq.dec",unit="deg",width="12",precision="12"),
                      Field(votable,ID="firstsideprob", name="Side lobe probability", datatype="float", ucd="stat.probability"),
                      Field(votable,ID="firstpeakflux", name="1.4GHz peak flux", datatype="float", unit="mJy/beam", ucd="phot.flux;em.radio.750-1500;stat.max"),
                      Field(votable,ID="firstintflux", name="1.4GHz total flux density", datatype="float", unit="mJy", ucd="phot.flux.density;em.radio.750-1500MHz"),
                      Field(votable,ID="firstrms", name="Local rms noise", datatype="float", unit="mJy/beam", ucd="stat.stdev;instr.det.noise"),
                      Field(votable,ID="firstmaj", name="Deconvolved major axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="firstmin", name="Deconvolved minor axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="firstPA", name="Deconvolved position angle", datatype="float",unit="arcsec", ucd="pos.posAng"),
                      Field(votable,ID="firstfmaj", name="Fitted major axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize;meta.modelled"),
                      Field(votable,ID="firstfmin", name="Fitted minor axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize;meta.modelled"),
                      Field(votable,ID="firstfPA", name="Fitted position angle", datatype="float", unit="arcsec", ucd="pos.posAng;meta.modelled"),
                      Field(votable,ID="firstimname", name="Field name", datatype="char", arraysize='*', ucd="meta.id;obs.field"),
                      Field(votable,ID="firstnumsdss", name="Number of SDSS matches within 8 arcsec", datatype="int", ucd="meta.number"),
                      Field(votable,ID="firstnearsdsssep", name="Closest SDSS separation", datatype="float", unit="arcsec", ucd="pos.angDistance"),
                      Field(votable,ID="firstsdssi", name="i mag of nearest SDSS source", datatype="char", arraysize='*', unit="mag.", ucd="phot.mag;em.opt.I"),
                      Field(votable,ID="firstsdssclass", name="Classification of nearest SDSS source", arraysize='*', datatype="char", ucd="meta.code.class;src.class"),
                      Field(votable,ID="firstnum2mass", name="Number of 2MASS matches within 8 arcsec", datatype="float", ucd="meta.number"),
                      Field(votable,ID="firstnear2masssep", name="Closest 2MASS separation", datatype="float", unit="arcsec", ucd="pos.angDistance"),
                      Field(votable,ID="first2massk", name="K-band magnitude of nearest 2MASS source", datatype="char", arraysize='*', unit="mag.", ucd="phot.mag;em.IR.K")]

            formatdata=[]
            for line in rawdata:
                radd=hmsra2ra(float(line[1]),float(line[2]),float(line[3]))
                decdd=dmsdec2dec(line[4][0],float(line[4]),float(line[5]),float(line[6]))
                formatdata=formatdata+[(float(line[0]),radd,decdd,float(line[7]),float(line[8]),float(line[9]),
                                        float(line[10]),float(line[11]),float(line[12]),float(line[13]),
                                        float(line[14]),float(line[15]),float(line[16]),line[17],int(line[18]),
                                        float(line[19]),line[20],line[21],int(line[22]),float(line[23]),line[24])]
                             
        elif self.id=='NVSS':
            vofields=[Field(votable,ID="nvssoffset", name="Offset from field center", datatype="float", ucd="pos.angDistance", unit="arcsec"),
                      Field(votable,ID="nvssradeg", name="Right Ascension",  datatype="double", ref="J2000", ucd="pos.eq.ra",unit="deg",width="12",precision="12"),
                      Field(votable,ID="nvssdecdeg", name="Declination", datatype="double", ref="J2000", ucd="pos.eq.dec",unit="deg",width="12",precision="12"),
                      Field(votable,ID="nvssintflux", name="1.4GHz total flux density", datatype="float", unit="mJy", ucd="phot.flux.density;em.radio.750-1500MHz"),                
                      Field(votable,ID="nvssmaj", name="Deconvolved major axis FWHM", datatype="char", arraysize='*', unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="nvssmin", name="Deconvolved minor axis FWHM", datatype="char", arraysize='*', unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="nvssPA", name="Deconvolved position angle", datatype="char",unit="arcsec", ucd="pos.posAng"),
                      Field(votable,ID="nvssimname", name="Field name", datatype="char", arraysize='*', ucd="meta.id;obs.field"),
                      Field(votable,ID="nvssxpix", name="X Pixel", datatype="float", ucd="pos.cartesian.x"),
                      Field(votable,ID="nvssypix", name="Y Pixel", datatype="float", ucd="pos.cartesian.y")]
            
            formatdata=[]
            for line in rawdata:
                radd=hmsra2ra(float(line[1]),float(line[2]),float(line[3]))
                decdd=dmsdec2dec(line[4][0],float(line[4]),float(line[5]),float(line[6]))
                formatdata=formatdata+[(float(line[0]),radd,decdd,float(line[7]),line[9],line[10],
                                        line[11],line[15],float(line[16]),float(line[17]))]

        elif self.id=='WENSS':

            vofields=[Field(votable,ID="wensssname", name="Source Name", datatype="char", arraysize='*'),
                      Field(votable,ID="wenssradeg", name="Right Ascension",  datatype="double", ref="J2000", ucd="pos.eq.ra",unit="deg",width="12",precision="12"),
                      Field(votable,ID="wenssdecdeg", name="Declination", datatype="double", ref="J2000", ucd="pos.eq.dec",unit="deg",width="12",precision="12"),
                      Field(votable,ID="wensspeakflux", name="325MHz peak flux", datatype="float", unit="mJy/beam", ucd="phot.flux;em.radio.200-400MHz;stat.max"),
                      Field(votable,ID="wenssintflux", name="325MHz total flux density", datatype="float", unit="mJy", ucd="phot.flux.density;em.radio.200-400MHz"),
                      Field(votable,ID="wenssmaj", name="Deconvolved major axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="wenssmin", name="Deconvolved minor axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="wenssPA", name="Deconvolved position angle", datatype="float",unit="arcsec", ucd="pos.posAng"),
                      Field(votable,ID="wenssimname", name="Field name", datatype="char", arraysize='*', ucd="meta.id;obs.field")]
            
            formatdata=[]
            for line in rawdata:
                radd=hmsra2ra(float(line[1]),float(line[2]),float(line[3]))
                decdd=dmsdec2dec(line[4][0],float(line[4]),float(line[5]),float(line[6]))
                formatdata=formatdata+[(line[0],radd,decdd,float(line[10]),float(line[11]),float(line[12]),
                                        float(line[13]),float(line[14]),line[16])]

        elif self.id=='SUMSS':

            vofields=[Field(votable,ID="sumsssname", name="Source Name", datatype="char", arraysize='*'),
                      Field(votable,ID="sumssradeg", name="Right Ascension",  datatype="double", ref="J2000", ucd="pos.eq.ra",unit="deg",width="12",precision="12"),
                      Field(votable,ID="sumssdecdeg", name="Declination", datatype="double", ref="J2000", ucd="pos.eq.dec",unit="deg",width="12",precision="12"),
                      Field(votable,ID="sumssraerr", name="Right Ascension Error", datatype="float", ucd="stat.error;pos.eq.ra",unit="arcsec"),
                      Field(votable,ID="sumssdecerr", name="Declination Error", datatype="float", ucd="stat.error;pos.eq.dec",unit="deg"),
                      Field(votable,ID="sumsspeakflux", name="843MHz peak flux", datatype="float", unit="mJy/beam", ucd="phot.flux;em.radio.750-1500;stat.max"),
                      Field(votable,ID="sumsspfuxerror", name="843MHz peak flux error", datatype="float", unit="mJy/beam", ucd="stat.error;phot.flux;em.radio.750-1500;stat.max"),
                      Field(votable,ID="sumssintflux", name="843MHz total flux density", datatype="float", unit="mJy", ucd="phot.flux.density;em.radio.750-1500MHz"),
                      Field(votable,ID="sumssintfluxerror", name="843MHz total flux density error", datatype="float", unit="mJy", ucd="stat.error;phot.flux.density;em.radio.750-1500"),
                      Field(votable,ID="sumssfmaj", name="Fitted major axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize;meta.modelled"),
                      Field(votable,ID="sumssfmin", name="Fitted minor axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize;meta.modelled"),
                      Field(votable,ID="sumssfPA", name="Fitted position angle", datatype="float", unit="arcsec", ucd="pos.posAng;meta.modelled"),
                      Field(votable,ID="sumssmaj", name="Deconvolved major axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="sumssmin", name="Deconvolved minor axis FWHM", datatype="float", unit="arcsec", ucd="phys.angSize"),
                      Field(votable,ID="sumssPA", name="Deconvolved position angle", datatype="char",arraysize='*',unit="deg", ucd="pos.posAng"),
                      Field(votable,ID="sumssimname", name="Mosaic name", datatype="char", arraysize='*',ucd="meta.id;obs.field"),
                      Field(votable,ID="sumssnummosaic", name="Number of included mosaics", datatype="int", ucd="meta.number"),
                      Field(votable,ID="sumssxpix", name="X-Pixel position of source", datatype="float", ucd="pos.cartesian.x"),
                      Field(votable,ID="sumssypix", name="Y-Pixel position of source", datatype="float", ucd="pos.cartesian.y")]
            
            formatdata=[]
            for line in rawdata:
                radd=hmsra2ra(float(line[1]),float(line[2]),float(line[3]))
                decdd=dmsdec2dec(line[4][0],float(line[4]),float(line[5]),float(line[6]))
                formatdata=formatdata+[(line[0],radd,decdd,float(line[8]),float(line[9]),float(line[10]),
                                        float(line[11]),float(line[12]),float(line[13]),float(line[14]),
                                        float(line[15]),float(line[16]),float(line[17]),float(line[18]),
                                        line[19],line[20],int(line[21]),float(line[22]),float(line[23]))]
                
        else:
            #Own table - need to deal with this case
            vofields=[]
            rawdata=[]
      
        return(vofields,formatdata)

    #def area(self):
    #Not yet implemented...


class position:

    def __init__(self,ra,dec,type='hms'):

        if type=='hms':
            self.ra_dd=hmsra2ra(ra[0],ra[1],ra[2])
            self.dec_dd=dmsdec2dec(dec[0],dec[1],dec[2],dec[3])          
        elif type=='dd':
            if (ra<0):
                ra=ra%360
            self.ra_dd=float(ra)
            self.dec_dd=float(dec)
        
        self.ra_hms=deg2sex(self.ra_dd,self.dec_dd)[:3]
        self.dec_hms=deg2sex(self.ra_dd,self.dec_dd)[3:]

    #
    # Return the offset in arcseconds between the postion and
    # and input position.
    #
    def separation(self,position):

        ra1=self.ra_dd
        ra2=position.ra_dd
        dec1=self.dec_dd
        dec2=position.dec_dd

        if ra1==ra2 and dec1==dec2:
            return [0.0,0.0,0.0]
        if ra1>250 and ra2<110:
            ra2=ra2+360
        if ra1<110 and ra2>250:
            ra1=ra1+360

        a=cosd(dec1)*cosd(ra1)*cosd(dec2)*cosd(ra2)
        b=cosd(dec1)*sind(ra1)*cosd(dec2)*sind(ra2)
        c=sind(dec1)*sind(dec2)
        separation=3600.0*degrees(acos(a+b+c))
        
        decoff=(dec2-dec1)*3600
        raoff=(ra2-ra1)*cosd((dec1+dec2)/2)*3600
        
        return separation
    #
    # Return the RA offset between this position and imput position (arcsec).
    #
    def raoff(self,position):

        avdec = (self.dec_dd+position.dec_dd)/2.0
        return (self.ra_dd-position.ra_dd)*cosd(avdec)*3600.0

    #
    # Return the Dec. offset between this postion and input position (arcsec).
    #
    def decoff(self,position):

        return (self.dec_dd-position.dec_dd)*3600.0

    #possibilities:
    #deg B1950,galcoord


class radio_source:

    def __init__(self,id,pos,pflux,flux,size):
        self.id=id
        self.position=pos
        self.pflux=pflux
        self.flux=flux
        self.size_is_dec=size[3]
        self.maj=size[0]
        self.min=size[1]
        self.pa=size[2]


    #
    # Check if the source is isolated with respect
    # to the input array of 'radio_source' within a given radius
    # above an optional given flux. Return the number of source within the radius.
    #
    def isisolated(self,sourcearray,radius,flux=0.0):
        
        numisolated=0
        for source in sourcearray:
            if source.flux>flux:
                if self.position.separation(source.position)<radius:
                    numisolated=numisolated+1
                    
        return numisolated


    #
    # Check if this source is resolved (>2.33sigma beam; 98% confidence)
    #
    def isresolved(self,survey,sigma):
        dmaj=False
        dmin=False
        if self.size_is_dec>0:
            if (self.maj>0):
                dmaj=True
            if (self.min>0):
                dmin=True
        else:
            beam=survey.resolution(self.position)
            sizeerror=self.sizeerrorcondon(sigma,survey)
            majcompare=beam[0]+(2.33*sizeerror[0])
            mincompare=beam[1]+(2.33*sizeerror[1])
            dmajcomp,dmincomp,dpacomp=deconv([self.maj,self.min,self.pa],[majcompare,mincompare,beam[2]])
            if (dmajcomp>0):
                dmaj=True
            if (dmincomp>0):
                dmin=True
        return(dmaj,dmin)

    #
    # Implement errors for elliptical gaussians in the presence of correlated noise
    # as per Condon (1998), MNRAS.
    #
    def rhosqt1(self,localrms,ncorr):
        return ((self.maj*self.min)/(4.0*ncorr)*(self.pflux**2.0)/(localrms**2.0))
    
    def rhosqt2(self,ncorr):
        return (1.0 + (ncorr)/(self.maj**2.0))

    def rhosqt3(self,ncorr):
        return (1.0 + (ncorr)/(self.min**2.0))

    def sizeerrorcondon(self,localrms,survey):

        beam_maj=survey.resolution(self.position)[0]
        beam_min=survey.resolution(self.position)[1]
        ncorr=beam_maj*beam_min
        majerr = sqrt((((self.maj**2.0)*2.0)/((self.rhosqt1(localrms,ncorr))*(self.rhosqt2(ncorr)**2.5)*(self.rhosqt3(ncorr)**0.5))))
        minerr = sqrt((((self.min**2.0)*2.0)/((self.rhosqt1(localrms,ncorr))*(self.rhosqt2(ncorr)**0.5)*(self.rhosqt3(ncorr)**2.5))))
        if self.maj==self.min:
            majcorr=self.maj+0.0001
        else:
            majcorr=self.maj
        paerr = degrees(sqrt((4.0/((self.rhosqt1(localrms,ncorr))*(self.rhosqt2(ncorr)**2.5)*(self.rhosqt3(ncorr)**0.5)))*
                             (((majcorr*self.min)/((majcorr**2.0)-(self.min**2.0)))**2.0)))
        if paerr>90:
            paerr=90
        return math.sqrt(majerr**2.0),math.sqrt(minerr**2.0),int(paerr)
        
    def fluxerrorcondon(self,localrms,survey):
        beam_maj=survey.resolution(self.position)[0]
        beam_min=survey.resolution(self.position)[1]
        sizeerr=self.sizeerrorcondon(localrms,survey)

        #Noise correlation length scale
        ncorr=beam_maj*beam_min

        #peak error
        #print self.pflux
        perr=sqrt((0.03*self.pflux)**2.0 + (((self.pflux**2.0)*2.0)/((self.rhosqt1(localrms,ncorr))*(self.rhosqt2(ncorr)**1.5)*(self.rhosqt3(ncorr)**1.5))))
        #print perr
        #print sizeerr[0]**2.0
        #print self.maj**2.0
        #integrated flux error
        interr=sqrt((((perr**2.0)/(self.pflux**2.0))+
                     ((ncorr)/(self.maj*self.min)*(((sizeerr[0]**2.0)/(self.maj**2.0))+
                                                   ((sizeerr[1]**2.0)/(self.min**2.0)))))*self.flux**2.0)
        print (sizeerr[0]**2.0)/(self.maj**2.0), (sizeerr[1]**2.0)/(self.min**2.0), (ncorr)/(self.maj*self.min), (perr**2.0)/(self.pflux**2.0), self.flux
        return(perr,interr)
        
    def posnerrorcondon(self,localrms,survey):
        beam_maj=survey.resolution(self.position)[0]
        beam_min=survey.resolution(self.position)[1]
        ncorr=beam_maj*beam_min
        sigmaxsq = (2.0*(self.maj**2.0))/(8.0*math.log(2.0)*
                                          (self.rhosqt1(localrms,ncorr)*(self.rhosqt2(ncorr)**2.5)*(self.rhosqt3(ncorr)**0.5)))
        sigmaysq = (2.0*(self.min**2.0))/(8.0*math.log(2.0)*
                                          (self.rhosqt1(localrms,ncorr)*(self.rhosqt2(ncorr)**0.5)*(self.rhosqt3(ncorr)**2.5)))
        raerr = sqrt(sigmaxsq*(math.sin(radians(self.pa))**2.0) + sigmaysq*(math.cos(radians(self.pa))**2.0))
        decerr = sqrt(sigmaxsq*(math.cos(radians(self.pa))**2.0) + sigmaysq*(math.sin(radians(self.pa))**2.0))
        return(raerr,decerr)           

#
# Class to deal with external comparison catalogues.
# Needs a survey object (to define which external survey we are using)
# and an imagedata_cat object (to define the data we are comparing with).
#
class external_cat:

    def __init__(self,survey,saddata,compareinfo=[0,0,0]):

        self.saddata=saddata
        self.survey=survey
        self.compareinfo=compareinfo
        self.rawcatdata=survey.getdata(saddata)
        self.catdata=survey.makearray(self.rawcatdata)
        self.doneflux=False
        self.doneposn=False
        self.doneann=False
        self.donetwofreq=False
        print '**** '+str(len(self.rawcatdata))+' sources found in '+survey.name
        

    colours=iter(['RED','GREEN','CYAN','ORANGE','YELLOW','WHITE'])
    #
    # For each radio_source object in self.catdata, return all the radio_source objects in imagedata
    # that are within threshold (arcsec) of the position. Return a list of dictionarys of the form:
    # 'matches': list of radio_source matches within offset, 'source': copy of the matchind source.
    # 
    def getmatches_external(self,offset):

        matchdata={'sourcesurvey': self.survey,
                   'matchsurvey': self.saddata.survey,
                   'matches': []}
        for sourcenum,source in enumerate(self.catdata):
            matchdata['matches'].append({'source':source, 'match':[]})
            for source2 in self.saddata.sourcedata:
                thissep=source.position.separation(source2.position)
                if thissep<offset:
                    matchdata['matches'][sourcenum]['match'].append(source2)

        return(matchdata)

    #
    # For each radio_source object in imagedata, return all the radio_source objects in catdata
    # This basically does the reverse of getmatches_external.
    #
    def getmatches_imagedata(self,offset):

        matchdata={'sourcesurvey': self.saddata.survey,
                   'matchsurvey': self.survey,
                   'matches': []}
        for sourcenum,source in enumerate(self.saddata.sourcedata):

            matchdata['matches'].append({'source': source, 'match': []})
            for source2 in self.catdata:
                thissep=source.position.separation(source2.position)
                if thissep<offset:
                    matchdata['matches'][sourcenum]['match'].append(source2)

        return(matchdata)
    
    #
    # Given a vo object, return a votable object from self.survey containing the external catalogue. 
    #
    def makevo(self,votable,description='Unknown catalogue'):

        table=Table(votable,ID=self.survey.name,name=self.survey.name+' Data')
        table.description=description
        table.format='tabledata'

        metadata,formatdata=self.survey.vo_metadata(self.rawcatdata,votable)

        table.fields.extend(metadata)
        table.create_arrays(len(formatdata))
        for num,thisitem in enumerate(formatdata):
            table.array[num]=tuple(thisitem)
        return table

    #
    # Make a kvis annotation file to plot fitted ellipses from the catalogue.
    # outname=output filename- will have '.ann' extension appended.
    # colour=colour for the contours.
    #
    def makeann(self,outname,colour='BLUE'):

        file=open(outname+'.ann','w')
        file.write('COLOUR '+colour+'\n')
        file.write('PA SKY\n')

        for source in self.catdata:
            ra=source.position.ra_dd
            dec=source.position.dec_dd
            if source.size_is_dec>0:
                if source.min>0:
                    fmaj,fmin,fpa=elconv([source.maj,source.min,source.pa],
                                         self.survey.resolution(source.position))
                else:
                    beam=self.survey.resolution(source.position)
                    fmaj,fmin,fpa=beam[0],beam[1],beam[2]
            else:
                fmaj=source.maj
                fmin=source.min
                fpa=source.pa
                
            file.write('ELLIPSE %f %f %f %f %f\n'%(ra,dec,fmaj/3600.0,fmin/3600.0,fpa))

        file.close()
        self.doneann=True
        return([outname+'.ann'])

    # Do flux density related plots.
    # 1. Plot the flux from sad vs. flux from survey.
    # 2. Plot a spectral index histogram.
    # 3+4. Plot the flux offset vs distance from image pointing center.
    # 5. Plot the flux vs sad vs. flux from suvrey with image source numbers.
    def plotflux(self,outname,doplots=1):

        from makeplots import plotfluxes,plotfluxeserror,plotsihisto,plotfluxpboff,plotfluxpbrad,plotfluxsnums

        plotsmade=[]

        sad_res=self.saddata.survey.resolution(self.saddata.pos)
        survey_res=self.survey.resolution(self.saddata.pos)

        # Which survey has the highest resolution??
        if (sad_res[0]*sad_res[1])>(survey_res[0]*survey_res[1]):
            radius=sad_res[0]/2.0
            matchdata=self.getmatches_imagedata(radius)
        else:
            radius=survey_res[0]/2.0
            matchdata=self.getmatches_external(radius)

        self.matchdata=matchdata
        # Make separate arrays for extended sources and point sources
        ps_matches,ps_flux_source,ps_flux_match,ps_beamdists,ps_raoffs,ps_decoffs,ps_ids=[],[],[],[],[],[],[]
        es_matches,es_flux_source,es_flux_match,es_beamdists,es_raoffs,es_decoffs,es_ids=[],[],[],[],[],[],[]
        es_fluxerror_source,es_fluxerror_match,ps_fluxerror_source,ps_fluxerror_match=[],[],[],[]
        es_ramatch,es_decmatch,ps_ramatch,ps_decmatch=[],[],[],[]
        for matches in matchdata['matches']:
            if len(matches['match'])==1 and not numpy.array(matches['source'].isresolved(matchdata['sourcesurvey'],matchdata['sourcesurvey'].sigma)).all() and not numpy.array(matches['match'][0].isresolved(matchdata['matchsurvey'],matchdata['matchsurvey'].sigma)).all():
                #Point source
                ps_matches.append(matches)
                ps_flux_source.append(matches['source'].flux)
                ps_flux_match.append(matches['match'][0].flux)
                ps_beamdists.append(matches['source'].position.separation(self.saddata.pos)/60.0)
                ps_raoffs.append(matches['source'].position.raoff(self.saddata.pos)/60.0)
                ps_decoffs.append(matches['source'].position.decoff(self.saddata.pos)/60.0)
                if matchdata['matchsurvey'].id==self.saddata.survey.id:
                    ps_ids.append(str(matches['match'][0].id))
                else:
                    ps_ids.append(str(matches['source'].id))
                ps_ramatch.append(matches['source'].position.ra_hms)
                ps_decmatch.append(matches['source'].position.dec_hms)
                ps_fluxerror_source.append(matches['source'].fluxerrorcondon(matchdata['sourcesurvey'].sigma,matchdata['sourcesurvey'])[1])
                ps_fluxerror_match.append(matches['match'][0].fluxerrorcondon(matchdata['matchsurvey'].sigma,matchdata['matchsurvey'])[1])
            elif len(matches['match'])>0:
                #Extended source
                es_matches.append(matches)
                es_flux_source.append(matches['source'].flux)
                es_fluxerror_source.append(matches['source'].fluxerrorcondon(matchdata['sourcesurvey'].sigma,matchdata['sourcesurvey'])[1])
                sum=0.0
                error=0.0
                for match in matches['match']:
                    sum=sum+match.flux
                    error=error+(match.fluxerrorcondon(matchdata['matchsurvey'].sigma,matchdata['matchsurvey'])[1])**2.0
                es_fluxerror_match.append(math.sqrt(error))
                es_flux_match.append(sum)
                es_beamdists.append(matches['source'].position.separation(self.saddata.pos)/60.0)
                es_raoffs.append(matches['source'].position.raoff(self.saddata.pos)/60.0)
                es_decoffs.append(matches['source'].position.decoff(self.saddata.pos)/60.0)
                label=''
                if matchdata['matchsurvey'].id==self.saddata.survey.id:
                    for match in matches['match']:
                        label=label+str(match.id)+':'
                    label=label[:-1]
                    es_ids.append(label)
                else:
                    es_ids.append(str(matches['source'].id))
                es_ramatch.append(matches['source'].position.ra_hms)
                es_decmatch.append(matches['source'].position.dec_hms)
        ps_flux_source=numpy.array(ps_flux_source)
        ps_flux_match=numpy.array(ps_flux_match)
        ps_beamdists=numpy.array(ps_beamdists)
        ps_raoffs=numpy.array(ps_raoffs)
        ps_decoffs=numpy.array(ps_decoffs)
        ps_fluxerror_source=numpy.array(ps_fluxerror_source)
        ps_fluxerror_match=numpy.array(ps_fluxerror_match)
        es_flux_source=numpy.array(es_flux_source)
        es_flux_match=numpy.array(es_flux_match)
        es_beamdists=numpy.array(es_beamdists)
        es_raoffs=numpy.array(es_raoffs)
        es_decoffs=numpy.array(es_decoffs)
        es_fluxerror_source=numpy.array(es_fluxerror_source)
        es_fluxerror_match=numpy.array(es_fluxerror_match)
        
        predictfactor=numpy.log10(matchdata['sourcesurvey'].frequency/matchdata['matchsurvey'].frequency) * (-0.7)
        if matchdata['sourcesurvey'].id==self.saddata.survey.id:
            # Predict flux of other catalogue using alpha=-0.7
            ps_flux_predict = numpy.power(10,predictfactor + numpy.log10(ps_flux_match))
            es_flux_predict = numpy.power(10,predictfactor + numpy.log10(es_flux_match))
            ps_flux_orig    = ps_flux_match
            es_flux_orig    = es_flux_match
            ps_fluxerror_predict = ps_fluxerror_match*ps_flux_predict/ps_flux_match
            es_fluxerror_predict = es_fluxerror_match*es_flux_predict/es_flux_match
            ps_flux         = ps_flux_source
            es_flux         = es_flux_source
            ps_fluxerror      = ps_fluxerror_source
            es_fluxerror      = es_fluxerror_source
            ps_flux_offset  = numpy.log10(ps_flux_predict) - numpy.log10(ps_flux_source)
            es_flux_offset  = numpy.log10(es_flux_predict) - numpy.log10(es_flux_source)
        else:
            ps_flux_predict = numpy.power(10,-predictfactor + numpy.log10(ps_flux_source))
            es_flux_predict = numpy.power(10,-predictfactor + numpy.log10(es_flux_source))
            ps_flux_orig    = ps_flux_source
            es_flux_orig    = es_flux_source
            ps_fluxerror_predict = ps_fluxerror_source*ps_flux_predict/ps_flux_source
            es_fluxerror_predict = es_fluxerror_source*es_flux_predict/es_flux_source
            ps_flux         = ps_flux_match
            es_flux         = es_flux_match
            ps_fluxerror      = ps_fluxerror_match
            es_fluxerror      = es_fluxerror_match
            ps_flux_offset  = numpy.log10(ps_flux_match) - numpy.log10(ps_flux_predict)
            es_flux_offset  = numpy.log10(es_flux_match) - numpy.log10(es_flux_predict)

        #Spectral indices
        ps_si = numpy.log10(ps_flux_source/ps_flux_match)/numpy.log10(matchdata['sourcesurvey'].frequency/matchdata['matchsurvey'].frequency)
        es_si = numpy.log10(es_flux_source/es_flux_match)/numpy.log10(matchdata['sourcesurvey'].frequency/matchdata['matchsurvey'].frequency)
            
        # Combine point and extended source data
        flux_source  = numpy.append(ps_flux_source,es_flux_source)
        flux_match   = numpy.append(ps_flux_match,es_flux_match)
        flux_offset  = numpy.append(ps_flux_offset,es_flux_offset)
        flux         = numpy.append(ps_flux,es_flux)
        flux_predict = numpy.append(ps_flux_predict,es_flux_predict)
        si           = numpy.append(ps_si,es_si)
        raoffs       = numpy.append(ps_raoffs,es_raoffs)
        decoffs      = numpy.append(ps_decoffs,es_decoffs)
        fluxerror_source = numpy.append(ps_fluxerror_source,es_fluxerror_source)
        fluxerror_match  = numpy.append(ps_fluxerror_match,es_fluxerror_match)
        fluxerror    = numpy.append(ps_fluxerror,es_fluxerror)
        fluxerror_predict = numpy.append(ps_fluxerror_predict,es_fluxerror_predict)
        ids          = ps_ids+es_ids

 
        self.ps_raoffs=ps_raoffs
        self.es_raoffs=es_raoffs
        self.ps_decoffs=ps_decoffs
        self.es_decoffs=es_decoffs
        self.ps_beamdists=ps_beamdists
        self.es_beamdists=es_beamdists
        self.ids=ids
        self.fluxoff=scipy.mean(flux_offset)
        self.psfluxoff=scipy.mean(ps_flux_offset)
        self.flux=flux
        self.ps_flux=ps_flux
        self.es_flux=es_flux
        self.flux_predict=flux_predict
        self.ps_flux_predict=ps_flux_predict
        self.es_flux_predict=es_flux_predict
        self.es_matches=es_matches
        self.ps_matches=ps_matches
        self.matches=numpy.append(ps_matches,es_matches)
        self.ps_flux_orig=ps_flux_orig
        self.es_flux_orig=es_flux_orig
        self.flux_orig=numpy.append(ps_flux_orig,es_flux_orig)
        self.ps_ramatch=ps_ramatch
        self.ps_decmatch=ps_decmatch
        self.es_ramatch=es_ramatch
        self.es_decmatch=es_decmatch
        
        
        if len(flux)>0 and doplots:

            plotsmade = plotsmade+[plotfluxes(ps_flux_predict,ps_flux,es_flux_predict,es_flux,self.psfluxoff,self.fluxoff,
                                              self.saddata.survey.telescope,self.survey.name,outname)]
            plotsmade = plotsmade+[plotfluxeserror(flux_predict,flux,fluxerror_predict,fluxerror,self.psfluxoff,self.fluxoff,
                                                   self.saddata.survey.telescope,self.survey.name,outname)]
            if (len(flux)>15):
                plotsmade = plotsmade+[plotsihisto(flux_source,flux_match,matchdata['sourcesurvey'].frequency,
                                                   matchdata['matchsurvey'].frequency,outname)]
            plotsmade = plotsmade+[plotfluxpboff(raoffs,decoffs,flux_offset,makeaipsimage(self.saddata.aipsdata),outname)]
            plotsmade = plotsmade+[plotfluxpbrad(ps_beamdists,ps_flux_offset,es_beamdists,es_flux_offset,outname)]
            plotsmade = plotsmade+[plotfluxsnums(flux,flux_predict,ids,self.saddata.survey.telescope,self.survey.name,outname)]
            print '**** Made '+self.survey.name+' flux comparison plots: '+str(plotsmade)
            self.doneflux=True
        elif len(flux)==0:
            print '**** No vaid matches found with '+self.survey.name+' so cannot do flux comparison.'
        return plotsmade
    
    # Do position offset plots.
    # 1. Plot the position offset as a function of distance from the image pointing center.
    # 2. Plot the ra offset vs the dec offset.
    # This will generate self.raoffsets and self.decoffsets at the 10sigma level from the chosen survey.
    def plotposn(self,outname):

        from CATparams import localsource,posnnumsigma

        plotsmade=[]

        sad_res=self.saddata.survey.resolution(self.saddata.pos)
        survey_res=self.survey.resolution(self.saddata.pos)

        # Which survey has the highest resolution??
        if (sad_res[0]*sad_res[1])>(survey_res[0]*survey_res[1]):
            radius=(sad_res[0]*localsource)/2.0
            matchdata=self.getmatches_imagedata(radius)           
        else:
            radius=(survey_res[0]*localsource)/2.0
            matchdata=self.getmatches_external(radius)
        print matchdata
        offsets,beamdists,raoffs,decoffs,raerrors,decerrors=[],[],[],[],[],[]
        for matches in matchdata['matches']:
            # Only single matches wanted.
            if len(matches['match'])==1:
                # Point sources.
                if not numpy.array(matches['source'].isresolved(matchdata['sourcesurvey'],matchdata['sourcesurvey'].sigma)).all():
                    if not numpy.array(matches['match'][0].isresolved(matchdata['matchsurvey'],matchdata['matchsurvey'].sigma)).all():
                        # >15 sigma sources (more position accuracy)
                        print posnnumsigma
                        if matches['source'].pflux>(posnnumsigma*matchdata['sourcesurvey'].sigma):
                            if matches['match'][0].pflux>(posnnumsigma*matchdata['matchsurvey'].sigma):
                                
                                # Include this match in plots
                                # Get distances from image pointings, offsets, raoffsets, decoffsets
                                offsets.append(matches['source'].position.separation(matches['match'][0].position))
                                raoffs.append(matches['source'].position.raoff(matches['match'][0].position))
                                decoffs.append(matches['source'].position.decoff(matches['match'][0].position))
                                beamdists.append(matches['source'].position.separation(self.saddata.pos)/60.0)
                                sourceerror=matches['source'].posnerrorcondon(matchdata['sourcesurvey'].sigma,matchdata['sourcesurvey'])
                                matcherror=matches['match'][0].posnerrorcondon(matchdata['matchsurvey'].sigma,matchdata['matchsurvey'])
                                raerrors.append(sqrt((sourceerror[0]**2.0) + (matcherror[0]**2.0)))
                                decerrors.append(sqrt((sourceerror[1]**2.0) + (matcherror[1]**2.0)))
                                
        self.surveyposnmatch=matchdata

        if len(offsets)>0:
            # Now lets plot plot 1
            pylab.plot(beamdists,offsets,'ro')
            pylab.xlabel('Position offset from phase center (arcmin)')
            pylab.ylabel('Position offset between '+matchdata['sourcesurvey'].name+' and '+matchdata['matchsurvey'].name+' (arcsec)')
            pylab.savefig(outname+'_BEAMDIST.png',format='png')
            pylab.close('all')
            plotsmade=plotsmade+[outname+'_BEAMDIST.png']
            
            # Now lets plot plot 2
            self.ramed=numpy.median(raoffs)
            self.rastd=scipy.std(raoffs,ddof=1)
            self.decmed=numpy.median(decoffs)
            self.decstd=scipy.std(decoffs,ddof=1)
            self.posnum=len(offsets)
            radius=survey_res[0]/1.5
            strra='$<RA>=${0:3.1f} arcsec\n$\sigma(RA)=${1:4.1f} arcsec'.format(self.ramed,self.rastd)
            strdec='$<Dec.>$={0:3.1f} arcsec\n$\sigma(Dec.)=${1:4.1f} arcsec'.format(self.decmed,self.decstd)
            sad_elx,sad_ely=mkellipse(0,0,sad_res[0]/2.0,sad_res[1]/2.0,sad_res[2]+90.0)
            survey_elx,survey_ely=mkellipse(0,0,survey_res[0]/2.0,survey_res[1]/2.0,survey_res[2]+90.0)
            pylab.plot(raoffs,decoffs,'ro')
            pylab.plot(sad_elx,sad_ely,'black')
            pylab.plot(survey_elx,survey_ely,'b--')
            pylab.xlim(-radius,radius)
            pylab.ylim(-radius,radius)
            pylab.axvline(color='black')
            pylab.axhline(color='black')
            pylab.text(-0.9*radius,0.8*radius,strdec)
            pylab.text(0.4*radius,-0.9*radius,strra)
            pylab.text(0.2*radius,0.9*radius,self.saddata.survey.name+' beam (black ellipse)',size=10)
            pylab.text(0.2*radius,0.8*radius,self.survey.name+' beam (blue ellipse)',size=10)
            ylab='Declination Offset ('+matchdata['sourcesurvey'].name+'-'+matchdata['matchsurvey'].name+') (arcsec)'
            xlab='RA Offset ('+matchdata['sourcesurvey'].name+'-'+matchdata['matchsurvey'].name+') (arcsec)'
            pylab.ylabel(ylab)
            pylab.xlabel(xlab)
            pylab.savefig(outname+'_XYOFF.png',format='png')
            pylab.close('all')
            plotsmade=plotsmade+[outname+'_XYOFF.png']

            # Now lets plot plot 3
            self.ramed=numpy.median(raoffs)
            self.rastd=scipy.std(raoffs)
            self.decmed=numpy.median(decoffs)
            self.decstd=scipy.std(decoffs)
            self.posnum=len(offsets)
            strra='$<RA>=${0:3.1f}\n$\sigma(RA)=${1:4.2f}'.format(self.ramed,self.rastd)
            strdec='$<Dec.>$={0:3.1f}\n$\sigma(Dec.)=${1:4.2f}'.format(self.decmed,self.decstd)
            sad_elx,sad_ely=mkellipse(0,0,sad_res[0],sad_res[1],sad_res[2]+90.0)
            survey_elx,survey_ely=mkellipse(0,0,survey_res[0],survey_res[1],survey_res[2]+90.0)
            pylab.errorbar(raoffs,decoffs,yerr=decerrors,xerr=raerrors,fmt=None,ecolor='black')
            pylab.plot(sad_elx,sad_ely,'black')
            pylab.plot(survey_elx,survey_ely,'b--')
            pylab.xlim(-1.0*radius,radius)
            pylab.ylim(-1.0*radius,radius)
            pylab.axvline(color='black')
            pylab.axhline(color='black')
            pylab.text(-0.9*radius,0.8*radius,strdec)
            pylab.text(0.5*radius,-0.9*radius,strra)
            pylab.text(0.2*radius,0.9*radius,self.saddata.survey.name+' beam (black ellipse)',size=10)
            pylab.text(0.2*radius,0.8*radius,self.survey.name+' beam (blue ellipse)',size=10)
            ylab='Declination Offset ('+matchdata['sourcesurvey'].name+'-'+matchdata['matchsurvey'].name+') (arcsec)'
            xlab='RA Offset ('+matchdata['sourcesurvey'].name+'-'+matchdata['matchsurvey'].name+') (arcsec)'
            pylab.ylabel(ylab)
            pylab.xlabel(xlab)
            pylab.savefig(outname+'_XYOFF_ERRORS.png',format='png')
            pylab.close('all')
            plotsmade=plotsmade+[outname+'_XYOFF_ERRORS.png']

            print '**** Made '+self.survey.name+' position comparison plots: '+str(plotsmade)
            self.doneposn=True
        else:
            print '**** No valid matches found with '+self.survey.name+'so cannot do position comparison.'

        return(plotsmade)


    #
    # Look in the array self.compareinfo to decide which plots to make
    # compareinfo[0]= Do position comparison
    # compareinfo[1]= Do flux density comparison
    # compareinfo[2]= Make annotation file
    #
    def docomparison(self,outputroot):


        plotsmade=[]
        
        if (numpy.array(self.compareinfo)>0).any():
            #Position comparison
            if self.compareinfo[0]==1:
                outname=outputroot+self.survey.name+'_POSN'
                plotsmade=plotsmade+self.plotposn(outname)
            #Flux comparison
            if self.compareinfo[1]==1:
                outname=outputroot+self.survey.name+'_FLUX'
                plotsmade=plotsmade+self.plotflux(outname)
            #Annotation file
            if self.compareinfo[2]==1:
                outname=outputroot+self.survey.name
                plotsmade=plotsmade+self.makeann(outname,external_cat.colours.next())
                
  
        return plotsmade




UNITSRE = re.compile('^Fluxes expressed in units of *([A-Za-z/]+)$')
DATARE = re.compile(' *[0-9]')
SEPARATOR = ' '
PLACEHOLDER = '?'
ALLFIELDS = ['ra', 'dec', 'num', 'rah', 'ram', 'ras',
             'decd', 'decm', 'decs',
             'converge', 'peak', 'peakerr', 'flux', 'fluxerr', 'dx', 'dy',
             'maj', 'min', 'pa',
             'dmaj', 'dmin', 'dpa',
             'majfit', 'minfit', 'pafit',
             'majdec', 'mindec', 'padec',
             'majlow', 'minlow', 'palow',
             'majhi', 'minhi', 'pahi',
             'xpix', 'ypix', 'maxresid', 'mosaic']


ALLFORMATS = ['%09.5f', '%+10.5f', '%5d', '%02d', '%02d', '%5.2f', '%3s', '%02d', '%4.1f',
             '%1s', '%10.5f', '%9.5f', '%10.5f', '%9.5f', '%4.2f', '%4.2f',
             '%5.2f', '%5.2f', '%6.2f',
             '%4.2f', '%4.2f', '%2d',
             '%5.2f', '%5.2f', '%6.2f',
             '%5.2f', '%5.2f', '%6.2f',
             '%5.2f', '%5.2f', '%6.2f',
             '%5.2f', '%5.2f', '%6.2f',
             '%6.1f', '%6.1f', '%6.4f', '%s']

ALLNULLS =  [None, None, None, None, None, None, None, None, None,
             '%c', None, '%c', None, '%c', '%c', '%c',
             None, None, None,
             '%c', '%c', None,
             None, None, None,
             '%c', '%c', '%c',
             '%c', '%c', '%c',
             '%c', '%c', '%c',
             None, None, '%c', None]


#Given two external survey objects, find the matches <offset
#in data2 for each source in data1.
def getmatches(data1,data2,offset):
  matchdata={'sourcesurvey': data1.survey,
             'matchsurvey': data2.survey,
             'matches': []}
  for sourcenum,source in enumerate(data1.catdata):
    matchdata['matches'].append({'source':source, 'match':[]})
    for source2 in data2.catdata:
      thissep=source.position.separation(source2.position)
      if thissep<offset:
        matchdata['matches'][sourcenum]['match'].append(source2)

  return(matchdata)

  

def readunits(file):
  line = file.readline()
  while len(line) != 0:
    unitsmatch = UNITSRE.match(line)
    if unitsmatch != None:
      return unitsmatch.group(1)
    line = file.readline()
  if len(line) == 0:
    raise 'Could not find flux units line'

def convertunits(str):
  if str == 'milliJY/BEAM':
    return 1.0
  elif str == 'microJY/BEAM':
    return 1.0e-3
  elif str == 'JY/BEAM':
    return 1.0e3
  elif str == 'RATIO':
    return 1.0
  elif str == 'UNDEFINE':
    return 1.0
  else:
    raise 'Unknown flux units of "%s"' % str

class source:
  def __init__(self,parent):
    self.parent = parent
    self.mosaic = self.parent.MOSAIC
  def __firsthalf__(self, line):
    self.num = int(string.lstrip(line[0:5]))
    if line[6] != ' ':
      self.converge  = line[6]
    else:
      self.converge = None
    self.peak = float(string.lstrip(line[7:16]))
    self.peak = self.peak*self.parent.UNITSFACTOR
    self.peakerr = float(string.lstrip(line[17:24]))
    self.peakerr = self.peakerr*self.parent.UNITSFACTOR
    self.flux = float(string.lstrip(line[25:34]))
    self.flux = self.flux*self.parent.UNITSFACTOR
    if line[36] != '*':
      self.fluxerr = float(string.lstrip(line[35:42]))
      self.fluxerr = self.fluxerr*self.parent.UNITSFACTOR
    else:
      self.fluxerr = None
# Next 5 lines Added by tmauch on 21-2-2001
# Will multiply flux by scaling factor if
# VSAD requires it.    
    if self.converge == '*':
      self.peak = self.peak*1000
      self.peakerr = self.peakerr*1000
      self.flux = self.flux*1000
      try:
        self.fluxerr = self.fluxerr*1000
      except:
        self.fluxerr = None
    self.rah = int(string.lstrip(line[43:46]))
    self.ram = int(string.lstrip(line[46:49]))
    self.ras = float(string.lstrip(line[49:57]))
    self.ra = hmsra2ra(self.rah, self.ram, self.ras)
    self.decd = string.lstrip(line[57:61])
    self.decsign = self.decd[0]
    self.decm = int(string.lstrip(line[61:64]))
    self.decs = float(string.lstrip(line[64:71]))
    self.dec = dmsdec2dec(self.decsign, float(self.decd), self.decm, self.decs)
    self.racos = ra2racos(self.ra, self.dec)
    try:
      self.dx = float(string.lstrip(line[72:80]))
      self.dy = float(string.lstrip(line[80:88]))
    except:
      self.dx = None
      self.dy = None
    self.maj = float(string.lstrip(line[89:98]))
    self.min = float(string.lstrip(line[98:106]))
    self.pa = float(string.lstrip(line[106:112]))
    try:
      self.dmaj = float(string.lstrip(line[113:120]))
      self.dmin = float(string.lstrip(line[120:127]))
    except:
      self.dmaj = None
      self.dmin = None
    try:
      self.dpa = float(string.lstrip(line[127:131]))
    except:
      self.dpa = None
  def __secondhalf__(self, line):
    if int(string.lstrip(line[0:5])) != self.num:
      raise 'Second half does not match first half'
    self.majfit = float(string.lstrip(line[5:14]))
    self.minfit = float(string.lstrip(line[14:22]))
    self.pafit = float(string.lstrip(line[22:29]))
    try:
      self.majdec = float(string.lstrip(line[29:40]))
      self.mindec = float(string.lstrip(line[40:48]))
      self.padec = float(string.lstrip(line[48:55]))
    except:
      self.majdec = None
      self.mindec = None
      self.padec = None
    try:
      self.majlow = float(string.lstrip(line[55:66]))
      self.minlow = float(string.lstrip(line[66:74]))
      self.palow = float(string.lstrip(line[74:81]))
    except:
      self.majlow = None
      self.minlow = None
      self.palow = None
    try:
      self.majhi = float(string.lstrip(line[81:92]))
      self.minhi = float(string.lstrip(line[92:100]))
      self.pahi = float(string.lstrip(line[100:107]))
    except:
      self.majhi = None
      self.minhi = None
      self.pahi = None
    self.xpix = float(string.lstrip(line[107:114]))
    self.ypix = float(string.lstrip(line[114:121]))
    try:
        self.maxresid = float(string.lstrip(line[122:131]))
    except:
        self.maxresid = None
  def __reprattr__(self, i):
    value = getattr(self, self.parent.fields[i])
    if value == None:
      return self.parent.nulls[i] % self.parent.placeholder
    else:
      return self.parent.formats[i] % value
  def __repr__(self):
    result = self.__reprattr__(0)
    for i in range(1, len(self.parent.fields)):
      result = result + self.parent.separator + self.__reprattr__(i)
    return result

class file:
  fields = ALLFIELDS
  formats = ALLFORMATS
  nulls = ALLNULLS
  separator = SEPARATOR
  placeholder = PLACEHOLDER
  def __init__(self, file):
    self.FILEPATH = file.name
    self.FILENAME = os.path.split(self.FILEPATH)[1]
    (self.MOSAIC, self.FILEEXT) = os.path.splitext(self.FILENAME)
    self.UNITS = readunits(file)
    self.UNITSFACTOR = convertunits(self.UNITS)
    self.data = []
    line = file.readline()
    half = 1
    while len(line) != 0:
      if DATARE.match(line) != None:
        if line[0:5] == '    1':
          if half == 2:
            break
          else:
            half = 2
        s = source(self)
        s.__firsthalf__(line)
        self.data.append(s)
      line = file.readline()
    i = 0
    while len(line) != 0:
      if DATARE.match(line) != None:
        s = self.data[i]
        s.__secondhalf__(line)
        i = i + 1
      line = file.readline()


class sad_catalogue:

    def __init__(self,fitsfilename,aipsimagedata,rmsimagedata,meanimagedata,saddata,maxwidth,imrms,fluxlimit=0.0):

      self.fitsfile=fitsfilename
      self.aipsdata=aipsimagedata
      self.saddata=saddata
      self.rmsdata=rmsimagedata
      self.meandata=meanimagedata
      self.rawdata=[]
      image=makeaipsimage(aipsimagedata)
      imhead=image.header()
      # store the survey info
      
      self.survey = survey(imhead['teles'],imhead['beamMaj']*3600.0, imhead['beamMin']*3600.0, imhead['beamPA'],
                           in_sigma=imrms*1000.0,in_frequency=imhead['crval'][2]/1e6 if imhead['ctype'][2] in ['FREQ    ','SPECLNMF '] else imhead['crval'][3]/1e6,
                           in_telescope=imhead['teles'],in_id='SAD')

      # Store the image data (pointing center, image radius)
      self.pos = position(imhead['crval'][0],imhead['crval'][1],'dd')
      sizedata = imrange(image)
      rasize = fabs(sizedata[0]-sizedata[1])
      decsize = fabs(sizedata[2]-sizedata[3])
      self.imsize = max(rasize,decsize)*3600.0/2

      
      for line in saddata.data:
        source=str(line).split()
        discard=False
        # Discard sources that are too wide.
        if (float(source[16])>self.survey.resolution(self.pos)[0]*maxwidth) or (float(source[17])>self.survey.resolution(self.pos)[1]*maxwidth):
          discard=True
        # Check that fitted position is 'CATparams.edgepix' pixels inside the image.
        elif (float(source[34])+CATparams.edgepix>float(imhead['inaxes'][0])) or (float(source[35])+CATparams.edgepix>float(imhead['inaxes'][1])):
          discard=True
        elif (float(source[34])-CATparams.edgepix<0.0) or (float(source[35])-CATparams.edgepix<0.0):
          discard=True
        # Make sure that interpolation works on this source.
        elif isnan(pixinterp(float(source[34]),float(source[35]),rmsimagedata)):
          discard=True
        # Add this to the raw data array if it passes the tests.
        if discard==False:
          self.rawdata.append(str(line).split())

      self.numsources=len(self.rawdata)
      print '****  Found '+str(self.numsources)+' sources.'
      # Check we actually have some sources to store.
      if self.numsources==0:
        print '**** WARNING: No useable sources in SAD output!' 
        self.sourcedata=[]
        return
      
      # Now set up an array of radio_source objects- which will be the
      # main storage for the catalogue data.
      tempsourcedata=self.survey.makearray(self)

      # Check which sources are above the flux limit specified in CATparams
      # and discard those that arent.
      self.sourcedata=[]
      rawarray=self.rawdata[:]
      self.rawdata=[]
      for sourcenum,source in enumerate(tempsourcedata):
          if source.flux>=fluxlimit:
              self.sourcedata.append(source)
              self.rawdata.append(rawarray[sourcenum])
      self.numfluxlimit=len(self.sourcedata)
    #
    # Print out the raw catalogue to disk in 'filename'.
    #
    def writeraw(self,filename):
      rawfile = open(filename,mode='w')
      for i in self.saddata.data:
        rawfile.write(str(i))
        rawfile.write('\n')
      rawfile.close()
      
      return
      
    #
    # Print out a formatted catalogue to disk in 'filename'
    #
    def writecat(self,filename):
      #catfilename = random_filename(length=8,suffix='.cat')
      catfile = open(filename,mode='w')
      for num,source in enumerate(self.rawdata):
        # Get the local rms noise.
        localrms = pixinterp(float(source[34]),float(source[35]),self.rmsdata)*1000.0
        localmean = pixinterp(float(source[34]),float(source[35]),self.meandata)*1000.0
        sourcedata = self.sourcedata[num]

        # Sad Source number
        catfile.write('{0:5d}  '.format(int(source[2])))

        # Right Ascension
        catfile.write('{0} {1} {2} '
                      .format(*sourcedata.position.ra_hms))

        # Declination
        catfile.write('{0:>3s} {1} {2}  '
                      .format(*sourcedata.position.dec_hms))

        # RA and Dec. Error
        catfile.write('{0:05.2f} {1:05.2f}   '
                      .format(*sourcedata.posnerrorcondon(localrms,self.survey))) 

        peakerr,interr = sourcedata.fluxerrorcondon(localrms,self.survey)

        # Peak Flux and its error
        catfile.write('{0:10.5f} {1:9.5f}  '
                      .format(sourcedata.pflux,peakerr))
        
        # Total Flux and its error
        print sourcedata.flux, interr, localrms
        catfile.write('{0:10.5f} {1:9.5f}    '
                      .format(sourcedata.flux,interr))

        # Fitted source size
        catfile.write('{0:6.2f} {1:6.2f} {2:6.2f}  '
                      .format(sourcedata.maj,sourcedata.min,sourcedata.pa))
        
        # Fitted source size errors
        catfile.write('{0:5.2f} {1:5.2f} {2:2d}    '
                      .format(*sourcedata.sizeerrorcondon(localrms,self.survey)))

        # Pixel Position in image
        catfile.write('{0:6.1f} {1:6.1f}  '
                      .format(float(source[34]),float(source[35])))

        # Image name
        catfile.write(self.fitsfile)

        # Local mean and local rms
        catfile.write('{0:9.4f} {1:7.4f}\n'
                      .format(localrms,localmean))

                            
      catfile.close()
      #os.rename(catfilename,filename)

    #
    # Print out a kvis annotation file to disk in 'filename'
    #
    def makeann(self,filename,colour='BLUE'):
      annfile = open(filename,mode='w')

      annfile.write('COLOUR '+colour+'\n')
      annfile.write('PA SKY\n')
      annfile.write('FONT hershey10\n')

      for source in self.sourcedata:
        ra = source.position.ra_dd
        dec = source.position.dec_dd
        maj = source.maj/3600.
        min = source.min/3600.
        num = source.id
        annfile.write('ELLIPSE {0} {1} {2} {3} {4}\n'.format(ra,dec,maj,min,source.pa))
        annfile.write('TEXT {0} {1} {2}\n'.format(ra,dec,num))
      annfile.close()

    #
    # Sort out a vo table from the catalogue
    #
    def makevo(self,votable,description='Unknown catalogue'):

        table=Table(votable,ID=self.survey.name,name=self.survey.name+' Data')
        table.description=description
        table.format='tabledata'

        metadata,formatdata=self.survey.vo_metadata(self,votable)
        
        table.fields.extend(metadata)
        table.create_arrays(len(formatdata))
        for num,thisitem in enumerate(formatdata):
            table.array[num]=tuple(thisitem)
        return table
          
    #
    # Print a catalogue header to 'filename' in 'outputdir'.
    #
    def printheader(self,filename):

        tempfile = open(filename, 'w')
        tempfile.write('Column description of final formatted catalogue (.cat & .vo files):\n')
        tempfile.write('Column 1: Source number.\n')
        tempfile.write('Columns 2+3+4: J2000 Right Ascension (hh mm ss.ss).\n')
        tempfile.write('Columns 5+6+7: J2000 Declination (dd \'\' \"\".\").\n')
        tempfile.write('Column 8: Error in Right Ascension (acrsec).\n')
        tempfile.write('Column 9: Error in Declination (arcsec).\n')
        tempfile.write('Column 10: Fitted peak flux (mJy/beam).\n')
        tempfile.write('Column 11: Statistical error of fitted peak flux (mJy/beam).\n')
        tempfile.write('Column 12: Total flux density (mJy).\n')
        tempfile.write('Column 13: Statistical error of total flux density (mJy).\n')
        tempfile.write('Column 14: Fitted major axis of elliptical Gaussian (arcsec).\n')
        tempfile.write('Column 15: Fitted minor axis of elliptical Gaussian (arcsec).\n')
        tempfile.write('Column 16: Fitted position angle of elliptical Gaussian (deg. E. of N.).\n')
        tempfile.write('Column 17: Statistical error of fitted major axis (arcsec).\n')
        tempfile.write('Column 18: Statistical error of fitted minor axis (arcsec).\n')
        tempfile.write('Column 19: Statistical error of fitted position angle (deg.).\n')
        tempfile.write('Column 20: X-pixel position in image of fitted position.\n')
        tempfile.write('Column 21: Y-pixel position in image of fitted position.\n')
        tempfile.write('Column 22: Image name.\n')
        tempfile.write('Column 23: Local rms in image at source position (mJy/beam).\n')
        tempfile.write('Column 24: Local mean in image at source psotion (mJy/beam).\n')
        tempfile.close()

        return(filename)

    #
    #Print some image statistics to and the
    #results of comparison surveys in 'externalsurveys'
    #to 'filename' in 'outputdir'
    #
    def getimstats(self,filename,externalsurveys=[],imagelocaldr=[0,0]):

        tempfile=open(filename,'w')
      	tempfile.write('\n')
        tempfile.write('Catalogue made from '+self.survey.telescope+' image on '+time.strftime('%c',time.localtime())+'\n')
        tempfile.write('Image Name: '+self.fitsfile+'\n')
	tempfile.write('Observed RA: ' + self.pos.ra_hms[0] + ' ' + self.pos.ra_hms[1] + ' ' + self.pos.ra_hms[2]+ '\n')
	tempfile.write('Observed Dec: ' + self.pos.dec_hms[0] + ' ' + self.pos.dec_hms[1] + ' ' + self.pos.dec_hms[2] + '\n')
        tempfile.write('Observing Frequency: '+str(self.survey.frequency)+' MHz\n')
        tempfile.write('Beam Size: {0:3.1f}x{1:3.1f} arcsec at PA: {2:3.0f} deg.\n'.format(*self.survey.resolution(self.pos)))
        tempfile.write('Image rms: {0:10.7f} mJy/beam\n'.format(self.survey.sigma))
        tempfile.write('Maximum pixel in image: {0:12.7f} mJy/beam\n'.format(maxpix(self.aipsdata,pb=CATparams.pbcorrect)*1000.0))
	tempfile.write('Global image dynamic range is {0:6.1f}\n'.format(maxpix(self.aipsdata,pb=CATparams.pbcorrect)*1000.0/self.survey.sigma))
        if (imagelocaldr[1]>0):
            # hrk changed 04.11.2011
            tempfile.write('Local dynamic range around '+str(imagelocaldr[1])+' bright sources is {0:6.1f}\n'.format(imagelocaldr[0]))
        else:
            tempfile.write('No sources to determine local dynamic range. Used global DR for rms model.\n')
        tempfile.write('Found '+str(self.numsources)+' sources down to '+str(CATparams.cutoff)+' sigma.\n')
	tempfile.write('Final catalogue contains '+str(self.numfluxlimit)+' sources down to '+str(CATparams.fluxcutoff)+' mJy.\n')
        tempfile.write('\n ############################ \n\n')

        for i in externalsurveys:
          if i.donetwofreq:
            tempfile.write(i.twofreqnames + ' COMPARISON:\n\n')
            tempfile.write(str(i.twofreqnum) + ' sources found in image area. \n')
            tempfile.write(i.twofreqnames + ' log average flux offset: {0:5.3f}\n'.format(i.twofreqfluxoffall))
            tempfile.write(i.twofreqnames + ' log average flux offset for point sources: {0:5.3f}\n'.format(i.twofreqps_fluxoff))
            tempfile.write('\n')
            tempfile.write('\n')
        
        imhead=makeaipsimage(self.aipsdata).header()
        rapix=fabs(imhead['cdelt'][0]*3600.0)
        decpix=fabs(imhead['cdelt'][1]*3600.0)
        for i in externalsurveys:
          thisname=i.survey.name
          tempfile.write(thisname+' COMPARISON:\n\n')
          tempfile.write(str(len(i.catdata))+' sources found in image area from '+thisname+'\n')
          tempfile.write(thisname+' beam Size: {0:3.1f}x{1:3.1f} arcsec at PA: {2:3.0f} deg.\n\n'.format(*i.survey.resolution(self.pos)))
          if i.doneflux:
            tempfile.write(thisname+' log average flux offset: {0:5.3f}\n'.format(i.fluxoff))
            tempfile.write(thisname+' log average flux offset for point sources: {0:5.3f}\n'.format(i.psfluxoff))
            tempfile.write('\n')
          if i.doneposn:
            tempfile.write(thisname+' position offset statistics for '+str(i.posnum)+' '+ str(CATparams.posnnumsigma) +'sigma point sources:\n')
            tempfile.write(i.surveyposnmatch['sourcesurvey'].name+' - ' +i.surveyposnmatch['matchsurvey'].name +' median RA offset: {0:5.2f} rms: {1:5.3f} (arcsec), offset: {2:5.2f} rms: {3:5.3f} (pixels)\n'.format(i.ramed,i.rastd,i.ramed/rapix,i.rastd/rapix))
            tempfile.write(i.surveyposnmatch['sourcesurvey'].name+' - ' +i.surveyposnmatch['matchsurvey'].name +' median Dec offset: {0:5.2f} rms: {1:5.3f} (arcsec), offset: {2:5.2f} rms: {3:5.3f} (pixels)\n'.format(i.decmed,i.decstd,i.decmed/decpix,i.decstd/decpix))
            tempfile.write('\n')
          tempfile.write('\n')

        tempfile.close()
        return(filename)


    #
    #Write a brief description of the filenames produced.
    #
    def printfiles(self,filename,externalsurveys,sad=0,cat=0,raw=0,ann=0,vo=0):
      
      tempfile=open(filename,'w')

      inname=self.fitsfile.rstrip('.FITS')

      tempfile.write('Output Files:\n')
      tempfile.write('The following catalogue files have been produced:\n')

      if sad>0:
        tempfile.write(inname+'.sad: Contains the raw output from the AIPS SAD task\n')
      if raw>0:
        tempfile.write(inname+'.raw: Contains the raw output fromt he AIPS SAD task formatted such that each source fits on one line\n')
      if cat>0:
        tempfile.write(inname+'.cat: Contains the final formatted catalogue in as a text file (column descriptions are below)\n')
      if vo>0:
        tempfile.write(inname+'.vo: Contains the final formatted catalogue as a vo table (can be input into TOPCAT).\n')
      if ann>0:
        tempfile.write(inname+'.ann: Contains the final catalogue as a kvis annotation file. \n')
        
      for i in externalsurveys:
        if i.doneflux:
          thisname=i.survey.name+'_FLUX'
          tempfile.write('\nFiles beginning in '+thisname+':\n')
          tempfile.write(thisname+'_BEAMDIST: Shows the offset in flux density between '+thisname+' and '+self.survey.name+' as a function of distance from the phase center of th input image.\n')
          tempfile.write(thisname+'_COMP: Shows a comparison between the predicted '+i.survey.name+' flux desnity at '+str(self.survey.frequency)+'MHz, using alpha=-0.7.\n')
          tempfile.write(thisname+'_SNUMS: The same as '+thisname+'_COMP, but with the fitted source numbers at the datapoints.\n')
          tempfile.write(thisname+'_PBEAM: Shows the x,y position of the sources in the image with the flux offset in a colour scale.\n')
          tempfile.write(thisname+'_SIHISTO: Shows the spectral index distribution between '+i.survey.name+' and '+self.survey.name+'.\n')
        if i.doneposn:
          thisname=i.survey.name+'_POSN'
          tempfile.write('\nFiles beginning in '+thisname+':\n')
          tempfile.write(thisname+'_BEAMDIST: Shows the position offsets between '+i.survey.name+' and '+self.survey.name+' as a function of distance from the phase center of the input image.\n')
          tempfile.write(thisname+'_XYOFF: Shows the RA and Declination offsets between the two surveys, the beams of the two surveys are also overlaid.\n')

      tempfile.close()
      return(filename)

#
#makefluxfile:
#given a list of fluxes, offsets, and predicted fluxes
#and sprectral indices, produce a flux comparison 
#text file.
#
def makefluxcompfile(ras,decs,beamdists,radists,decdists,fluxorig,freqorig,fluxgmrt,freqgmrt,fluxpredict,alphas,outfile):

  def plotheader(file):
    
    file.write('# The following file contains the raw data from matching the catalogue\n')
    file.write('# to other data for flux density comparison.\n')
    file.write('# Columns are:\n')
    file.write('# 1-3: Right Ascension (J2000 hh mm ss.sss)\n')
    file.write('# 4-6: Declination (J2000 dd mm ss.ss)\n')
    file.write('# 7: Offset of source from image center (arcminutes)\n')
    file.write('# 8: RA offset of source from image center (arcminutes)\n')
    file.write('# 9: Dec. offset of source from image center (arcminutes)\n')
    file.write('# 10: E=Extended Source P=Point Source\n')
    file.write('# 11: Flux density of source in comparison catalogue (mJy)\n')
    file.write('# 12: Frequency of comparison catalogue (MHz)\n')
    file.write('# 13: Flux density of source from catalogue (sum of # of components in Col. 10) (mJy)\n')
    file.write('# 14: Frequency of catalogue\n')
    file.write('# 15: Flux predicted from Col. 11, at catalogue frequency (Col. 14), using SI from Col. 16\n')
    file.write('# 16: Spectral Index (alpha) used to obtain predicted flux in Col. 15 (S~nu^alpha)\n')
    file.write('#\n')

  textout=open(outfile,'w')
  plotheader(textout)
  #Point sources
  psdata=zip(ras[0],decs[0],beamdists[0],radists[0],decdists[0],fluxorig[0],fluxgmrt[0],fluxpredict[0],alphas[0])
  esdata=zip(ras[1],decs[1],beamdists[1],radists[1],decdists[1],fluxorig[1],fluxgmrt[1],fluxpredict[1],alphas[1])
  for data in psdata:
    outstr= '{0:2s} {1:2s} {2:6s} {3:>3s} {4:2s} {5:5s} {6:8.2f} {7:8.2f} {8:8.2f} {9:1s} {10:8.3f} {11:10.1f} {12:8.3f} {13:10.1f} {14:8.3f} {15:5.2f} \n'.format(data[0][0],data[0][1],data[0][2],data[1][0],data[1][1],data[1][2],data[2],data[3],data[4],'P',data[5],freqorig,data[6],freqgmrt,data[7],data[8])
    textout.write(outstr)
  for data in esdata:
    outstr= '{0:2s} {1:2s} {2:6s} {3:>3s} {4:2s} {5:5s} {6:8.2f} {7:8.2f} {8:8.2f} {9:1s} {10:8.3f} {11:10.1f} {12:8.3f} {13:10.1f} {14:8.3f} {15:5.2f} \n'.format(data[0][0],data[0][1],data[0][2],data[1][0],data[1][1],data[1][2],data[2],data[3],data[4],'E',data[5],freqorig,data[6],freqgmrt,data[7],data[8])
    textout.write(outstr)
  textout.close()

  return(outfile+'_FLUXCOMP.txt')

      
#
#fluxcompare:
#given a list of data from external surveys, determine which
#will give the best predicted flux at the obs freq. Then
#if a 2 frequency solution can be found- produce a flux
#comparison using this 2 frequency solution. Also produce
#a text file containing the data from the best flux comparison.
#
def fluxcompare(externallist,outputdir):

  from makeplots import plot2freqpredict

  outfiles=[]

  #Get a list of surveynames
  surveynames=[survey.survey.id for survey in externallist]
  
  #Start with the NVSS this is the basis of the best flux comparison:
  if ('NVSS' in surveynames) and (externallist[surveynames.index('NVSS')].doneflux):
    primarydata=externallist[surveynames.index('NVSS')]
    if 'WENSS' in surveynames:
      extradata=externallist[surveynames.index('WENSS')]
    elif 'SUMSS' in surveynames:
      extradata=externallist[surveynames.index('SUMSS')]
    else:
      extradata=None
  #Otherwise check with WENSS and SUMSS
  elif ('SUMSS' in surveynames) and (externallist[surveynames.index('SUMSS')].doneflux): 
    primarydata=externallist[surveynames.index('SUMSS')]
    if 'FIRST' in surveynames:
      extradata=externallist[surveynames.index('FIRST')]
    else:
      extradata=None
      
  elif ('WENSS' in surveynames) and (externallist[surveynames.index('WENSS')].doneflux):
    primarydata=externallist[surveynames.index('WENSS')]
    if 'FIRST' in surveynames:
      extradata=externallist[surveynames.index('FIRST')]
    else:
      extradata=None

  #Finally check first
  elif ('FIRST' in surveynames) and (externallist[surveynames.index('FIRST')].doneflux):
    primarydata=externallist[surveynames.index('FIRST')]
    extradata=None
  
  #Otherwise nothing to produce!!!
  else:
    print '**** No flux comparison available to produce flux output file!'
    return []

  #We now have primarydata and extradata arrays which we can use to get good predicted fluxes.
  if (not extradata==None):
    # Crossmatch and remake flux comparison plots with 2 comparison surveys.
    outfiles=outfiles+plot2freqpredict(primarydata,extradata,outputdir)
    return outfiles
  else:
    # Just print the predicted fluxes from the main survey
    
    outfiles.append(makefluxcompfile([primarydata.ps_ramatch,primarydata.es_ramatch],[primarydata.ps_decmatch,primarydata.es_decmatch],
                                     [primarydata.ps_beamdists,primarydata.es_beamdists],[primarydata.ps_raoffs,primarydata.es_raoffs],
                                     [primarydata.ps_decoffs,primarydata.es_decoffs],[primarydata.ps_flux_orig,primarydata.es_flux_orig],
                                     primarydata.survey.frequency,[primarydata.ps_flux,primarydata.es_flux],primarydata.saddata.survey.frequency,
                                     [primarydata.ps_flux_predict,primarydata.es_flux_predict],
                                     [[-0.7]*len(primarydata.ps_flux),[-0.7]*len(primarydata.es_flux)],outputdir+primarydata.survey.name))
    
  return outfiles



      
