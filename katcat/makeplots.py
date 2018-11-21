import matplotlib
matplotlib.use('Agg')
from math import *
import numpy
import pylab
import scipy
from astrop import *
from catalogue import *

# Do a flux vs flux plot for extended and pointsources between two surveys.
def plotfluxes(ps_predict,ps_flux,es_predict,es_flux,ps_meanoff,meanoff,name1,name2,outname):

    allflux=numpy.append(ps_flux,es_flux)
    allfluxpredict=numpy.append(ps_predict,es_predict)

    pylab.plot([0.1,1,10000000],[10**(-1-ps_meanoff),10**(0-ps_meanoff),10**(7-ps_meanoff)],'r-')
    pylab.plot([0.1,1,10000000],[10**(-1-meanoff),10**(0-meanoff),10**(7-meanoff)],'k-')
    pylab.plot([0.1,1,10000000],[0.1,1,10000000],'b-')
    pylab.plot(es_predict,es_flux,'go')
    pylab.plot(ps_predict,ps_flux,'ro')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlim(numpy.amin(numpy.append(allfluxpredict,allflux))/2.0,numpy.amax(numpy.append(allfluxpredict,allflux))*2.0)
    pylab.ylim(pylab.xlim())
    strflux='Log mean flux density offset is {0:5.3f} (black line)\nPoint source only offset is {1:5.3f} (red line)\n{2:d} point sources (red points)\n{3:d} other sources (green points)'.format(meanoff,ps_meanoff,len(ps_flux),len(es_flux))
    pylab.text(10**(log10(pylab.xlim()[0])+(0.05*(log10(pylab.xlim()[1])-log10(pylab.xlim()[0])))),10**(log10(pylab.ylim()[1])-0.16*(log10(pylab.ylim()[1])-log10(pylab.ylim()[0]))),strflux)
    pylab.ylabel('Flux density '+name1+' (mJy)')
    pylab.xlabel('Predicted flux density from '+name2+' (mJy)')
    pylab.savefig(outname+'_COMP.png',format='png')
    pylab.close('all')
    return outname+'_COMP.png'

# Exactly the same as plotfluxes above- but with errorbars nstead of points
def plotfluxeserror(predict,flux,errpredict,errflux,ps_meanoff,meanoff,name1,name2,outname):

    pylab.plot([0.1,1,10000000],[10**(-1-ps_meanoff),10**(0-ps_meanoff),10**(7-ps_meanoff)],'r-')
    pylab.plot([0.1,1,10000000],[10**(-1-meanoff),10**(0-meanoff),10**(7-meanoff)],'k-')
    pylab.plot([0.1,1,10000000],[0.1,1,10000000],'b-')
    pylab.errorbar(predict,flux,xerr=errpredict,yerr=errflux,fmt=None,ecolor='black')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlim(numpy.amin(numpy.append(predict,flux))/2.0,numpy.amax(numpy.append(predict,flux))*2.0)
    pylab.ylim(pylab.xlim())
    strflux='Log mean flux density offset is {0:5.3f} (black line)\nPoint source only offset is {1:5.3f} (red line)\n'.format(meanoff,ps_meanoff)
    pylab.text(10**(log10(pylab.xlim()[0])+(0.05*(log10(pylab.xlim()[1])-log10(pylab.xlim()[0])))),10**(log10(pylab.ylim()[1])-0.16*(log10(pylab.ylim()[1])-log10(pylab.ylim()[0]))),strflux)
    pylab.ylabel('Flux density '+name1+' (mJy)')
    pylab.xlabel('Predicted flux density from '+name2+' (mJy)')
    pylab.savefig(outname+'_ERRORS.png',format='png')
    pylab.close('all')
    return outname+'_ERRORS.png'



# Plot a spectral index histogram between flux1 and flux2 at freq1 and freq2.
def plotsihisto(flux1,flux2,freq1,freq2,outname):

    si=numpy.log(flux1/flux2)/numpy.log(freq1/freq2)
    pylab.hist(si,int(len(si)/5.0))
    pylab.xlim(-3,3)
    pylab.ylim(ymin=0)
    pylab.text(-2.9,1,'Mean alpha is {0:4.2f}'.format(scipy.mean(si)))
    pylab.xlabel('Spectral Index between '+str(freq1)+'MHz and '+str(freq2)+'MHz.')
    pylab.ylabel('Number')
    pylab.savefig(outname+'_SIHISTO.png',format='png')
    pylab.close('all')
    return outname+'_SIHISTO.png'

# Plot the flux offset as a function of ra,dec distance from the primary beam.
def plotfluxpboff(raoffs,decoffs,fluxoff,aipsimage,outname):

     #image=makeaipsimage(aipsimage)
     width=(imrange(aipsimage)[3]-imrange(aipsimage)[2])/2.0
     widthamin=width*60.0

     pylab.xlim(-widthamin,widthamin)
     pylab.ylim(-widthamin,widthamin)
     plt=pylab.scatter(raoffs,decoffs,c=fluxoff,marker='o')
     pylab.xlabel('RA offset (arcmin.)')
     pylab.ylabel('Declination offset (arcmin.)')
     cbar=pylab.colorbar(plt)
     cbar.set_label('log flux offset')
     pylab.savefig(outname+'_PBEAM.png',format='png')
     pylab.close('all')
     return outname+'_PBEAM.png'

# Plot the flux offset as a function of radial distance from the primary beam.
def plotfluxpbrad(ps_offsets,ps_fluxoff,es_offsets,es_fluxoff,outname):

    figure=pylab.figure()
    ax=figure.add_subplot(111)
    ax.plot(ps_offsets,ps_fluxoff,'ro')
    ax.plot(es_offsets,es_fluxoff,'go')
    ax.axhline(-0.7*numpy.log10(1285./843.),color='r',linestyle='--')
    ax.text(0.01,0.95,'Point sources (red circles)',transform = ax.transAxes)
    ax.text(0.01,0.91,'Other sources (green circles)',transform = ax.transAxes)
    ax.set_xlabel('Position offset from phase center (arcmin)')
    ax.set_ylabel('log flux offset')
    pylab.savefig(outname+'_BEAMDIST.png',format='png')
    pylab.close(figure)
    return outname+'_BEAMDIST.png'

# Plot the flux vs predicted flux with source numbers from the catalogue on the plot.
def plotfluxsnums(flux,flux_predict,ids,name1,name2,outname):

    pylab.plot([0.1,1,1000000],[0.1,1,1000000],'b-')
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.xlim(numpy.amin(numpy.append(flux_predict,flux))/2.0,numpy.amax(numpy.append(flux_predict,flux))*2.0)
    pylab.ylim(pylab.xlim())
    for x,y,lab in zip(flux_predict,flux,ids):
        pylab.text(x,y,lab,size=5,alpha=0.4)
    pylab.ylabel('Flux density '+name1+' (mJy)')
    pylab.xlabel('Predicted flux density from '+name2+' (mJy)')
        
    pylab.savefig(outname+'_SNUMS.png',format='png')
    pylab.close('all')
    return outname+'_SNUMS.png'

     
def plot2freqpredict(primarydata,extradata,outputdir):

    if primarydata.matchdata['sourcesurvey'].id=='SAD':
        comparesurveys=primarydata.matchdata['matchsurvey'].name+'+'+extradata.survey.name
        sadname=primarydata.matchdata['sourcesurvey'].name
    else:
        comparesurveys=primarydata.matchdata['sourcesurvey'].name+'+'+extradata.survey.name
        sadname=primarydata.matchdata['matchsurvey'].name
        
    prim_sourcesurvey=primarydata.matchdata['sourcesurvey']

    # Some simplifying variables
    if prim_sourcesurvey.id=='SAD':
        gmrt_freq=prim_sourcesurvey.frequency
        freq1=primarydata.matchdata['matchsurvey'].frequency
        freq2=extradata.survey.frequency
    else:
        gmrt_freq=primarydata.matchdata['matchsurvey'].frequency
        freq1=prim_sourcesurvey.frequency
        freq2=extradata.survey.frequency

    # Match the extra data to the 'source' array in primarydata
    # Which has the worst resolution?
    prim_res=prim_sourcesurvey.resolution(primarydata.saddata.pos)
    extr_res=extradata.survey.resolution(extradata.saddata.pos)

    if (prim_res[0]*prim_res[1])>(extr_res[0]*extr_res[1]):
        radius=prim_res[0]/2.0
        matchdata=getmatches(primarydata,extradata,radius)
    else:
        radius=extr_res[0]/2.0
        matchdata=getmatches(extradata,primarydata,radius)

    prim_ps_source=[ps['source'] for ps in primarydata.ps_matches]
    prim_all_source=[all['source'] for all in primarydata.matches]
    ps_gmrtflux,es_gmrtflux,ps_predictflux,es_predictflux=[],[],[],[]
    ps_beamdists,es_beamdists,ps_raoffs,ps_decoffs,es_raoffs,es_decoffs=[],[],[],[],[],[]
    ps_ids,es_ids,ps_si,es_si=[],[],[],[]
    ps_ra,ps_dec,es_ra,es_dec=[],[],[],[]
    ps_origflux,es_origflux=[],[]
    # Now get the flux arrays sorted out
    if matchdata['sourcesurvey']==prim_sourcesurvey:
        # primarydata is the 'source' array.
        for matches in matchdata['matches']:
            #Point source in all 3 surveys??
            if len(matches['match'])==1 and not numpy.array(matches['match'][0].isresolved(matchdata['matchsurvey'],matchdata['matchsurvey'].sigma)).all() and matches['source'] in prim_ps_source:
                #YES this is a point source!!
                if prim_sourcesurvey.id=='SAD':
                    ps_gmrtflux.append(matches['source'].flux)
                    ps_gmrterror.append(matches['source'].fluxerrorcondon(matchdata['sourcesurvey'].sigma,matchdata['sourcesurvey']))
                    ps_predictflux.append(
                        fluxpredict(matches['match'][0].flux,primarydata.ps_flux_orig[prim_ps_source.index(matches['source'])],
                                    freq2,freq1,gmrt_freq))
                    ps_si.append(getsi(matches['match'][0].flux,primarydata.ps_flux_orig[prim_ps_source.index(matches['source'])],
                                 freq2,freq1))
                    ps_origflux.append(primarydata.ps_flux_orig[prim_ps_source.index(matches['source'])])
                else:
                    ps_gmrtflux.append(primarydata.ps_flux[prim_ps_source.index(matches['source'])])
                    ps_predictflux.append(
                        fluxpredict(matches['source'].flux,matches['match'][0].flux,freq1,freq2,gmrt_freq))
                    ps_origflux.append(matches['source'].flux)
                    ps_si.append(getsi(matches['source'].flux,matches['match'][0].flux,freq1,freq2))
                ps_ra.append(matches['source'].position.ra_hms)
                ps_dec.append(matches['source'].position.dec_hms)
                ps_beamdists.append(matches['source'].position.separation(primarydata.saddata.pos)/60.0)
                ps_raoffs.append(matches['source'].position.raoff(primarydata.saddata.pos)/60.0)
                ps_decoffs.append(matches['source'].position.decoff(primarydata.saddata.pos)/60)
                ps_ids.append(primarydata.ids[prim_all_source.index(matches['source'])])

            elif len(matches['match'])>0 and matches['source'] in prim_all_source:
                #NO this is not a point source!!
                thisflux=0.0
                for thismatch in matches['match']:
                    thisflux=thisflux+thismatch.flux
                if prim_sourcesurvey.id=='SAD':
                    es_gmrtflux.append(matches['source'].flux)
                    es_predictflux.append(
                        fluxpredict(thisflux,primarydata.flux_orig[prim_all_source.index(matches['source'])],
                                    freq2,freq1,gmrt_freq))
                    es_origflux.append(primarydata.flux_orig[prim_all_source.index(matches['source'])])
                    es_si.append(getsi(thisflux,primarydata.flux_orig[prim_all_source.index(matches['source'])],
                                       freq2,freq1))
                else:
                    es_gmrtflux.append(primarydata.flux[prim_all_source.index(matches['source'])])
                    es_predictflux.append(fluxpredict(matches['source'].flux,thisflux,freq1,freq2,gmrt_freq))
                    es_origflux.append(matches['source'].flux)
                    es_si.append(getsi(matches['source'].flux,thisflux,freq1,freq2))
                es_ra.append(matches['source'].position.ra_hms)
                es_dec.append(matches['source'].position.dec_hms)
                es_beamdists.append(matches['source'].position.separation(primarydata.saddata.pos)/60.0)
                es_raoffs.append(matches['source'].position.raoff(primarydata.saddata.pos)/60.0)
                es_decoffs.append(matches['source'].position.decoff(primarydata.saddata.pos)/60.0)
                es_ids.append(primarydata.ids[prim_all_source.index(matches['source'])])
    else:
        # primarydata is the 'match' array.
        for matches in matchdata['matches']:
            #Point source in all 3 surveys??
            if len(matches['match'])==1 and not numpy.array(matches['source'].isresolved(matchdata['sourcesurvey'],matchdata['sourcesurvey'].sigma)).all() and matches['match'][0] in prim_ps_source:
                #YES we have a point source!!
                if prim_sourcesurvey.id=='SAD':
                    ps_gmrtflux.append(matches['match'][0].flux)
                    ps_predictflux.append(
                        fluxpredict(matches['source'].flux,primarydata.ps_flux_orig[prim_ps_source.index(matches['match'][0])],freq2,freq1,gmrt_freq))
                    ps_origflux.append(primarydata.ps_flux_orig[prim_ps_source.index(matches['match'][0])])
                    ps_si.append(getsi(matches['source'].flux,primarydata.ps_flux_orig[prim_ps_source.index(matches['match'][0])],freq2,freq1))
                else:
                    ps_gmrtflux.append(primarydata.ps_flux[prim_ps_source.index(matches['match'][0])])
                    ps_predictflux.append(fluxpredict(matches['source'].flux,matches['match'][0].flux,freq2,freq1,gmrt_freq))
                    ps_origflux.append(matches['match'][0].flux)
                    ps_si.append(getsi(matches['source'].flux,matches['match'][0].flux,freq2,freq1))
                ps_ra.append(matches['source'].position.ra_hms)
                ps_dec.append(matches['source'].position.dec_hms)
                ps_beamdists.append(matches['source'].position.separation(primarydata.saddata.pos)/60.0)
                ps_raoffs.append(matches['source'].position.raoff(primarydata.saddata.pos)/60.0)
                ps_decoffs.append(matches['source'].position.decoff(primarydata.saddata.pos)/60)
                ps_ids.append(primarydata.ids[prim_all_source.index(matches['match'][0])])
            elif len(matches['match'])>0:
                #No not a point source!!
                thisflux=0.0
                this_prim_match_flux=0.0
                threematch=True
                label=''
                for thismatch in matches['match']:
                    if thismatch in prim_all_source:
                        thisflux=thisflux+thismatch.flux
                        this_prim_match_flux=this_prim_match_flux+primarydata.flux[prim_all_source.index(thismatch)]
                        label=label+primarydata.ids[prim_all_source.index(thismatch)]+':'
                        label=label[:-1]
                    else:
                        threematch=False
                if threematch:
                    if prim_sourcesurvey.id=='SAD':
                        es_gmrtflux.append(thisflux)
                        es_predictflux.append(fluxpredict(this_prim_match_flux,matches['source'].flux,freq1,freq2,gmrt_freq))
                        es_origflux.append(this_prim_match_flux)
                        es_si.append(getsi(this_prim_match_flux,matches['source'].flux,freq1,freq2))
                    else:
                        es_gmrtflux.append(this_prim_match_flux)
                        es_predictflux.append(fluxpredict(thisflux,matches['source'].flux,freq1,freq2,gmrt_freq))
                        es_origflux.append(thisflux)
                        es_si.append(getsi(thisflux,matches['source'].flux,freq1,freq2))
                    es_ra.append(matches['source'].position.ra_hms)
                    es_dec.append(matches['source'].position.dec_hms)
                    es_beamdists.append(matches['source'].position.separation(primarydata.saddata.pos)/60.0)
                    es_raoffs.append(matches['source'].position.raoff(primarydata.saddata.pos)/60.0)
                    es_decoffs.append(matches['source'].position.decoff(primarydata.saddata.pos)/60)
                    es_ids.append(label)
    
    #Now we have our arrays- so we can make 2frequency versions of the plots.
    ps_gmrtflux=numpy.array(ps_gmrtflux)
    es_gmrtflux=numpy.array(es_gmrtflux)
    ps_predictflux=numpy.array(ps_predictflux)
    es_predictflux=numpy.array(es_predictflux)

    gmrtflux=numpy.append(ps_gmrtflux,es_gmrtflux)
    predictflux=numpy.append(ps_predictflux,es_predictflux)

    ids=ps_ids+es_ids
    raoffs=ps_raoffs+es_raoffs
    decoffs=ps_decoffs+es_decoffs

    ps_fluxoff=numpy.log10(ps_predictflux)-numpy.log10(ps_gmrtflux)
    es_fluxoff=numpy.log10(es_predictflux)-numpy.log10(es_gmrtflux)
    fluxoff=numpy.log10(predictflux)-numpy.log10(gmrtflux)
    ps_meanoff=scipy.mean(ps_fluxoff)
    meanoff=scipy.mean(fluxoff)

    plotsmade=[]

    if len(gmrtflux>0):
        plotsmade=plotsmade+[plotfluxes(ps_predictflux,ps_gmrtflux,es_predictflux,es_gmrtflux,ps_meanoff,meanoff,sadname,comparesurveys,outputdir+comparesurveys)]
        plotsmade=plotsmade+[plotfluxsnums(gmrtflux,predictflux,ids,sadname,comparesurveys,outputdir+comparesurveys)]
        plotsmade=plotsmade+[plotfluxpbrad(ps_beamdists,ps_fluxoff,es_beamdists,es_fluxoff,outputdir+comparesurveys)]
        plotsmade=plotsmade+[plotfluxpboff(raoffs,decoffs,fluxoff,outputdir+comparesurveys)]
        plotsmade=plotsmade+[makefluxcompfile([ps_ra,es_ra],[ps_dec,es_dec],[ps_beamdists,es_beamdists],[ps_raoffs,es_raoffs],
                                          [ps_decoffs,es_decoffs],[ps_origflux,es_origflux],freq1,
                                          [ps_gmrtflux,es_gmrtflux],gmrt_freq,[ps_predictflux,es_predictflux],[ps_si,es_si],outputdir+comparesurveys)]
        print '**** Found '+comparesurveys+' data. Doing 2 frequency spectral index comparison.'
        print '**** Made '+comparesurveys+' flux comparison plots: '+str(plotsmade)
        primarydata.donetwofreq=True
        primarydata.twofreqnames=comparesurveys
        primarydata.twofreqfluxoffall=meanoff
        primarydata.twofreqps_fluxoff=ps_meanoff
        primarydata.twofreqnum=len(gmrtflux)

    return plotsmade
    


def plotpcoff(posndetcat,imagedata,pcoffplot,posnfilename,telescope='GMRT'):

    from searchcats import *
    import matplotlib
    matplotlib.use('Agg')
    import pylab
    from HRKCatalib import *

    # Get phase center of image 

    rahead=getheader(imagedata,'RA---SIN')
    dechead=getheader(imagedata,'DEC--SIN')
    radms=deg2radec(rahead[0],dechead[0],1)[0:3]
    decdms=deg2radec(rahead[0],dechead[0],1)[3:6]
    
    # Get offsets from phase center
    offsets=[]
    beamdists=[]
    for source in posndetcat:
        if len(source[0])>0:
            minoffset=source[0][0][0]
            for offsetdata in source[0]:
                if (offsetdata[0]<minoffset):
                    minoffset=offsetdata[0]
            offsets.append(minoffset)

            thisra=source[1][1:4]
            thisdec=source[1][4:7]
            beamdists.append(separation(radms,decdms,thisra,thisdec)[0]/60)
            
    # now plot
    pylab.xlim(xmin=0)
    pylab.plot(beamdists,offsets,'ro')
    pylab.xlabel('Offset from phase center (arcmin)')
    pylab.ylabel('Offset between '+telescope+' and '+posnfilename+' (arcsec)')
    pylab.savefig(pcoffplot,format='png')
    pylab.close('all')


def plotxy(firstdetcat,offset,beam,posnplot,posnfilename,survey,telescope='GMRT'):

    import matplotlib
    matplotlib.use('Agg')

    import pylab
    import scipy

    def mkellipse(X,Y,major,minor,PA,dxdy=100):

        from numpy import linspace
        from numpy import pi
        from numpy import cos,sin

        an = linspace(0,2*pi,dxdy)
        
        ex = X + (major * cos(an) * cos(PA*pi/180.) - minor * sin(an) * sin(PA*pi/180.));
        
        ey = Y + (major * cos(an) * sin(PA*pi/180.) + minor * sin(an) * cos(PA*pi/180.));
        return(ex,ey)


    #Construct a list of xoffsets in arcseconds and yoffsets in arcseconds
    raoffs=[]
    decoffs=[]
    sourcecount=0
    for source in firstdetcat:
        if len(source[0])>0:
            minoffset=source[0][0][0]
            minraoff=source[0][0][1]
            mindecoff=source[0][0][2]
	    if survey == 'FIRST':
                minfirstpflux=float(source[2][8])
                minfirstmaj=float(source[2][11])
	    else:
                minfirstpflux=4.0
                minfirstmaj=2.0
            for offsetindex in range(len(source[0])):
                if (source[0][offsetindex][0]<minoffset):
                    minoffset=source[0][offsetindex][0]
                    minraoff=source[0][offsetindex][1]
                    mindecoff=source[0][offsetindex][2]
                    if survey == 'FIRST':
		         minfirstpflux=float(source[2+offsetindex][8])
                         minfirstmaj=float(source[2+offsetindex][11])
            if (minfirstpflux > 3.0 and minfirstmaj < 3.0):
	        raoffs.append(-1.0*minraoff)
                decoffs.append(-1.0*mindecoff)
    #Get means of SD's of raoffs and decoffs.
    ramean=0
    rastd=0
    decmean=0
    decstd=0


    ramean=scipy.mean(raoffs)
    rastd=scipy.std(raoffs)
    decmean=scipy.mean(decoffs)
    decstd=scipy.std(decoffs)

    strra='$<RA>=${0:3.1f}\n$\sigma(RA)=${1:4.2f}'.format(ramean,rastd)
    strdec='$<Dec.>$={0:3.1f}\n$\sigma(Dec.)=${1:4.2f}'.format(decmean,decstd)
    elx,ely=mkellipse(0,0,beam[0]/2,beam[1]/2,beam[2])
    pylab.plot(raoffs,decoffs,'ro')
    pylab.plot(elx,ely,'black')
    pylab.xlim(-1.0*offset,offset)
    pylab.ylim(-1.0*offset,offset)
    
    pylab.axvline()
    pylab.axhline()

    pylab.text(-0.9*offset,0.8*offset,strdec)
    pylab.text(0.5*offset,-0.9*offset,strra)
    ylab='Declination Offset ('+telescope+'-'+posnfilename+') (arcsec)'
    xlab='RA Offset ('+telescope+'-'+posnfilename+') (arcsec)'
    pylab.ylabel(ylab)
    pylab.xlabel(xlab)
    pylab.savefig(posnplot,format='png')
    pylab.close('all')
    return ramean,rastd,decmean,decstd

