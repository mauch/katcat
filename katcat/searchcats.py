import string
import os
import urllib2
import urllib
import gzip
import time
import ObitTalkUtil

from math import *
from astrop import *

FIRSTCATURL='http://sundog.stsci.edu/cgi-bin/searchfirst'
NVSSCATURL='http://www.cv.nrao.edu/cgi-bin/NVSSlist.pl'
SUMSSURL='http://www.astrop.physics.usyd.edu.au/sumsscat/sumsscat.Mar-11-2008'
WENSS='WENSS.Z'
SUMSS='SUMSS.Z'

# Get the WENSS (catid=1) or SUMSS (catid=2) catalogue (should be on the disk) at a given position.
def getincdata(ra,dec,offset,catid):

    def parsesumsscat(line):
        outdata=[]

        outdata.append('SUMSS'+line[0:2].strip()+line[3:5].strip()+line[6:8].strip()+line[13:16].strip()+line[17:19].strip()+line[20:22].strip())
        outdata.append(line[0:2].strip()) #rah
        outdata.append(line[3:5].strip()) #ram
        outdata.append(line[6:11].strip()) #ras
        outdata.append(line[13:16].strip()) #decd
        outdata.append(line[17:19].strip()) #decm
        outdata.append(line[20:24].strip()) #decs
        outdata.append(line[53:62].strip()) #intflux
        outdata.append(line[25:30].strip()) #raerr
        outdata.append(line[31:35].strip()) #decerr
        outdata.append(line[36:45].strip()) #pflux
        outdata.append(line[46:52].strip()) #perr
        outdata.append(line[53:62].strip()) #intflux
        outdata.append(line[63:69].strip()) #interr
        outdata.append(line[70:77].strip()) #fitmaj
        outdata.append(line[78:83].strip()) #fitmin
        outdata.append(line[84:89].strip()) #fitpa
        outdata.append(line[90:97].strip()) #decmaj
        outdata.append(line[98:103].strip()) #decmin
        outdata.append(line[104:109].strip()) #decpa
        outdata.append(line[110:119].strip()) #mosaicname
        outdata.append(line[120:121].strip()) #nummosaics
        outdata.append(line[122:128].strip()) #xpix
        outdata.append(line[129:135].strip()) #ypix

        return outdata

    def parsewensscat(line):
        outdata=[]
        
        outdata.append(line[0:16]) #Source name
        outdata.append(line[42:44].strip()) #rah
        outdata.append(line[45:47].strip()) #ram
        outdata.append(line[48:53].strip()) #ras
        outdata.append(line[54:57].strip()) #decd
        outdata.append(line[58:60].strip()) #decm
        outdata.append(line[61:65].strip()) #decs
        outdata.append(line[79:86].strip()) #Integrated flux
        outdata.append(line[67:68]) #source type
        outdata.append(line[69:70]) #fit flag
        
        outdata.append(line[72:78].strip()) #Peak flux
        outdata.append(line[79:86].strip()) #Integrated flux
        
        outdata.append(line[87:91].strip()) #Major axis length
        outdata.append(line[92:95].strip()) #Minor axis length
        outdata.append(line[96:99].strip()) #PA
        
        outdata.append(line[100:104].strip()) #local rms
        outdata.append(line[105:114]) #Frame name of source

        return outdata

    if catid==1:

        wenssoutput=[]
        # ASSUME FITS DISK IS disk=0 here!
        WENSSCAT = ObitTalkUtil.FITSDir.FITSdisks[0] + '/' + WENSS
        if not os.path.exists(WENSSCAT):
            print 'WENSS catalogue not found!'
            print 'Please ensure that '+WENSSCAT+' is present'
        else:
            wenssfile = gzip.open(WENSSCAT,'r')
            for source in wenssfile:
                rawenss=source[42:54].split()
                decwenss=source[54:66].split()
            
                if separation(rawenss,decwenss,ra,dec)[0]<offset:
                     #We want this source
                     wenssoutput = wenssoutput + [parsewensscat(source)]

                     
        return wenssoutput

    elif catid==2:

        sumssfile=False
        sumssoutput=[]
        SUMSSCAT = ObitTalkUtil.FITSDir.FITSdisks[0] + '/' + SUMSS
        if not os.path.exists(SUMSSCAT):
            print 'SUMSS catalogue not found!'
            print 'I will try and download it.'
            sumssfile=geturl(SUMSSURL,'',30)
            if not sumssfile:
                print 'Cannot reach online SUMSS catalogue either.'
                print 'No SUMSS data available.'
        else:
            sumssfile = gzip.open(SUMSSCAT,'r')

        if sumssfile:
            for source in sumssfile:
                rasumss=source[0:11].split()
                decsumss=source[13:24].split()

                if separation(rasumss,decsumss,ra,dec)[0]<offset:
                    sumssoutput = sumssoutput + [parsesumsscat(source)]
        return sumssoutput
            
# Attempt to obtain open a url 'tries' times with timeout 'timeout' with a 5 second wait between each try.
def geturl(url,params,timeout,tries=5):

    #print 'Attempting to open URL:'+url

    for i in range(1,tries+1):
        #print 'Attempt '+str(i)+' of '+str(tries)
        try:
            urlsocket=urllib2.urlopen(url,params,timeout)
            urldata=urlsocket.read()
        except:
            #print "FAIL"
            time.sleep(5)
        else:
            break

    if i==tries:
        print "Cannot access "+url
        return False
    else:
        #print 'SUCCESS'
        return urldata

        
    
# At a given ra,dec and required search radius, search the FIRST catalogue and
# produce a list of all of the FIRST data at given position. An empty array is returned
# if there is nothing. Else the reutrned array contains [[source data 1],[source data 2],...]
def getfirstdata(ra,dec,offset,posnfilename):

    firsturl=False
    firstdata=[]
    if os.path.exists(posnfilename):
        #If file is available we will use that.
        firsturl=gzip.open(posnfilename,'r')
        try:
            firsturl.readline()
        except:
            firsturl=open(posnfilename,'r')
        firsturl.seek(0)
        for source in firsturl:
            rafiledd=source.split()[1]
            decfiledd=source.split()[2]
            radecsex=deg2radec(float(rafiledd),float(decfiledd))

            rafile=radecsex[0:3]
            decfile=radecsex[3:6]
            if separation(ra,dec,rafile,decfile)[0]<offset:
                outdata=source.split()[0:1] + list(rafile) + list(decfile) + source.split()[3:]
                firstdata=firstdata+[outdata]
    else:
        #Otherwise we'll use the FIRST website.
        firstparams = urllib.urlencode({'RA': ' '.join(ra+dec),'Radius': offset,'Text':1, 'Equinox': 'J2000'})
        firsturl = geturl(FIRSTCATURL,firstparams,10)

        if not firsturl:
            print "**** No FIRST data available."
            return firstdata
           
        urldata=firsturl.split('\n')
        numsources=0
        for tempstr in urldata:
            if len(tempstr)>0 and tempstr.split()[len(tempstr.split())-1]=='arcsec':
                if tempstr.split()[1]=='No':
                    print "**** No FIRST data available."
                    return firstdata
                else:
                    numsources=float(tempstr.split()[1])
                    break

#FIRST cat only returns 500 sources at a time so we need to iterate every 500 sources.
        numloops=int(ceil(numsources/500))

        for x in range(numloops):
            firstparams = urllib.urlencode({'RA': ' '.join(ra+dec),'Radius': offset,'Text':1, 'Equinox': 'J2000', 'PStart':x*500, '.cgifields': 'Text'})
            firsturl = geturl(FIRSTCATURL,firstparams,10)
            offsetcounter=0
            urldata=firsturl.split('\n')

            for tempstr in urldata:
                if len(tempstr)>0 and tempstr[0]=='#':
                    offsetcounter+=1
            del(urldata[0:offsetcounter])

            for lines in urldata:
                if len(lines)>6:
                    firstdata=firstdata+[lines.split()]

    return firstdata

#Function to take a catalogue and match each source from cat1 in cat2.
#Assumed FIRST input is from the output of getfirstdata function
def getmatches(cat1,cat2,offset):

    outcat=[]

    for source in cat1:
        ra=source[1:4]
        dec=source[4:7]
        offsetarray=[]
        detsarray=[]
        for source2 in cat2:
            ra2=source2[1:4]
            dec2=source2[4:7]

            thisoffset=separation(ra,dec,ra2,dec2)

            if (thisoffset[0]<offset):
                offsetarray=offsetarray + [thisoffset]
                detsarray=detsarray+[source2]

        outcat=outcat+[[offsetarray]+[source]+detsarray]
    
    return outcat

# Does the same as the above 'getfirstdata' function but for the NVSS catalogue instead.
def getnvssdata(ra,dec,offset):

    def parsenvssline1(line):

        outdata=[]
        outdata.append(line[0:2].strip()) #rah
        outdata.append(line[3:5].strip()) #ram
        outdata.append(line[6:11].strip()) #ras
        outdata.append(line[12:15].strip()) #decd
        outdata.append(line[16:18].strip()) #decm
        outdata.append(line[19:23].strip()) #decs
        outdata.append(line[30:37].strip()) #flux need to ensure its position 7 in output
        outdata.append(line[24:30].strip()) #offset

        outdata.append(line[38:43].strip()) #maj
        outdata.append(line[44:49].strip()) #min
        outdata.append(line[50:55].strip()) #pa

        outdata.append(line[56:58].strip()) #residual

        outdata.append(line[59:65].strip()) #pflux
        outdata.append(line[66:71].strip()) #pangle

        outdata.append(line[72:80].strip()) #mosaic

        outdata.append(line[81:88].strip()) #xpix
        outdata.append(line[89:96].strip()) #ypix

        return(outdata)
        

    def parsenvssline2(line):

	outdata=[]

        outdata.append(line[0:11].strip()) #raerror
        outdata.append(line[12:23].strip()) #decerror

        outdata.append(line[24:29].strip()) #offseterr

        outdata.append(line[29:36].strip()) #fluxerr

        outdata.append(line[36:42].strip()) #majerr
        outdata.append(line[42:48].strip()) #minerr
        outdata.append(line[48:54].strip()) #paerr

        outdata.append(line[54:58].strip()) #residualerr

        outdata.append(line[58:64].strip()) #pfluxerr
        if (len(line)>65):
            outdata.append(line[64:70].strip()) #pangleerr
        else:
            outdata.append('')

        return outdata


    def parsenvssfile(nvssurl,nvssfiturl):

        nvssdata=[]
        tempdata=iter(nvssurl)
        count=0;
        templine=tempdata.next()
        while True:
            if templine and templine[0].isdigit():
                line=parsenvssline1(templine)
                tempdata.next()
                fmaj,fmin,fpa=-1.,-1.,-1.

                for tempfitline in nvssfiturl:
                    if tempfitline and tempfitline[0].isdigit():
                        fitline=parsenvssline1(tempfitline)
                        if separation(fitline[0:3],fitline[3:6],line[0:3],line[3:6])[0]<0.01:
                            fmaj,fmin,fpa=fitline[8],fitline[9],fitline[10]
                if fmaj>-1.0:
                    line=line+[fmaj,fmin,fpa]
                    thisoffset=separation(ra,dec,line[0:3],line[3:6])[0]
                    line=[str(thisoffset)]+line
                    if thisoffset<offset:
                        nvssdata=nvssdata+[line]
            try:
                templine=tempdata.next()
            except:
                break
            
        return nvssdata


    if os.path.exists('./NVSS.z'):
        #Read a file.
        nvssurl = gzip.open('./NVSS.z','r')
        nvssurldata=nvssurl.readlines()
        nvssdata=parsenvssfile(nvssurldata)
    else:
        #Get it from the web!
        nvssparams = urllib.urlencode({'Equinox': 3,'DecFit': 0,'FluxDensity': 0,'PolFluxDensity': 0,'RA': ' '.join(ra),'Dec': ' '.join(dec) ,'searchrad':offset})
        nvssurl = geturl(NVSSCATURL,nvssparams,10)
        if not nvssurl:
            print '**** No NVSS data available'
            return []
        else:
            nvssparams = urllib.urlencode({'Equinox': 3,'DecFit': 1,'FluxDensity': 0,'PolFluxDensity': 0,'RA': ' '.join(ra),'Dec': ' '.join(dec) ,'searchrad':offset})
            nvssfiturl = geturl(NVSSCATURL,nvssparams,10)
            if not nvssfiturl:
                return []
            else:
                nvssfiturldata = nvssfiturl.split('\n')
                nvssurldata = nvssurl.split('\n')
                nvssdata=parsenvssfile(nvssurldata,nvssfiturldata)
    return nvssdata
        
def getowndata(ra,dec,offset,posnfilename):

    outdata=[]

    if os.path.exists(posnfilename):
        print 'Using catalogue in '+posnfilename
        posnfile=gzip.open(posnfilename)
        try:
            posnfile.readline()
        except:
            posnfile=open(posnfilename,'r')
        posnfile.seek(0)
        for source in posnfile:
            rafiledd=source.split()[1]
            decfiledd=source.split()[2]
            radecsex=deg2radec(float(rafiledd),float(decfiledd))
            rafile=radecsex[0:3]
            decfile=radecsex[3:6]
            if separation(ra,dec,rafile,decfile)[0]<offset:
                outdata=outdata+[source.split()[0:1]+list(rafile)+list(decfile)+source.split()[3:]]
    else:
        print 'FILE: '+posnfile+' not found'

    return(outdata)

#If Survey=1:
#Go through the entries in GMRT cat and find all of the FIRST associations within a certain offset.
#If Survey=2:
#Go through the entries in the NVSS cat covering the region of the GMRT image and find the GMRT 
#associations within each source.
#If survey=3 do WENSS (requires WENSS.z to be present.
#If survey=4 do SUMSS.
#This will return a list of GMRT matched with a nested list of first matches after each one
#in the format [[[offsetarray],[GMRT data 1],[FIRST data 1],[FIRST data 2],...],[[offsetarray],[GMRT data 2], [FIRST data 1], ...],...]
#survey=1 =>FIRST, survey=2 =>NVSS, survey=3=>WENSS, survey=4=>SUMSS, survey=5=>OWN CATALOGUE
def findmatches(gmrtcat,gmrtimage,offset,survey,posnfilename):
    
    rahead=getheader(gmrtimage,'RA---SIN')
    dechead=getheader(gmrtimage,'DEC--SIN')
    radms=deg2radec(rahead[0],dechead[0],1)[0:3]
    decdms=deg2radec(rahead[0],dechead[0],1)[3:6]
    sizedata=imrange(gmrtimage)
    rasize=fabs(sizedata[0]-sizedata[1])
    decsize=fabs(sizedata[2]-sizedata[3])
    imsize=ceil(max(rasize,decsize))*3600/2.
    if survey==1:
        firstcatdata=getfirstdata(radms,decdms,imsize,posnfilename)
        if firstcatdata:
            matchdata=getmatches(gmrtcat,firstcatdata,offset)
        else:
            matchdata=False
    elif survey==2:
        nvsscatdata=getnvssdata(radms,decdms,imsize)
        if nvsscatdata:
            matchdata=getmatches(nvsscatdata,gmrtcat,offset)	
        else:
            matchdata=False
    elif survey==3:
        wensscatdata=getincdata(radms,decdms,imsize,1)
        if wensscatdata:
            matchdata=getmatches(wensscatdata,gmrtcat,offset)
        else:
            matchdata=False
    elif survey==4:
        sumsscatdata=getincdata(radms,decdms,imsize,2)
        if sumsscatdata:
            matchdata=getmatches(sumsscatdata,gmrtcat,offset)
        else:
            matchdata=False
    elif survey==5:
        # User provided catalogue in 'posnfilename'
        # Read the data and check that it lies within the area of the GMRT image. 
        data=getowndata(radms,decdms,imsize,posnfilename)
        if data:
            matchdata=getmatches(data,gmrtcat,offset)
        else:
            matchdata=False

    return matchdata

