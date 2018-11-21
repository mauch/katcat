from math import *
#from HRKCatalib import *

def cosd(angle):
    return(cos(radians(angle)))

def sind(angle):
    return(sin(radians(angle)))

def separation(ra1s,dec1s,ra2s,dec2s):

    ra1,dec1 = hmsra2ra(float(ra1s[0]),float(ra1s[1]),float(ra1s[2])),dmsdec2dec(dec1s[0][0],float(dec1s[0]),float(dec1s[1]),float(dec1s[2]))
    ra2,dec2 = hmsra2ra(float(ra2s[0]),float(ra2s[1]),float(ra2s[2])),dmsdec2dec(dec2s[0][0],float(dec2s[0]),float(dec2s[1]),float(dec2s[2]))

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

    return [separation,raoff,decoff]


#####
#Return x,y coordinates of an ellipse
#with axes 'major' and 'minor' at
#position angle 'PA' with center X,Y
def mkellipse(X,Y,major,minor,PA,dxdy=100):

    from numpy import linspace
    from numpy import pi
    from numpy import cos,sin
    
    an = linspace(0,2*pi,dxdy)
    
    ex = X + (major * cos(an) * cos(PA*pi/180.) - minor * sin(an) * sin(PA*pi/180.));
    
    ey = Y + (major * cos(an) * sin(PA*pi/180.) + minor * sin(an) * cos(PA*pi/180.));
    return(ex,ey)

######
#Deconvolve elliptical gaussian g1[maj1,min1,pa1]
#from elliptical gaussian g2[maj2,min2,pa2]
#Equations taken from Wild, 1970, AuJPh
def deconv(g1,g2):

    if g1[0]==g1[1]:
        maj1corr=g1[0]+0.0001
    else:
        maj1corr=g1[0]

    if g2[0]==g2[1]:
        maj2corr=g2[0]+0.0001
    else:
        maj2corr=g2[0]

    g1sd=maj1corr*maj1corr-g1[1]*g1[1]
    g2sd=maj2corr*maj2corr-g2[1]*g2[1]
    g1ss=maj1corr*maj1corr+g1[1]*g1[1]
    g2ss=maj2corr*maj2corr+g2[1]*g2[1]
    
    pad=g1[2]-g2[2]

    betasq=g1sd*g1sd + g2sd*g2sd - 2.0*g1sd*g2sd*cosd(2.0*pad)
    beta=sqrt(betasq)

    #Can we deconvolve??
    if (g1ss-g2ss+beta)<0:
        decmaj=-1.0
    else:
        decmaj=sqrt(g1ss-g2ss+beta)/2.0

    if (g1ss-g2ss-beta)<0:
        decmin=-1.0
    else:
        decmin=sqrt(g1ss-g2ss-beta)/2.0

    if (decmaj>0.0) or (decmin>0.0):
        tan2decpa=(g1sd*sind(2.0*g1[2]) - g2sd*sind(2.0*g2[2]))/(g1sd*cosd(2.0*g1[2]) - g2sd*cosd(2.0*g2[2]))
        decpa=degrees(atan(tan2decpa))/2.0
    else:
        decpa=0.0
    
    return(decmaj,decmin,decpa)

######
#Convolve an elliptical gaussian g1[maj1,min1,pa1]
#with elliptial gaussian g2[maj2,min2,pa2]
#Equations taken from Wild, 1970, AuJPh
#
def elconv(g1,g2):

    g1sd=g1[0]*g1[0]-g1[1]*g1[1]
    g2sd=g2[0]*g2[0]-g2[1]*g2[1]
    g1ss=g1[0]*g1[0]+g1[1]*g1[1]
    g2ss=g2[0]*g2[0]+g2[1]*g2[1]
    
    pad=g1[2]-g2[2]

    betasq=g1sd*g1sd + g2sd*g2sd + 2.0*g1sd*g2sd*cosd(2.0*pad)
    beta=sqrt(betasq)

    cmaj=sqrt(g1ss+g2ss+beta)/2.0
    cmin=sqrt(g1ss+g2ss-beta)/2.0

    cos2cpa=(g1sd*cosd(2.0*g1[2]) + g2sd*cosd(2.0*g2[2]))/beta
    # Fix rounding error
    if cos2cpa>=1.0:
        cos2cpa = 0.999999999
    cpa=degrees(acos(cos2cpa))/2.0

    return(cmaj,cmin,cpa)

def hmsra2ra(h,m,s):
  return (h + m/60.0 + s/3600.0) * 15.0

def ra2racos(ra, dec):
  # Convert RA to RA/sin(decradians)
  return ra*cos(dec/180*pi)

def dmsdec2dec(ds,d,m,s):
  if ds != "-":
    return d + m/60.0 + s/3600.0
  else:
    return d - m/60.0 - s/3600.0

def stringlist2ra(stringlist):
  # Convert RA string list to decimal degrees
  return hmsra2ra(float(stringlist[0]), float(stringlist[1]), float(stringlist[2]))

def stringlist2dec(stringlist):
  # Convert RA string list to decimal degrees
  return dmsdec2dec(float(stringlist[0]), float(stringlist[1]), float(stringlist[2]))

def deltara(ra1, ra2):
  # Wrap RA to give shortest distance
  dRA = ra1 - ra2
  if dRA < 0.0:
    dRA = dRA * -1
  if dRA > 180.0:
    dRA = 360.0 - dRA
  return dRA

def deg2sex(ra,dec):
  # Convert decimal degrees into hms,dms

  arcsec = ra/15.0
  ra1    = int(arcsec)
  ra2    = (arcsec - int(arcsec))*60.0
  ra3    = int(ra2)
  ra4    = (ra2 - int(ra2))*60.0

  decs=''
  dec1   = int(dec)
  dec2   = (dec - int(dec))*60.0
  dec3   = int(dec2)
  dec4   = (dec2 - int(dec2))*60.0
  if (dec<0):
    decs='-'
    dec1=int(fabs(dec1))
    dec3=int(fabs(dec2))
    dec4=float(fabs(dec4))

  return '{0:02d}'.format(ra1),'{0:02d}'.format(ra3),'{0:05.2f}'.format(ra4),decs+'{0:02d}'.format(dec1),'{0:02d}'.format(dec3),'{0:05.2f}'.format(dec4)

#given two fluxes at freq1 and freq2 compute alpha
#(the power-law spectral index between them)
def getsi(flux1,flux2,freq1,freq2):

  return(log(flux1/flux2)/log(freq1/freq2))

#Given two fluxes (flux 1 and 2) at (freq1 and freq2)
#predict the flux at freq3 assuming a power-law
def fluxpredict(flux1,flux2,freq1,freq2,freq3):

  alpha=log(flux1/flux2)/log(freq1/freq2)

  flux3=alpha*(log(freq3/freq1))+log(flux1)

  return(exp(flux3))


#Find the offset between two picel positions using pythogras
def pixoffset(pixx1,pixy1,pixx2,pixy2):

    xoff=pixx1-pixx2
    yoff=pixy1-pixy2
    offset=sqrt(xoff**2.0 + yoff**2.0)

    return(offset)
