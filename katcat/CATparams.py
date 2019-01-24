#################
#
# GMRT Cataloguer
#
# Parameters to be chaged by user
#
#################

#################

# Aips Usernumber

aipsno = 100


#################

# Cutoff in sigma for the catalogue (5sigma is normal)

cutoff = 5

#################

# Flux density cutoff (in milliJy) for the catalogue (this is applied as well as the sigma cut above)

fluxcutoff = 0.0


#################

# Is the input image primary beam corrected? The code will attempt to remove any primary
# beam correction from the image if it has been applied to create a bright source model.
# IF pbcorrect=1 (DEFAULT) then input image is assumed to be corrected for the primary beam
# and this correction will be removed to create the bright source model.
# IF pbcorrect=0 then input image is assumed to NOT be corrected for the primary beam
# and no primary beam correction will be removed from the image.
# If your image is not primary beam corrected or you are using a mosaic image then
# you must set this parameter to 0!!

pbcorrect = 0


#################

# Should the output catalogue be primary beam corrected or not? The default is
# to produce a catalogue that includes the primary beam correction. Setting output_pbcorrect
# to 0 will remove the primary beam correction in the output catalogue if pbcorrect=1 . If pbcorrect
# is 0 then this parameter will have no effect.
#

output_pbcorrect = 1

#################

# Cutoff in sigma for the detection of bright sources in the image
# around which the local rms will be increased artifically to remove
# artefacts. 100 sigma seems to work well here. Decrease this if you
# are seeing artefacts around fainter sources.

highcut = 999


#################

# Set the beam paramaters in the image manually with bmaj(arcsec),bmin(arcsec),bpa(deg.).
# Use this if your image doesn't have beam parameters in its header (they will be added by the code)
# BE VERY CAREFUL- as this will override any (probably correct) beam parameters that are in the
# FITS header!!!!

user_bmaj = 0.0
user_bmin = 0.0
user_bpa = 0


################

# The following parameters define the output files that will be returned by the code.
# The default is to make an annotation file, a formatted catalogue and a vo table
# The user can switch these off- or add the raw data, or the raw sad output as well.

keepsad = 0
keepraw = 0
keepcat = 1
keepann = 1
keepvo  = 1


################

# Set up a more streamlined way of orgaising output.
# First for each survey (FIRST,NVSS,WENSS,SUMSS- define a 4 variable array with entries:
# DO Position?, Do Flux?, Do kvis Annotation?
# Extra variable for 

FIRSTcompare = [1,1,1]
NVSScompare  = [0,0,0]
WENSScompare = [1,1,1]
SUMSScompare = [1,1,1]

#OWNCATcompare = [['MYCAT','MYCATFILE',1400,5,[1,0,0,0]],['MYCAT2','MYCAT2FILE',....],....]


#################

# These are advanced parameters that you shouldnt change if you don;t know what they are.
# Some set the inputs to the AIPS task 'sad'.
#

dorms = 1  #If =1 then output the rms model image.

rmsbox = 100 #The size of the box (in beams) used to comput the rms noise per pixel.

gain = 0.05  #For sad

ifact = 0.0 #(ifact*cutoff)=icut (for sad)

maxfitwidth = 5 #Largest source to fit in beams.

mindr = 50 #The minimum acceptable local dynamic range (sources above this next to bright sources will end up in the catalogue)
 
localwidth = 0.05 #Fraction of the image size to accept as a source 'local' to a bright source.

noisefact = 1.2 #Factor by which to increase the model noise derived from the local DR. This is to catch stray noise peaks close to bright sources.

edgepix = 20 #Number of pixels next to the edge of the image to reject fits (spurious fits can occur close to the edge of the images)

localsource = 0.5 #Number of beams around a source to accept a match as a 'connected' source.

posnnumsigma = 20
