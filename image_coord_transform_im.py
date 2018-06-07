import numpy as np
from astropy.io import ascii
from astropy.io import fits
import catalog_match_im_trans as cmti
from pyraf import iraf
import os

'''
Run with import image_coord_transform_im as pcti

#this matches the input and reference catalog and runs geomap and
#geoxytran on these catalogs
pcti.cat_match_im(septol, fullcat_trans = 1, im_trans = 1)

'''

def cat_match_im(septol, **kwargs):

    '''Written by Gregory Rudnick 7 June 2018

    PURPOSE:

    Read in an input image, an input catalog with x-y coordinates, a
    reference catalog with x-y coordinates, and a list of x-y
    coordinates for three matched stars.  Then match the catalogs
    within some tolerance using an initial transformation derived from
    the three matched stars.  Use the transformation to transform the
    input image and the input catalog.

    Plot the differences in each coordinate.

    This version differs from the original GOGREEN file in that I've
    taken out all of the references to GOGREEN paths so it should be generic.

    INPUT PARAMETERS:

    septol: the maximum separation allowed for a match in pixels

    NEEDED INPUT FILES (enter in cat_read() subroutine) :

    input catalog with x and y coordinates

    reference catalog

    file with coordinates of three objects in the input and reference
    image.  This contains the x-y coordinates of 1-3 reference list
    tie points in the first line, followed by the x-y coordinates of
    the 1-3 corresponding input tie point in succeeding lines.



    OPTIONAL KEYWORD PARAMETERS

    im_trans: if set==1 then transform image using geotran

    '''

    #read in the z-band GOGREEN catalog and the reference catalog.
    #@@change@@ need to put in whatever cluster name you use.
    clustname = 'CL1216-1201'
    (ref_dat, gg_catname, incat_dat, incatname, initcoordfile) \
        = cat_read()

    #rename inputs to make code more readable
    xref = np.array(ref_dat['X_IMAGE'])
    yref = np.array(ref_dat['Y_IMAGE'])
    xin = np.array(incat_dat['X_IMAGE'])
    yin = np.array(incat_dat['Y_IMAGE'])

    #match catalogs against each other
    #return matched values
    #catpath = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SPT2106'
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_im.' + clustname + '.in'
    (xrefm,yrefm,xinm,yinm, lims) \
        = cmti.cat_im_match(xref, yref, xin, yin, septol, \
                            icfile = initcoordfile, matchfile = geomap_infile)

    print("limits of coordinates are",lims)

    #make a plot of the residuals
    pretrans_plotfile = catpath + '/pretrans.' + clustname + '_coordiff_im.pdf'
    cmti.match_diff_im_plot(xrefm,yrefm,xinm,yinm, plotfile = pretrans_plotfile)

    #run geomap to compute the transformation and geoxytran to
    #transform the input coordinates using that solution.  This also
    #plots the transformed coordinates.
    (dbfile, geomap_infile) =  georun_im(clustname, lims)
                
    #transform image using geotran
    if 'im_trans' in kwargs.keys():
        if kwargs['im_trans'] == 1:
            (impath, refimpath) = preimage_read()
            imtrans(impath, refimpath, dbfile,geomap_infile,lims)
            
def imtrans(impath, refimpath, dbfile, geomap_infile, lims):

    '''Run geotran on an image

    INPUT:

    impath : the full path to the image to be transformed

    refimpath : the full path to the reference image

    dbfile: the geomap database file
    
    geomap_infile: the name of the input coordinate file for the
    geomap transformation, which is used as the record identifier for
    the transfor

    lims: a dictionary that defines the minimum and maximum 

    '''

    #find the x and y limits of the reference image
    refim = fits.open(refimpath)
    xmax = refim[0].header['NAXIS1']
    ymax = refim[0].header['NAXIS2']
    xmin = 1
    ymin = 1
    refim.close()

    outpath = impath.replace('.fits[1]','.trans.fits')

    iraf.geotran(impath, outpath, dbfile, geomap_infile, xmin = xmin, xmax = xmax, \
                 ymin = ymin, ymax = ymax, xscale = 1.0, yscale = 1.0)
            
def georun_im(clustname, lims):

    #run geomap and geoxytran.  I have to run this within this module
    #as geomap crashes when run from an external module.

    #in this version 
    
    catpath = '.'
    geomap_infile = catpath + '/geomap_coords_im.' + clustname + '.in'
    dbfile = catpath + '/' + clustname + '_geomap_im.db'
    iraf.geomap.xmin=lims['xmin']
    iraf.geomap.xmax=lims['xmax']
    iraf.geomap.ymin=lims['ymin']
    iraf.geomap.ymax=lims['ymax']
    iraf.geomap.maxiter=3
    iraf.geomap.reject=3.0
    iraf.geomap.fitgeometry='rxyscale'
       
    iraf.geomap(geomap_infile, dbfile)
    
    geotran_outfile = catpath + '/geomap_coords_im.' + clustname + '.out'

    #remove geotran output file if it already exists
    if os.path.isfile(geotran_outfile) is True:
        cmdstr = 'rm ' + geotran_outfile
        os.system(cmdstr)

    #transform the input coordinates in the geomap input file.  
    iraf.geoxytran(geomap_infile, geotran_outfile, dbfile, geomap_infile, direction="backward",xcolumn=3,ycolumn=4)
    
    #read in the geomap input catalog that was transformed using geoxytran.
    trans_ggdat = ascii.read(geotran_outfile)
    xtrans = trans_ggdat['xin']
    ytrans = trans_ggdat['yin']
    xref = trans_ggdat['xref']
    yref = trans_ggdat['yref']

    #make a plot of the new residuals
    posttrans_plotfile = catpath + '/posttrans.' + clustname + '_coordiff_im.pdf'
    cmti.match_diff_im_plot(xref,yref,xtrans,ytrans, plotfile = posttrans_plotfile)

    return dbfile, geomap_infile


def cat_trans_im(incat, dbfile, geomap_infile, refcat, clustname, septol,  **kwargs):

    '''This routine reads in a FITS catalog, writes a temporary output
    ASCII catalog, then transforms this catalog using geoxytran.  It
    then appends these new coordinates as new columns to the original catalog.

    It is usually run by hand after the geomap solution has been derived.

    INPUT

    incat: the input SExtractor catalog

    dbfile: the geomap database file

    geomap_infile: used as the database record

    refcat: the original astrometric reference catalog

    clustname: the name of the cluster

    septol: the maximum separation allowed for a match in pixels

    OPTIONAL KEYWORDS

    xmin, xmax, ymin, ymax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the ra and decdiff plots.
    If one is given, all must be given.

    '''

    #read in GMOS sextractor photometry catalog with (x,y) coordinates
    cat_dat = ascii.read(incat)

    #read in original astrometric catalog
    ref_dat = ascii.read(refcat)
    xref = np.array(ref_dat['X_IMAGE'])
    yref = np.array(ref_dat['Y_IMAGE'])

    #tmp coordinate files for geoxytran
    tmpin = 'tmp_geoxytran_im_in'
    tmpout = 'tmp_geoxytran_im_out'

    #make a new name for the transformed file
    newcat = incat.replace('.sexcat','.trans.cat')

    if os.path.isfile(tmpout) is True:
        cmdstr = 'rm ' + tmpout
        os.system(cmdstr)
        cmdstr = 'rm ' + newcat
        os.system(cmdstr)

    #output temporary ASCII file with coordinates
    fo = open(tmpin, "w")
    fo.write("# x y\n")
    for i,val in enumerate(cat_dat['X_IMAGE']):
        fo.write('{} {}\n'.format(cat_dat['X_IMAGE'][i],cat_dat['Y_IMAGE'][i]))
    fo.close()

    iraf.geoxytran(tmpin, tmpout, dbfile, geomap_infile, direction="backward",\
                   xcolumn=1, ycolumn = 2)

    #read in ascii output file
    trans_dat = ascii.read(tmpout)
    xtrans = np.array(trans_dat['x'])
    ytrans = np.array(trans_dat['y'])

    #replace the old RAs and DECs with new RAs and DECs in place    
    tcat_dat = cat_dat
    tcat_dat['X_IMAGE'] = xtrans
    tcat_dat['Y_IMAGE'] = ytrans

    #write the new catalog file with the transformed coordinates
    ascii.write(tcat_dat, newcat, format='commented_header')
    
    #match new catalog against original catalog
    mfile = "allcat_im_match.txt"
    (xrefm,yrefm,xtransm,ytransm, translims) = cmti.cat_im_match(xref, yref, \
                                                                 xtrans, ytrans, septol, \
                                                                 matchfile = mfile)

    #make a plot of the residuals
    allcattrans_plotfile = 'allcat_trans.' + clustname + '_coordiff_im.pdf'
    #passes ra and dec limits if they are defined to find source
    #outside of ra and dec lims.  Assumes that if one keyword is givenSPT2106_GMOS_z.v0.sexcat'
    #that all are given
    
    if 'xmin' in kwargs.keys():
        cmti.match_diff_im_plot(xrefm,yrefm,xtransm,ytransm, plotfile = allcattrans_plotfile, \
                            xmin = kwargs['xmin'], xmax = kwargs['xmax'], \
                            ymin = kwargs['ymin'], ymax = kwargs['ymax'])
    else:
        cmti.match_diff_im_plot(xrefm,yrefm,xtransm,ytransm, plotfile = allcattrans_plotfile)

    
def cat_read():

    #read in catalogs

    #read in input catalog with x and y coordinates
    #@@change to location of your input catalog@@
    incat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SPT2106/SPT2106_GMOS_z.v0.sexcat'
    incat_dat = ascii.read(incat)
    
    #reference catalog
    #@@change to location of the reference catalog
    refcat = '/Users/grudnick/Work/GOGREEN/Catalogs/Astrometric/SPT2106/SPT2106_J.v0.sexcat'

    ref_dat = ascii.read(refcat)

    #file with coordinates of three objects in the input and reference
    #image
    #@@change to name of file with three coordinates.@@
    #This contains the x-y coordinates of 1-3 reference list tie
    #points in the first line, followed by the x-y coordinates of the
    #1-3 corresponding input tie point in succeeding lines.
    initcoordfile = 'SPT2106_initref.coord.txt'


    return gg_dat, ref_dat, refcat, incat_dat, incat, initcoordfile;

def preimage_read():

    '''
    return the paths to the reference and input image

    INPUT PARAMETERS:

    OUTPUT

    The paths to both images
    
    '''

    #full path to input image
    #@@change to input image path@@
    impath = '/Users/grudnick/Work/GOGREEN/Data/Preimaging/SPT2106/GMOS/Z/mfrgS20150411S0405_add.fits'

    #full path to reference image
    #@@change to reference image path@@
    refimpath = '/Users/grudnick/Work/GOGREEN/Data/Imaging/Fourstar/Reduced/Sep2016/SPT2106_J_v02_ipe.fits'
    
    return impath,refimpath
