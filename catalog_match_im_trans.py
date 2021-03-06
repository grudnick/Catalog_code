import numpy as np
from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from pyraf import iraf 
import os
   

def cat_im_match(xref, yref, xin, yin, septol, **kwargs):

    '''Written by Gregory Rudnick 9 January 2018

    PURPOSE:

    Take two lists of x and y coordinates in the same units and match
    them within some tolerance.

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    xref, yref: the reference coordinates in pixels.  Numpy arrays.

    xin, yin: the input coordinates in pixels.  If performing
    coordinate transforms, these would be the ones to be transformed.
    Numpy arrays

    septol: the maximum separation allowed for a match.  Units are
    pixels

    OPTIONAL KEYWORD PARAMETERS

    icfile: The filename with the input reference coordinates
    for xyxymatch.  These need to be created by hand from the images

    matchfile: a name of the file that will contain the reference and
    input coordinates.  Suitable for geomap input.  Default is 'xyxy_out.txt'

    OUTPUT

    arrays containing row-matched (x,y) coordinates of the reference and input
    coordinates that fall with septol

    '''

    # #loop through all coordiantes and find closest match.  This is an
    # #N^2 process.  I can make it faster later.
    # #initialize distance coordinate
    # dist = np.full(xref.size, 1.e6)
    # #initialize array with indices of closest match
    # imatch = np.zeros(xref.size, dtype=np.int8)
    # for iref,valref in enumerate(xref):
    #     for iin,valin in enumerate(xin):
    #         #print("hi",iref,iin,xref[iref],xin[iin],yref[iref],yin[iin])
    #         disttest = np.sqrt( (xref[iref] - xin[iin])**2 +(yref[iref] - yin[iin])**2)

    #         #find the closest reference point to each input point
    #         if disttest < dist[iref]:
    #             #print(disttest,iref,iin)
    #             dist[iref] = disttest
    #             imatch[iref] = iin

            
    # #select close matches, where I assume septol is in pixels, just like the input catalog
    # iclose = np.where(dist < septol)
    # #convert tuple to pure array
    # iclose=iclose[0]

    #generate input files for xyxymatch from reference coordinates
    xyxy_refin = 'xyxymatch_refcoords.txt'
    fo = open(xyxy_refin, "w")
    for i,val in enumerate(xref):
        #print(i,iclose[i])
        #print(raref[iclose[i]],decref[iclose[i]],rain[idx[iclose[i]]],decin[idx[iclose[i]]])
        fo.write('{} {}\n'.format(xref[i],yref[i]))
    fo.close()

    xyxy_inin = 'xyxymatch_incoords.txt'
    fo = open(xyxy_inin, "w")
    for i,val in enumerate(xin):
        fo.write('{} {}\n'.format(xin[i],yin[i]))
    fo.close()

    #set xyxymatch output name if keyword is given
    xyxy_out = 'xyxy_match.txt'
    #remove geotran output file if it already exists
    if os.path.isfile(xyxy_out) is True:
        cmdstr = 'rm ' + xyxy_out
        os.system(cmdstr)

        
    #run xyxymatch and include a set of reference points if given as input
    if 'icfile' in kwargs.keys():
        iraf.xyxymatch(xyxy_inin, xyxy_refin, xyxy_out, septol, refpoints=kwargs['icfile'], \
                       xcolumn=1,ycolumn=2,xrcolumn=1,yrcolumn=2,matching="tolerance")
    else: 
        iraf.xyxymatch(xyxy_inin, xyxy_refin, xyxy_out, septol, refpoints="", \
                       xcolumn=1,ycolumn=2,xrcolumn=1,yrcolumn=2,matching="tolerance")

    #read in xyxymatch output to get coordinate limits and to return
    #ordered matched coordinates
    xyxy_cat = ascii.read(xyxy_out)
    xref_m = np.array(xyxy_cat['col1'])
    yref_m = np.array(xyxy_cat['col2'])
    xin_m = np.array(xyxy_cat['col3'])
    yin_m = np.array(xyxy_cat['col4'])
    
    #write files only if "matchfile" keyword is set
    #this reformat of the xyxymatch code is just so that my other routines work
    keys = sorted(kwargs.keys())
    for kw in keys:
        if kw == 'matchfile':
            #print(kwargs[kw])
            #open file for writing and write a header
            fo = open(kwargs[kw], "w")
            fo.write("# xref yref xin yin\n")
            for i,val in enumerate(xref_m):
                #print(i,iclose[i])
                #print(raref[iclose[i]],decref[iclose[i]],rain[idx[iclose[i]]],decin[idx[iclose[i]]])
                fo.write('{} {} {} {}\n'.format(xref_m[i],yref_m[i],xin_m[i],yin_m[i]))
            fo.close()

    #store the limits of the coordinates
    lims = {'xmax' : np.amax(xref_m), 'xmin' : np.amin(xref_m), 'ymax' : np.amax(yref_m), 'ymin' : np.amin(yref_m)}
    
    #return all matches within the tolerance
    return xref_m, yref_m, xin_m, yin_m, lims



def match_diff_im_plot(xrefm, yrefm, xinm, yinm, **kwargs):

    #put here because this has to be imported after georun is executed
    #to prevent the code from crashing
    import matplotlib.pyplot as plt
    '''PURPOSE: Plot panels of differences in x and y for a set of
    matched catalogs.

    INPUT PARAMETERS:

    xrefm, yrefm, xinm, yinm: A set of coordinates for matched
    objects.  These arrays must all be the same length.

    OPTIONAL KEWORD PARAMETERS

    plotfile: the name of the file containing the plot

    xmin, xmax, ymin, ymax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the x and ydiff plots.
    If one is given, all must be given.

    OUTPUT:

    a plot of the differences between the x and y  coordinates

    '''

    #set LaTeX fonts for labels
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    
    lims = {'xmax' : np.amax(xrefm), 'xmin' : np.amin(xrefm), 'ymax' : np.amax(yrefm), 'ymin' : np.amin(yrefm)}
    print("plotlims are ",lims)

    padfac = 0.03
    xdiff = (xrefm - xinm)
    ydiff = (yrefm - yinm)
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    xlims = [lims['xmin'] * (1. - padfac), lims['xmax'] * (1. + padfac)]

    #handles how to do padding if y-values are positive or negative
    if lims['ymin'] < 0:
        ylims = [lims['ymin'] * (1. + padfac/10.), lims['ymax'] * (1. - padfac/10.)]
    else:
        ylims = [lims['ymin'] * (1. - padfac/10.), lims['ymax'] * (1. + padfac/10.)]
    yline = [0,0]
    yline1 = [-3.125, -3.125]
    yline2 = [3.125, 3.125]

    #finds source outside of x and y lims.  Assumes that if one
    #keyword is given that all are given
    if 'xmin' in kwargs.keys():
        iin = np.where((xrefm >= kwargs['xmin']) & (xrefm <= kwargs['xmax']) & \
                       (yrefm >= kwargs['ymin']) & (yrefm <= kwargs['ymax']))
        iout = np.where(((xrefm < kwargs['xmin']) | (xrefm > kwargs['xmax'])) | \
                        ((yrefm < kwargs['ymin']) | (yrefm > kwargs['ymax'])))
    
    plt.clf()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    #color code the points in and outside the area where transform is
    #valid, if the x and y limits are given.
    if 'xmin' in kwargs.keys():
        ax1.scatter(xrefm[iout], xdiff[iout], color='y', edgecolors='k')
        ax1.scatter(xrefm[iin], xdiff[iin], color='b', edgecolors='k')
        medxdiff = np.median(xdiff[iin])
    else:
        ax1.scatter(xrefm, xdiff, color='b', edgecolors='k')
        medxdiff = np.median(xdiff)

    ax1.text(1.002 * xlims[0], 2.5, medxdiff, color='r')
    ax1.plot(xlims, yline, color='r')
    ax1.plot(xlims, yline1, color='r', linestyle = ':')
    ax1.plot(xlims, yline2, color='r', linestyle = ':')
    ax1.set_xlim(xlims)
    ax1.set_ylim([-19., 19.])
    ax1.set_ylabel(r'xref - xin')

    #color code the points in and outside the area where transform is
    #valid, if the x and y limits are given.
    if 'xmin' in kwargs.keys():
        ax2.scatter(yrefm[iout], xdiff[iout], color='y', edgecolors='k')
        ax2.scatter(yrefm[iin], xdiff[iin], color='b', edgecolors='k')
    else:
        ax2.scatter(yrefm, xdiff, color='b', edgecolors='k')

    #ax2.scatter(decrefm, radiff)
    ax2.plot(ylims, yline, color='r')
    ax2.plot(ylims, yline1, color='r', linestyle = ':')
    ax2.plot(ylims, yline2, color='r', linestyle = ':')
    #ax2.set_xticks(np.arange(min(decrefm), max(decrefm)+0.04, 0.04))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax2.set_xlim(ylims)
    ax2.set_ylim([-19.,19.])

    #color code the points in and outside the area where transform is
    #valid, if the x and y limits are given.
    if 'xmin' in kwargs.keys():
        ax3.scatter(xrefm[iout], ydiff[iout], color='y', edgecolors='k')
        ax3.scatter(xrefm[iin], ydiff[iin], color='b', edgecolors='k')
        medydiff = np.median(ydiff[iin])
    else:
        ax3.scatter(xrefm, ydiff, color='b', edgecolors='k')
        medydiff = np.median(ydiff)

    ax3.text(1.002 * xlims[0], 2.5, medydiff, color='r')
    #ax3.scatter(xrefm, ydiff)
    ax3.plot(xlims, yline, color='r')
    ax3.plot(xlims, yline1, color='r', linestyle = ':')
    ax3.plot(xlims, yline2, color='r', linestyle = ':')
    ax3.set_xlim(xlims)
    ax3.set_ylim([-19.,19.])
    ax3.set_xlabel(r'x')
    ax3.set_ylabel(r'yref - yin')
    
    #color code the points in and outside the area where transform is
    #valid, if the x and y limits are given.
    if 'xmin' in kwargs.keys():
        ax4.scatter(yrefm[iout], ydiff[iout], color='y', edgecolors='k')
        ax4.scatter(yrefm[iin], ydiff[iin], color='b', edgecolors='k')
    else:
        ax4.scatter(yrefm, ydiff, color='b', edgecolors='k')

    #    ax4.scatter(yrefm, ydiff)
    ax4.plot(ylims, yline, color='r')
    ax4.plot(ylims, yline1, color='r', linestyle = ':')
    ax4.plot(ylims, yline2, color='r', linestyle = ':')
    #ax4.set_xticks(np.arange(min(yrefm), max(yrefm)+0.04, 0.04))
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax4.set_xlim(ylims)
    ax4.set_ylim([-19.,19.])
    ax4.set_xlabel(r'y')

    keys = sorted(kwargs.keys())
    for kw in keys:
        if kw == 'plotfile':
            plt.savefig(kwargs[kw])
    #plt.show()
    #print("done with plot")
    #plt.close()
    
