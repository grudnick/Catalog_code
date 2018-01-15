import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from pyraf import iraf 
   

def cat_sky_match(raref, decref, rain, decin, septol, **kwargs):

    '''Written by Gregory Rudnick 9 January 2018

    PURPOSE:

    Take two lists of ra and dec coordinates in the same units and match
    them within some tolerance.

    Plot the differences in each coordinate.

    INPUT PARAMETERS:

    raref, decref: the reference coordinates in degrees.  Numpy arrays.

    rain, decin: the input coordinates in degrees.  If performing
    coordinate transforms, these would be the ones to be transformed.
    Numpy arrays

    septol: the maximum separation allowed for a match.  Units are
    arcseconds

    OPTIONAL KEWYORD PARAMETERS

    matchfile: a name of the file that will contain the reference and
    input coordinates.  Suitable for geomap input.

    OUTPUT

    arrays containing row-matched ras and decs of the reference and input
    coordinates that fall with septol

    '''

    refcat = SkyCoord(ra=raref*u.degree, dec=decref*u.degree)
    incat = SkyCoord(ra=rain*u.degree, dec=decin*u.degree)
    
    #match catalogs
    (idx, d2d, d3d) = refcat.match_to_catalog_sky(incat)
    #convert distance to arcseconds
    #d2d = d2d*3600

    #print(idx,d2d)
    
    #select close matches
    iclose = np.where(d2d < septol/3600.*u.deg)
    #convert tuple to pure array
    iclose=iclose[0]

    #print(d2d[indx])
    #print(np.column_stack((raref[iclose],rain[idx[iclose]],d2d[iclose]/u.deg*3600.)))


    #write files
    keys = sorted(kwargs.keys())
    for kw in keys:
        if kw == 'matchfile':
            #print(kwargs[kw])
            #open file for writing and write a header
            fo = open(kwargs[kw], "w")
            fo.write("# raref decref rain decin\n")
            for i,val in enumerate(iclose):
                #print(i,iclose[i])
                #print(raref[iclose[i]],decref[iclose[i]],rain[idx[iclose[i]]],decin[idx[iclose[i]]])
                fo.write('{} {} {} {}\n'.format(raref[iclose[i]],decref[iclose[i]],rain[idx[iclose[i]]],decin[idx[iclose[i]]]))
            fo.close()

    #store the limits of the coordinates
    lims = {'ramax' : np.amax(raref[iclose]), 'ramin' : np.amin(raref[iclose]), 'decmax' : np.amax(decref[iclose]), 'decmin' : np.amin(decref[iclose])}
    
    #return all matches within the tolerance
    return raref[iclose], decref[iclose], rain[idx[iclose]], decin[idx[iclose]], lims



def match_diff_plot(rarefm, decrefm, rainm, decinm, **kwargs):

    '''PURPOSE: Plot panels of differences in ra and dec for a set of
    matched catalogs.

    INPUT PARAMETERS:

    rarefm, decrefm, rainm, decinm: A set of coordinates for matched
    objects.  These arrays must all be the same length.

    OPTIONAL KEWORD PARAMETERS

    plotfile: the name of the file containing the plot

    OUTPUT:

    a plot of the differences between the RA and DEC coordinates

    '''

    #set LaTeX fonts for labels
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)

    
    lims = {'ramax' : np.amax(rarefm), 'ramin' : np.amin(rarefm), 'decmax' : np.amax(decrefm), 'decmin' : np.amin(decrefm)}
    print("plotlims are ",lims)

    padfac = 0.001
    radiff = (rarefm - rainm)*3600.
    decdiff = (decrefm - decinm)*3600.
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    ralims = [lims['ramin'] * (1. - padfac), lims['ramax'] * (1. + padfac)]

    #handles how to do padding if DECs are positive or negative
    if lims['decmin'] < 0:
        declims = [lims['decmin'] * (1. + padfac/10.), lims['decmax'] * (1. - padfac/10.)]
    else:
        declims = [lims['decmin'] * (1. - padfac/10.), lims['decmax'] * (1. + padfac/10.)]
    yline = [0,0]
    yline1 = [-0.5, -0.5]
    yline2 = [0.5, 0.5]

    plt.clf()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    ax1.scatter(rarefm, radiff)
    ax1.plot(ralims, yline, color='r')
    ax1.plot(ralims, yline1, color='r', linestyle = ':')
    ax1.plot(ralims, yline2, color='r', linestyle = ':')
    ax1.set_xlim(ralims)
    ax1.set_ylim([-3.,3.])
    ax1.set_ylabel(r'\Delta RA')

    ax2.scatter(decrefm, radiff)
    ax2.plot(declims, yline, color='r')
    ax2.plot(declims, yline1, color='r', linestyle = ':')
    ax2.plot(declims, yline2, color='r', linestyle = ':')
    #ax2.set_xticks(np.arange(min(decrefm), max(decrefm)+0.04, 0.04))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax2.set_xlim(declims)
    ax2.set_ylim([-3.,3.])

    ax3.scatter(rarefm, decdiff)
    ax3.plot(ralims, yline, color='r')
    ax3.plot(ralims, yline1, color='r', linestyle = ':')
    ax3.plot(ralims, yline2, color='r', linestyle = ':')
    ax3.set_xlim(ralims)
    ax3.set_ylim([-3.,3.])
    ax3.set_xlabel(r'RA')
    ax3.set_ylabel(r'\Delta Dec')
    
    ax4.scatter(decrefm, decdiff)
    ax4.plot(declims, yline, color='r')
    ax4.plot(declims, yline1, color='r', linestyle = ':')
    ax4.plot(declims, yline2, color='r', linestyle = ':')
    #ax4.set_xticks(np.arange(min(decrefm), max(decrefm)+0.04, 0.04))
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax4.set_xlim(declims)
    ax4.set_ylim([-3.,3.])
    ax4.set_xlabel(r'Dec')

    keys = sorted(kwargs.keys())
    for kw in keys:
        if kw == 'plotfile':
            plt.savefig(kwargs[kw])
    #plt.show()
    #print("done with plot")
    #plt.close()
    
