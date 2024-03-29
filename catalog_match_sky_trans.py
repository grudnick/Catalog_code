import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import median_abs_deviation as mad
from scipy.optimize import curve_fit
  

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

    OPTIONAL KEYWORD PARAMETERS

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


    #write files only if "matchfile" keyword is set
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

    print("RA for matching pairs in spectro catalog")
    print(raref[iclose])
    print("DEC for matching pairs in spectro catalog")
    print(decref[iclose])
    #store the limits of the coordinates
    lims = {'ramax' : np.amax(raref[iclose]), 'ramin' : np.amin(raref[iclose]), 'decmax' : np.amax(decref[iclose]), 'decmin' : np.amin(decref[iclose])}
    
    #return all matches within the tolerance
    return raref[iclose], decref[iclose], rain[idx[iclose]], decin[idx[iclose]], lims



def match_diff_sky_plot(rarefm, decrefm, rainm, decinm, detrend = False, **kwargs):

    import matplotlib.pyplot as plt

    '''PURPOSE: Plot panels of differences in ra and dec for a set of
    matched catalogs.

    INPUT PARAMETERS:

    rarefm, decrefm, rainm, decinm: A set of coordinates for matched
    objects.  These arrays must all be the same length.

    OPTIONAL KEWORD PARAMETERS

    plotfile: the name of the file containing the plot

    ramin, ramax, decmin, decmax.  These are the limits over which the
    transform was originally computed.  If these are given then it
    uses those limits to color the points in the ra and decdiff plots.
    If one is given, all must be given.

    detrend: Default - False.  If True, subtract the best fit line
    from the data and recompute statistics

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
    ralims = [lims['ramin'] - 2./60., lims['ramax'] + 2./60.]
    #ralims = [lims['ramin'] * (1. - padfac), lims['ramax'] * (1. + padfac)]

    #handles how to do padding if DECs are positive or negative
    if lims['decmin'] < 0:
        declims = [lims['decmin'] * (1. + padfac/10.), lims['decmax'] * (1. - padfac/10.)]
    else:
        declims = [lims['decmin'] * (1. - padfac/10.), lims['decmax'] * (1. + padfac/10.)]
    yline = [0,0]
    yline1 = [-0.1, -0.1]
    yline2 = [0.1, 0.1]

    #finds source outside of ra and dec lims.  Assumes that if one
    #keyword is given that all are given
    if 'ramin' in kwargs.keys():
        iin = np.where((rarefm >= kwargs['ramin']) & (rarefm <= kwargs['ramax']) & \
                       (decrefm >= kwargs['decmin']) & (decrefm <= kwargs['decmax']))
        iout = np.where(((rarefm < kwargs['ramin']) | (rarefm > kwargs['ramax'])) | \
                        ((decrefm < kwargs['decmin']) | (decrefm > kwargs['decmax'])))
    
    plt.clf()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10,10))

    #plot limits
    ypmax = 1.0
    ypmin = -1.0
    nsig_rej = 3.0
    
    #color code the points in and outside the area where transform is
    #valid, if the ra and dec limits are given.
    if 'ramin' in kwargs.keys():
        medradiff = np.median(radiff[iin])
        sigradiff = 1.4826 * mad(radiff[iin])
        #only take points within 3 sigma of the median
        ifit = iin & (abs(radiff - medradiff) < nsig_rej * sigradiff)
        inofit = iin & (abs(radiff - medradiff) >= nsig_rej * sigradiff)

        ax1.scatter(rarefm[iout], radiff[iout], color='y', edgecolors='k')
        ax1.scatter(rarefm[iin], radiff[iin], color='c', edgecolors='k')
        ax1.scatter(rarefm[ifit], radiff[ifit], color='b', edgecolors='k')

        ##fit line to data
        popt,pcov = curve_fit(flin, rarefm[ifit], radiff[ifit])
    else:
        medradiff = np.median(radiff)
        sigradiff = 1.4826 * mad(radiff)

        #only take points within 3 sigma of the median
        ifit = (abs(radiff - medradiff) < nsig_rej * sigradiff)
        inofit = (abs(radiff - medradiff) >= nsig_rej * sigradiff)

        ##fit line to data
        popt,pcov = curve_fit(flin, rarefm[ifit], radiff[ifit])
        ax1.scatter(rarefm, radiff, color='c', edgecolors='k')
        ax1.scatter(rarefm[ifit], radiff[ifit], color='b', edgecolors='k')

    #plot line fit
    slope= popt[0]
    yint = popt[1]
    print('slope = ',slope, 'yint = ', yint)
    xfit = np.array([min(rarefm),max(rarefm)])
    yfit = xfit * slope + yint
    ax1.plot(xfit,yfit,'g--')

        #print(1.002 * ralims[0],medradiff)    
    #ax1.text(1.002 * ralims[0], 2.5, medradiff, color='r')
    labsize=14
    ax1.text(np.median(rarefm) - 3./60., 0.85*ypmax,
                 r'median($\Delta RA$) = ' + str(round(medradiff,2)), color='r',fontsize=labsize)
    ax1.text(np.median(rarefm) - 3./60., 0.7*ypmax,
                 r'$\sigma(\Delta RA$) = ' + str(round(sigradiff,2)), color='r',fontsize=labsize)
    #compute change in delta RA over range in RA
    ax1.text(np.median(rarefm) - 3./60., 0.55*ypmax,
                 r'$\Delta_{tot} RA$ = ' + str(round((max(rarefm) - min(rarefm)) * slope,2)), color='r',fontsize=labsize)

    ax1.plot(ralims, yline, color='r')
    ax1.plot(ralims, yline1, color='r', linestyle = ':')
    ax1.plot(ralims, yline2, color='r', linestyle = ':')
    ax1.set_xlim(ralims)
    ax1.set_ylim([ypmin,ypmax])
    ax1.set_ylabel(r'\Delta RA')

    #detrend data if keyword is set
    if detrend:
        radiffcor = radiff - (rarefm * slope + yint)
        ax1.scatter(rarefm[ifit],radiffcor[ifit], color='',edgecolor='m')
        medradiffcor = np.median(radiffcor[ifit])
        sigradiffcor = 1.4826 * mad(radiffcor)
        ax1.text(np.median(rarefm) - 3./60., 0.7*ypmin,
                r'median($\Delta_{cor} RA$) = ' + str(round(medradiffcor,2)), color='m',fontsize=labsize)
        ax1.text(np.median(rarefm) - 3./60., 0.85*ypmin,
                r'$\sigma(\Delta_{cor} RA$) = ' + str(round(sigradiffcor,2)), color='m',fontsize=labsize)

       #rarefm, radiff

    #color code the points in and outside the area where transform is
    #valid, if the ra and dec limits are given.
    if 'ramin' in kwargs.keys():
        #only take points within 3 sigma of the median
        ifit = iin & (abs(radiff - medradiff) < nsig_rej * sigradiff)
        inofit = iin & (abs(radiff - medradiff) >= nsig_rej * sigradiff)

        ax2.scatter(decrefm[iout], radiff[iout], color='y', edgecolors='k')
        ax2.scatter(decrefm[iin], radiff[iin], color='c', edgecolors='k')
        ax2.scatter(decrefm[ifit], radiff[ifit], color='b', edgecolors='k')
        popt,pcov = curve_fit(flin, decrefm[ifit], radiff[ifit])       
    else:
        #only take points within 3 sigma of the median
        ifit = (abs(radiff - medradiff) < nsig_rej * sigradiff)
        inofit = (abs(radiff - medradiff) >= nsig_rej * sigradiff)

        ax2.scatter(decrefm, radiff, color='c', edgecolors='k')
        ax2.scatter(decrefm[ifit], radiff[ifit], color='b', edgecolors='k')
        popt,pcov = curve_fit(flin, decrefm[ifit], radiff[ifit])

    #plot line fit
    slope= popt[0]
    yint = popt[1]
    print('slope = ',slope, 'yint = ', yint)
    xfit = np.array([min(decrefm),max(decrefm)])
    yfit = xfit * slope + yint
    ax2.plot(xfit,yfit,'g--')

    #compute change in delta RA over range in DEC
    ax2.text(np.median(decrefm) - 3./60., 0.55*ypmax,
                 r'$\Delta_{tot} RA$ = ' + str(round((max(decrefm) - min(decrefm)) * slope,2)), color='r',fontsize=labsize)

    #ax2.scatter(decrefm, radiff)
    ax2.plot(declims, yline, color='r')
    ax2.plot(declims, yline1, color='r', linestyle = ':')
    ax2.plot(declims, yline2, color='r', linestyle = ':')
    #ax2.set_xticks(np.arange(min(decrefm), max(decrefm)+0.04, 0.04))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax2.set_xlim(declims)
    ax2.set_ylim([ypmin,ypmax])


    #detrend data if keyword is set
    if detrend:
        radiffcor = radiff - (decrefm * slope + yint)
        ax2.scatter(decrefm[ifit],radiffcor[ifit], color='',edgecolor='m')
        medradiffcor = np.median(radiffcor[ifit])
        sigradiffcor = 1.4826 * mad(radiffcor)
        ax2.text(np.median(decrefm) - 3./60., 0.7*ypmin,
                r'median($\Delta_{cor} RA$) = ' + str(round(medradiffcor,2)), color='m',fontsize=labsize)
        ax2.text(np.median(decrefm) - 3./60., 0.85*ypmin,
                r'$\sigma(\Delta_{cor} RA$) = ' + str(round(sigradiffcor,2)), color='m',fontsize=labsize)

    #color code the points in and outside the area where transform is
    #valid, if the ra and dec limits are given.
    if 'ramin' in kwargs.keys():
        meddecdiff = np.median(decdiff[iin])
        sigdecdiff = 1.4826 * mad(decdiff[iin])

        #only take points within 3 sigma of the median
        ifit = iin & (abs(decdiff - meddecdiff) < nsig_rej * sigdecdiff)
        inofit = iin & (abs(decdiff - meddecdiff) >= nsig_rej * sigdecdiff)

        ax3.scatter(rarefm[iout], decdiff[iout], color='y', edgecolors='k')
        ax3.scatter(rarefm[iin], decdiff[iin], color='c', edgecolors='k')
        ax3.scatter(rarefm[ifit], decdiff[ifit], color='b', edgecolors='k')
        ##fit line to data
        popt,pcov = curve_fit(flin, rarefm[ifit], decdiff[ifit])

    else:
        meddecdiff = np.median(decdiff)
        sigdecdiff = 1.4826 * mad(decdiff)

        #only take points within 3 sigma of the median
        ifit = (abs(decdiff - meddecdiff) < nsig_rej * sigdecdiff)
        inofit = (abs(decdiff - meddecdiff) >= nsig_rej * sigdecdiff)

        ax3.scatter(rarefm, decdiff, color='c', edgecolors='k')
        ax3.scatter(rarefm[ifit], decdiff[ifit], color='b', edgecolors='k')

        ##fit line to data
        popt,pcov = curve_fit(flin, rarefm[ifit], decdiff[ifit])


    #ax3.text(1.002 * ralims[0], 2.5, meddecdiff, color='r')
    ax3.text(np.median(rarefm) - 3./60., 0.85*ypmax,
                 r'median($\Delta Dec$) = ' + str(round(meddecdiff,2)), color='r',fontsize=labsize)
    ax3.text(np.median(rarefm) - 3./60., 0.7*ypmax,
                 r'$\sigma(\Delta Dec$) = ' + str(round(sigdecdiff,2)), color='r',fontsize=labsize)

    #plot line fit
    slope= popt[0]
    yint = popt[1]
    print('slope = ',slope, 'yint = ', yint)
    xfit = np.array([min(rarefm),max(rarefm)])
    yfit = xfit * slope + yint
    ax3.plot(xfit,yfit,'g--')

    #compute change in delta DEC over range in RA
    ax3.text(np.median(rarefm) - 3./60., 0.55*ypmax,
                 r'$\Delta_{tot} Dec$ = ' + str(round((max(rarefm) - min(rarefm)) * slope,2)), color='r',fontsize=labsize)
    #ax3.scatter(rarefm, decdiff)
    ax3.plot(ralims, yline, color='r')
    ax3.plot(ralims, yline1, color='r', linestyle = ':')
    ax3.plot(ralims, yline2, color='r', linestyle = ':')
    ax3.set_xlim(ralims)
    ax3.set_ylim([ypmin,ypmax])
    ax3.set_xlabel(r'RA')
    ax3.set_ylabel(r'\Delta Dec')
    
    if detrend:
        decdiffcor = decdiff - (rarefm * slope + yint)
        ax3.scatter(rarefm[ifit],decdiffcor[ifit], color='',edgecolor='m')
        meddecdiffcor = np.median(decdiffcor[ifit])
        sigdecdiffcor = 1.4826 * mad(decdiffcor)
        ax3.text(np.median(rarefm) - 3./60., 0.7*ypmin,
                r'median($\Delta_{cor} Dec$) = ' + str(round(meddecdiffcor,2)), color='m',fontsize=labsize)
        ax3.text(np.median(rarefm) - 3./60., 0.85*ypmin,
                r'$\sigma(\Delta_{cor} Dec$) = ' + str(round(sigdecdiffcor,2)), color='m',fontsize=labsize)

    #color code the points in and outside the area where transform is
    #valid, if the ra and dec limits are given.
    if 'ramin' in kwargs.keys():
        #only take points within 3 sigma of the median
        ifit = iin & (abs(decdiff - meddecdiff) < nsig_rej * sigdecdiff)
        inofit = iin & (abs(decdiff - meddecdiff) >= nsig_rej * sigdecdiff)

        ax4.scatter(decrefm[iout], decdiff[iout], color='y', edgecolors='k')
        ax4.scatter(decrefm[iin], decdiff[iin], color='c', edgecolors='k')
        ax4.scatter(decrefm[ifit], decdiff[ifit], color='b', edgecolors='k')
        ##fit line to data
        popt,pcov = curve_fit(flin, decrefm[ifit], decdiff[ifit])

    else:
        #only take points within 3 sigma of the median
        ifit = (abs(decdiff - meddecdiff) < nsig_rej * sigdecdiff)
        inofit = (abs(decdiff - meddecdiff) >= nsig_rej * sigdecdiff)

        ax4.scatter(decrefm, decdiff, color='c', edgecolors='k')
        ax4.scatter(decrefm[ifit], decdiff[ifit], color='b', edgecolors='k')
        ##fit line to data
        popt,pcov = curve_fit(flin, decrefm[ifit], decdiff[ifit])

    #plot line fit
    slope= popt[0]
    yint = popt[1]
    print('slope = ',slope, 'yint = ', yint)
    xfit = np.array([min(decrefm),max(decrefm)])
    yfit = xfit * slope + yint
    ax4.plot(xfit,yfit,'g--')


    #compute change in delta DEC over range in DEC
    ax4.text(np.median(decrefm) - 3./60., 0.55*ypmax,
                 r'$\Delta_{tot} RA$ = ' + str(round((max(decrefm) - min(decrefm)) * slope,2)), color='r',fontsize=labsize)
    #    ax4.scatter(decrefm, decdiff)
    ax4.plot(declims, yline, color='r')
    ax4.plot(declims, yline1, color='r', linestyle = ':')
    ax4.plot(declims, yline2, color='r', linestyle = ':')
    #ax4.set_xticks(np.arange(min(decrefm), max(decrefm)+0.04, 0.04))
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax4.set_xlim(declims)
    ax4.set_ylim([ypmin,ypmax])
    ax4.set_xlabel(r'Dec')

    if detrend:
        decdiffcor = decdiff - (decrefm * slope + yint)
        ax4.scatter(decrefm[ifit],decdiffcor[ifit], color='',edgecolor='m')
        meddecdiffcor = np.median(decdiffcor[ifit])
        sigdecdiffcor = 1.4826 * mad(decdiffcor)
        ax4.text(np.median(decrefm) - 3./60., 0.7*ypmin,
                r'median($\Delta_{cor} Dec$) = ' + str(round(meddecdiffcor,2)), color='m',fontsize=labsize)
        ax4.text(np.median(decrefm) - 3./60., 0.85*ypmin,
                r'$\sigma(\Delta_{cor} Dec$) = ' + str(round(sigdecdiffcor,2)), color='m',fontsize=labsize)

    keys = sorted(kwargs.keys())
    for kw in keys:
        if kw == 'plotfile':
            plt.savefig(kwargs[kw])
            plt.show()
    #print("done with plot")
    #plt.close()
    
    return medradiff, meddecdiff

#function for a straight line fit
def flin(x,A,B):
    return A*x + B
