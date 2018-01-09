import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
   

def cat_sky_match(raref, decref, rain, decin):

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

    septol: the maximum separation allowed for a match in the native
    units of the catalog

    '''

    refcat = SkyCoord(ra=raref*u.degree, dec=decref*u.degree)
    incat = SkyCoord(ra=rain*u.degree, dec=decin*u.degree)
    
    #match catalogs
    idx, d2d, d3d = refcat.match_to_catalog_sky(catalog)

    print(idx)

