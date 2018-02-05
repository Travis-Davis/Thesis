# DEMSI_grid_gen.py - Generate Grid for use in DEMSI
#
# [X,Y,mask]=DEMSI_grid_gen(minlat,maxlat,minlon,maxlon,res)
#
# EXPLAIN WHAT THE FUNCTION DOES!
#
# Input:
# minlat  - min. latitude (degrees North)
# minlon  - min. longitude (degrees East)
# maxlat  - max. latitude (degrees North)
# maxlon  - max. longitude (degrees East)
# res     - grid resolution (km)
#
# Output:
# X    - polar stereographic X array (km)
# Y    - polar stereographic Y array (km)
# mask - Array of dimension [X,Y] with boolean value for land presence (1: z > 0)
#
# Where (X,Y)=(0,0) is at the pole.
#
# See also: demsi_xy2geodetic
#           demsi_geodetic2xy
#
# Created by Travis Davis January 6 2018
# Department of Oceanography, Naval Postgraduate School

import numpy as np
import scipy
from demsi_xy2geodetic import demsi_xy2geodetic
from demsi_geodetic2xy import demsi_geodetic2xy

def demsi_gridgen(minlat, maxlat, minlon, maxlon, res):

    # set defaults
    E2 = 0.006693883
    E = np.sqrt(E2) # eccentricity of Hughes ellipsoid
    RE = 6378.273   # earth radius in km
    SLAT = 70.      # latitude of true distance
    PI = np.pi      # PI
    CDR = 180./PI   # Conversion constant from degrees to radians

    # convert latitude and longitude to radians
    lat = np.divide(np.abs(lat),CDR)
    lon = np.divide(lon,CDR)

    T = np.divide(np.tan(PI/4.-lat/2.),np.power(np.divide(1.-E*np.sin(lat), \
        1.+E*np.sin(lat)),E/2.))
    M = np.divide(np.cos(lat),np.sqrt(1.-E2*np.power(np.sin(lat),2.)))

    SL = SLAT/CDR
    TC = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL),1.+ \
         E*np.sin(SL)),E/2.))
    MC = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))

    # special case of true latitude being at the pole
    if np.abs(90.-SLAT)<1.E-5:
        RHO = np.array(np.dot(2.*RE,np.divide(T,np.sqrt(np.dot(np.power(1.+E,1.+E), \
              np.power(1.-E,1.-E))))))
    else:
        RHO = np.array(np.dot(RE,np.divide(np.dot(MC,T),TC)))
        Y = np.dot(-SGN,np.dot(RHO,np.cos(np.dot(SGN,lon))))
        X = np.dot(SGN,np.dot(RHO,np.sin(np.dot(SGN,lon))))

    # special case of being at the pole
    if np.abs(lat) >= PI/2.:
        X = 0.
        Y = 0.

    if RHO == 0.:   # Calculate scale at the pole
        K = np.dot(0.5,np.divide(MC,np.dot(TC,np.sqrt(np.dot(np.power(1.+E,1.+E), \
            np.power(1.-E,1.-E))))))
    else:           # Or elsewhere
        K = np.divide(RHO,np.dot(RE,M))

    return (X, Y, SLAT, K)
