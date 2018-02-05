import numpy as npimport scipydef ncgeodetictoxy(lat, lon, SGN):    # Local Variables: M, E, RHO, MC, RE, debug, SLAT, SL, lon, TC, lat, T, CDR, SGN, Y, X, PI, K, E2    # Function calls: disp, cos, sqrt, mfilename, abs, tan, pi, sin, ncgeodetictoxy    # NCGEODETICTOXY - Convert lat-long to X and Y coordinates on a stereographic grid    #    # function [X,Y,SLAT,K]=ncgeodetictoxy(lat,lon,SGN)    #    # This function converts latitudes and longitudes to x and y coordinates    # on a polar stereographic projection, and is essentially a JPL algorithm    # rewritten to work in Matlab.    #    # Input:    # lat  - latitude (degrees North)    # lon  - longitude (degrees East)    # SGN  - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere    #    # Output:    # X   - polar stereographic X coordinate (km)    # Y   - polar stereographic Y coordinate (km)    # SLAT - latitude of true distance (default is 70 degrees if omitted)    # K   - length scale on the stereographic projection    #    # Where (X,Y)=(0,0) is at the pole.    #    # *-------------------------------------------------------------------------*    # * This subroutine converts from geodetic latitude and longitude to Polar  *    # * Stereographic (X,Y) coordinates for the polar regions.  The equations   *    # * are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *    # * Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *    # * Printing Office.  See JPL Technical Memorandum 3349-8401 for further    *    # * details.                     			                                *    # * CODE CONVERTED FROM A JPL SCRIPT BY C. S. Morris and  V. J. Troisi      *    # *-------------------------------------------------------------------------*    #    # In addition to Morris and Troisi routine, a calculation of scale has also been    # made using equations on pp 161 and 315 of Snyder, J.P., 1987, "Map Projections:    # A Working Manual", USGS. Note that the equation for calculating K at the pole    # on page 161 (equation 21-35) incorrectly includes the radius of the Earth.    #    # See also: ncxytogeodetic    #    # Converted by Andrew Roberts 2012    # Department of Oceanography, Naval Postgraduate School    #    # $Id$
    # set defaults    E2 = 0.006693883    E = np.sqrt(E2) # eccentricity of Hughes ellipsoid    RE = 6378.273   # earth radius in km    SLAT = 70.      # latitude of true distance    PI = np.pi      # PI    CDR = 180./PI   # Conversion constant from degrees to radians
    
    # convert latitude and longitude to radians
    #lat = np.mat(lat)
    #lon = np.mat(lon)    lat = np.divide(np.abs(lat),CDR)    lon = np.divide(lon,CDR)
    
    #lat = np.array(lat)
    #lon = np.array(lon)
    
    T = np.divide(np.tan(PI/4.-lat/2.),np.power(np.divide(1.-E*np.sin(lat),1.+E*np.sin(lat)),E/2.))    M = np.divide(np.cos(lat),np.sqrt(1.-E2*np.power(np.sin(lat),2.)))
    
    SL = SLAT/CDR
    TC = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL),1.+E*np.sin(SL)),E/2.))
    MC = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))
    if np.abs(90.-SLAT)<1.E-5:
        RHO = np.array(np.dot(2.*RE,np.divide(T,np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1.-E))))))    else:        RHO = np.array(np.dot(RE,np.divide(np.dot(MC,T),TC)))            Y = np.dot(-SGN,np.dot(RHO,np.cos(np.dot(SGN,lon))))
    X = np.dot(SGN,np.dot(RHO,np.sin(np.dot(SGN,lon))))
    latcheck = np.abs(lat)
    loncheck = np.abs(lon)

    X[latcheck >= PI/2.] = 0.
    Y[loncheck >= PI/2.] = 0.
    
    K = np.divide(RHO,RHO)
    if RHO.any == 0.:   # Calculate scale at the pole
        K[RHO == 0.] = np.dot(0.5,np.divide(MC,np.dot(TC,np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1.-E))))))
    
    # Or elsewhere
    K[RHO != 0.] = np.divide(RHO,np.dot(RE,M))        return (X, Y, SLAT, K)