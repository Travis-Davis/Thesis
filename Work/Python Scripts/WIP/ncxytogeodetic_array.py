import numpy as npimport scipydef ncxytogeodetic(X, Y, SGN):    # Local Variables: E, RHO, CM, RE, debug, SLAT, SL, lon, CHI, Y, ALAT, T, CDR, SGN, lat, X, ALONG, PI, E2    # Function calls: disp, cos, atan, sqrt, mfilename, wrapTo180, atan2, abs, ncxytogeodetic, tan, pi, sin    # NCXYTOGEODETIC - Convert X and Y coordinates on stereographic grid to lat-long    #    # function [lat,lon,SLAT]=ncxytogeodetic(X,Y,SGN)    #    # This function converts x and y coordinates on a polar stereographic grid    # to latitudes and longitudes, and is essentially a JPL algorithm rewritten    # to work in Matlab.    #    # Input:    # X   - polar stereographic X coordinate (km)    # Y   - polar stereographic Y coordinate (km)    # SGN - Hemisphere specifier: 1=Northern Hemisphere, -1=Southern Hemisphere    #    # Where (X,Y)=(0,0) is at the pole.    #    # Output:    # lat   - Latitude    # lon   - Longitude    # SLAT  - Latitude of true distance    #    # *-------------------------------------------------------------------------*    # * This subroutine converts from Polar Stereographic (X,Y) coordinates     *    # * to geodetic latitude and longitude for the polar regions. The equations *    # * are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *    # * Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *    # * Printing Office.  See JPL Technical Memorandum 3349-8401 for further    *    # * details.                                                                *    # * CODE CONVERTED FROM A JPL SCRIPT BY C. S. Morris and  V. J. Troisi      *    # * With added area calculations					                        *    # *-------------------------------------------------------------------------*    #    # See also: ncgeodetictoxy    #    # Converted by Andrew Roberts 2012    # Department of Oceanography, Naval Postgraduate School    #    # $Id$        # set defaults    E2 = 0.006693883    E = np.sqrt(E2)     # eccentricity of Hughes ellipsoid    RE = 6378.273       # earth radius in km    SLAT = 70.          # latitude of true distance    PI = np.pi          # PI    CDR = 180./np.pi    # Conversion constant from degrees to radians
    
    X = np.array(X)
    Y = np.array(Y)

    SL  = SLAT*PI/180.
    RHO = np.sqrt(np.power(X,2.)+np.power(Y,2.))
    CM  = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))
    T   = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL),1.+E*np.sin(SL)),E/2.))
    
    if np.abs((SLAT-90.))<1.E-5:        T = np.dot(RHO,np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1-E))))/(2*RE)    else:        T = np.divide(np.dot(RHO,T),RE*CM)                CHI = PI/2.-2.*np.arctan(T)
        ALAT = CHI+E2/2.+5.0/24.*np.power(E2,2.)+np.power(E2,3.)/12.*np.sin(2.*CHI)+7.0/48.*np.power(E2,2)+29.0/240.*np.power(E2,3)*np.sin(4*CHI)+7.0/120.*np.power(E2,3.)*np.sin(6.*CHI)    ALAT = SGN*ALAT*180./PI    ALONG = np.arctan2(SGN*X,-SGN*Y)    ALONG = SGN*ALONG*180./PI

    if RHO.any<0.1:
        ALAT[RHO < .1] = PI*SGN/2
        ALONG[RHO < .1] = 0.0    #-------------------------------------------------------------------------*

    lon = ALONG
    lat = ALAT        return (lat, lon)