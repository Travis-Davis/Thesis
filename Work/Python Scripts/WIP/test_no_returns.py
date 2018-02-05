# ncgeodetictoxy test script
import numpy as np
import scipy

# Test Case #1
lat1 = 80.
lon1 = 0.1
SGN1 = 1

E2 = 0.006693883
E = np.sqrt(E2) # eccentricity of Hughes ellipsoid
RE = 6378.273   # earth radius in km
SLAT = 70.      # latitude of true distance
PI = np.pi      # PI
CDR = 180./PI   # Conversion constant from degrees to radians
    
# convert latitude and longitude to radians
lat = lat1
lon = lon1
SGN = SGN1

lat = np.divide(np.abs(lat),CDR)
lon = np.divide(lon,CDR)
    
T = np.divide(np.tan(PI/4.-lat/2.),np.power(np.divide(1.-E*np.sin(lat),1.+E*np.sin(lat)),E/2.))
M = np.divide(np.cos(lat),np.sqrt(1.-E2*np.power(np.sin(lat),2.)))
    
SL = SLAT/CDR
TC = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL),1.+E*np.sin(SL)),E/2.))
MC = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))
    
if np.abs(90.-SLAT)<1.E-5:
    RHO = np.array(np.dot(2.*RE,np.divide(T,np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1.-E))))))
else:
    RHO = np.array(np.dot(RE,np.divide(np.dot(MC,T),TC)))

Y = np.dot(-SGN,np.dot(RHO,np.cos(np.dot(SGN,lon))))
X = np.dot(SGN,np.dot(RHO,np.sin(np.dot(SGN,lon))))

if np.abs(lat) >= PI/2.:
    X = 0.
    Y = 0.
    
if RHO == 0.:   # Calculate scale at the pole
    K = np.dot(0.5,np.divide(MC,np.dot(TC,np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1.-E))))))
else:           # Or elsewhere
    K = np.divide(RHO,np.dot(RE,M))

SL  = SLAT*PI/180.
RHO = np.sqrt(np.power(X,2.)+np.power(Y,2.))
CM  = np.divide(np.cos(SL),np.sqrt(1.-E2*np.power(np.sin(SL),2.)))
T   = np.divide(np.tan(PI/4.-SL/2.),np.power(np.divide(1.-E*np.sin(SL),1.+E*np.sin(SL)),E/2.))
    
if np.abs((SLAT-90.))<1.E-5:
    T = RHO*np.sqrt(np.dot(np.power(1.+E,1.+E),np.power(1.-E,1-E)))/(2*RE)
else:
    T = np.divide(np.dot(RHO,T),RE*CM)

CHI = PI/2.-2.*np.arctan(T)
    
ALAT  = CHI+(E2/2.+5.0/24.*np.power(E2,2.)+np.power(E2,3.)/12.)*np.sin(2.*CHI)+(7.0/48.*np.power(E2,2.)+29.0/240.*np.power(E2,3.))*np.sin(4.*CHI)+7.0/120.*np.power(E2,3.)*np.sin(6.*CHI)
ALAT  = SGN*ALAT*180./PI
ALONG = np.arctan2(SGN*X,-SGN*Y)
ALONG = SGN*ALONG*180./PI
    
if RHO < 0.1:
    ALAT  = PI*SGN/2
    ALONG = 0.0
    
#-------------------------------------------------------------------------*
    
lon1r = ALONG
lat1r = ALAT

print lat1r/lat1,lon1r/lon1r
