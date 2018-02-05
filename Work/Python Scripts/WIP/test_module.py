# ncgeodetictoxy test script
import numpy as np
import scipy
from ncgeodetictoxy_single import ncgeodetictoxy
from ncxytogeodetic_single import ncxytogeodetic

# Test Case #1
lat1 = 80.
lon1 = 0.1
SGN1 = 1

# Test Case #2
lat2 = 60.
lon2 = 45.
SGN2 = 1

# Test Case #3
lat3 = -60.
lon3 = 45.
SGN3 = -1

# Test Evaluation
(X1,Y1,SLAT1,K1) = ncgeodetictoxy(lat1,lon1,SGN1)
(X2,Y2,SLAT2,K2) = ncgeodetictoxy(lat2,lon2,SGN2)
(X3,Y3,SLAT3,K3) = ncgeodetictoxy(lat3,lon3,SGN3)

(lat1r,lon1r)    = ncxytogeodetic(X1,Y1,SGN1)
(lat2r,lon2r)    = ncxytogeodetic(X2,Y2,SGN2)
(lat3r,lon3r)    = ncxytogeodetic(X3,Y3,SGN3)

print lat1r/lat1,lon1r/lon1r
print lat2r/lat2,lon2r/lon2r
print lat3r/lat3,lon3r/lon3r
