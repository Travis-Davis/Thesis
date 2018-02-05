# ncgeodetictoxy test script
import numpy as np
import scipy
from ncgeodetictoxy import ncgeodetictoxy

# Test Case #1
lat1 = 80.
lon1 = 0.
SGN1 = 1

# Test Case #2
lat2 = 60.
lon2 = 45.
SGN2 = 1

# Test Case #3
lat3 = 60.
lon3 = 45.
SGN3 = -1

# Test 1 Evaluation
(X1,Y1,SLAT1,K1) = ncgeodetictoxy(lat1,lon1,SGN1)
print 'Test Case #1'
print 'lat = ',lat1,'| lon = ',lon1,'| SGN = ',SGN1
print 'X = ',X1
print 'Y = ',Y1

# Test 2 Evaluation
(X2,Y2,SLAT2,K2) = ncgeodetictoxy(lat2,lon2,SGN2)
print 'Test Case #2'
print 'lat = ',lat2,'| lon = ',lon2,'| SGN = ',SGN2
print 'X = ',X2
print 'Y = ',Y2

# Test 3 Evaluation
(X3,Y3,SLAT3,K3) = ncgeodetictoxy(lat3,lon3,SGN3)
print 'Test Case #3'
print 'lat = ',lat3,'| lon = ',lon3,'| SGN = ',SGN3
print 'X = ',X3
print 'Y = ',Y3
