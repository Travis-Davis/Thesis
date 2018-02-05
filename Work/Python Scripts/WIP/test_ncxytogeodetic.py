# ncgeodetictoxy test script
import numpy as np
import scipy
from ncxytogeodetic import ncxytogeodetic

# Test Case #1
X1 = 0.
Y1 = -1085.94318794
SGN1 = 1.

# Test Case #2
X2 = 2349.87883566
Y2 = -2349.8788566
SGN2 = 1.

# Test Case #3
X3 = 2349.87883566
Y3 = 2349.87883566
SGN3 = -1.

# Test 1 Evaluation
(lat1,lon1) = ncxytogeodetic(X1,Y1,SGN1)
print 'Test Case #1'
print 'X = ',X1,'| Y = ',Y1,'| SGN = ',SGN1
print 'lat = ',lat1
print 'lon = ',lon1

# Test 2 Evaluation
(lat2,lon2) = ncxytogeodetic(X2,Y2,SGN2)
print 'Test Case #2'
print 'X = ',X2,'| Y = ',Y2,'| SGN = ',SGN2
print 'lat = ',lat2
print 'lon = ',lon2

# Test 3 Evaluation
(lat3,lon3) = ncxytogeodetic(X3,Y3,SGN3)
print 'Test Case #3'
print 'X = ',X3,'| Y = ',Y3,'| SGN = ',SGN3
print 'lat = ',lat3
print 'lon = ',lon3
