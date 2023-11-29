# Radial Collimator Designer
#
# Authors: S. Kowarik, L. Pithan, A. Zykov,  H. Wilming, L. Bogula,
#   S. Boitano, P. Schaefer, H. Pithan, and F. Carla
#
# Python version by Jonathan Rawle, April 2020

import numpy as np
import scipy.optimize
from scipy.special import erfc

## Geometrical parameters of collimator structure
detDis = 500    # Sample detector distance in mm
xdis1 = 100      # Distance from sample to front side on x axis in mm
xdis2 = 356     # Distance from sample to back side on x mm axis
detYsize = 260  # Horizontal detector size in mm
detZsize = 290   # Vertical detector size in mm
detYminus = -detYsize/2 # Offset to define direct beam position in y - direction 
detZminus = -15 # Offset to define direct beam position in z - direction 
oAng = 1  # Opening angle for one indivitual channel in degree 
ijrange = 30 # max number of grid points in each directions ... must be chosen large enough to cover the whole collimator area; depends on detector size and cannel opening angle
rot = -15

# Parameters that should usually stay the scame for all calculations
margin = np.radians( 14 * oAng) # adapt for width of soller channels !!!
R1 = 0.8*xdis1 # defines the radius of the front side calculation; will be cropped later 
R2 = 1.5*xdis2 # defines the radius of the back side calculation; will be cropped later 
outerStructurePath = "RadialCollimatorBox_py.json" # Path to save dimensions of outer box for CAD package
innerStructurePath =  "RadialCollimatorGrid_py.json" # Path to save honeycomb structure for CAD package

## Transformation functions and constants

# Basis vectors
a1 = np.array([np.sin(np.pi/3), np.cos(np.pi/3)])
a2 = np.array([np.sin(np.pi/3), -np.cos(np.pi/3)])

# Scale factor
alphafactor = 0.0027779

# Transformations
def SphereToCartesian(*args):
    r, theta, phi =  np.concatenate([np.array([a]).flatten() for a in args])
    return [r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)]

def SolveSphereToCartesian(r, theta, phi, K):
    f = SphereToCartesian(r, theta, phi)
    return abs(f[0] - K)

def ScaledSphereToCartesian(thetaphi, abu):
    return SphereToCartesian([1, np.pi * np.arctan(abu[0] * thetaphi[0]**abu[1]), thetaphi[1]])

def ScaledSphereToCartesian2(thetaphi, abu, r):
    return SphereToCartesian([1, np.pi * np.arctan(r * abu[0] * thetaphi[0]**abu[1]), thetaphi[1]])

def ToSpherical(xyz):
    return (np.arccos(xyz[2] / np.sqrt(xyz[0]**2 + xyz[1]**2 + xyz[2]**2)), np.arctan2(xyz[1], xyz[0]))

# Rotation matrices

def Ry():
    return np.array([[np.cos(np.pi/2), 0, np.sin(np.pi/2)],
                     [        0,       1,       0],
                     [-np.sin(np.pi/2), 0, np.cos(np.pi/2)]])

def Rx(phi):
    global rot
    a = np.array([[np.cos(np.radians(-rot)), np.sin(np.radians(-rot)), 0],
                  [-np.sin(np.radians(-rot)), np.cos(np.radians(-rot)), 0],
                  [             0,                     0,             1]])
    b = np.array([[1,     0,          0],
                  [0, np.cos(phi), np.sin(phi)],
                  [0, -np.sin(phi), np.cos(phi)]])
    return(np.dot(a, b))

# Transfer functions from lattice points to theta, phi couples

def NonLinScale(x):
    return np.exp(0.045 * x) + 0.5 * erfc((x - 4) * 0.1) + 2/(x + 0.1)

def ToPolarCoordinates(v):
    return [np.sqrt(v[0]**2 + v[1]**2), np.arctan2(v[1], v[0])]

def PointOnUnitSphere(i, j, a):
    if (i==0 and j==0):
        return np.linalg.multi_dot([ Ry(), Rx(a[2]), np.array([0, 0, 1]) ])
    else:
        return np.linalg.multi_dot([ Ry(), Rx(a[2]),
                                     ScaledSphereToCartesian2(ToPolarCoordinates(i*a1 + j*a2), a,
                                                             NonLinScale(np.linalg.norm(i*a1 + j*a2)))
                                     ])

def ToSpherical_PointOnUnitSphere(i, j, a):
    return [x for x in (np.concatenate(([i, j], ToSpherical(PointOnUnitSphere(i, j, a)))))]
    
def Normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return (v / norm)

class infDict(dict):
    def __getitem__(self, key):
        return dict.get(self, key, np.inf)

def getWall(latticepointSp, cen, r1, r2):
    a = latticepointSp[(cen[0]-1, cen[1])]
    b = latticepointSp[(cen[0], cen[1]-1)]
    c = latticepointSp[(cen[0]+1, cen[1]-1)]
    d = latticepointSp[(cen[0]+1, cen[1])]
    tmp = []
    if(np.max(a) != np.inf and np.max(a) != np.inf):
        tmp.append([SphereToCartesian(r1, a), SphereToCartesian(r2, a), SphereToCartesian(r2, b), SphereToCartesian(r1, b)])
    if(np.max(b) != np.inf and np.max(c) != np.inf):
        tmp.append([SphereToCartesian(r1, b), SphereToCartesian(r2, b), SphereToCartesian(r2, c), SphereToCartesian(r1, c)])
    if(np.max(c) != np.inf and np.max(d) != np.inf):
        tmp.append([SphereToCartesian(r1, c), SphereToCartesian(r2, c), SphereToCartesian(r2, d), SphereToCartesian(r1, d)])
    return tmp

def soller(detDis, detYminus, detZminus, detYsize, detZsize, xdis1, xdis2, oAng, alphafactor, margin, ijrange, rot, R1, R2, format='json'):

    # Calculation of outside box    
    theta1, phi1 = ToSpherical([detDis, detYminus, detZminus])
    theta2, phi2 = ToSpherical([detDis, detYminus, detZminus + detZsize])
    theta3, phi3 = ToSpherical([detDis, detYminus + detYsize, detZminus + detZsize])
    theta4, phi4 = ToSpherical([detDis, detYminus + detYsize, detZminus])

    c1ThetaMinPhiMin = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta2, phi2, xdis1)).x, theta2, phi2), 1)
    c1ThetaMaxPhiMin = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta1, phi1, xdis1)).x, theta1, phi1), 1)
    c1ThetaMinPhiMax = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta3, phi3, xdis1)).x, theta3, phi3), 1)
    c1ThetaMaxPhiMax = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta4, phi4, xdis1)).x, theta4, phi4), 1)

    c2ThetaMinPhiMin = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta2, phi2, xdis2)).x, theta2, phi2), 1)
    c2ThetaMaxPhiMin = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta1, phi1, xdis2)).x, theta1, phi1), 1)
    c2ThetaMinPhiMax = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta3, phi3, xdis2)).x, theta3, phi3), 1)
    c2ThetaMaxPhiMax = np.round(SphereToCartesian(scipy.optimize.minimize_scalar(SolveSphereToCartesian,  args=(theta4, phi4, xdis2)).x, theta4, phi4), 1)

    # For further calculation
    thetaMin = np.min([theta1, theta2, theta3, theta4])
    thetaMax = np.max([theta1, theta2, theta3, theta4])
    phiMin = np.min([phi1, phi2, phi3, phi4])
    phiMax = np.max([phi1, phi2, phi3, phi4])

    # Calculate the faces / corners of the box
    faceFront = [c1ThetaMinPhiMin, c1ThetaMaxPhiMin, c1ThetaMaxPhiMax, c1ThetaMinPhiMax]
    faceBack = [c2ThetaMinPhiMin, c2ThetaMaxPhiMin, c2ThetaMaxPhiMax, c2ThetaMinPhiMax]
    fNormFront = Normalize(np.cross(c1ThetaMaxPhiMax - c1ThetaMaxPhiMin,
                                    c1ThetaMinPhiMin - c1ThetaMaxPhiMax))
    fNormBack = -Normalize(np.cross(c2ThetaMaxPhiMax - c2ThetaMaxPhiMin, 
                                    c2ThetaMinPhiMin - c2ThetaMaxPhiMax))

    faceLeft = [c1ThetaMinPhiMin, c2ThetaMinPhiMin, c2ThetaMaxPhiMin, c1ThetaMaxPhiMin]
    faceRight = [c1ThetaMinPhiMax, c2ThetaMinPhiMax, c2ThetaMaxPhiMax, c1ThetaMaxPhiMax]
    fNormLeft = Normalize(np.cross(c1ThetaMinPhiMin - c2ThetaMinPhiMin, 
                                   c1ThetaMinPhiMin - c2ThetaMaxPhiMin))
    fNormRight = -Normalize(np.cross(c1ThetaMinPhiMax - c2ThetaMinPhiMax,
                                     c1ThetaMinPhiMax - c2ThetaMaxPhiMax))
    
    faceTop = [c1ThetaMaxPhiMin, c2ThetaMaxPhiMin, c2ThetaMaxPhiMax, c1ThetaMaxPhiMax]
    faceBottom = [c1ThetaMinPhiMin, c2ThetaMinPhiMin, c2ThetaMinPhiMax, c1ThetaMinPhiMax]
    fNormTop = Normalize(np.cross(c1ThetaMaxPhiMin - c2ThetaMaxPhiMin,
                                  c1ThetaMaxPhiMin - c2ThetaMaxPhiMax))
    fNormBottom = -Normalize(np.cross(c1ThetaMinPhiMin - c2ThetaMinPhiMin, 
                                      c1ThetaMinPhiMin - c2ThetaMinPhiMax))

    box = [[faceTop, fNormTop], [faceBottom, fNormBottom], [faceLeft, fNormLeft],
           [faceRight, fNormRight], [faceFront, fNormFront], [faceBack, fNormBack]]

    cornersExport = ('"c1ThetaMinPhiMin"=' + str(c1ThetaMinPhiMin.tolist()) +
                                      ';\n"c1ThetaMinPhiMax"=' + str(c1ThetaMinPhiMax.tolist()) +
                                      ';\n"c1ThetaMaxPhiMin"=' + str(c1ThetaMaxPhiMin.tolist()) +
                                      ';\n"c1ThetaMaxPhiMax"=' + str(c1ThetaMaxPhiMax.tolist()) +
                                      ';\n"c2ThetaMinPhiMin"=' + str(c2ThetaMinPhiMin.tolist()) +
                                      ';\n"c2ThetaMinPhiMax"=' + str(c2ThetaMinPhiMax.tolist()) +
                                      ';\n"c2ThetaMaxPhiMin"=' + str(c2ThetaMaxPhiMin.tolist()) +
                                      ';\n"c2ThetaMaxPhiMax"=' + str(c2ThetaMaxPhiMax.tolist())).replace(' ','')
    if format == 'scad':
        boxExport = 'mybox=' + str(box).replace('array(','').replace(')','').replace(' ','') + ';'
        cornersExport = cornersExport.replace('"', '') + ';\n'    
    else:
        boxExport = ('{"mybox":' + str(box).replace('array(','').replace(')','').replace(' ','') + ',').replace('.,','.0,').replace('.]','.0]')
        cornersExport = (cornersExport.replace('=', ':').replace(';', ',') + '}\n').replace('.,','.0,').replace('.]','.0]')
    f = open(outerStructurePath, 'w')
    f.write(boxExport + '\n' + cornersExport)
    f.close()


    # Calculation of inner comb structure
    alpha = oAng * alphafactor
    params = [alpha, 1.000, np.radians(0)]
    # 3rd parameter: rotation of grid around beam axis - superseded by rot parameter?

    thetaMinG = thetaMin - margin
    thetaMaxG = thetaMax + margin
    phiMinG = phiMin - margin
    phiMaxG = phiMax + margin
    
    # Generating list of lattice points that represent the center of the combs
    centerIndex = np.fromfunction(lambda i,j: np.array([i-j, i+2*j-3*ijrange]), (ijrange*2+1,ijrange*2+1)).reshape(2,-1).transpose()
    centers = np.array(list(map(ToSpherical_PointOnUnitSphere, centerIndex[:,0], centerIndex[:,1], np.tile(params, (len(centerIndex[:,0]), 1)))))
    centersSpRed = centerIndex[np.where(np.logical_and(np.logical_and(centers[:,2] >= thetaMinG, centers[:,2] <= thetaMaxG), np.logical_and(centers[:,3] >= phiMinG, centers[:,3] <= phiMaxG)))]

    # Generating list of lattice points that produce the honey-comb structure
    latIndex = np.concatenate((np.add(centerIndex, np.tile([1,0], (len(centerIndex), 1))), np.add(centersSpRed, np.tile([0,1], (len(centersSpRed), 1)))))
    latticeSpRed =  np.array(list(map(ToSpherical_PointOnUnitSphere, latIndex[:,0], latIndex[:,1], np.tile(params, (len(latIndex[:,0]), 1)))))
    latticepointSp = infDict(map(lambda i: [([i][0][0],[i][0][1]),[[i][0][2],[i][0][3]]], latticeSpRed))
    myWalls = []
    for c in centersSpRed:
        myWalls.extend(getWall(latticepointSp, c, R1, R2))

    if format == 'scad':
        wallExport = 'mywalls=' + str(myWalls).replace(' ','') + ';'
    else:
        wallExport = ('{"mywalls":' + str(myWalls).replace(' ','') + '}').replace('.,','.0,').replace('.]','.0]')
    f = open(innerStructurePath, 'w')
    f.write(wallExport)
    f.close()
    
###################################    
soller(detDis, detYminus, detZminus, detYsize, detZsize, xdis1, xdis2, oAng, alphafactor, margin, ijrange, rot, R1, R2, 'json')
