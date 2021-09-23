#!/usr/bin/env python 3.6.8 Anaconda 64-bit
# -*- coding: utf-8 -*-
#Created on Mon Sep 20 03:30:52 2021
#@author: frankgrijalva

#imports and customs check
import numpy as np

#functions
def Lambdanator (m, grateSpacing, distToTarget, diffracMax0Dist, diffracLocalMaxDist):

     DM = grateSpacing / m
     theta0 = diffracMax0Dist / distToTarget

     lamb = ((DM) * (np.cos(np.arctan(theta0)/2) - np.cos(-1*(np.arctan(theta0)/2) + np.arctan(diffracLocalMaxDist/distToTarget))))
     return lamb

#vars
distToBoard = 8.02 #±.001m
distFromBeamToRef = 0.85 #±.0005m
angleOfRef = np.arctan(distFromBeamToRef/distToBoard)/2
diffGratingSpacing = 1E-3 #≠0, "known" number in m

#array of local maximum spacing data
mDistDillArr = np.asarray([0.85, 0.92, 0.99, 1.05, 1.1, 1.15, 1.2, 1.24, 1.28, 1.32, 1.36, 1.39, 1.43])


#array for ans storage
mLambDillArr = np.zeros_like(mDistDillArr)

#looping over m
for i in range(len(mDistDillArr)):
     #feeding the func "1" for the 0th max so it doesn't divide by 0
     if i == 0:
          mLambDillArr[0] = Lambdanator((1), diffGratingSpacing, distToBoard, distFromBeamToRef, mDistDillArr[1])
     else:
          mLambDillArr[i] = Lambdanator((i), diffGratingSpacing, distToBoard, distFromBeamToRef, mDistDillArr[i])
     print("\n M"+str(i), "dist. from M0:", mDistDillArr[i], "\n calculated lambda:", mLambDillArr[i])


#checking results
meanWavelength = np.mean(mLambDillArr)*1E9
print("\n\n\nMean wavelength:", meanWavelength, "nm")



#error calcs
#more imports
import sympy as sp
#sympy symbol declaration
m, d, y, x, z = sp.symbols("m d y x z")

#length uncertainties, x is central max, y is target dist, z is local max.
dx = 0.0005 #m
dy = 0.001 #m
dz = 0.0005 #m

#function to diff
wavEQ = (d/m) * (sp.cos(sp.atan(x/y)/2) - sp.cos(-.5*(sp.atan(x/y)/2) + sp.atan(z/y)))

#finding partial derivatives
lWRTx = sp.diff(wavEQ, x)
lWRTy = sp.diff(wavEQ, y)
lWRTz = sp.diff(wavEQ, z)


#3.47 from Taylor - err for a function of many vars
dq = sp.sqrt((lWRTx * dx)**2 + (lWRTy * dy)**2 + (lWRTz * dz)**2)
print("\npartial diff for error:",dq)

#checking the results
dqNum = dq.subs([(m,1), (d,1E-3), (x, distFromBeamToRef), (y,distToBoard), (z,mDistDillArr[1])]) * 1E9
print("\nErr:", dqNum, "nm\n")

#outputting final ans stuffs
print("Wavelength", meanWavelength, "±", dqNum, "nm")
