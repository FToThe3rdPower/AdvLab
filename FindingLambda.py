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
distToBoard = 8.02 #±.02m
distFromBeamToRef = 0.85 #±.01m
angleOfRef = np.arctan(distFromBeamToRef/distToBoard)/2
diffGratingSpacing = 1E-3 #mm

#arrays of data
mDistDillArr = np.asarray([distFromBeamToRef, distFromBeamToRef + .07, distFromBeamToRef + .14, distFromBeamToRef + .2, distFromBeamToRef + .25, distFromBeamToRef + .30, distFromBeamToRef + .35,
               distFromBeamToRef + .39, distFromBeamToRef + .43, distFromBeamToRef + .47, distFromBeamToRef + .51, distFromBeamToRef + .54, distFromBeamToRef + .58])


#array for ans storage
mLambDillArr = np.zeros_like(mDistDillArr)

#looping over m
for i in range(len(mDistDillArr)):
     #feeding the func "1' for the 0th max so it doesn't divide by 0
     if i == 0:
          mLambDillArr[0] = Lambdanator((1), diffGratingSpacing, distToBoard, distFromBeamToRef, mDistDillArr[1])
     else:
          mLambDillArr[i] = Lambdanator((i), diffGratingSpacing, distToBoard, distFromBeamToRef, mDistDillArr[i])
     print("\n M"+str(i), "dist. from M0:", mDistDillArr[i], "\n calculated lambda:", mLambDillArr[i])


#checking results
print("\n\n\n\n", np.mean(mLambDillArr))