#!/usr/bin/env python 3.6.8 Anaconda 64-bit
# -*- coding: utf-8 -*-
#Created on Tue Sep  7 11:53:07 2021
#@author: frankgrijalva


#Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spy


#Funcs for later use
def Curvy(x, a, b, background):
     return a*np.exp(b*x) + background

def XtraCurvy(x, a, b, c, d, background):
     return a*np.exp(b*x) + c*np.exp(d*x) + background

def XXtraCurvy(x, a, b, c, d, f, g, background):
     return a*np.exp(b*x) + c*np.exp(d*x) + f*np.exp(g*x) + background

def RMSE(resi):
     return np.sqrt(np.mean((resi)**2))

def DataCFPlotty(XData, YData, curveFitFunc, fitParams, figTitle, title, xLabel, yLabel, dataStyleStr, curveStyleStr, plotStyle):
     plt.figure(figTitle)
     plt.title(title)

     if (plotStyle == "semilogy"):
          plt.semilogy(XData, YData, dataStyleStr, label="Data")
          plt.semilogy(XData, curveFitFunc(XData, *fitParams), curveStyleStr, label="Curve Fit")

     else:
          plt.plot(XData, YData, dataStyleStr, label="Data")
          plt.plot(XData, curveFitFunc(XData, *fitParams), curveStyleStr, label="Curve Fit")

     plt.xlabel(xLabel)
     plt.ylabel(yLabel)
     plt.legend()
     return plt.show()


#Data handling
dataArray = np.loadtxt("/Users/frankgrijalva/Desktop/School/427/Labs/First 3/Fake rad src/fakeRadSrc.txt", dtype=float, skiprows=3)
##second data collection run because the first one might be bad
dataArray1 = np.loadtxt("/Users/frankgrijalva/Desktop/School/427/Labs/First 3/Fake rad src/fakeRadSrc.txt", dtype=float, skiprows=3)
##third data capture because the first 2 were kinda bad
estebanData = np.loadtxt("/Users/frankgrijalva/Desktop/School/427/Labs/First 3/Fake rad src/half_life_data_1000a.txt", dtype=float, skiprows=3)

#counter normalization to 0 and sample interval correction
xData = (dataArray[:,0] - dataArray[0,0])/5000
x1Data = (dataArray1[:,0] - dataArray1[0,0])/5000
xEstebansData1k = (estebanData[:,0] - estebanData[0,0])/1000

#naming the slices made it easier to read
smpl1cap1 = dataArray[:,1]
smpl1cap2 = dataArray1[:,1]
# smpl2cap1 = dataArray[:80,2]
smpl2capEst1k = estebanData[:,2]
smpl3cap1 = dataArray[:,3]
smpl3capEst1k = estebanData[:,3]
smpl4cap1 = dataArray[:,4]
smpl4capEst1k = estebanData[:,4]



#curve fitting
doubleFit, dubCovMat = spy.curve_fit(XtraCurvy, xData, smpl3cap1, p0=(1,-1,1,-1,0))
fit4, covMat4 = spy.curve_fit(XXtraCurvy, xData, dataArray[:,4], p0=(1,-1,1,-1,1,-1,0))

#second dataset fits for comparison
fit11, covMat11 = spy.curve_fit(Curvy, x1Data, dataArray1[:,1], p0=(1,-1,0))
fit44, covMat44 = spy.curve_fit(XXtraCurvy, x1Data, dataArray1[:,4], p0=(1,-1,1,-1,1,-1,0))


#esteban dataset fit for comparison
fit2E, covMat2E = spy.curve_fit(Curvy, xEstebansData1k, smpl2capEst1k, p0=(1,-1,0), method="trf")
fit4E, covMat4E = spy.curve_fit(XtraCurvy, xEstebansData1k, smpl4capEst1k, p0=(1,-1,1,-1,0))


#calculating the parameters we're after from the fit functions
#Calculating the decay constant, tau
tau1 = 1/fit11[1]
tau2E1k = 1/fit2E[1]
tau3a = 1/doubleFit[1]
tau3b = 1/doubleFit[3]
tau4a = 1/fit4[1]
tau4b = 1/fit4[3]
tau4c = 1/fit4[5]

#Calculating half-life, correcting for the negative in the variable
smpl1HL = (-np.log(2) * tau1)
smpl2E1kHL = (-np.log(2) * tau2E1k)
smpl3HLa = (-np.log(2) * tau3a)
smpl3HLb = (-np.log(2) * tau3b)
smpl4HLa = (-np.log(2) * tau4a)
smpl4HLb = (-np.log(2) * tau4b)
smpl4HLc = (-np.log(2) * tau4c)


#calculating initial abundance
A0_smpl1 = fit11[0]
A0_smpl2 = fit2E[0]
A0_smpl3a = doubleFit[0]
A0_smpl3b = doubleFit[2]
A0_smpl4a = fit4[0]
A0_smpl4b = fit4[2]
A0_smpl4c = fit4[4]


#calculating the background radiation in each sample
bckgrnd1 = fit11[2]
bckgrnd2 = fit2E[2]
bckgrnd3 = doubleFit[4]
bckgrnd4 = fit4[6]


#calculating the errors
smpl1StdDev = np.sqrt(np.diag(covMat11))
smpl2StdDev = np.sqrt(np.diag(covMat2E))
smpl3StdDev = np.sqrt(np.diag(dubCovMat))
smpl4StdDev = np.sqrt(np.diag(covMat4))

#naming them in "english"
smpl1A0Err = np.round(smpl1StdDev[0], 2)
smpl1TauErr = np.round(smpl1StdDev[1], 3)
smpl1BckErr = np.round(smpl1StdDev[2], 2)

smpl2A0err = np.round(smpl2StdDev[0], 2)
smpl2TauErr = np.round(smpl2StdDev[1], 3)
smpl2BckErr = np.round(smpl2StdDev[2], 2)

smpl3A0ErrA = np.round(smpl3StdDev[0], 2)
smpl3TauErrA = np.round(smpl3StdDev[1], 3)
smpl3A0ErrB = np.round(smpl3StdDev[2], 2)
smpl3TauErrB = np.round(smpl3StdDev[3], 3)
smpl3BckErr = np.round(smpl3StdDev[4], 2)

smpl4A0ErrA = np.round(smpl4StdDev[0], 2)
smpl4TauErrA = np.round(smpl4StdDev[1], 3)
smpl4A0ErrB = np.round(smpl4StdDev[2], 2)
smpl4TauErrB = np.round(smpl4StdDev[3], 3)
smpl4A0ErrC = np.round(smpl4StdDev[4], 2)
smpl4TauErrC = np.round(smpl4StdDev[5], 3)
smpl4BckErr = np.round(smpl4StdDev[6], 2)





#outputting the half life
print("\n\nElement 1 was initially present at",A0_smpl1, "±", smpl1A0Err, "units with \u03C41/2 =", smpl1HL,"±", smpl1TauErr, "(sec) with", bckgrnd1,"±", smpl1BckErr, "units of background radiation")
print("\n\nElement 2 from Esteban's data was initially present at",A0_smpl2, "±", smpl2A0err, "units with \u03C41/2 =", smpl2E1kHL,"±", smpl2TauErr, "(sec) with", bckgrnd1,"±", smpl1BckErr, "units of background radiation")
print("\n\nElement 3a was initially present at", A0_smpl3a, "±", smpl3A0ErrA, "units with \u03C41/2_a:", smpl3HLa, "±", smpl3TauErrA,"(sec)", "\nelement 3b was initially present at",A0_smpl3b, "±", smpl3A0ErrB, "units with \u03C41/2-b:", smpl3HLb, "±", smpl3TauErrB, "(sec) with", bckgrnd3,"±", smpl3BckErr, "units of background radiation")
print("\n\nElement 4a was initially present at", A0_smpl4a, "±", smpl4A0ErrA, "units with \u03C41/2_a:", smpl4HLa, "±", smpl4TauErrA,"(sec)", "\nelement 4b was initially present at",A0_smpl4b, "±", smpl4A0ErrB,"units with \u03C41/2_b:", smpl4HLb, "±", smpl4TauErrB,"(sec)", "\nelement 4c was initially present at",A0_smpl4c, "±", smpl4A0ErrC, "units with \u03C41/2_c:", smpl4HLc, "±", smpl4TauErrC,"(sec) with", bckgrnd4,"±", smpl4BckErr, "units of background radiation")


#residuals
res11 = (Curvy(x1Data, *fit11) - smpl1cap2) / Curvy(x1Data, *fit11)

res2E =(Curvy(xEstebansData1k, *fit2E) - smpl2capEst1k)/Curvy(xEstebansData1k, *fit2E)

res3 = (XtraCurvy(xData, *doubleFit) - smpl3cap1) / XtraCurvy(xData, *doubleFit)

res4 = (XXtraCurvy(xData, *fit4) - smpl4cap1) / XXtraCurvy(xData, *fit4)



#RMSE vals
print("\n\nRMSE 1:", RMSE(res11))
print("\nRMSE Esteban 1k ch 2", RMSE(res2E))
print("\nRMSE 3:", RMSE(res3))
print("\nRMSE 4:", RMSE(res4), "\n\n\n\n")



#Plotting
plt.close("all")

DataCFPlotty(x1Data, smpl1cap2, Curvy, fit11, "First sample", "Sample 1 data", "S", "Count", ".r", "-b", "reg")

DataCFPlotty(x1Data, smpl1cap2, Curvy, fit11, "First sample log", "Sample 1 data log", "S", "Count", ".r", "-b", "semilogy")

DataCFPlotty(xData, smpl3cap1, XtraCurvy, doubleFit, "Third sample" , "Sample 3 data", "S", "Count", ".k", "-g", 'reg')

DataCFPlotty(xData, smpl3cap1, XtraCurvy, doubleFit, "Third sample log" , "Sample 3 data log", "S", "Count", ".k", "-g", 'semilogy')

DataCFPlotty(xData, smpl4cap1, XXtraCurvy, fit4, "Fourth sample", "Sample 4 data", "S", "Count", ".c", "-k", "reg")


#Checking esteban's data
DataCFPlotty(xEstebansData1k, smpl2capEst1k, Curvy, fit2E, "Esteban 1s sample 2", "Esteban 1k sample 2", "S", "Count", ".k", "-r", "reg")

DataCFPlotty(xEstebansData1k, smpl2capEst1k, Curvy, fit2E, "Esteban 1s sample 2 log", "Esteban 1k sample 2 log", "S", "Count", ".k", "-r", "semilogy")




#residual figs
plt.figure("smpl1 resi")
plt.title("Sample 1 residuals")
plt.plot(res11, ".r", label='run 1')
plt.xlabel("Data point")
plt.ylabel("Difference\n(CF-Data)/(CF)/(CF)")
plt.legend()


plt.figure("smpl 2 estb 1k resi")
plt.title("Sample 2 Esteban 1s data residuals")
plt.plot(res2E, ".r", label="run 2 Estb 1k")
plt.xlabel("Data point")
plt.ylabel("Difference\n(CF-Data)/(CF)/(CF)")
plt.legend()

plt.figure("smpl 3 resi")
plt.title("Sample 3 residuals")
plt.plot(res3, ".r", label="run 1")
plt.xlabel("Data point")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()

plt.figure("smpl 4 residuals")
plt.title("Sample 4 residuals")
plt.plot(res4, ".r", label="3 exp fit")
plt.xlabel("Data point")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()


plt.show()