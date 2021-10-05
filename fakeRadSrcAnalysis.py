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
dataArray1 = np.loadtxt("/Users/frankgrijalva/Desktop/School/427/Labs/First 3/Fake rad src/fakeRadSrc1.txt", dtype=float, skiprows=3)
##third data capture because the first 2 were kinda bad
estebanData = np.loadtxt("/Users/frankgrijalva/Desktop/School/427/Labs/First 3/Fake rad src/half_life_data_1000a.txt", dtype=float, skiprows=3)

#counter normalization to 0 and sample interval correction to sec
xData = (dataArray[:,0] - dataArray[0,0])/5000
x1Data = (dataArray1[:,0] - dataArray1[0,0])/5000
xEstebansData1k = (estebanData[:,0] - estebanData[0,0])/1000

#naming the slices made it easier to read
smpl1cap1 = dataArray[:,1]
smpl1cap2 = dataArray1[:,1]
# smpl2cap1 = dataArray[:80,2]
smpl2capEst1k = estebanData[:35,2]
smpl3cap1 = dataArray[:,3]
smpl3capEst1k = estebanData[:,3]
smpl4cap1 = dataArray[:,4]
smpl4capEst1k = estebanData[:,4]



#curve fitting
doubleFit, dubCovMat = spy.curve_fit(XtraCurvy, xEstebansData1k, smpl3capEst1k, p0=(1,-1,1,-1,0))
fit4, covMat4 = spy.curve_fit(XtraCurvy, xData, dataArray[:,4], p0=(750,-5,1060,-14,20))

#second dataset fits for comparison
fit11, covMat11 = spy.curve_fit(Curvy, xData, dataArray[:,1], p0=(1,-1,0))
fit44, covMat44 = spy.curve_fit(XXtraCurvy, x1Data, dataArray1[:,4], p0=(1,-1,2,-2,5,-5,30))


#esteban dataset fit for comparison
fit2E, covMat2E = spy.curve_fit(Curvy, xEstebansData1k[:35], smpl2capEst1k, p0=(1,-1,0), method="trf")
fit4E, covMat4E = spy.curve_fit(XtraCurvy, xEstebansData1k, smpl4capEst1k, p0=(1,-1,1,-1,0))


#calculating the parameters we're after from the fit functions
#Calculating the decay constant, tau
tau1 = np.round(1/fit11[1], 3)
tau2E1k = np.round(1/fit2E[1], 3)
tau3a =np.round( 1/doubleFit[1], 3)
tau3b = np.round(1/doubleFit[3], 3)
tau4a = np.round(1/fit4[1], 3)
tau4b = np.round(1/fit4[3], 3)


#Calculating half-life, correcting for the negative in the variable
smpl1HL = np.round((-np.log(2) * tau1), 3)
smpl2E1kHL = np.round((-np.log(2) * tau2E1k), 3)
smpl3HLa = np.round((-np.log(2) * tau3a), 3)
smpl3HLb = np.round((-np.log(2) * tau3b), 3)
smpl4HLa = np.round((-np.log(2) * tau4a), 3)
smpl4HLb = np.round((-np.log(2) * tau4b), 3)


#calculating initial abundance
A0_smpl1 = np.round(fit11[0], 2)
A0_smpl2 = np.round(fit2E[0], 2)
A0_smpl3a = np.round(doubleFit[0], 2)
A0_smpl3b = np.round(doubleFit[2], 2)
A0_smpl4a = np.round(fit4[0], 2)
A0_smpl4b = np.round(fit4[2], 2)


#calculating the background radiation in each sample
bckgrnd1 = fit11[2]
bckgrnd2 = fit2E[2]
bckgrnd3 = doubleFit[4]
bckgrnd4 = fit4[4]


#calculating the errors
smpl1StdDev = np.sqrt(np.diag(covMat11))
smpl2StdDev = np.sqrt(np.diag(covMat2E))
smpl3StdDev = np.sqrt(np.diag(dubCovMat))
smpl4StdDev = np.sqrt(np.diag(covMat4))

#naming them in "english"
smpl1A0Err = np.round(smpl1StdDev[0], 2)
smpl1TauErr = np.round(smpl1StdDev[1] * smpl1HL, 3)
smpl1BckErr = np.round(smpl1StdDev[2], 2)

smpl2A0err = np.round(smpl2StdDev[0], 2)
smpl2TauErr = np.round(smpl2StdDev[1] * smpl2E1kHL, 3)
smpl2BckErr = np.round(smpl2StdDev[2], 2)

smpl3A0ErrA = np.round(smpl3StdDev[0], 2)
smpl3TauErrA = np.round(smpl3StdDev[1] * smpl3HLa, 3)
smpl3A0ErrB = np.round(smpl3StdDev[2], 2)
smpl3TauErrB = np.round(smpl3StdDev[3] * smpl3HLb, 3)
smpl3BckErr = np.round(smpl3StdDev[4], 2)

smpl4A0ErrA = np.round(smpl4StdDev[0], 2)
smpl4TauErrA = np.round(smpl4StdDev[1] * smpl4HLa, 3)
smpl4A0ErrB = np.round(smpl4StdDev[2], 2)
smpl4TauErrB = np.round(smpl4StdDev[3] * smpl4HLb, 3)
smpl4BckErr = np.round(smpl4StdDev[4], 2)


#calculating relative abundance
relA0_smpl3a = str(np.round(((A0_smpl3a / (A0_smpl3a + A0_smpl3b)) * 100),2))+"±"+str(np.round(smpl3A0ErrA*(A0_smpl3a / (A0_smpl3a + A0_smpl3b))))
relA0_smpl3b = str(np.round(((A0_smpl3b / (A0_smpl4a + A0_smpl3b)) * 100),2))+"±"+str(np.round(smpl3A0ErrB*(A0_smpl3b / (A0_smpl3a + A0_smpl3b))))
relA0_smpl4a = str(np.round(((A0_smpl4a / (A0_smpl4a + A0_smpl4b)) * 100),2))+"±"+str(np.round(smpl4A0ErrA*(A0_smpl4a / (A0_smpl4a + A0_smpl4b))))
relA0_smpl4b = str(np.round(((A0_smpl4b / (A0_smpl4a + A0_smpl4b)) * 100),2))+"±"+str(np.round(smpl4A0ErrB*(A0_smpl4b / (A0_smpl4a + A0_smpl4b))))



#outputting the half life
print("\n\n\nElement 1 was initially",A0_smpl1, "±", smpl1A0Err, "units with a \u03C41/2 =", smpl1HL,"±", smpl1TauErr, "(sec) with", bckgrnd1,"±", smpl1BckErr, "units of background radiation")
print("\n\n\nElement 2 from Esteban's data was initially",A0_smpl2, "±", smpl2A0err, "units with a \u03C41/2 =", smpl2E1kHL,"±", smpl2TauErr, "(sec) with", bckgrnd2,"±", smpl1BckErr, "units of background radiation")
print("\n\n\nElement 3a was initially present at", A0_smpl3a, "±", smpl3A0ErrA, "units, or", relA0_smpl3a + "% of Sample 3","with a \u03C41/2_a:", smpl3HLa, "±", smpl3TauErrA,"(sec)", "\n\nElement 3b was initially present at",A0_smpl3b, "±", smpl3A0ErrB, "units", relA0_smpl3b + "% of Sample 3","with a \u03C41/2-b:", smpl3HLb, "±", smpl3TauErrB, "(sec) with", bckgrnd3,"±", smpl3BckErr, "units of background radiation")
print("\n\n\nElement 4a was initially", A0_smpl4a, "±", smpl4A0ErrA, "units or", relA0_smpl4a + "% of Sample 4 with a \u03C41/2_a:", smpl4HLa, "±", smpl4TauErrA,"(sec)", "\n\nElement 4b was initially present at",A0_smpl4b, "±", smpl4A0ErrB,"units, or", relA0_smpl4b + "% of Sample 4 with a \u03C41/2_b:", smpl4HLb, "±", smpl4TauErrB,"(sec) with", bckgrnd4,"±", smpl4BckErr, "units of background radiation")


#residuals
res11 = (Curvy(xData, *fit11) - smpl1cap1) / Curvy(xData, *fit11)

res2E =(Curvy(xEstebansData1k[:35], *fit2E) - smpl2capEst1k)/Curvy(xEstebansData1k[:35], *fit2E)

res3 = (XtraCurvy(xEstebansData1k, *doubleFit) - smpl3capEst1k) / XtraCurvy(xEstebansData1k, *doubleFit)

res4 = (XtraCurvy(xData, *fit4) - smpl4cap1) / XtraCurvy(xData, *fit4)



#RMSE vals
print("\n\n\nRMSE 1:", RMSE(res11))
print("\nRMSE Esteban 1k ch 2:", RMSE(res2E))
print("\nRMSE Esteban 1k ch 3:", RMSE(res3))
print("\nRMSE 4:", RMSE(res4), "\n\n\n\n")



#Plotting
plt.close("all")

DataCFPlotty(xData, smpl1cap1, Curvy, fit11, "First sample log", "Sample 1 data log", "Second", "Count", ".r", "-b", "semilogy")

DataCFPlotty(xEstebansData1k[:35], smpl2capEst1k, Curvy, fit2E, "Esteban 1k sample 2 log", "Sample 2 log", "Second", "Count", ".k", "-r", "semilogy")

DataCFPlotty(xEstebansData1k, smpl3capEst1k, XtraCurvy, doubleFit, "Third sample log" , "Sample 3 data log", "Second", "Count", ".k", "-g", 'semilogy')

DataCFPlotty(xData, smpl4cap1, XtraCurvy, fit4, "Fourth sample log", "Sample 4 data log", "Second", "Count", ".c", "-k", "semilogy")



#residual figs
plt.figure("smpl1 data / sec resi")
plt.title("Sample 1 residuals")
plt.plot(xData, res11, ".r", label='run 1')
plt.xlabel("Second")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()

plt.figure("smpl 2 estb 1k data / sec resi")
plt.title("Sample 2 residuals")
plt.plot(xEstebansData1k[:35], res2E, ".r", label="run 2 Estb 1k")
plt.xlabel("Second")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()

plt.figure("smpl 3 ET 1k data / sec resi")
plt.title("Sample 3 residuals")
plt.plot(xEstebansData1k, res3, ".r", label="run 1")
plt.xlabel("Second")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()

plt.figure("smpl 4 data / sec resi")
plt.title("Sample 4 residuals")
plt.plot(xData, res4, ".r", label="run 1")
plt.xlabel("Second")
plt.ylabel("Difference\n(CF-Data)/(CF)")
plt.legend()


plt.show()
