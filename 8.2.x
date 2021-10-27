#Err analysis by Taylor 8.x, Chi^2 fitting homework


#imports
import numpy as np
import matplotlib.pyplot as plt


#data (arrays must be same size, but they were manually entered for each problem)
x = np.arange(200, 1000, 100)
y = np.asarray([5.1, 5.5, 5.9, 6.8, 7.4, 7.5, 8.6, 9.4])
n = len(x)


#funcAlicious
def aSum(arrayToSum):
     return np.sum(arrayToSum)

def aSqSum(arrayToSum):
     return np.sum(np.power(arrayToSum, 2))

def aTimesbSum(aArr, bArr):
     return np.sum(aArr*bArr)

def lineFunc(m, c, xArr):
     y = m*xArr+c
     return y


delta = (n * aSqSum(x)) - (aSum(x))**2

a = ((aSqSum(x) * aSum(y)) - (aSum(x) * aTimesbSum(x, y))) / delta

b = (n * aTimesbSum(x,y) - (aSum(x) * aSum(y))) / delta


#err stuff
siggy = np.ones(n, dtype=float)
for i in range(n):
     siggy[i] = (y[i]-a-b*x[i])**2

stdDevY = np.sqrt((1/(n-2)) * np.sum(siggy))


#plotting
plt.close("all")
plt.figure('#2')
plt.errorbar(x, y, yerr=stdDevY, xerr=None, fmt=".k", ecolor="coral", label="data")
plt.plot(x, lineFunc(b, a, x), "g", label="χ^2 linear curve fit") #yes, this is the correct order. b=slope, a=intercept
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
