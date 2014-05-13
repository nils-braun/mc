#!/usr/bin/env python3
# Exercise 2 from MC-Ereignisgeneratoren
# May 2014

import pylab
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

G_F = 10e-5
M_W = 80
Gamma_W = 2

def calcSumM(Theta, sHat):
    return (2*G_F * M_W**8 * (1 + np.cos(Theta))**2)/((sHat - M_W)**2 + M_W**2 * Gamma_W**2)


def pdf(x, sHat):
    return x**(0.2 - 0.3*np.log(sHat)) * (1 - x)**0.1


if __name__ == "__main__":
    xValues = np.arange(0.1, 1.0, 0.01)
    sHatValues = np.arange(50, 100, 1)
    thetaValues = np.arange(-1.0, 1.0, 0.1)

    xValuesGrid, sHatValuesGrid = np.meshgrid(xValues, sHatValues)
    pdfValues = pdf(xValuesGrid, sHatValuesGrid)

    thetaValuesGrid, sHatValuesGrid = np.meshgrid(thetaValues, sHatValues)
    sumMValues = calcSumM(thetaValuesGrid, sHatValuesGrid)

    fig = pylab.figure()
    ax = fig.gca(projection='3d')

    ax.plot_surface(thetaValuesGrid, sHatValuesGrid, sumMValues)

    pylab.show()

