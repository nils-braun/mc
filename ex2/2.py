# -*- coding: utf-8 -*-
"""
Created on Fri May 09 10:07:35 2014

@author: s0cke
"""

import numpy as np
import matplotlib.pyplot as plt

sqrtS = 2000
MassW = 80
GammaW = 2
GF = 10**(-5)

#Partonverteilung
def f(x, s):
    return x**(0.2-0.3*np.log(s))*(1-x)**0.1
    
def plotF(log=True):
    t = np.arange(0., 1., 0.001)
    plt.plot(t, f(t, 2000))
    if log == True:
        plt.yscale('log')
        

def getEvent():
    #Tau,y,p2e,p3e = np.random.random(4,)
    #y = y*10 - 5
    #x2 = np.sqrt(Tau/np.exp(2*np.abs(y)))
    #x1 = Tau/x2

    x1,x2,p2e,p3e = np.random.random(4,)
    y = 0.5*np.log(x1/x2)
    Tau = x1*x2
    sqrtSHut = np.sqrt(Tau)*sqrtS
    p2e = p2e*sqrtSHut - sqrtSHut/2.0
    p3e = p3e*sqrtSHut - sqrtSHut/2.0
    p2n = -p2e
    p3n = sqrtS/2. * (x1-x2)-p3e
    
    if (sqrtS**2*(-p2e**2*(x1+x2)**2 + x1*(-2.0*p3e+sqrtS*x1)*x2*(2*p3e+sqrtS*x2)) < 0):
        return getEvent()
    else:
        p1e = np.random.choice([1,-1]) * np.sqrt(sqrtS**2*(-p2e**2*(x1+x2)**2 + x1*(-2.0*p3e+sqrtS*x1)*x2*(2*p3e+sqrtS*x2))) / np.sqrt(sqrtS**2*(x1+x2)**2)
        p1n = -p1e
        
        p0e = np.sqrt(p1e**2+p2e**2+p3e**2)
        p0n = np.sqrt(p1n**2+p2n**2+p3n**2)
        #sqrtSHut = np.sqrt(2*(p0e*p0n-p1e*p1n-p2e*p2n-p3e*p3n))
        
        return p0e, p1e, p2e, p3e, p0n, p1n, p2n, p3n, x1, x2, Tau, y, sqrtSHut
        
        
def cuts(p0e, p1e, p2e, p3e, p0n, p1n, p2n, p3n, x1, x2, Tau, y, sqrtSHut):
    if ((sqrtSHut < 50) or (sqrtSHut > 100)):
        return False
    if np.abs(y) >= 5:
        return False
    return True
    

def getDSigma():
    x = getEvent()
    while (not cuts(*x)):
        x = getEvent()
    #print x
    cosTheta = x[3] / np.sqrt(x[1]**2 + x[2]**2 + x[3]**2)
    sum = 2*GF**2*MassW**8*(1+cosTheta)**2 / ((x[12]-MassW**2)**2 + MassW**2*GammaW**2)
    #print sum
    
    dSigma = f(x[8],x[12]) * f(x[9],x[12]) * sum * 1/(2*x[12]) * 1/(2*np.pi)**2 * 1/(4*x[0]*x[4])
    #print dSigma
    return dSigma, x
    

def getMaxDSigma():
    maxDSigma = -1
    for i in range(10000):
        maxDSigma = np.maximum(getDSigma()[0], maxDSigma)
    print maxDSigma
    
    
def plot(events):
    #plt.hist(np.sqrt(events[:,1]**2 + events[:,2]**2 + events[:,3]**2), bins=40)
    plt.hist(np.sqrt(events[:,5]**2 + events[:,6]**2), bins=40)
        

if __name__ == '__main__':
    maxDSigma = 4.16315676394e-06
    sigma = 0
    counter = 0
    #events = list()
    #for i in range(10000):
    #    dSigma,x = getDSigma()
    #    if np.random.random() < dSigma/maxDSigma:
    #        sigma += dSigma
    #        counter += 1.0
    #        events.append(np.array(x))
    #print sigma/counter
    #print sigma/counter*0.3894/1000.0
    #print counter
    
    events = np.load('eventsData.npy')
    plot(events)
    
    #np.save("eventsData", np.array(events))