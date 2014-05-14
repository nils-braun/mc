# -*- coding: utf-8 -*-
"""
Created on Fri May 09 10:07:35 2014

@author: s0cke, nilpferd1991
"""

import numpy as np
import matplotlib.pyplot as plt

# Globale Variablendefintionen
sqrtS = 2000
MassW = 80
GammaW = 2
GF = 10**(-5)


class Event:
    """ Klasse zur Speicherung und Erzeugung eines Events """

    def __init__(self):
        self.p0e, self.p1e, self.p2e, self.p3e = 0, 0, 0, 0
        self.p0n, self.p1n, self.p2n, self.p3n = 0, 0, 0, 0
        self.x1, self.x2, self.tau, self.y, self.sqrt_s_hut = 0, 0, 0, 0, 0
        self.cos_theta, self.phi = 0, 0

        self.get_event()

    @staticmethod
    def random(size=-1):
        """ Gibt eine Anzahl an zufälligen Elementen zurück
        :rtype : np.array
        :param size: Die Anzahl der zufälligen Elemente
        TODO: Eigenen Zufallszahlengenerator implementieren
        """

        if size == -1:
            return np.random.random()
        else:
            return np.random.random(size, )

    def _f(self, x):
        """ Berechnet die Parton-pdf

        :param x: der Partonenbruchteil
        :rtype : double
        """
        return x**(0.2-0.3*np.log(self.sqrt_s_hut**2))*(1-x)**0.1

    def get_event(self):
        """ Gibt ein zufällig gewürfeltes Event zurück
        Ändert:
        p0e, p1e, p2e, p3e, p0n, p1n, p2n, p3n, x1, x2, Tau, y, sqrtSHut

        TODO: Andere Variablen verwenden!
        :rtype : none
        """

        # Berechnen der gleichverteilten Zufallsvariablen
        random = self.random(4)
        # tau zwischen 50/SqrtS und 100/SqrtS
        self.tau = 50/sqrtS * random[0] + 50/sqrtS
        # y zwischen -5 und 5
        self.y = 10*random[1] - 5
        # cos Theta zwischen -1 und 1
        self.cos_theta = 2*random[2] - 1
        # phi zwischen 0 und 2 pi
        self.phi = 2 * np.pi * random[3]

        # Berechnen der restlichen Parameter

        

    def get_sum_m(self):
        return 2*GF**2*MassW**8*(1+self.cos_theta)**2 / ((self.sqrt_s_hut**2-MassW**2)**2 + MassW**2*GammaW**2)

    def get_dsigma(self):
        """ Berechnet dSigma für das gegebene Event
        TODO: Jakobideterminante implementieren!
        """
        return self._f(self.x1) * self._f(self.x2) * \
               self.get_sum_m() * 1/(2*self.sqrt_s_hut**2) * 1/(2*np.pi)**2 * 1/(4*self.p0e*self.p0n)
        
    # TODO: Sollte in get_event implementiert werdn!
    #def cuts(p0e, p1e, p2e, p3e, p0n, p1n, p2n, p3n, x1, x2, Tau, y, sqrtSHut):
    #    if ((sqrtSHut < 50) or (sqrtSHut > 100)):
    #        return False
    #    if np.abs(y) >= 5:
    #        return False
    #    return True


def get_max_dsigma():
    max_dsigma = -1
    for i in range(10000):
        event = Event()
        max_dsigma = np.maximum(event.get_dsigma(), max_dsigma)
    return max_dsigma
    
# TODO: implementieren
#def plot(events):
    #plt.hist(np.sqrt(events[:,1]**2 + events[:,2]**2 + events[:,3]**2), bins=40)
    #plt.hist(np.sqrt(events[:,5]**2 + events[:,6]**2), bins=40)
        

if __name__ == '__main__':
    #maxDSigma = get_max_dsigma()
    maxDSigma = 4.16315676394e-06
    sigma = 0
    counter = 0
    events = list()
    for i in range(10000):
        x = Event()
        dSigma = x.get_dsigma()
        if Event.random() < dSigma/maxDSigma:
            sigma += dSigma
            counter += 1.0
            events.append(x)
    print(sigma/counter, sigma/counter*0.3894/1000.0, counter)
    
    #events = np.load('eventsData.npy')
    #plot(events)
    
    #np.save("eventsData", np.array(events))