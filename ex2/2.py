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
        :rtype : none
        """

        # Berechnen der gleichverteilten Zufallsvariablen
        random = self.random(4)
        # tau zwischen 50/SqrtS und 100/SqrtS
        self.tau = 50**2/sqrtS**2 * random[0] + 50**2/sqrtS**2
        # abs(y) kleiner log(tau)/2
        self.y = np.log(self.tau)*random[1] - np.log(self.tau)/2.0
        # cos Theta zwischen -1 und 1
        self.cos_theta = 2*random[2] - 1
        # phi zwischen 0 und 2 pi
        self.phi = 2 * np.pi * random[3]

        # Berechnen der restlichen Parameter
        self.x2 = np.sqrt(self.tau/np.exp(2*self.y))
        self.x1 = self.tau / self.x2
        self.sqrt_s_hut = np.sqrt(self.tau) * sqrtS

        # TODO: Berechnung von pe und pn!
        # Aus y = 1/2 log(1 + beta / 1 - beta)
        beta = (np.exp(2 * self.y) - 1)/(np.exp(2 * self.y) + 1)
        gamma = 1/np.sqrt(1 - beta**2)
        sin_theta = np.sqrt(1 - self.cos_theta**2)
        # _star bezeichnet die Größen im Schwerpunktsystem
        p1e_star = self.sqrt_s_hut/2 * sin_theta * np.cos(self.phi)
        p1n_star = -p1e_star
        p2e_star = self.sqrt_s_hut/2 * sin_theta * np.sin(self.phi)
        p2n_star = -p2e_star
        p3e_star = self.sqrt_s_hut/2 * self.cos_theta
        p3n_star = -p3e_star
        p0n_star = p0e_star = np.sqrt(p1e_star**2 + p2e_star**2 + p3e_star**2)

        # Boost ins Laborsystem
        self.p0e = gamma * (p0e_star - beta * p3e_star)
        self.p0n = gamma * (p0n_star - beta * p3n_star)
        self.p1e = p1e_star
        self.p1n = p1n_star
        self.p2e = p2e_star
        self.p2n = p2n_star
        self.p3e = gamma*(p3e_star - beta * p0e_star)
        self.p3n = gamma*(p3n_star - beta * p0n_star)

    def get_sum_m(self):
        return 2*GF**2*MassW**8*(1+self.cos_theta)**2 / ((self.sqrt_s_hut**2-MassW**2)**2 + MassW**2*GammaW**2)

    def get_dsigma(self):
        """ Berechnet dSigma für das gegebene Event
        """
        return self._f(self.x1) * self._f(self.x2) * 1/(32*np.pi**2) * \
               self.get_sum_m() * 1/(2*self.sqrt_s_hut**2) * 1/(2*np.pi)**2 * 1/(4*self.p0e*self.p0n)


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
        #dSigma = x.get_dsigma()
        #if Event.random() < dSigma/maxDSigma:
        #    sigma += dSigma
        #    counter += 1.0
        #    events.append(x)
    #print(sigma/counter, sigma/counter*0.3894/1000.0, counter)
    
    #events = np.load('eventsData.npy')
    #plot(events)
    
    #np.save("eventsData", np.array(events))