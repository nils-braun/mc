# -*- coding: utf-8 -*-
"""
Created on Fri May 09 10:07:35 2014

@author: s0cke, nilpferd1991
"""

# TODO: Unterschied von sigma bei crude/importance sampling

import numpy as np
import matplotlib.pyplot as plt
import os.path
from matplotlib.backends.backend_pdf import PdfPages

# Globale Variablendefintionen
sqrtS = 2000
MassW = 80
GammaW = 2
GF = 10**(-5)

# primitive Integration?
crude_mc = 0


class Event:
    """ Klasse zur Speicherung und Erzeugung eines Events """

    outputFile = "eventsData_" + str(sqrtS/1000) + "TeV_crude" + str(crude_mc) + ".npy"

    def __init__(self):
        self.p0e, self.p1e, self.p2e, self.p3e = 0, 0, 0, 0
        self.p0n, self.p1n, self.p2n, self.p3n = 0, 0, 0, 0
        self.x1, self.x2, self.tau, self.y, self.sqrt_s_hut = 0, 0, 0, 0, 0
        self.cos_theta, self.phi = 0, 0

        self.dsigma = np.nan

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
        # cos Theta zwischen -1 und 1
        self.cos_theta = 2*random[2] - 1
        # phi zwischen 0 und 2 pi
        self.phi = 2 * np.pi * random[3]

        if crude_mc == 0:
            # importance sampling mit Breit Wigner
            # rho zwischen 0 und 1
            self.rho = random[0]
            self.sqrt_s_hut = np.sqrt(MassW**2 + MassW*GammaW * np.tan(np.arctan((50**2-MassW**2)/MassW/GammaW)*self.rho+np.arctan((100**2-MassW**2)/MassW/GammaW)*(1-self.rho)))
            self.tau = self.sqrt_s_hut**2 / sqrtS**2
            
        else:
            # primitive Integration (crude mc)
            # tau zwischen 50/SqrtS und 100/SqrtS
            self.tau = (100**2 - 50**2)/sqrtS**2 * random[0] + 50**2/sqrtS**2
            self.sqrt_s_hut = np.sqrt(self.tau) * sqrtS
            
        # abs(y) kleiner log(tau)/2
        self.y = np.log(self.tau)*random[1] - np.log(self.tau)/2.0
        self.x2 = np.sqrt(self.tau/np.exp(2*self.y))
        self.x1 = self.tau / self.x2

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
        """ Gibt die Summe über M wie auf dem Blatt angeben zurück
        """
        return 2*GF**2*MassW**8*(1+ np.random.choice([1,-1])* self.cos_theta)**2 / ((self.sqrt_s_hut**2-MassW**2)**2 + MassW**2*GammaW**2)

    def get_dsigma(self):
        """ Berechnet dSigma für das gegebene Event, falls dies nocht schon vorher geschiehen ist und abgespeichert wurde
        """
        if self.dsigma is np.nan:
            if crude_mc == 0:
                # Jakobideterminante für importance sampling berücksichtigen
                self.dsigma = self._f(self.x1) * self._f(self.x2) * 1/(32*np.pi**2) * self.get_sum_m() * 1/(2*self.sqrt_s_hut**2) * 1/(2*np.pi)**2 * 1/(4*self.p0e*self.p0n) / (sqrtS**2 * MassW * GammaW) * ((self.sqrt_s_hut**2 - MassW**2)**2+MassW**2*GammaW**2)
            else:
                self.dsigma = self._f(self.x1) * self._f(self.x2) * 1/(32*np.pi**2) * self.get_sum_m() * 1/(2*self.sqrt_s_hut**2) * 1/(2*np.pi)**2 * 1/(4*self.p0e*self.p0n)
            return self.dsigma
        else:
            return self.dsigma

    def __getitem__(self, item):
        """ Eine Hilfsfunktion zur Ausgabe der einzelnen Variablen als Liste
        """
        return (self.p0e, self.p1e, self.p2e, self.p3e,  # 0, 1, 2, 3
                self.p0n, self.p1n, self.p2n, self.p3n,  # 4, 5, 6, 7
                self.x1, self.x2,                        # 8, 9
                self.tau, self.y,                        # 10, 11
                self.sqrt_s_hut,                         # 12
                self.cos_theta, self.phi,                # 13, 14
                self.get_dsigma())[item]                 # 15


def plot(plot_events):
    good_events = np.array([list(x) for x in plot_events])
    pp = PdfPages('multipage.pdf')
    # Transversalimpuls des Elektrons
    plt.hist(np.sqrt(good_events[:, 1]**2 + good_events[:, 2]**2), bins=40)
    plt.title('pT Elektron')
    plt.xlabel('pT in GeV')
    pp.savefig()
    plt.close()
    # Transversalimpuls des Neutrinos
    plt.hist(np.sqrt(good_events[:, 5]**2 + good_events[:, 6]**2), bins=40)
    plt.title('pT Neutrino')
    plt.xlabel('pT in GeV')
    pp.savefig()
    plt.close()
    # invariante Masse der Partonen
    plt.hist(good_events[:, 12], bins=40)
    plt.title('invariante Masse der Partonen')
    plt.xlabel('^s in GeV')
    pp.savefig()
    plt.close()
    # cos Theta
    plt.hist(good_events[:, 13], bins=40)
    plt.title('cos Theta')
    pp.savefig()
    plt.close()
    # Rapidität
    plt.hist(good_events[:, 11], bins=40)
    plt.title('Rapidität')
    pp.savefig()
    plt.close()
    # x1
    plt.hist(good_events[:, 8], bins=40)
    plt.title('x1')
    pp.savefig()
    plt.close()
    # x2
    plt.hist(good_events[:, 8], bins=40)
    plt.title('x2')
    pp.savefig()
    plt.close()
    pp.close()
        

if __name__ == '__main__':

    debug = 0

    events = list()
    goodEvents = list()

    if debug == 0:
        print("calculating new points...")
        maxDSigma = 0
        for i in range(100000):
            x = Event()
            dSigma = x.get_dsigma()
            events.append(x)
            maxDSigma = np.maximum(maxDSigma, dSigma)

        for x in events:
            if Event.random() <= x.get_dsigma()/maxDSigma:
                goodEvents.append(x)

    if os.path.isfile(Event.outputFile):
        print("reading old points...")
        loadedEvents = np.load(Event.outputFile)
        goodEvents.extend(loadedEvents)

    np.save(Event.outputFile, goodEvents)

    sigma = sum([x.get_dsigma()/len(goodEvents) for x in goodEvents])
    print("sigma =", sigma*0.3894/1000.0, ", count =", len(goodEvents))

    plot(goodEvents)