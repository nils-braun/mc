__author__ = 'mapr'

import numpy as np


class LorentzVector:
    def __init__(self, e, px, py, pz):
    #def __init__(self, vector) vector[e, px, py, pz]
        self.e = e
        self.px = px
        self.py = py
        self.pz = pz

    def invariant_mass(self):
        return self.e**2 - self.px**2 - self.py**2 - self.pz**2

    def boost(self, _, vx, vy, vz):
        vsquared = vx**2 + vy**2 + vz**2
        gamma = 1./np.sqrt(1-vsquared)
