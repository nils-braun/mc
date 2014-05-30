__author__ = 'mapr'

import Particle
import numpy as np
from vector import LorentzVector
import pylab as plt


def kaellem(a, b, c):
    return a**2 + b**2 + c**2 - 2*a*b - 2*a*c - 2*b*c


def generator(n, s):
    m1squared = np.random.random() * s
    m2squared = np.random.random() * (s + m1squared - np.sqrt(4 * s * m1squared))
    cosTheta = np.random.random() * 2 - 1
    phi = np.random.random() * 2 * np.pi

    p = np.sqrt(1./(4*s) * kaellem(s, m1squared, m2squared))

    e1 = (s + m1squared - m2squared) / (2*np.sqrt(s))
    e2 = (s + m2squared - m1squared) / (2*np.sqrt(s))
    p1 = LorentzVector(p * cosTheta * np.cos(phi), p * cosTheta * np.sin(phi), p * np.sqrt(1-cosTheta**2), e1)
    p2 = LorentzVector(-p * cosTheta * np.cos(phi), -p * cosTheta * np.sin(phi), -p * np.sqrt(1-cosTheta**2), e2)

    factor = (s + m1squared - np.sqrt(4 * s * m1squared)) / 8. / np.pi * np.sqrt(kaellem(s, m1squared, m2squared))

    momenta = []
    if n == 2:
        momenta.append(p1)
        momenta.append(p2)
    else:
        momenta.append(p1)
        factor2, momenta2 = generator(n-1, m2squared)
        factor *= factor2
        for momentum in momenta2:
            vabssquared = 1 - m2squared / e2**2
            vx = p2.x / np.sqrt(m2squared) * np.sqrt(1-vabssquared)
            vy = p2.y / np.sqrt(m2squared) * np.sqrt(1-vabssquared)
            vz = p2.z / np.sqrt(m2squared) * np.sqrt(1-vabssquared)
            momentum.Boost(vx, vy, vz)

        momenta.extend(momenta2)

    return factor, momenta



def main():
    s = 90**2  # GeV**2
    factor, momenta = generator(2, s)

    energyConservation = 0
    momentumConservation = [0, 0, 0]
    finalStateMass = []
    for lorentzvector in momenta:
        energyConservation += lorentzvector.t
        momentumConservation[0] += lorentzvector.x
        momentumConservation[1] += lorentzvector.y
        momentumConservation[2] += lorentzvector.z
        finalStateMass.append(lorentzvector.M2())

    return momenta[0].M2(), momenta[1].M2()


if __name__ == '__main__':
    m1List = list()
    m2List = list()
    for i in xrange(0, int(10e3)):
        m1, m2 = main()
        m1List.append(m1)
        m2List.append(m2)

    plt.hist(m1List)
    plt.show()
    plt.hist(m2List)
    plt.show()