#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def t1(t, W):
	random = np.random.rand()
	return random**(1/W) * t

def returnPartonList():
	# Startwert
	t = 1
	# Kritische Untere Grenze
	tc = t/1000.0
	# W-Parameter
	W = 1

	partonList = list()

	while t > tc:
		t = t1(t, W)
		partonList.append(t)

	return partonList

lenList = list()
energyList = list()

for i in xrange(10000):
	tmp = returnPartonList()
	energyList.extend(tmp)
	lenList.append(len(tmp))

plt.hist(lenList)
plt.show()

plt.hist(energyList, bins=40)
plt.show()
