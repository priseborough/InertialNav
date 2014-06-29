#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('RngFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'IRNG', 'VRNG'])

figTas = plt.figure()

ax1 = figTas.add_subplot(111)
SRNG = pow(data['VRNG'],0.5)
ax1.set_title("Range Finder Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('range (m)')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['IRNG'], color='b', label='Range')
ax1.plot(data['time'], SRNG, color='r')
ax1.plot(data['time'], -SRNG, color='r')

plt.show()
