#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('OptFlowFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'ILOSX', 'VLOSX', 'ILOSY', 'VLOSY'])

figOptFlow = plt.figure()

ax1 = figOptFlow.add_subplot(211)
SLOSX = pow(data['VLOSX'],0.5)
ax1.set_title("Optical Flow Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X LOS rate (rad/sec)')
ax1.plot(data['time'], data['ILOSX'], color='b')
ax1.plot(data['time'], SLOSX, color='r')
ax1.plot(data['time'], -SLOSX, color='r')

ax2 = figOptFlow.add_subplot(212)
SLOSY = pow(data['VLOSY'],0.5)
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Y LOS rate (rad/sec)')
ax2.plot(data['time'], data['ILOSY'], color='b')
ax2.plot(data['time'], SLOSY, color='r')
ax2.plot(data['time'], -SLOSY, color='r')

plt.show()
