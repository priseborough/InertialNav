#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('RngFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'IRNG', 'VRNG', 'EST', 'DIST', 'BARO', 'ESTHAGL1', 'ESTHAGL2'])

figTas = plt.figure()

ax1 = figTas.add_subplot(111)
SRNG = pow(data['VRNG'],0.5)
ax1.set_title("Range Finder Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('range (m)')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['IRNG'], color='b', label='Range Innovations')
ax1.plot(data['time'], SRNG, color='r')
ax1.plot(data['time'], -SRNG, color='r')
ax1.plot(data['time'], data['EST'], color='g', label='Terrain offset estimate')
ax1.plot(data['time'], data['DIST'], color='m', label='Raw distance measurement')
ax1.plot(data['time'], data['BARO'], color='k', label='Barometric altitude')
ax1.plot(data['time'], data['ESTHAGL1'], color='c', label='HAGL estimate from 2-state EKF')
ax1.plot(data['time'], data['ESTHAGL2'], color='y', label='HAGL estimate from main EKF')
ax1.legend()

plt.show()
