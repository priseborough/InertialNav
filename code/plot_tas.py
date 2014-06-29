#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('TasFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'ITAS', 'VTAS'])

figTas = plt.figure()

ax1 = figTas.add_subplot(111)
STAS = pow(data['VTAS'],0.5)
ax1.set_title("True Airspeed Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['ITAS'], color='b', label='Vtas')
ax1.plot(data['time'], STAS, color='r')
ax1.plot(data['time'], -STAS, color='r')

plt.show()
