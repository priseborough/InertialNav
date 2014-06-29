#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('MagFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'IMX', 'VMX', 'IMY', 'VMY', 'IMZ', 'VMZ'])

fig = plt.figure()

ax1 = fig.add_subplot(311)
SMX = pow(data['VMX'],0.5)
ax1.set_title("Magnetometer Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['IMX'], color='b', label='Mag X')
ax1.plot(data['time'], SMX, color='r')
ax1.plot(data['time'], -SMX, color='r')

ax2 = fig.add_subplot(312)
SMY = pow(data['VMY'],0.5)
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax2.set_ylim([-0.0025,0.0025])
ax2.plot(data['time'], data['IMY'], color='b', label='Mag Y')
ax2.plot(data['time'], SMY, color='r')
ax2.plot(data['time'], -SMY, color='r')

ax3 = fig.add_subplot(313)
SMZ = pow(data['VMZ'],0.5)
ax3.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax3.set_ylim([-0.0025,0.0025])
ax3.plot(data['time'], data['IMZ'], color='b', label='Mag Z')
ax3.plot(data['time'], SMZ, color='r')
ax3.plot(data['time'], -SMZ, color='r')

plt.show()
