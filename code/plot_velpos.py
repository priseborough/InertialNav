#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('VelPosFuse.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'IVX', 'VVX', 'IVY', 'VVY', 'IVZ', 'VVZ', 'IPX', 'VPX', 'IPY', 'VPY', 'IPZ', 'VPZ'])

figVel = plt.figure()

ax1 = figVel.add_subplot(311)
SVX = pow(data['VVX'],0.5)
#ax1.set_title("Velocity Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['IVX'], color='b', label='Velocity Innovations')
ax1.plot(data['time'], SVX, color='r', label='sigma X')
ax1.plot(data['time'], -SVX, color='r')
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

ax2 = figVel.add_subplot(312)
SVY = pow(data['VVY'],0.5)
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax2.set_ylim([-0.0025,0.0025])
ax2.plot(data['time'], data['IVY'], color='b', label='Vel Y')
ax2.plot(data['time'], SVY, color='r', label='sigma Y')
ax2.plot(data['time'], -SVY, color='r')

ax3 = figVel.add_subplot(313)
SVZ = pow(data['VVZ'],0.5)
ax3.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax3.set_ylim([-0.0025,0.0025])
ax3.plot(data['time'], data['IVZ'], color='b', label='Vel Z')
ax3.plot(data['time'], SVZ, color='r', label='sigma Z')
ax3.plot(data['time'], -SVZ, color='r')
ax3.set_ylabel('Z')
plt.savefig('velInnov.png', bbox_inches='tight')

figPos = plt.figure()

ax1 = figPos.add_subplot(311)
SPX = pow(data['VPX'],0.5)
#ax1.set_title("Position Innovations")
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X')
#ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['IPX'], color='b', label='Position Innovations')
ax1.plot(data['time'], SPX, color='r', label='sigma')
ax1.plot(data['time'], -SPX, color='r')
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

ax2 = figPos.add_subplot(312)
SPY = pow(data['VPY'],0.5)
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax2.set_ylim([-0.0025,0.0025])
ax2.plot(data['time'], data['IPY'], color='b', label='Pos Y')
ax2.plot(data['time'], SPY, color='r')
ax2.plot(data['time'], -SPY, color='r')

ax3 = figPos.add_subplot(313)
SPZ = pow(data['VPZ'],0.5)
ax3.set_xlabel('time (s)')
ax2.set_ylabel('Y')
#ax3.set_ylim([-0.0025,0.0025])
ax3.plot(data['time'], data['IPZ'], color='b', label='Pos Z')
ax3.plot(data['time'], SPZ, color='r')
ax3.plot(data['time'], -SPZ, color='r')
ax3.set_ylabel('Z')
plt.savefig('posInnov.png', bbox_inches='tight')

plt.show()
