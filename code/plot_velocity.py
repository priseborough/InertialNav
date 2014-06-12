#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

adata = np.genfromtxt('NTUN.txt', skip_header=1,
	skip_footer=1, names=['time', 'time2', 'empty1', 'empty2', 'empty3', 'empty4', 'empty5', 'airspeed', 'empty6'], autostrip=True, filling_values=0)

data = np.genfromtxt('EulDataOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'roll', 'roll_onb', 'pitch', 'pitch_onb', 'yaw', 'yaw_onb', 'empty1', 'empty2'])

vdata = np.genfromtxt('ValidationOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'roll_int', 'pitch_int', 'yaw_int', 'north', 'east', 'down'])


fig = plt.figure()

ax1 = fig.add_subplot(311)

ax1.set_title("Airspeed / Pitch Estimated / Pitch Onboard")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Airspeed (m/s)')
# ax1.set_ylim([-0.0025,0.0025])
ax1.plot(adata['time'], adata['airspeed'], color='r', label='Pn')

ax2 = fig.add_subplot(312)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Pitch estimated (deg)')
# ax2.set_ylim([-0.0025,0.0025])
ax2.plot(data['time'], data['pitch'], color='g', label='pitch')
ax2.plot(data['time'], data['pitch_onb'], color='r', label='pitch_onb')

ax3 = fig.add_subplot(313)
 
ax3.set_xlabel('time (s)')
ax3.set_ylabel('Pitch onboard (deg)')
# ax3.set_ylim([-0.0025,0.0025])
#ax3.plot(data['time'], data['pitch'], color='g', label='pitch')
#ax3.plot(data['time'], data['pitch_onb'], color='b', label='pitch_onb')
ax3.plot(vdata['time'], vdata['pitch_int'], color='r', label='pitch integral')

plt.show()