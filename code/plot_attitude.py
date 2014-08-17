#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math
data = np.genfromtxt('EulDataOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'roll', 'roll_onb', 'pitch', 'pitch_onb', 'yaw', 'yaw_onb', 'empty1', 'empty2'])
gdata = np.genfromtxt('GPSrawOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'lat', 'lon', 'alt', 'veln', 'vele', 'veld'])
mdata = np.genfromtxt('MAG.txt', skip_header=1,
	skip_footer=1, names=['time', 'time2', 'magx', 'magy', 'magz', 'offx', 'offy', 'offz'], autostrip=True, filling_values=0)

fig = plt.figure()

ax1 = fig.add_subplot(311)

ax1.set_title("Attitude estimate")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('angle (degrees)')

ax2 = fig.add_subplot(312)

ax2.set_title("Attitude estimate")    
ax2.set_xlabel('time (s)')
ax2.set_ylabel('angle (degrees)')

ax3 = fig.add_subplot(313)

ax3.set_title("Attitude estimate")    
ax3.set_xlabel('time (s)')
ax3.set_ylabel('angle (degrees)')

data['roll'] = np.multiply(data['roll'], 180 / math.pi)
data['pitch'] = np.multiply(data['pitch'], 180 / math.pi)
data['yaw'] = np.multiply(data['yaw'], 180 / math.pi)

data['roll_onb'] = np.multiply(data['roll_onb'], 180 / math.pi)
data['pitch_onb'] = np.multiply(data['pitch_onb'], 180 / math.pi)
data['yaw_onb'] = np.multiply(data['yaw_onb'], 180 / math.pi)

ax1.plot(data['time'], data['roll'], color='r', label='roll estimate')
ax2.plot(data['time'], data['pitch'], color='g', label='pitch estimate')
ax3.plot(data['time'], data['yaw'], color='b', label='yaw estimate')

ax1.plot(data['time'], data['roll_onb'], color='m', label='roll onboard')
ax2.plot(data['time'], data['pitch_onb'], color='c', label='pitch onboard')
ax3.plot(data['time'], data['yaw_onb'], color='k', label='yaw onboard')
# ax3.plot(gdata['time'], np.math.atan2(gdata['vele'], gdata['veln']))
# ax3.plot(mdata['time'], np.math.atan2(mdata['magy'], mdata['magx']))

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc=2)
handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=2)

plt.show()
