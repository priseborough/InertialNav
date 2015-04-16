#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('GlitchOffsetOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'pNorth', 'pEast', 'vNorth', 'vEast'])

fig = plt.figure()

ax1 = fig.add_subplot(211)

ax1.set_title("GPS glitch offset")    
#ax1.set_xlabel('time (s)')
ax1.set_ylabel('position offset (m)')
ax1.plot(data['time'], data['pNorth'], color='b', label='pNorth')
ax1.plot(data['time'], data['pEast'], color='r', label='pEast')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)


ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('velocity offset (m/sec)')
ax2.plot(data['time'], data['vNorth'], color='b', label='vNorth')
ax2.plot(data['time'], data['vEast'], color='r', label='vEast')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc=2)

plt.show()
