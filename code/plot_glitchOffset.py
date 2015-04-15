#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

data = np.genfromtxt('GlitchOffsetOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'north', 'east'])

fig = plt.figure()

ax1 = fig.add_subplot(211)

ax1.set_title("GPS glitch offset")    
#ax1.set_xlabel('time (s)')
ax1.set_ylabel('north position (m)')
ax1.plot(data['time'], data['north'], color='b', label='north')


ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('east position (m)')
ax2.plot(data['time'], data['east'], color='b', label='east')

plt.show()
