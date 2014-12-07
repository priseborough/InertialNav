#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

	# State vector:
	# 0-3: quaternions (q0, q1, q2, q3)
	# 4-6: Velocity - m/sec (North, East, Down)
	# 7-9: Position - m (North, East, Down)
	# 10-12: Delta Angle bias - rad (X,Y,Z)
	# 13-14: Wind Vector  - m/sec (North,East)
	# 15-17: Earth Magnetic Field Vector - milligauss (North, East, Down)
	# 18-20: Body Magnetic Field Vector - milligauss (X,Y,Z)

data = np.genfromtxt('FlowRawOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'flowRadX', 'flowRadY', 'gyroX', 'gyroY', 'flowMsX', 'flowMsY', 'GpsVn', 'GpsVe', 'distance'])

fig = plt.figure()

ax1 = fig.add_subplot(211)

ax1.set_title("Flow angular rate estimate")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('flow X (rad)')
ax1.plot(data['time'], data['flowRadX'], color='r', label='flow x')
#ax1.plot(data['time'], data['gyroX'], color='b', label='angular rate x')
#ax1.plot(data['time'], data['predFlowX'], color='g', label='pred flow x')

ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('flow Y (rad)')
ax2.plot(data['time'], data['flowRadY'], color='r', label='flow y')
#ax2.plot(data['time'], data['gyroY'], color='b', label='angular rate y')
#ax2.plot(data['time'], data['predFlowY'], color='g', label='pred flow y')

figvel = plt.figure()

av1 = figvel.add_subplot(311)

av1.set_title("Flow Velocity estimate")    
av1.set_xlabel('time (s)')
av1.set_ylabel('north velocity (m/s)')
av1.plot(data['time'], data['flowMsX'], color='r', label='Flow Vn')
av1.plot(data['time'], data['GpsVn'], color='g', label='GPS Vn')

av2 = figvel.add_subplot(312)
  
av2.set_xlabel('time (s)')
av2.set_ylabel('east velocity (m/s)')
av2.plot(data['time'], data['flowMsY'], color='r', label='Flow Ve')
av2.plot(data['time'], data['GpsVe'], color='g', label='GPS Ve')

av3 = figvel.add_subplot(313)
 
av3.set_xlabel('time (s)')
av3.set_ylabel('ground distance (m)')
av3.plot(data['time'], data['distance'], color='b', label='distance')

plt.show()
