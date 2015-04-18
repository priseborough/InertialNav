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

data = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
	'Bx', 'By', 'Bz', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd'])

data2 = np.genfromtxt('RefVelPosDataOut.txt', delimiter=' ', skip_header=1,
        skip_footer=1, names=['time', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd'])

fig = plt.figure()

ax1 = fig.add_subplot(311)

ax1.set_title("Position estimate")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('north position (m)')
ax1.plot(data['time'], data['Pn'], color='b', label='Pn')
ax1.plot(data2['time'], data2['Pn'], color='r', label='Pn')
handles, labels = ax1.get_legend_handles_labels()
labels = ['estimate', 'onboard']
ax1.legend(handles, labels, loc=2)


ax2 = fig.add_subplot(312)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('east position (m)')
ax2.plot(data['time'], data['Pe'], color='b', label='Pe')
ax2.plot(data2['time'], data2['Pe'], color='r', label='Pe')
handles, labels = ax2.get_legend_handles_labels()
labels = ['estimate', 'onboard']
ax2.legend(handles, labels, loc=2)

ax3 = fig.add_subplot(313)
 
ax3.set_xlabel('time (s)')
ax3.set_ylabel('down position (m)')
ax3.plot(data['time'], data['Pd'], color='b', label='Pd')
ax3.plot(data2['time'], data2['Pd'], color='r', label='Pd')
handles, labels = ax3.get_legend_handles_labels()
labels = ['estimate', 'onboard']
ax3.legend(handles, labels, loc=2)

figvel = plt.figure()

av1 = figvel.add_subplot(311)

av1.set_title("Velocity estimate")    
av1.set_xlabel('time (s)')
av1.set_ylabel('north velocity (m/s)')
av1.plot(data['time'], data['Vn'], color='b', label='Vn')
av1.plot(data2['time'], data2['Vn'], color='r', label='Vn')
handles, labels = av1.get_legend_handles_labels()
labels = ['estimate', 'onboard']
av1.legend(handles, labels, loc=2)

av2 = figvel.add_subplot(312)
  
av2.set_xlabel('time (s)')
av2.set_ylabel('east velocity (m/s)')
av2.plot(data['time'], data['Ve'], color='b', label='Ve')
av2.plot(data2['time'], data2['Ve'], color='r', label='Ve')
handles, labels = av2.get_legend_handles_labels()
labels = ['estimate', 'onboard']
av2.legend(handles, labels, loc=2)

av3 = figvel.add_subplot(313)
 
av3.set_xlabel('time (s)')
av3.set_ylabel('down velocity (m/s)')
av3.plot(data['time'], data['Vd'], color='b', label='Vd')
av3.plot(data['time'], data2['Vd'], color='r', label='Vd')
handles, labels = av3.get_legend_handles_labels()
labels = ['estimate', 'onboard']
av3.legend(handles, labels, loc=2)

plt.show()
