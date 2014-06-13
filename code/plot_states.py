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
	# 13: Accel offset
	# 14-15: Wind Vector  - m/sec (North,East)
	# 16-18: Earth Magnetic Field Vector - milligauss (North, East, Down)
	# 19-21: Body Magnetic Field Vector - milligauss (X,Y,Z)
	# 22: Terrain
try:
	data = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
		skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
		'Bx', 'By', 'Bz', 'Aoff', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd', 'dist'])
except ValueError:
	try:
		data = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
		skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
		'Bx', 'By', 'Bz', 'Aoff', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd'])
	except ValueError:
		data = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
			skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
			'Bx', 'By', 'Bz', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd'])

fig = plt.figure()

ax1 = fig.add_subplot(611)

ax1.set_title("Offsets")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('X gyro offset')
ax1.set_ylim([-0.0025,0.0025])
ax1.plot(data['time'], data['Bx'], color='r', label='Pn')

ax2 = fig.add_subplot(612)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Y gyro offset')
ax2.set_ylim([-0.0025,0.0025])
ax2.plot(data['time'], data['By'], color='g', label='Pe')

ax3 = fig.add_subplot(613)
 
ax3.set_xlabel('time (s)')
ax3.set_ylabel('Z gyro offset')
ax3.set_ylim([-0.0025,0.0025])
ax3.plot(data['time'], data['Bz'], color='b', label='Pd')

ax4 = fig.add_subplot(614)
 
ax4.set_xlabel('time (s)')
ax4.set_ylabel('Mag body offset N')
ax4.set_ylim([-0.4,0.4])
ax4.plot(data['time'], data['Mbn'], color='b', label='Pd')

ax5 = fig.add_subplot(615)
 
ax5.set_xlabel('time (s)')
ax5.set_ylabel('Mag body offset E')
ax5.set_ylim([-0.4,0.4])
ax5.plot(data['time'], data['Mbe'], color='b', label='Pd')

ax6 = fig.add_subplot(616)
 
ax6.set_xlabel('time (s)')
ax6.set_ylabel('Mag body offset D')
ax6.set_ylim([-0.4,0.4])
ax6.plot(data['time'], data['Mbd'], color='b', label='Pd')

plt.show()