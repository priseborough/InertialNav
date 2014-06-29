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

ax1 = fig.add_subplot(211)

ax1.set_title("Wind Velocity")    
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Wind North')
ax1.plot(data['time'], data['Wn'], color='r', label='Wind N')

ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('Wind East')
ax2.plot(data['time'], data['We'], color='g', label='Wind E')

plt.show()