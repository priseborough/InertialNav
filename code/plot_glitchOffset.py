#!/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import numpy as np
import math

# timestamp, velN, velE, velD, posN, posE, posOffN, posOffE, velOffN, velOffE
data = np.genfromtxt('GlitchOffsetOut.txt', delimiter=' ', skip_header=1,
	skip_footer=1, names=['time', 'velN', 'velE', 'velD', 'posN', 'posE', 'posD', 'posOffN', 'posOffE', 'velOffN', 'velOffE'])

try:
	data2 = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
		skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
		'Bx', 'By', 'Bz', 'Aoff', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd', 'dist'])
except ValueError:
	try:
		data2 = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
		skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
		'Bx', 'By', 'Bz', 'Aoff', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd'])
	except ValueError:
		data2 = np.genfromtxt('StateDataOut.txt', delimiter=' ', skip_header=1,
			skip_footer=1, names=['time', 'q1', 'q2', 'q3', 'q4', 'Vn', 'Ve', 'Vd', 'Pn', 'Pe', 'Pd',
			'Bx', 'By', 'Bz', 'Wn', 'We', 'Mn', 'Me', 'Md', 'Mbn', 'Mbe', 'Mbd'])

fig = plt.figure(10)

ax1 = fig.add_subplot(211)

ax1.set_title("GPS glitch offset")    
ax1.set_ylabel('offset (m)')
ax1.plot(data['time'], data['posOffN'], color='b', label='posOffN')
ax1.plot(data['time'], data['posOffE'], color='r', label='posOffE')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)


ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('offset (m/sec)')
ax2.plot(data['time'], data['velOffN'], color='b', label='velOffN')
ax2.plot(data['time'], data['velOffE'], color='r', label='velOffE')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc=2)

fig = plt.figure(11)

ax1 = fig.add_subplot(211)

ax1.set_title("GPS data")    
ax1.set_ylabel('(m/sec)')
ax1.plot(data['time'], data['velN'], color='r', label='velN')
ax1.plot(data['time'], data['velE'], color='g', label='velE')
ax1.plot(data['time'], data['velD'], color='b', label='velD')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)

ax2 = fig.add_subplot(212)
  
ax2.set_xlabel('time (s)')
ax2.set_ylabel('(m)')
ax2.plot(data['time'], data['posN'], color='r', label='posN')
ax2.plot(data['time'], data['posE'], color='g', label='posE')
ax2.plot(data['time'], data['posD'], color='b', label='posD')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc=2)

fig = plt.figure(12)

ax1 = fig.add_subplot(311)

ax1.set_title("GPS velNED vs. EKF vel")    
ax1.set_ylabel('(m/sec)')
ax1.plot(data['time'], data['velN'], color='r', label='velN')
ax1.plot(data2['time'], data2['Vn'], color='b', label='Vn')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc=2)

ax2 = fig.add_subplot(312)
  
ax2.set_ylabel('(m/sec)')
ax2.plot(data['time'], data['velE'], color='r', label='velE')
ax2.plot(data2['time'], data2['Ve'], color='b', label='Ve')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc=2)

ax3 = fig.add_subplot(313)
  
ax3.set_xlabel('time (s)')
ax3.set_ylabel('(m/sec)')
ax3.plot(data['time'], data['velD'], color='r', label='velD')
ax3.plot(data2['time'], data2['Vd'], color='b', label='Vd')
handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels, loc=2)


plt.show()
