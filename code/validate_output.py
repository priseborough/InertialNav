#!/usr/bin/env python

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

# XXX load GPS and other data as well and run a line-by-line (or time-by-time) comparison with thresholds
print("Filter output validation successful - because we don't validate yet!")

exit(0)
