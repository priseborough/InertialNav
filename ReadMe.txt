Files for prototype 23 and 21 state Extended Kalman filters designed for APMPlane implementation
Author: Paul Riseborough

This is an implementation of a strapdown inertial navigation system with an Extended Kalman Filter algorithm used 
to provide aiding using the following data sources:

3D velocity measurements (taken from UBlox6 or equivalent GPS)
2D position (taken from UBlox6 or equivalent GPS)
Height measurement (could be GPS or barometric height or some combination)
3-axis body fixed magnetic flux measurements
True airspeed measurement
Laser range finder measurement (assumed to be fixed in body axes and pointing downwards with a specified pitch offset)

The IMU delta angles and delta velocities are assumed to be simple integrals of the corresponding angular rates 
and axial accelerations, with no coning or sculling compensation applied to them before input to the filter.

The filter estimates the following states:
4 quaternion parameters
3 North,East,Down velocity components
3 North,East,Down  position components
3 IMU delta angle bias components
1 IMU delta velocity bias in Z (not included in the 21-state filter)
2 North,East wind velocity components
3 North,East,Down  earth magnetic flux components
3 X,Y,Z body fixed magnetic flux components (these are opposite sign to the compass offsets used by APM)
1 Offset of terrain along down axis (not included in 21 sate filter)

The filter is designed to run in parallel with the existing APM AHRS complementary filter, firstly to provide
a bootstrap for initial alignment, and secondly to provide a watchdog reference to detect filter divergence.

No optimisation of tuning parameters has been performed, other than to roughly set them based on innovation sequence
RMS errors obtained from the test data set.

The three test data sets were obtained from a PX4 FMU and digital airspeed sensor, installed in a Skywalker X-8 
airframe. The 2nd and 3rd data sets include aerobatic manoeuvres (loops and axial rolls), but no spinning due to
airframe limitations.

No range finder measurements have been available, so range finder data has been synthesised from the baro height data
to test the filter maths.

========================================================================================================================
Instructions To Run Simulink Model:

You will need Matlab + Simulink 2012a or later to run this model

Instructions to run:

1) Add 'plots', 'scripts' and 'TestData' directories to your Matlab path
2) Make 'models' your working directory
3) Run the RunNavFilterTestHarness24.m script file to run the 24-state filter, or RunNavFilterTestHarness21.m to run
   the 21-state filter
4) Test data will be loaded and the model will be built, run and plots generated.

You can load other test data by modifying the file load command at the top of the LoadNavFilterTestData.m script file.

========================================================================================================================
Instructions To Run C++ code test harness

1) Build the main - closed loop test harness - single precision.cpp file
2) Make sure the ATT, GPS, IMU, MAG and NTUN.txt files are in the directory with the executable
3) Run and exit when the Euler difference display stops updating
4) The program will put the results into the folliwing space deliminted data files:

CovDataOut.txt
EulDataOut.txt
MagFuse.txt
RefVelPosDataOut.txt
StateDataOut.txt
TasFuse.txt
VelPosFuse.txt

5) These text files can be plotted using the Matlab script PlotCcodeOutput.m , if you don't have Matlab a the
   script file will show you what data is in each column so you can use your plotting tool of choice
