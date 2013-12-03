Files for a prototype 24 state Extended Kalman filter designed for APMPlane implementation
Author: Paul Riseborough

This is an implementation of a strapdown inertial navigation system with an Extended Kalman Filter algorithm used 
to provide aiding using the following data sources:

3D velocity measurements (taken from UBlox6 or equivalent GPS)
2D position (taken from UBlox6 or equivalent GPS)
Height measurement (could be GPS or barometric height or some combination)
3-axis body fixed magnetic flux measurements
True airspeed measurement

The IMU delta angles and delta velocities are assumed to be simple integrals of the corresponding angular rates 
and axial accelerations, with no coning or sculling compensation applied to them before input to the filter.

The filter estimates the following states:
4 quaternion parameters
3 North,East,Down velocity components
3 North,East,Down  position components
3 IMU delta angle bias components
3 IMU delta velocity bias components
2 North,East wind velocity components
3 North,East,Down  earth magnetic flux components
3 X,Y,Z body fixed magnetic flux components (these are opposite sign to the compass offsets used by APM)

The filter is designed to run in parallel with the existing APM AHRS complementary filter, firstly to provide
a bootstrap for initial alignment, and secondly to provide a watchdog reference to detect filter divergence.

No optimisation of tuning parameters has been performed, other than to roughly set them based on innovation sequence
RMS errors obtained from the test data set.

The three test data sets were obtained from a PX4 FMU and digital airspeed sensor, installed in a Skywalker X-8 
airframe. The 2nd and 3rd data sets include aerobatic manoeuvres (loops and axial rolls), but no spinning due to
airframe limitations.

You will need Matlab + Simulink 2012a or later to run this model

Instructions to run:

1) Add 'plots', 'scripts' and 'TestData' directories to your Matlab path
2) Make 'models' your working directory
3) Run the RunNavFilterTestHarness.m script file
4) Test data will be loaded and the model will be built, run and plots generated.

You can load other test data by modifying the file load command at the top of the LoadNavFilterTestData.m script file.

