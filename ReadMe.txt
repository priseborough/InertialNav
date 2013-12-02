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

You will need Matlab + Simulink 2012a or later to run this model

Instructions to run:

1) Add models, plots, scripts and TestData directories to your Matlab path
2) Make models your working directory
3) Run the RunNavFilterTestHarness.m script file
4) Test data will be loaded and the model will be built, run and plots generated.

You can load other test data by modifying the file load command at the top of the LoadNavFilterTestData.m script file.

