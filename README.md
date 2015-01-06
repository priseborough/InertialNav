## Inertial Navigation Estimation Library ##

[![Build Status](https://travis-ci.org/priseborough/InertialNav.svg?branch=master)](https://travis-ci.org/priseborough/InertialNav)

Files for prototype 21, 22, 23 and 24 state Extended Kalman filters designed for APMPlane implementation
Author: Paul Riseborough

This is an implementation of a strapdown inertial navigation system with an Extended Kalman Filter algorithm used 
to provide aiding using the following data sources (depending on filter variant):

3D velocity measurements (taken from UBlox6 or equivalent GPS)
2D position (taken from UBlox6 or equivalent GPS)
Height measurement (could be GPS or barometric height or some combination)
3-axis body fixed magnetic flux measurements
True airspeed measurement
Terrain laser range finder measurements (aligned with vehicle Z axis)
Ground optical flow measurements (flow angular rates about vehicle X and Y axes)

The IMU delta angles and delta velocities are assumed to be simple integrals of the corresponding angular rates 
and axial accelerations, with no coning or sculling compensation applied to them before input to the filter.

The filter estimates the following states:
4 quaternion parameters
3 North,East,Down velocity components
3 North,East,Down  position components
3 IMU delta angle bias components
2 IMU delta velocity bias in X and Y (only included in 24 state derivation)
1 IMU delta velocity bias in Z (only included in the 22 and 23 state derivation)
2 North,East wind velocity components
3 North,East,Down  earth magnetic flux components
3 X,Y,Z body fixed magnetic flux components (these are opposite sign to the compass offsets used by APM)
1 Offset of terrain along down axis (only included in the 23 state filter derivation)

The filter is designed to run in parallel with the existing APM AHRS complementary filter, firstly to provide
a bootstrap for initial alignment, and secondly to provide a watchdog reference to detect filter divergence.

No optimisation of tuning parameters has been performed, other than to roughly set them based on innovation sequence
RMS errors obtained from the test data set.

The three test data sets were obtained from a PX4 FMU and digital airspeed sensor, installed in a Skywalker X-8 
airframe. The 2nd and 3rd data sets include aerobatic manoeuvres (loops and axial rolls), but no spinning due to
airframe limitations.

Some additional data sets incorporating optical flow and range finder measurements are included.


### Instructions To Run Simulink Model ###

Note : Simulink models are only available for 21 and 24 state architecture, and do not include range finder or optical
flow measurements.

You will need Matlab + Simulink 2012a or later to run this model

Instructions to run:

1. Add 'plots', 'scripts' and 'TestData' directories to your Matlab path
2. Make 'models' your working directory
3. Run the RunNavFilterTestHarness24.m script file to run the 24-state filter, or RunNavFilterTestHarness21.m to run
   the 21-state filter
4. Test data will be loaded and the model will be built, run and plots generated.

You can load other test data by modifying the file load command at the top of the LoadNavFilterTestData.m script file.


### Instructions To Run C++ code test harness ###

1. Go to the code directory and run make
2. Make sure the ATT, GPS, IMU, MAG and NTUN.txt files are in the directory with the executable (unzip one of the data ZIP files)
3. Run: ./estimator_closed_loop_test
4. The program will put the results into the following space deliminted data files
5. Run make plots to generate Python plots

 * CovDataOut.txt
 * EulDataOut.txt
 * MagFuse.txt
 * RefVelPosDataOut.txt
 * StateDataOut.txt
 * TasFuse.txt
 * VelPosFuse.txt
