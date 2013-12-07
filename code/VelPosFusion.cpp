/// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-
/*
  Fusion of Position & Velocity Measurements for a 24 State Extended Kalman Navigation Filter
  Based on https://github.com/priseborough/InertialNav
  Converted to C++ by Paul Riseborough

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Kalman24.h"

// some macros to make the conversion from MATLab easier
#define single(x) float(x)
#define pi M_PI

/*
  run one state and covariance update step using a combination of NED Vel and Pos measurements
 */
void Kalman24::FuseVelPosNED(float accNavMag,
									bool FuseGPSData,
									Vector3f VelNED,
									Vector3f PosNE,
									float StatesAtGpsTime[24],
									bool FuseHgtData,
									float HgtMea,
									float StatesAtHgtTime[24])
{

	// declare fail counters
	static uint32_t velFailCount;
	static uint32_t posFailCount;
	static uint32_t hgtFailCount;
	
	// declare health flags
	bool velHealth = false;
	bool posHealth = false;
	bool hgtHealth = false;
	
	// declare variables used to check measurement errors
	float varTotal;
	float varLimit;
	float velInnov[3];
	float posInnov[2];
	float hgtInnov;
	
	// index used to access filter states
	uint8_t stateIndex;

	// Form the observation vector
	float observation[6];
    for (uint8_t i = 0; i<=2; i++) observation[i] = VelNED[i];
    for (uint8_t i = 3; i<=4; i++) observation[i] = PosNE[i];
    observation[5] = -HgtMea;
	
    // Estimate the GPS Velocity, GPS horiz position and height measurement variances.
	float velErr = 0.15*accNavMag; // additional error in GPS velocities caused by manoeuvring
	float posErr = 0.15*accNavMag; // additional error in GPS position caused by manoeuvring
	float R_OBS[6];
    for (uint8_t i= 1; i<=3; i++) R_OBS[i-1] = 0.01 + velErr*velErr;
    for (uint8_t i= 4; i<=5; i++) R_OBS[i-1] = 4.0 + posErr*posErr;
	R_OBS[5] = 4.0;
	
	// Specify the count before forcing use of GPS or height data after invalid
	// data. We need to force GPS again sooner without
	// airspeed data as the nav velocities will be unconstrained.
	const float gpsRetryTime = 30.0;
	const float gpsRetryTimeNoTAS = 5.0;
	const float hgtRetryTime = 30.0;
	const float dt = 0.02;
	uint32_t maxVelFailCount;
	uint32_t maxVelFailCount;
	uint32_t maxHgtFailCount;
	if useAirspeed
	{
		maxVelFailCount = gpsRetryTime/dt;
		maxPosFailCount = maxVelFailCount;
	}
	 else
	{
		maxVelFailCount = gpsRetryTimeNoTAS/dt;
		maxPosFailCount = maxVelFailCount;
	}
	maxHgtFailCount = hgtRetryTime/dt;

	// Perform sequential fusion of GPS measurements. This assumes that the
	// errors in the different velocity and position components are
	// uncorrelated which is not true, however in the absence of covariance
	// data from the GPS receiver it is the only assumption we can make
	// so we might as well take advantage of the computational efficiencies
	// associated with sequential fusion

	if (FuseGPSData || FuseHgtData)
	{
		// calculate innovation variances
		for (uint8_t obsIndex = 0; i<=5; i++)
		{
			stateIndex = 4 + obsIndex;
			_varVelPosInnov[obsIndex] = P[stateIndex,stateIndex] + R_OBS[obsIndex];
		}
		// calculate innovations and check GPS data validity against limits using a 5-sigma test
		if FuseGPSData
		{
			for (uint8_t i = 0; i<=2; i++)
			{
				velInnov[i] = StatesAtGpsTime(i+4) - VelNED[i];
			}
			varTotal = velInnov[0]*velInnov[0] + velInnov[1]*velInnov[1] + velInnov[2]*velInnov[2];
			varLimit = 25.0*(varInnov[0] + varInnov[1] + varInnov[2]);
			if (varTotal < varlimit) || (velFailCount > maxVelFailCount)
			{
				velHealth = true;
				velFailCount = 0;
			}
			else
			{
				velHealth = false;
				velFailCount = velFailCount + 1;
			}
			posInnov[0] = StatesAtGpsTime(7) - PosNE[0];
			posInnov[1] = StatesAtGpsTime(8) - PosNE[1];
			varTotal = posInnov[0]*posInnov[0] + posInnov[1]*posInnov[1];
			varLimit = 25.0*(varInnov[3] + varInnov[4]);
			if (varTotal < varlimit) || (posFailCount > maxPosFailCount)
			{
				posHealth = true;
				posFailCount = 0;
			}
			else
			{
				posHealth = false;
				posFailCount = posFailCount + 1;
			}
		}
		// calculate innovations and check height data validity against limits using a 5-sigma test
		if FuseHgtData
		{
			hgtInnov = StatesAtHgtTime(9) + HgtMea;
			if (hgtInnov*hgtInnov) < 25.0*varInnov(5) || (hgtFailCount > maxHgtFailCount)
			{
				hgtHealth = true;
				hgtFailCount = 0;
			}
			else
			{
				hgtHealth = false;
				hgtFailCount = hgtFailCount + 1;
			}
		}
		// Set range for sequential fusion of velocity and position measurements depending
		// on which data is available
		uint8_t startIndex;
		if FuseGPSData
		{
			startIndex = 0;
		}
		else
		{
			startIndex = 5;
		}
		uint8_t endIndex;
		if FuseHgtData
		{
			endIndex = 5;
		}
		else
		{
			endIndex = 4;
		}
		// Fuse measurements sequentially
		for (uint8_t obsIndex= startIndex; i<=endIndex; i++)
		{
			// Apply data health checks
			if (velHealth && (obsIndex >= 0 && obsIndex <= 2)) || (posHealth && (obsIndex == 3 || obsIndex == 4)) || (hgtHealth && (obsIndex == 5))
			{	
				stateIndex = 4 + obsIndex;
				// Calculate the measurement innovation, using states from a
				// different time coordinate if fusing height data
				if (obsIndex == 5)
				{
					_velPosInnov[obsIndex] = StatesAtHgtTime[stateIndex] - observation[obsIndex];
				}
				else
				{
					_velPosInnov[obsIndex] = StatesAtGpsTime[stateIndex] - observation[obsIndex];
				}
				// Calculate the Kalman Gain
				// Calculate innovation variances - also used for data logging
				_varVelPosInnov[obsIndex] = P[stateIndex,stateIndex] + R_OBS[obsIndex];
				float SK = 1.0/_varVelPosInnov[obsIndex];
				// Check the innovation for consistency and don't fuse if > TBD Sigma
				// Currently disabled for testing
				for (uint8_t i= 0; i<=23; i++)
				{
					K[i] = P[i,stateIndex]*SK;
				}
				// Calculate state corrections
				for (uint8_t i = 0; i<=23; i++)
				{
					xk[i] = K[i] * _velPosInnov[obsIndex];
				}
				// Update the covariance - take advantage of direct observation of a
				// single state at index = stateIndex to reduce computations
				// Optimised implementation of standard equation P = (I - K*H)*P;
				for (uint8_t rowIndex= 0; i<=23; i++)
				{
					for (uint8_t colIndex= 0; i<=23; i++)
					{
						_KHP[rowIndex,colIndex] = K[rowIndex]*P[stateIndex,colIndex];
					}
				}
				for (uint8_t rowIndex= 0; i<=23; i++)
				{
					for (uint8_t colIndex= 0; i<=23; i++)
					{
						P[rowIndex,colIndex] = P[rowIndex,colIndex] - _KHP[rowIndex,colIndex];
					}
				}
				// Apply the state corrections and re-normalise the quaternions
				for (uint8_t rowIndex = 0; i<=23; i++)
				{
					states[i] = states[i] - xk[i];
				}
				// Normalise the quaternion states
				float quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
				if (quatMag > 1e-12)
				{
					for (uint8_t rowIndex = 0; i<=3; i++)
					{
						states[i] = states[i] / quatMag;
					}
				}
			}
		}
	// Force symmetry on the covariance matrix to prevent ill-conditioning
	// of the matrix which would cause the filter to blow-up
	for (uint8_t rowIndex = 1; i<=23; i++)
	{
		for (uint8_t colIndex = rowIndex-1; i<=22; i++)
		{
			float temp = 0.5*(P[rowIndex,colIndex] + P[colIndex,rowIndex]);
			P[rowIndex,colIndex] = temp;
			P[colindex,rowIndex] = temp;
		}
	}
	// Output corrected covariance matrix and state vector
	_nextP = P;
	_nextStates = states;
	}
}