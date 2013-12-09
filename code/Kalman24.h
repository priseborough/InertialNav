/// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-
/*
  24 state EKF based on work by Paul Riseborough

  Converted from Matlab to C++ by Andrew Tridgell

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

#include <AP_Math.h>

class Kalman24
{
public:
    // constructor
    Kalman24() {}

    void CovariancePrediction(Vector3f deltaAngle, 
                              Vector3f deltaVelocity,
                              float dt,
                              bool onGround);

private:
    // note that some of these are temporary variables, but we don't 
    // want huge spikes in stack usage, so better to make these class variables
	
	// Generic
    float _processNoise[24];
	float _KHP[24][24];
	float _Tnb[3][3];
    float _states[24];
    float _nextStates[24];
    float _P[24][24];
    float _nextP[24][24];
	
	// Covariance Predition
    float _SF[21];
    float _SG[8];
    float _SQ[11];
    float _SPP[13];
	

	// Velocity and Position Measurement Fusion
	float _velPosInnov[6];
	float _varVelPosInnov[6];
	
	// Magnetometer Measurement Fusion
	float _magInnov[3];
	float _varMagInnov[3];
	float _SH_MAG[9];
	float _H_MAG[1][24]
	float _SK_MX[6]
	float _SK_MY[6]
	float _SK_MZ[6]
	float _K_MAG[24][1]

    float P(uint8_t i, uint8_t j) { return _P[i-1][j-1]; }
    float &nextP(uint8_t i, uint8_t j) { return _predP[i-1][j-1]; }
    float &SF(uint8_t i) { return _SF[i-1]; }
    float &SG(uint8_t i) { return _SG[i-1]; }
    float &SQ(uint8_t i) { return _SQ[i-1]; }
    float &SPP(uint8_t i) { return _SPP[i-1]; }
    float &processNoise(uint8_t i) { return _processNoise[i-1]; }
    float &states(uint8_t i) { return _states[i-1]; }

    void predPzeroRows(uint8_t first, uint8_t last);
    void predPzeroCols(uint8_t first, uint8_t last);
};

