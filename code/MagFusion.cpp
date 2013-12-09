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
void Kalman24::FuseMagnetometer(float accNavMag,
									bool FuseGPSData,
									Vector3f MagData,
									Vector3f PosNE,
									float StatesAtMeasTime[24],
									bool useCompass)
{
static float q0 = 1.0;
static float q1 = 0.0;
static float q2 = 0.0;
static float q3 = 0.0;
static float magN = 0.0;
static float magE = 0.0;
static float magD = 0.0;
static uint8_t obsIndex = 1;
float Tnb[3][[3] = {
{1.0,0.0,0.0} ,
{0.0,1.0,0.0} ,
{0.0,0.0,1.0}
};
static float MagPred[3] = {0.0,0.0,0.0};

// define magnetometer measurement error variance (milligauss).
const float R_MAG = 2500.0;

// Perform sequential fusion of Magnetometer measurements.
// This assumes that the errors in the dfferent componenets are
// uncorrelated which is not true, however in the absence of covariance
// data fit is the only assumption we can make
// so we might as well take advantage of the computational efficiencies
// associated with sequential fusion
if (useCompass && FuseData) || (useCompass && (obsIndex == 2 || obsIndex == 3))
    
    // Sequential fusion of XYZ components to spread processing load across
    // three prediction time steps. 
    
    // Calculate observation jacobians and Kalman gains
    K_MAG = single(zeros(24,1));
    H_MAG = single(zeros(1,24));
    if FuseData
        // Copy required states to local variable names
        q0 = StatesAtMeasTime[0];
        q1 = StatesAtMeasTime[1];
        q2 = StatesAtMeasTime[2];
        q3 = StatesAtMeasTime[3];
        magN = StatesAtMeasTime[18);
        magE = StatesAtMeasTime[19];
        magD = StatesAtMeasTime[20];
        magXbias = StatesAtMeasTime[21];
        magYbias = StatesAtMeasTime[22];
        magZbias = StatesAtMeasTime[23];
        // rotate predicted earth components into body axes and calculate
        // predicted measurements
        Tnb = single([q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q0*q3) , 2*(q1*q3-q0*q2) ;...
            2*(q1*q2 - q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q0*q1);...
            2*(q1*q3 + q0*q2) , 2*(q2*q3 - q0*q1) , q0^2 - q1^2 - q2^2 + q3^2]);
		MagPred[0] = Tnb[0][0]*magN + Tnb[0][1]*magE + Tnb[0][2]*magD + magXbias;
		MagPred[1] = Tnb[1][0]*magN + Tnb[1][1]*magE + Tnb[1][2]*magD + magXbias;
		MagPred[2] = Tnb[2][0]*magN + Tnb[2][1]*magE + Tnb[2][2]*magD + magZbias;
		
        // Calculate observation jacobian common terms
        SH_MAG(1) = 2*magD*q3 + 2*magE*q2 + 2*magN*q1;
        SH_MAG(2) = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
        SH_MAG(3) = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
        SH_MAG(4) = q3^2;
        SH_MAG(5) = q2^2;
        SH_MAG(6) = q1^2;
        SH_MAG(7) = q0^2;
        SH_MAG(8) = 2*magN*q0;
        SH_MAG(9) = 2*magE*q3;
        // calculate observation jacobian
        H_MAG(1,1) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        H_MAG(1,2) = SH_MAG(1);
        H_MAG(1,3) = 2*magE*q1 - 2*magD*q0 - 2*magN*q2;
        H_MAG(1,4) = SH_MAG(3);
        H_MAG(1,19) = SH_MAG(6) - SH_MAG(5) - SH_MAG(4) + SH_MAG(7);
        H_MAG(1,20) = 2*q0*q3 + 2*q1*q2;
        H_MAG(1,21) = 2*q1*q3 - 2*q0*q2;
        H_MAG(1,22) = 1;
        // calculate Kalman gain
        SK_MX = single(zeros(1,6));
        SK_MX(1) = 1/(P(22,22) + R_MAG + P(2,22)*SH_MAG(1) + P(4,22)*SH_MAG(3) - P(19,22)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) - (2*magD*q0 - 2*magE*q1 + 2*magN*q2)*(P(22,3) + P(2,3)*SH_MAG(1) + P(4,3)*SH_MAG(3) - P(19,3)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,3)*(2*q0*q3 + 2*q1*q2) - P(21,3)*(2*q0*q2 - 2*q1*q3) - P(3,3)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(22,1) + P(2,1)*SH_MAG(1) + P(4,1)*SH_MAG(3) - P(19,1)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,1)*(2*q0*q3 + 2*q1*q2) - P(21,1)*(2*q0*q2 - 2*q1*q3) - P(3,1)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(1)*(P(22,2) + P(2,2)*SH_MAG(1) + P(4,2)*SH_MAG(3) - P(19,2)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,2)*(2*q0*q3 + 2*q1*q2) - P(21,2)*(2*q0*q2 - 2*q1*q3) - P(3,2)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(3)*(P(22,4) + P(2,4)*SH_MAG(1) + P(4,4)*SH_MAG(3) - P(19,4)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,4)*(2*q0*q3 + 2*q1*q2) - P(21,4)*(2*q0*q2 - 2*q1*q3) - P(3,4)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7))*(P(22,19) + P(2,19)*SH_MAG(1) + P(4,19)*SH_MAG(3) - P(19,19)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,19)*(2*q0*q3 + 2*q1*q2) - P(21,19)*(2*q0*q2 - 2*q1*q3) - P(3,19)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(20,22)*(2*q0*q3 + 2*q1*q2) - P(21,22)*(2*q0*q2 - 2*q1*q3) - P(3,22)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + (2*q0*q3 + 2*q1*q2)*(P(22,20) + P(2,20)*SH_MAG(1) + P(4,20)*SH_MAG(3) - P(19,20)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,20)*(2*q0*q3 + 2*q1*q2) - P(21,20)*(2*q0*q2 - 2*q1*q3) - P(3,20)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,20)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (2*q0*q2 - 2*q1*q3)*(P(22,21) + P(2,21)*SH_MAG(1) + P(4,21)*SH_MAG(3) - P(19,21)*(SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7)) + P(20,21)*(2*q0*q3 + 2*q1*q2) - P(21,21)*(2*q0*q2 - 2*q1*q3) - P(3,21)*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P(1,21)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(1,22)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MX(2) = SH_MAG(4) + SH_MAG(5) - SH_MAG(6) - SH_MAG(7);
        SK_MX(3) = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
        SK_MX(4) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MX(5) = 2*q0*q2 - 2*q1*q3;
        SK_MX(6) = 2*q0*q3 + 2*q1*q2;
        K_MAG(1,1) = SK_MX(1)*(P(1,22) + P(1,2)*SH_MAG(1) + P(1,4)*SH_MAG(3) + P(1,1)*SK_MX(4) - P(1,3)*SK_MX(3) - P(1,19)*SK_MX(2) + P(1,20)*SK_MX(6) - P(1,21)*SK_MX(5));
        K_MAG(2,1) = SK_MX(1)*(P(2,22) + P(2,2)*SH_MAG(1) + P(2,4)*SH_MAG(3) + P(2,1)*SK_MX(4) - P(2,3)*SK_MX(3) - P(2,19)*SK_MX(2) + P(2,20)*SK_MX(6) - P(2,21)*SK_MX(5));
        K_MAG(3,1) = SK_MX(1)*(P(3,22) + P(3,2)*SH_MAG(1) + P(3,4)*SH_MAG(3) + P(3,1)*SK_MX(4) - P(3,3)*SK_MX(3) - P(3,19)*SK_MX(2) + P(3,20)*SK_MX(6) - P(3,21)*SK_MX(5));
        K_MAG(4,1) = SK_MX(1)*(P(4,22) + P(4,2)*SH_MAG(1) + P(4,4)*SH_MAG(3) + P(4,1)*SK_MX(4) - P(4,3)*SK_MX(3) - P(4,19)*SK_MX(2) + P(4,20)*SK_MX(6) - P(4,21)*SK_MX(5));
        K_MAG(5,1) = SK_MX(1)*(P(5,22) + P(5,2)*SH_MAG(1) + P(5,4)*SH_MAG(3) + P(5,1)*SK_MX(4) - P(5,3)*SK_MX(3) - P(5,19)*SK_MX(2) + P(5,20)*SK_MX(6) - P(5,21)*SK_MX(5));
        K_MAG(6,1) = SK_MX(1)*(P(6,22) + P(6,2)*SH_MAG(1) + P(6,4)*SH_MAG(3) + P(6,1)*SK_MX(4) - P(6,3)*SK_MX(3) - P(6,19)*SK_MX(2) + P(6,20)*SK_MX(6) - P(6,21)*SK_MX(5));
        K_MAG(7,1) = SK_MX(1)*(P(7,22) + P(7,2)*SH_MAG(1) + P(7,4)*SH_MAG(3) + P(7,1)*SK_MX(4) - P(7,3)*SK_MX(3) - P(7,19)*SK_MX(2) + P(7,20)*SK_MX(6) - P(7,21)*SK_MX(5));
        K_MAG(8,1) = SK_MX(1)*(P(8,22) + P(8,2)*SH_MAG(1) + P(8,4)*SH_MAG(3) + P(8,1)*SK_MX(4) - P(8,3)*SK_MX(3) - P(8,19)*SK_MX(2) + P(8,20)*SK_MX(6) - P(8,21)*SK_MX(5));
        K_MAG(9,1) = SK_MX(1)*(P(9,22) + P(9,2)*SH_MAG(1) + P(9,4)*SH_MAG(3) + P(9,1)*SK_MX(4) - P(9,3)*SK_MX(3) - P(9,19)*SK_MX(2) + P(9,20)*SK_MX(6) - P(9,21)*SK_MX(5));
        K_MAG(10,1) = SK_MX(1)*(P(10,22) + P(10,2)*SH_MAG(1) + P(10,4)*SH_MAG(3) + P(10,1)*SK_MX(4) - P(10,3)*SK_MX(3) - P(10,19)*SK_MX(2) + P(10,20)*SK_MX(6) - P(10,21)*SK_MX(5));
        K_MAG(11,1) = SK_MX(1)*(P(11,22) + P(11,2)*SH_MAG(1) + P(11,4)*SH_MAG(3) + P(11,1)*SK_MX(4) - P(11,3)*SK_MX(3) - P(11,19)*SK_MX(2) + P(11,20)*SK_MX(6) - P(11,21)*SK_MX(5));
        K_MAG(12,1) = SK_MX(1)*(P(12,22) + P(12,2)*SH_MAG(1) + P(12,4)*SH_MAG(3) + P(12,1)*SK_MX(4) - P(12,3)*SK_MX(3) - P(12,19)*SK_MX(2) + P(12,20)*SK_MX(6) - P(12,21)*SK_MX(5));
        K_MAG(13,1) = SK_MX(1)*(P(13,22) + P(13,2)*SH_MAG(1) + P(13,4)*SH_MAG(3) + P(13,1)*SK_MX(4) - P(13,3)*SK_MX(3) - P(13,19)*SK_MX(2) + P(13,20)*SK_MX(6) - P(13,21)*SK_MX(5));
        K_MAG(14,1) = SK_MX(1)*(P(14,22) + P(14,2)*SH_MAG(1) + P(14,4)*SH_MAG(3) + P(14,1)*SK_MX(4) - P(14,3)*SK_MX(3) - P(14,19)*SK_MX(2) + P(14,20)*SK_MX(6) - P(14,21)*SK_MX(5));
        K_MAG(15,1) = SK_MX(1)*(P(15,22) + P(15,2)*SH_MAG(1) + P(15,4)*SH_MAG(3) + P(15,1)*SK_MX(4) - P(15,3)*SK_MX(3) - P(15,19)*SK_MX(2) + P(15,20)*SK_MX(6) - P(15,21)*SK_MX(5));
        K_MAG(16,1) = SK_MX(1)*(P(16,22) + P(16,2)*SH_MAG(1) + P(16,4)*SH_MAG(3) + P(16,1)*SK_MX(4) - P(16,3)*SK_MX(3) - P(16,19)*SK_MX(2) + P(16,20)*SK_MX(6) - P(16,21)*SK_MX(5));
        K_MAG(17,1) = SK_MX(1)*(P(17,22) + P(17,2)*SH_MAG(1) + P(17,4)*SH_MAG(3) + P(17,1)*SK_MX(4) - P(17,3)*SK_MX(3) - P(17,19)*SK_MX(2) + P(17,20)*SK_MX(6) - P(17,21)*SK_MX(5));
        K_MAG(18,1) = SK_MX(1)*(P(18,22) + P(18,2)*SH_MAG(1) + P(18,4)*SH_MAG(3) + P(18,1)*SK_MX(4) - P(18,3)*SK_MX(3) - P(18,19)*SK_MX(2) + P(18,20)*SK_MX(6) - P(18,21)*SK_MX(5));
        K_MAG(19,1) = SK_MX(1)*(P(19,22) + P(19,2)*SH_MAG(1) + P(19,4)*SH_MAG(3) + P(19,1)*SK_MX(4) - P(19,3)*SK_MX(3) - P(19,19)*SK_MX(2) + P(19,20)*SK_MX(6) - P(19,21)*SK_MX(5));
        K_MAG(20,1) = SK_MX(1)*(P(20,22) + P(20,2)*SH_MAG(1) + P(20,4)*SH_MAG(3) + P(20,1)*SK_MX(4) - P(20,3)*SK_MX(3) - P(20,19)*SK_MX(2) + P(20,20)*SK_MX(6) - P(20,21)*SK_MX(5));
        K_MAG(21,1) = SK_MX(1)*(P(21,22) + P(21,2)*SH_MAG(1) + P(21,4)*SH_MAG(3) + P(21,1)*SK_MX(4) - P(21,3)*SK_MX(3) - P(21,19)*SK_MX(2) + P(21,20)*SK_MX(6) - P(21,21)*SK_MX(5));
        K_MAG(22,1) = SK_MX(1)*(P(22,22) + P(22,2)*SH_MAG(1) + P(22,4)*SH_MAG(3) + P(22,1)*SK_MX(4) - P(22,3)*SK_MX(3) - P(22,19)*SK_MX(2) + P(22,20)*SK_MX(6) - P(22,21)*SK_MX(5));
        K_MAG(23,1) = SK_MX(1)*(P(23,22) + P(23,2)*SH_MAG(1) + P(23,4)*SH_MAG(3) + P(23,1)*SK_MX(4) - P(23,3)*SK_MX(3) - P(23,19)*SK_MX(2) + P(23,20)*SK_MX(6) - P(23,21)*SK_MX(5));
        K_MAG(24,1) = SK_MX(1)*(P(24,22) + P(24,2)*SH_MAG(1) + P(24,4)*SH_MAG(3) + P(24,1)*SK_MX(4) - P(24,3)*SK_MX(3) - P(24,19)*SK_MX(2) + P(24,20)*SK_MX(6) - P(24,21)*SK_MX(5));
        _varMagInnov(1) = 1/SK_MX(1);
        
        // reset the observation index to 1 (we start by fusing the X
        // measurement)
        obsIndex = 1;
    elseif obsIndex == 2 // we are now fusing the Y measurement
        H_MAG(1,1) = SH_MAG(3);
        H_MAG(1,2) = SH_MAG(2);
        H_MAG(1,3) = SH_MAG(1);
        H_MAG(1,4) = 2*magD*q2 - SH_MAG(9) - SH_MAG(8);
        H_MAG(1,19) = 2*q1*q2 - 2*q0*q3;
        H_MAG(1,20) = SH_MAG(5) - SH_MAG(4) - SH_MAG(6) + SH_MAG(7);
        H_MAG(1,21) = 2*q0*q1 + 2*q2*q3;
        H_MAG(1,23) = 1;
        // calculate the Kalman gain
        SK_MY = single(zeros(5,1));
        SK_MY(1) = 1/(P(23,23) + R_MAG + P(1,23)*SH_MAG(3) + P(2,23)*SH_MAG(2) + P(3,23)*SH_MAG(1) - P(20,23)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - (2*q0*q3 - 2*q1*q2)*(P(23,19) + P(1,19)*SH_MAG(3) + P(2,19)*SH_MAG(2) + P(3,19)*SH_MAG(1) - P(20,19)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,19)*(2*q0*q3 - 2*q1*q2) + P(21,19)*(2*q0*q1 + 2*q2*q3) - P(4,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (2*q0*q1 + 2*q2*q3)*(P(23,21) + P(1,21)*SH_MAG(3) + P(2,21)*SH_MAG(2) + P(3,21)*SH_MAG(1) - P(20,21)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,21)*(2*q0*q3 - 2*q1*q2) + P(21,21)*(2*q0*q1 + 2*q2*q3) - P(4,21)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(23,4) + P(1,4)*SH_MAG(3) + P(2,4)*SH_MAG(2) + P(3,4)*SH_MAG(1) - P(20,4)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,4)*(2*q0*q3 - 2*q1*q2) + P(21,4)*(2*q0*q1 + 2*q2*q3) - P(4,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - P(19,23)*(2*q0*q3 - 2*q1*q2) + P(21,23)*(2*q0*q1 + 2*q2*q3) + SH_MAG(3)*(P(23,1) + P(1,1)*SH_MAG(3) + P(2,1)*SH_MAG(2) + P(3,1)*SH_MAG(1) - P(20,1)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,1)*(2*q0*q3 - 2*q1*q2) + P(21,1)*(2*q0*q1 + 2*q2*q3) - P(4,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(2)*(P(23,2) + P(1,2)*SH_MAG(3) + P(2,2)*SH_MAG(2) + P(3,2)*SH_MAG(1) - P(20,2)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,2)*(2*q0*q3 - 2*q1*q2) + P(21,2)*(2*q0*q1 + 2*q2*q3) - P(4,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(1)*(P(23,3) + P(1,3)*SH_MAG(3) + P(2,3)*SH_MAG(2) + P(3,3)*SH_MAG(1) - P(20,3)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,3)*(2*q0*q3 - 2*q1*q2) + P(21,3)*(2*q0*q1 + 2*q2*q3) - P(4,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7))*(P(23,20) + P(1,20)*SH_MAG(3) + P(2,20)*SH_MAG(2) + P(3,20)*SH_MAG(1) - P(20,20)*(SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7)) - P(19,20)*(2*q0*q3 - 2*q1*q2) + P(21,20)*(2*q0*q1 + 2*q2*q3) - P(4,20)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - P(4,23)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MY(2) = SH_MAG(4) - SH_MAG(5) + SH_MAG(6) - SH_MAG(7);
        SK_MY(3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MY(4) = 2*q0*q3 - 2*q1*q2;
        SK_MY(5) = 2*q0*q1 + 2*q2*q3;
        K_MAG(1,1) = SK_MY(1)*(P(1,23) + P(1,1)*SH_MAG(3) + P(1,2)*SH_MAG(2) + P(1,3)*SH_MAG(1) - P(1,4)*SK_MY(3) - P(1,20)*SK_MY(2) - P(1,19)*SK_MY(4) + P(1,21)*SK_MY(5));
        K_MAG(2,1) = SK_MY(1)*(P(2,23) + P(2,1)*SH_MAG(3) + P(2,2)*SH_MAG(2) + P(2,3)*SH_MAG(1) - P(2,4)*SK_MY(3) - P(2,20)*SK_MY(2) - P(2,19)*SK_MY(4) + P(2,21)*SK_MY(5));
        K_MAG(3,1) = SK_MY(1)*(P(3,23) + P(3,1)*SH_MAG(3) + P(3,2)*SH_MAG(2) + P(3,3)*SH_MAG(1) - P(3,4)*SK_MY(3) - P(3,20)*SK_MY(2) - P(3,19)*SK_MY(4) + P(3,21)*SK_MY(5));
        K_MAG(4,1) = SK_MY(1)*(P(4,23) + P(4,1)*SH_MAG(3) + P(4,2)*SH_MAG(2) + P(4,3)*SH_MAG(1) - P(4,4)*SK_MY(3) - P(4,20)*SK_MY(2) - P(4,19)*SK_MY(4) + P(4,21)*SK_MY(5));
        K_MAG(5,1) = SK_MY(1)*(P(5,23) + P(5,1)*SH_MAG(3) + P(5,2)*SH_MAG(2) + P(5,3)*SH_MAG(1) - P(5,4)*SK_MY(3) - P(5,20)*SK_MY(2) - P(5,19)*SK_MY(4) + P(5,21)*SK_MY(5));
        K_MAG(6,1) = SK_MY(1)*(P(6,23) + P(6,1)*SH_MAG(3) + P(6,2)*SH_MAG(2) + P(6,3)*SH_MAG(1) - P(6,4)*SK_MY(3) - P(6,20)*SK_MY(2) - P(6,19)*SK_MY(4) + P(6,21)*SK_MY(5));
        K_MAG(7,1) = SK_MY(1)*(P(7,23) + P(7,1)*SH_MAG(3) + P(7,2)*SH_MAG(2) + P(7,3)*SH_MAG(1) - P(7,4)*SK_MY(3) - P(7,20)*SK_MY(2) - P(7,19)*SK_MY(4) + P(7,21)*SK_MY(5));
        K_MAG(8,1) = SK_MY(1)*(P(8,23) + P(8,1)*SH_MAG(3) + P(8,2)*SH_MAG(2) + P(8,3)*SH_MAG(1) - P(8,4)*SK_MY(3) - P(8,20)*SK_MY(2) - P(8,19)*SK_MY(4) + P(8,21)*SK_MY(5));
        K_MAG(9,1) = SK_MY(1)*(P(9,23) + P(9,1)*SH_MAG(3) + P(9,2)*SH_MAG(2) + P(9,3)*SH_MAG(1) - P(9,4)*SK_MY(3) - P(9,20)*SK_MY(2) - P(9,19)*SK_MY(4) + P(9,21)*SK_MY(5));
        K_MAG(10,1) = SK_MY(1)*(P(10,23) + P(10,1)*SH_MAG(3) + P(10,2)*SH_MAG(2) + P(10,3)*SH_MAG(1) - P(10,4)*SK_MY(3) - P(10,20)*SK_MY(2) - P(10,19)*SK_MY(4) + P(10,21)*SK_MY(5));
        K_MAG(11,1) = SK_MY(1)*(P(11,23) + P(11,1)*SH_MAG(3) + P(11,2)*SH_MAG(2) + P(11,3)*SH_MAG(1) - P(11,4)*SK_MY(3) - P(11,20)*SK_MY(2) - P(11,19)*SK_MY(4) + P(11,21)*SK_MY(5));
        K_MAG(12,1) = SK_MY(1)*(P(12,23) + P(12,1)*SH_MAG(3) + P(12,2)*SH_MAG(2) + P(12,3)*SH_MAG(1) - P(12,4)*SK_MY(3) - P(12,20)*SK_MY(2) - P(12,19)*SK_MY(4) + P(12,21)*SK_MY(5));
        K_MAG(13,1) = SK_MY(1)*(P(13,23) + P(13,1)*SH_MAG(3) + P(13,2)*SH_MAG(2) + P(13,3)*SH_MAG(1) - P(13,4)*SK_MY(3) - P(13,20)*SK_MY(2) - P(13,19)*SK_MY(4) + P(13,21)*SK_MY(5));
        K_MAG(14,1) = SK_MY(1)*(P(14,23) + P(14,1)*SH_MAG(3) + P(14,2)*SH_MAG(2) + P(14,3)*SH_MAG(1) - P(14,4)*SK_MY(3) - P(14,20)*SK_MY(2) - P(14,19)*SK_MY(4) + P(14,21)*SK_MY(5));
        K_MAG(15,1) = SK_MY(1)*(P(15,23) + P(15,1)*SH_MAG(3) + P(15,2)*SH_MAG(2) + P(15,3)*SH_MAG(1) - P(15,4)*SK_MY(3) - P(15,20)*SK_MY(2) - P(15,19)*SK_MY(4) + P(15,21)*SK_MY(5));
        K_MAG(16,1) = SK_MY(1)*(P(16,23) + P(16,1)*SH_MAG(3) + P(16,2)*SH_MAG(2) + P(16,3)*SH_MAG(1) - P(16,4)*SK_MY(3) - P(16,20)*SK_MY(2) - P(16,19)*SK_MY(4) + P(16,21)*SK_MY(5));
        K_MAG(17,1) = SK_MY(1)*(P(17,23) + P(17,1)*SH_MAG(3) + P(17,2)*SH_MAG(2) + P(17,3)*SH_MAG(1) - P(17,4)*SK_MY(3) - P(17,20)*SK_MY(2) - P(17,19)*SK_MY(4) + P(17,21)*SK_MY(5));
        K_MAG(18,1) = SK_MY(1)*(P(18,23) + P(18,1)*SH_MAG(3) + P(18,2)*SH_MAG(2) + P(18,3)*SH_MAG(1) - P(18,4)*SK_MY(3) - P(18,20)*SK_MY(2) - P(18,19)*SK_MY(4) + P(18,21)*SK_MY(5));
        K_MAG(19,1) = SK_MY(1)*(P(19,23) + P(19,1)*SH_MAG(3) + P(19,2)*SH_MAG(2) + P(19,3)*SH_MAG(1) - P(19,4)*SK_MY(3) - P(19,20)*SK_MY(2) - P(19,19)*SK_MY(4) + P(19,21)*SK_MY(5));
        K_MAG(20,1) = SK_MY(1)*(P(20,23) + P(20,1)*SH_MAG(3) + P(20,2)*SH_MAG(2) + P(20,3)*SH_MAG(1) - P(20,4)*SK_MY(3) - P(20,20)*SK_MY(2) - P(20,19)*SK_MY(4) + P(20,21)*SK_MY(5));
        K_MAG(21,1) = SK_MY(1)*(P(21,23) + P(21,1)*SH_MAG(3) + P(21,2)*SH_MAG(2) + P(21,3)*SH_MAG(1) - P(21,4)*SK_MY(3) - P(21,20)*SK_MY(2) - P(21,19)*SK_MY(4) + P(21,21)*SK_MY(5));
        K_MAG(22,1) = SK_MY(1)*(P(22,23) + P(22,1)*SH_MAG(3) + P(22,2)*SH_MAG(2) + P(22,3)*SH_MAG(1) - P(22,4)*SK_MY(3) - P(22,20)*SK_MY(2) - P(22,19)*SK_MY(4) + P(22,21)*SK_MY(5));
        K_MAG(23,1) = SK_MY(1)*(P(23,23) + P(23,1)*SH_MAG(3) + P(23,2)*SH_MAG(2) + P(23,3)*SH_MAG(1) - P(23,4)*SK_MY(3) - P(23,20)*SK_MY(2) - P(23,19)*SK_MY(4) + P(23,21)*SK_MY(5));
        K_MAG(24,1) = SK_MY(1)*(P(24,23) + P(24,1)*SH_MAG(3) + P(24,2)*SH_MAG(2) + P(24,3)*SH_MAG(1) - P(24,4)*SK_MY(3) - P(24,20)*SK_MY(2) - P(24,19)*SK_MY(4) + P(24,21)*SK_MY(5));
        _varMagInnov(2) = 1/SK_MY(1);
    elseif obsIndex == 3 // we are now fusing the Z measurement
        // calculate the observation jacobian
        H_MAG(1,1) = SH_MAG(2);
        H_MAG(1,2) = 2*magN*q3 - 2*magE*q0 - 2*magD*q1;
        H_MAG(1,3) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        H_MAG(1,4) = SH_MAG(1);
        H_MAG(1,19) = 2*q0*q2 + 2*q1*q3;
        H_MAG(1,20) = 2*q2*q3 - 2*q0*q1;
        H_MAG(1,21) = SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7);
        H_MAG(1,24) = 1;
        // calculate the Kalman gain
        SK_MZ = single(zeros(6,1));
        SK_MZ(1) = 1/(P(24,24) + R_MAG + P(1,24)*SH_MAG(2) + P(4,24)*SH_MAG(1) + P(21,24)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) - (2*magD*q1 + 2*magE*q0 - 2*magN*q3)*(P(24,2) + P(1,2)*SH_MAG(2) + P(4,2)*SH_MAG(1) + P(21,2)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,2)*(2*q0*q2 + 2*q1*q3) - P(20,2)*(2*q0*q1 - 2*q2*q3) - P(2,2)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,2)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(8) + SH_MAG(9) - 2*magD*q2)*(P(24,3) + P(1,3)*SH_MAG(2) + P(4,3)*SH_MAG(1) + P(21,3)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,3)*(2*q0*q2 + 2*q1*q3) - P(20,3)*(2*q0*q1 - 2*q2*q3) - P(2,3)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,3)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(2)*(P(24,1) + P(1,1)*SH_MAG(2) + P(4,1)*SH_MAG(1) + P(21,1)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,1)*(2*q0*q2 + 2*q1*q3) - P(20,1)*(2*q0*q1 - 2*q2*q3) - P(2,1)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,1)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + SH_MAG(1)*(P(24,4) + P(1,4)*SH_MAG(2) + P(4,4)*SH_MAG(1) + P(21,4)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,4)*(2*q0*q2 + 2*q1*q3) - P(20,4)*(2*q0*q1 - 2*q2*q3) - P(2,4)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,4)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + (SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7))*(P(24,21) + P(1,21)*SH_MAG(2) + P(4,21)*SH_MAG(1) + P(21,21)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,21)*(2*q0*q2 + 2*q1*q3) - P(20,21)*(2*q0*q1 - 2*q2*q3) - P(2,21)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,21)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(19,24)*(2*q0*q2 + 2*q1*q3) - P(20,24)*(2*q0*q1 - 2*q2*q3) - P(2,24)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + (2*q0*q2 + 2*q1*q3)*(P(24,19) + P(1,19)*SH_MAG(2) + P(4,19)*SH_MAG(1) + P(21,19)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,19)*(2*q0*q2 + 2*q1*q3) - P(20,19)*(2*q0*q1 - 2*q2*q3) - P(2,19)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,19)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) - (2*q0*q1 - 2*q2*q3)*(P(24,20) + P(1,20)*SH_MAG(2) + P(4,20)*SH_MAG(1) + P(21,20)*(SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7)) + P(19,20)*(2*q0*q2 + 2*q1*q3) - P(20,20)*(2*q0*q1 - 2*q2*q3) - P(2,20)*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P(3,20)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2)) + P(3,24)*(SH_MAG(8) + SH_MAG(9) - 2*magD*q2));
        SK_MZ(2) = SH_MAG(4) - SH_MAG(5) - SH_MAG(6) + SH_MAG(7);
        SK_MZ(3) = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
        SK_MZ(4) = SH_MAG(8) + SH_MAG(9) - 2*magD*q2;
        SK_MZ(5) = 2*q0*q1 - 2*q2*q3;
        SK_MZ(6) = 2*q0*q2 + 2*q1*q3;
        K_MAG(1,1) = SK_MZ(1)*(P(1,24) + P(1,1)*SH_MAG(2) + P(1,4)*SH_MAG(1) - P(1,2)*SK_MZ(3) + P(1,3)*SK_MZ(4) + P(1,21)*SK_MZ(2) + P(1,19)*SK_MZ(6) - P(1,20)*SK_MZ(5));
        K_MAG(2,1) = SK_MZ(1)*(P(2,24) + P(2,1)*SH_MAG(2) + P(2,4)*SH_MAG(1) - P(2,2)*SK_MZ(3) + P(2,3)*SK_MZ(4) + P(2,21)*SK_MZ(2) + P(2,19)*SK_MZ(6) - P(2,20)*SK_MZ(5));
        K_MAG(3,1) = SK_MZ(1)*(P(3,24) + P(3,1)*SH_MAG(2) + P(3,4)*SH_MAG(1) - P(3,2)*SK_MZ(3) + P(3,3)*SK_MZ(4) + P(3,21)*SK_MZ(2) + P(3,19)*SK_MZ(6) - P(3,20)*SK_MZ(5));
        K_MAG(4,1) = SK_MZ(1)*(P(4,24) + P(4,1)*SH_MAG(2) + P(4,4)*SH_MAG(1) - P(4,2)*SK_MZ(3) + P(4,3)*SK_MZ(4) + P(4,21)*SK_MZ(2) + P(4,19)*SK_MZ(6) - P(4,20)*SK_MZ(5));
        K_MAG(5,1) = SK_MZ(1)*(P(5,24) + P(5,1)*SH_MAG(2) + P(5,4)*SH_MAG(1) - P(5,2)*SK_MZ(3) + P(5,3)*SK_MZ(4) + P(5,21)*SK_MZ(2) + P(5,19)*SK_MZ(6) - P(5,20)*SK_MZ(5));
        K_MAG(6,1) = SK_MZ(1)*(P(6,24) + P(6,1)*SH_MAG(2) + P(6,4)*SH_MAG(1) - P(6,2)*SK_MZ(3) + P(6,3)*SK_MZ(4) + P(6,21)*SK_MZ(2) + P(6,19)*SK_MZ(6) - P(6,20)*SK_MZ(5));
        K_MAG(7,1) = SK_MZ(1)*(P(7,24) + P(7,1)*SH_MAG(2) + P(7,4)*SH_MAG(1) - P(7,2)*SK_MZ(3) + P(7,3)*SK_MZ(4) + P(7,21)*SK_MZ(2) + P(7,19)*SK_MZ(6) - P(7,20)*SK_MZ(5));
        K_MAG(8,1) = SK_MZ(1)*(P(8,24) + P(8,1)*SH_MAG(2) + P(8,4)*SH_MAG(1) - P(8,2)*SK_MZ(3) + P(8,3)*SK_MZ(4) + P(8,21)*SK_MZ(2) + P(8,19)*SK_MZ(6) - P(8,20)*SK_MZ(5));
        K_MAG(9,1) = SK_MZ(1)*(P(9,24) + P(9,1)*SH_MAG(2) + P(9,4)*SH_MAG(1) - P(9,2)*SK_MZ(3) + P(9,3)*SK_MZ(4) + P(9,21)*SK_MZ(2) + P(9,19)*SK_MZ(6) - P(9,20)*SK_MZ(5));
        K_MAG(10,1) = SK_MZ(1)*(P(10,24) + P(10,1)*SH_MAG(2) + P(10,4)*SH_MAG(1) - P(10,2)*SK_MZ(3) + P(10,3)*SK_MZ(4) + P(10,21)*SK_MZ(2) + P(10,19)*SK_MZ(6) - P(10,20)*SK_MZ(5));
        K_MAG(11,1) = SK_MZ(1)*(P(11,24) + P(11,1)*SH_MAG(2) + P(11,4)*SH_MAG(1) - P(11,2)*SK_MZ(3) + P(11,3)*SK_MZ(4) + P(11,21)*SK_MZ(2) + P(11,19)*SK_MZ(6) - P(11,20)*SK_MZ(5));
        K_MAG(12,1) = SK_MZ(1)*(P(12,24) + P(12,1)*SH_MAG(2) + P(12,4)*SH_MAG(1) - P(12,2)*SK_MZ(3) + P(12,3)*SK_MZ(4) + P(12,21)*SK_MZ(2) + P(12,19)*SK_MZ(6) - P(12,20)*SK_MZ(5));
        K_MAG(13,1) = SK_MZ(1)*(P(13,24) + P(13,1)*SH_MAG(2) + P(13,4)*SH_MAG(1) - P(13,2)*SK_MZ(3) + P(13,3)*SK_MZ(4) + P(13,21)*SK_MZ(2) + P(13,19)*SK_MZ(6) - P(13,20)*SK_MZ(5));
        K_MAG(14,1) = SK_MZ(1)*(P(14,24) + P(14,1)*SH_MAG(2) + P(14,4)*SH_MAG(1) - P(14,2)*SK_MZ(3) + P(14,3)*SK_MZ(4) + P(14,21)*SK_MZ(2) + P(14,19)*SK_MZ(6) - P(14,20)*SK_MZ(5));
        K_MAG(15,1) = SK_MZ(1)*(P(15,24) + P(15,1)*SH_MAG(2) + P(15,4)*SH_MAG(1) - P(15,2)*SK_MZ(3) + P(15,3)*SK_MZ(4) + P(15,21)*SK_MZ(2) + P(15,19)*SK_MZ(6) - P(15,20)*SK_MZ(5));
        K_MAG(16,1) = SK_MZ(1)*(P(16,24) + P(16,1)*SH_MAG(2) + P(16,4)*SH_MAG(1) - P(16,2)*SK_MZ(3) + P(16,3)*SK_MZ(4) + P(16,21)*SK_MZ(2) + P(16,19)*SK_MZ(6) - P(16,20)*SK_MZ(5));
        K_MAG(17,1) = SK_MZ(1)*(P(17,24) + P(17,1)*SH_MAG(2) + P(17,4)*SH_MAG(1) - P(17,2)*SK_MZ(3) + P(17,3)*SK_MZ(4) + P(17,21)*SK_MZ(2) + P(17,19)*SK_MZ(6) - P(17,20)*SK_MZ(5));
        K_MAG(18,1) = SK_MZ(1)*(P(18,24) + P(18,1)*SH_MAG(2) + P(18,4)*SH_MAG(1) - P(18,2)*SK_MZ(3) + P(18,3)*SK_MZ(4) + P(18,21)*SK_MZ(2) + P(18,19)*SK_MZ(6) - P(18,20)*SK_MZ(5));
        K_MAG(19,1) = SK_MZ(1)*(P(19,24) + P(19,1)*SH_MAG(2) + P(19,4)*SH_MAG(1) - P(19,2)*SK_MZ(3) + P(19,3)*SK_MZ(4) + P(19,21)*SK_MZ(2) + P(19,19)*SK_MZ(6) - P(19,20)*SK_MZ(5));
        K_MAG(20,1) = SK_MZ(1)*(P(20,24) + P(20,1)*SH_MAG(2) + P(20,4)*SH_MAG(1) - P(20,2)*SK_MZ(3) + P(20,3)*SK_MZ(4) + P(20,21)*SK_MZ(2) + P(20,19)*SK_MZ(6) - P(20,20)*SK_MZ(5));
        K_MAG(21,1) = SK_MZ(1)*(P(21,24) + P(21,1)*SH_MAG(2) + P(21,4)*SH_MAG(1) - P(21,2)*SK_MZ(3) + P(21,3)*SK_MZ(4) + P(21,21)*SK_MZ(2) + P(21,19)*SK_MZ(6) - P(21,20)*SK_MZ(5));
        K_MAG(22,1) = SK_MZ(1)*(P(22,24) + P(22,1)*SH_MAG(2) + P(22,4)*SH_MAG(1) - P(22,2)*SK_MZ(3) + P(22,3)*SK_MZ(4) + P(22,21)*SK_MZ(2) + P(22,19)*SK_MZ(6) - P(22,20)*SK_MZ(5));
        K_MAG(23,1) = SK_MZ(1)*(P(23,24) + P(23,1)*SH_MAG(2) + P(23,4)*SH_MAG(1) - P(23,2)*SK_MZ(3) + P(23,3)*SK_MZ(4) + P(23,21)*SK_MZ(2) + P(23,19)*SK_MZ(6) - P(23,20)*SK_MZ(5));
        K_MAG(24,1) = SK_MZ(1)*(P(24,24) + P(24,1)*SH_MAG(2) + P(24,4)*SH_MAG(1) - P(24,2)*SK_MZ(3) + P(24,3)*SK_MZ(4) + P(24,21)*SK_MZ(2) + P(24,19)*SK_MZ(6) - P(24,20)*SK_MZ(5));
        _varMagInnov(3) = 1/SK_MZ(1);
    end
    // Calculate the measurement innovation
    _magInnov(obsIndex) = MagPred(obsIndex) - MagData(obsIndex);
    // Check the innovation for consistencey and don't fuse if > 5Sigma
    if ((_magInnov(obsIndex)^2) / _varMagInnov(obsIndex)) < +inf
        xk = K_MAG * _magInnov(obsIndex);
        states = states - xk;
        // normalise the quaternion states
        quatMag = sqrt(states(1)^2 + states(2)^2 + states(3)^2 + states(4)^2);
        if (quatMag > 1e-12)
            states(1:4) = states(1:4) / quatMag;
        end
        // correct the covariance P = (I - K*H)*P
        // take advantage of the empty columns in KH to reduce the
        // number of operations
        KH = K_MAG * H_MAG;
        for i = 1:24
            for j = 1:24
                for k = 1:4
                    KHP(i,j) = KHP(i,j) + KH(i,k)*P(k,j);
                end
                for k = 19:24
                    KHP(i,j) = KHP(i,j) + KH(i,k)*P(k,j);
                end
            end
        end
        P = P - KHP;
    end
    obsIndex = obsIndex + 1;
end
// Force symmetry on the covariance matrix to prevent ill-conditioning
// of the matrix which would cause the filter to blow-up
P = 0.5*(P + transpose(P));
// Set default output for states and covariance
corrP = P;
corrStates = states;
}