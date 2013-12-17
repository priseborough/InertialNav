/// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-
/*
 * Fusion of Position & Velocity Measurements for a 24 State Extended Kalman Navigation Filter
 * Based on https://github.com/priseborough/InertialNav
 * Converted to C by Paul Riseborough
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// run one state and covariance update step using a combination of NED Vel and Pos measurements
void FuseVelPosNED(
        float nextStates[24], // state output
        float nextP[24][24],  // covariance output
        float innovation[6], // innovation output
        float varInnov[6], // innovation variance output
        float states[24], // state input
        float P[24][24], // covariance input
        float accNavMag, // magnitude of navigation accel (- used to adjust GPS obs variance (m/s^2)
        bool FuseGPSData, // this boolean causes the PosNE and VelNED obs to be fused
        float VelNED[3], // North, East, Down velocity obs (m/s)
        bool useVelD, // this boolean casues the D component of the VelNED vector to be used
        float PosNE[3], // North, East position obs (m)
        float StatesAtGpsTime[24], // States at the effective measurement time for PosNE and VelNED measurements
        bool FuseHgtData, // this boolean causes the HgtMea obs to be fused
        float HgtMea, //  measured height (m)
        float StatesAtHgtTime[24], // States at the effective measurement time for the HgtMea measurement
        bool useAirspeed) // this boolean indicates that airspeed measurements are also being used
{

// declare variables used by fault isolation logic
    float gpsRetryTime = 30.0;
    float gpsRetryTimeNoTAS = 5.0;
    float hgtRetryTime = 30.0;
    float dt = 0.02;
    unsigned int maxVelFailCount;
    unsigned int maxPosFailCount;
    unsigned int maxHgtFailCount;
    static unsigned int velFailCount;
    static unsigned int posFailCount;
    static unsigned int hgtFailCount;
    bool velHealth = false;
    bool posHealth = false;
    bool hgtHealth = false;
    bool quatHealth = true;

// declare variables used to check measurement errors
    float velInnov[3] = {0.0,0.0,0.0};
    float posInnov[2] = {0.0,0.0};
    float hgtInnov = 0.0;

// declare indices used to access arrays
    int stateIndex;
    int obsIndex;
    int i;
    int j;
    int iMax;

// declare variables used by state and covariance update calculations
    float velErr;
    float posErr;
    float R_OBS[6];
    float observation[6];
    float KHP[24][24];
    float SK;
    float K[24];
    float quatMag;
    int startIndex;
    int endIndex;

//Default action - pass through states and covariances
    for (i=0; i<=23; i++)
    {
        nextStates[i] = states[i];
        for (j=0; j<=23; j++)
        {
            nextP[i][j] = P[i][j];
        }
    }

// Form the observation vector
    for (i=0; i<=2; i++) observation[i] = VelNED[i];
    for (i=3; i<=4; i++) observation[i] = PosNE[i];
    observation[5] = -(HgtMea);

// Estimate the GPS Velocity, GPS horiz position and height measurement variances.
    velErr = 0.15*accNavMag; // additional error in GPS velocities caused by manoeuvring
    posErr = 0.15*accNavMag; // additional error in GPS position caused by manoeuvring
    for (i=1; i<=3; i++) R_OBS[i-1] = 0.01 + velErr*velErr;
    for (i=4; i<=5; i++) R_OBS[i-1] = 4.0 + posErr*posErr;
    R_OBS[5] = 4.0;

// Specify the count before forcing use of GPS or height data after invalid
// data. We need to force GPS again sooner without
// airspeed data as the nav velocities will be unconstrained.
    if (useAirspeed)
    {
        maxVelFailCount = ceil(gpsRetryTime/dt);
        maxPosFailCount = maxVelFailCount;
    }
    else
    {
        maxVelFailCount = ceil(gpsRetryTimeNoTAS/dt);
        maxPosFailCount = maxVelFailCount;
    }
    maxHgtFailCount = ceil(hgtRetryTime/dt);

// Perform sequential fusion of GPS measurements. This assumes that the
// errors in the different velocity and position components are
// uncorrelated which is not true, however in the absence of covariance
// data from the GPS receiver it is the only assumption we can make
// so we might as well take advantage of the computational efficiencies
// associated with sequential fusion
    if (FuseGPSData || FuseHgtData)
    {
        // Set innovation variances to zero default
        for (i = 0; i<=5; i++)
        {
             varInnov[i] = 0.0;
        }
        // calculate innovations and check GPS data validity against limits using a 5-sigma test
        if (FuseGPSData)
        {
            if (useVelD) iMax = 2; else iMax = 1;
            for (i = 0; i<=iMax; i++)
            {
                velInnov[i] = StatesAtGpsTime[i+4] - VelNED[i];
                stateIndex = 4 + i;
                varInnov[i] = P[stateIndex][stateIndex] + R_OBS[i];
            }
            if ((velInnov[0]*velInnov[0] + velInnov[1]*velInnov[1] + velInnov[2]*velInnov[2]) < 25.0*(varInnov[0] + varInnov[1] + varInnov[2]) || (velFailCount > maxVelFailCount))
            {
                velHealth = true;
                velFailCount = 0;
            }
            else
            {
                velHealth = false;
                velFailCount = velFailCount + 1;
            }
            posInnov[0] = StatesAtGpsTime[7] - PosNE[0];
            posInnov[1] = StatesAtGpsTime[8] - PosNE[1];
            varInnov[3] = P[7][7] + R_OBS[3];
            varInnov[4] = P[8][8] + R_OBS[4];
            if ((posInnov[0]*posInnov[0] + posInnov[1]*posInnov[1]) < 100.0*(varInnov[3] + varInnov[4]) || (posFailCount > maxPosFailCount))
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
        if (FuseHgtData)
        {
            hgtInnov = StatesAtHgtTime[9] + HgtMea;
            varInnov[5] = P[9][9] + R_OBS[5];
            if ((hgtInnov*hgtInnov) < 25.0*varInnov[5] || (hgtFailCount > maxHgtFailCount))
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
        if (FuseGPSData)
        {
            startIndex = 0;
        }
        else
        {
            startIndex = 5;
        }
        if (FuseHgtData)
        {
            endIndex = 5;
        }
        else
        {
            endIndex = 4;
        }
        // Fuse measurements sequentially
        // Set quaternion health to good as starting default
        quatHealth = true;
        for (obsIndex= startIndex; i<=endIndex; i++)
        {
            // Apply data health checks
            if ((velHealth && (obsIndex >= 0 && obsIndex <= 2)) || (posHealth && (obsIndex == 3 || obsIndex == 4)) || (hgtHealth && (obsIndex == 5)))
            {
                stateIndex = 4 + obsIndex;
                // Calculate the measurement innovation, using states from a
                // different time coordinate if fusing height data
                if (obsIndex == 5)
                {
                    innovation[obsIndex] = StatesAtHgtTime[stateIndex] - observation[obsIndex];
                }
                else
                {
                    innovation[obsIndex] = StatesAtGpsTime[stateIndex] - observation[obsIndex];
                }
                // Calculate the Kalman Gain
                // Calculate innovation variances - also used for data logging
                varInnov[obsIndex] = P[stateIndex][stateIndex] + R_OBS[obsIndex];
                SK = 1.0/varInnov[obsIndex];
                // Check the innovation for consistency and don't fuse if > TBD Sigma
                // Currently disabled for testing
                for (i= 0; i<=23; i++)
                {
                    K[i] = P[i][stateIndex]*SK;
                }
                // Calculate state corrections and re-normalise the quaternions
                for (i = 0; i<=23; i++)
                {
                    states[i] = states[i] - K[i] * innovation[obsIndex];
                }
                quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
                if (quatMag < 0.9) quatHealth = false; // normalisation will have introduced significant errors
                if (quatMag > 1e-12) // divide by  0 protection
                {
                    for (i = 0; i<=3; i++)
                    {
                        states[i] = states[i] / quatMag;
                    }
                }
                // Update the covariance - take advantage of direct observation of a
                // single state at index = stateIndex to reduce computations
                // Optimised implementation of standard equation P = (I - K*H)*P;
                for (i= 0; i<=23; i++)
                {
                    for (j= 0; j<=23; j++)
                    {
                        KHP[i][j] = K[i] * P[stateIndex][j];
                    }
                }
                for (i= 0; i<=23; i++)
                {
                    for (j= 0; j<=23; j++)
                    {
                        P[i][j] = P[i][j] - KHP[i][j];
                    }
                }
            }
        }
    }
// don't update states and covariances if the quaternion health is bad
    if (quatHealth)
    {
        for (i= 0; i<=23; i++)
        {
            nextStates[i] = states[i];
            for (j= 0; j<=23; j++)
            {
                nextP[i][j] = P[i][j];
            }
        }
    }
}
