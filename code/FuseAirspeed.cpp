/// -*- tab-width: 4; Mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*-
/*
 * Fusion of Position & Velocity Measurements for a 24 State Extended Kalman Navigation Filter
 * Based on https://github.com/priseborough/InertialNav
 * Converted to C++ by Paul Riseborough
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

// run state and covariance update using true airspeed measurement fusion
void FuseAirspeed(
    float nextStates[24], // state output
    float nextP[24][24], // covariance output
    float innovation, // innovation output
    float varInnov, // innovation variance output
    float states[24], // state input
    float P[24][24], // covariance input
    bool FuseData, // boolean true when airspeed data is to be fused
    float VtasMeas, // true airspeed measurement (m/s)
    float StatesAtMeasTime[24], // filter states at the effective measurement time
    bool useAirspeed) // boolean true if airspeed data is being used
{
    float vn;
    float ve;
    float vd;
    float vwn;
    float vwe;
    unsigned short int i;
    unsigned short int j;
    unsigned short int k;
    const float R_TAS = 2.0;
    float SH_TAS[3];
    float SK_TAS[2];
    float H_TAS[24];
    float K_TAS[24];
    float KH[24][24];
    float KHP[24][24];
    bool quatHealth = true;
    float VtasPred;
    float quatMag;

    // Copy required states to local variable names
    vn = StatesAtMeasTime[4];
    ve = StatesAtMeasTime[5];
    vd = StatesAtMeasTime[6];
    vwn = StatesAtMeasTime[16];
    vwe = StatesAtMeasTime[17];

    // Need to check that it is flying before fusing airspeed data
    // Calculate the predicted airspeed
    VtasPred = sqrt((ve - vwe)*(ve - vwe) + (vn - vwn)*(vn - vwn) + vd*vd);
    // Perform fusion of True Airspeed measurement
    if (useAirspeed && FuseData && (VtasPred > 1.0) && (VtasMeas > 8.0))
    {
        // Calculate observation jacobians
        SH_TAS[0] = 1/(sqrt(sq(ve - vwe) + sq(vn - vwn) + sq(vd)));
        SH_TAS[1] = (SH_TAS[0]*(2*ve - 2*vwe))*0.5;
        SH_TAS[2] = (SH_TAS[0]*(2*vn - 2*vwn))*0.5;
        H_TAS[4] = SH_TAS[2];
        H_TAS[5] = SH_TAS[1];
        H_TAS[6] = vd*SH_TAS[0];
        H_TAS[16] = -SH_TAS[2];
        H_TAS[17] = -SH_TAS[1];
        // Calculate Kalman gains
        SK_TAS[0] = 1/(R_TAS + SH_TAS[2]*(P[4][4]*SH_TAS[2] + P[5][4]*SH_TAS[1] - P[16][4]*SH_TAS[2] - P[17][4]*SH_TAS[1] + P[6][4]*vd*SH_TAS[0]) + SH_TAS[1]*(P[4][5]*SH_TAS[2] + P[5][5]*SH_TAS[1] - P[16][5]*SH_TAS[2] - P[17][5]*SH_TAS[1] + P[6][5]*vd*SH_TAS[0]) - SH_TAS[2]*(P[4][16]*SH_TAS[2] + P[5][16]*SH_TAS[1] - P[16][16]*SH_TAS[2] - P[17][16]*SH_TAS[1] + P[6][16]*vd*SH_TAS[0]) - SH_TAS[1]*(P[4][17]*SH_TAS[2] + P[5][17]*SH_TAS[1] - P[16][17]*SH_TAS[2] - P[17][17]*SH_TAS[1] + P[6][17]*vd*SH_TAS[0]) + vd*SH_TAS[0]*(P[4][6]*SH_TAS[2] + P[5][6]*SH_TAS[1] - P[16][6]*SH_TAS[2] - P[17][6]*SH_TAS[1] + P[6][6]*vd*SH_TAS[0]));
        SK_TAS[1] = SH_TAS[1];
        K_TAS[0] = SK_TAS[0]*(P[0][4]*SH_TAS[2] - P[0][16]*SH_TAS[2] + P[0][5]*SK_TAS[1] - P[0][17]*SK_TAS[1] + P[0][6]*vd*SH_TAS[0]);
        K_TAS[1] = SK_TAS[0]*(P[1][4]*SH_TAS[2] - P[1][16]*SH_TAS[2] + P[1][5]*SK_TAS[1] - P[1][17]*SK_TAS[1] + P[1][6]*vd*SH_TAS[0]);
        K_TAS[2] = SK_TAS[0]*(P[2][4]*SH_TAS[2] - P[2][16]*SH_TAS[2] + P[2][5]*SK_TAS[1] - P[2][17]*SK_TAS[1] + P[2][6]*vd*SH_TAS[0]);
        K_TAS[3] = SK_TAS[0]*(P[3][4]*SH_TAS[2] - P[3][16]*SH_TAS[2] + P[3][5]*SK_TAS[1] - P[3][17]*SK_TAS[1] + P[3][6]*vd*SH_TAS[0]);
        K_TAS[4] = SK_TAS[0]*(P[4][4]*SH_TAS[2] - P[4][16]*SH_TAS[2] + P[4][5]*SK_TAS[1] - P[4][17]*SK_TAS[1] + P[4][6]*vd*SH_TAS[0]);
        K_TAS[5] = SK_TAS[0]*(P[5][4]*SH_TAS[2] - P[5][16]*SH_TAS[2] + P[5][5]*SK_TAS[1] - P[5][17]*SK_TAS[1] + P[5][6]*vd*SH_TAS[0]);
        K_TAS[6] = SK_TAS[0]*(P[6][4]*SH_TAS[2] - P[6][16]*SH_TAS[2] + P[6][5]*SK_TAS[1] - P[6][17]*SK_TAS[1] + P[6][6]*vd*SH_TAS[0]);
        K_TAS[7] = SK_TAS[0]*(P[7][4]*SH_TAS[2] - P[7][16]*SH_TAS[2] + P[7][5]*SK_TAS[1] - P[7][17]*SK_TAS[1] + P[7][6]*vd*SH_TAS[0]);
        K_TAS[8] = SK_TAS[0]*(P[8][4]*SH_TAS[2] - P[8][16]*SH_TAS[2] + P[8][5]*SK_TAS[1] - P[8][17]*SK_TAS[1] + P[8][6]*vd*SH_TAS[0]);
        K_TAS[9] = SK_TAS[0]*(P[9][4]*SH_TAS[2] - P[9][16]*SH_TAS[2] + P[9][5]*SK_TAS[1] - P[9][17]*SK_TAS[1] + P[9][6]*vd*SH_TAS[0]);
        K_TAS[10] = SK_TAS[0]*(P[10][4]*SH_TAS[2] - P[10][16]*SH_TAS[2] + P[10][5]*SK_TAS[1] - P[10][17]*SK_TAS[1] + P[10][6]*vd*SH_TAS[0]);
        K_TAS[11] = SK_TAS[0]*(P[11][4]*SH_TAS[2] - P[11][16]*SH_TAS[2] + P[11][5]*SK_TAS[1] - P[11][17]*SK_TAS[1] + P[11][6]*vd*SH_TAS[0]);
        K_TAS[12] = SK_TAS[0]*(P[12][4]*SH_TAS[2] - P[12][16]*SH_TAS[2] + P[12][5]*SK_TAS[1] - P[12][17]*SK_TAS[1] + P[12][6]*vd*SH_TAS[0]);
        K_TAS[13] = SK_TAS[0]*(P[13][4]*SH_TAS[2] - P[13][16]*SH_TAS[2] + P[13][5]*SK_TAS[1] - P[13][17]*SK_TAS[1] + P[13][6]*vd*SH_TAS[0]);
        K_TAS[14] = SK_TAS[0]*(P[14][4]*SH_TAS[2] - P[14][16]*SH_TAS[2] + P[14][5]*SK_TAS[1] - P[14][17]*SK_TAS[1] + P[14][6]*vd*SH_TAS[0]);
        K_TAS[15] = SK_TAS[0]*(P[15][4]*SH_TAS[2] - P[15][16]*SH_TAS[2] + P[15][5]*SK_TAS[1] - P[15][17]*SK_TAS[1] + P[15][6]*vd*SH_TAS[0]);
        K_TAS[16] = SK_TAS[0]*(P[16][4]*SH_TAS[2] - P[16][16]*SH_TAS[2] + P[16][5]*SK_TAS[1] - P[16][17]*SK_TAS[1] + P[16][6]*vd*SH_TAS[0]);
        K_TAS[17] = SK_TAS[0]*(P[17][4]*SH_TAS[2] - P[17][16]*SH_TAS[2] + P[17][5]*SK_TAS[1] - P[17][17]*SK_TAS[1] + P[17][6]*vd*SH_TAS[0]);
        K_TAS[18] = SK_TAS[0]*(P[18][4]*SH_TAS[2] - P[18][16]*SH_TAS[2] + P[18][5]*SK_TAS[1] - P[18][17]*SK_TAS[1] + P[18][6]*vd*SH_TAS[0]);
        K_TAS[19] = SK_TAS[0]*(P[19][4]*SH_TAS[2] - P[19][16]*SH_TAS[2] + P[19][5]*SK_TAS[1] - P[19][17]*SK_TAS[1] + P[19][6]*vd*SH_TAS[0]);
        K_TAS[20] = SK_TAS[0]*(P[20][4]*SH_TAS[2] - P[20][16]*SH_TAS[2] + P[20][5]*SK_TAS[1] - P[20][17]*SK_TAS[1] + P[20][6]*vd*SH_TAS[0]);
        K_TAS[21] = SK_TAS[0]*(P[21][4]*SH_TAS[2] - P[21][16]*SH_TAS[2] + P[21][5]*SK_TAS[1] - P[21][17]*SK_TAS[1] + P[21][6]*vd*SH_TAS[0]);
        K_TAS[22] = SK_TAS[0]*(P[22][4]*SH_TAS[2] - P[22][16]*SH_TAS[2] + P[22][5]*SK_TAS[1] - P[22][17]*SK_TAS[1] + P[22][6]*vd*SH_TAS[0]);
        K_TAS[23] = SK_TAS[0]*(P[23][4]*SH_TAS[2] - P[23][16]*SH_TAS[2] + P[23][5]*SK_TAS[1] - P[23][17]*SK_TAS[1] + P[23][6]*vd*SH_TAS[0]);
        varInnov = 1.0/SK_TAS[0];
        // Calculate the measurement innovation
        innovation = VtasPred - VtasMeas;
        // Check the innovation for consistency and don't fuse if > 5Sigma
        if ((innovation*innovation*SK_TAS[0]) < 25.0)
        {
            // correct the state vector
            for (j= 0; j<=23; j++)
            {
                states[j] = states[j] - K_TAS[j] * innovation;
            }
            // normalise the quaternion states
            quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
            if (quatMag < 0.9) quatHealth = false;
            if (quatMag > 1e-12)
            {
                for (j= 0; j<=3; j++)
                {
                    float quatMagInv = 1.0/quatMag;
                    states[j] = states[j] * quatMagInv;
                }
            }
            // correct the covariance P = (I - K*H)*P
            // take advantage of the empty columns in H to reduce the
            // number of operations
            for (i = 0; i<=23; i++)
            {
                for (j = 4; j<=6; j++)
                {
                    KH[i][j] = K_TAS[i] * H_TAS[j];
                }
                for (j = 16; j<=17; j++)
                {
                    KH[i][j] = K_TAS[i] * H_TAS[j];
                }
            }
            for (i = 0; i<=23; i++)
            {
                for (j = 0; j<=23; j++)
                {
                    for (k = 4; k<=6; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][j] * P[k][j];
                    }
                    for (k = 16; k<=17; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][j] * P[k][j];
                    }
                }
            }
            for (i = 0; i<=23; i++)
            {
                for (j = 0; j<=23; j++)
                {
                    P[i][j] = P[i][j] - KHP[i][j];
                }
            }
        }
        // only update states and convariance if quaternion health OK
        if (quatHealth)
        {
            for (i = 0; i<=23; i++)
                nextStates[i] = states[i];
            {
                for (j = 0; j<=23; j++)
                {
                    nextP[i][j] = P[i][j];
                }
            }
        }
    }
}
