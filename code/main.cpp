#include <math.h>
#include <stdint.h>
#include <stdio.h>

class Vector3f
{
private:
public:
    float x;
    float y;
    float z;
};

class Mat3f
{
private:
public:
    Vector3f x;
    Vector3f y;
    Vector3f z;
};

// overload + operator to provide a vector addition
Vector3f operator+( Vector3f vecIn1, Vector3f vecIn2)
{
    Vector3f vecOut;
    vecOut.x = vecIn1.x + vecIn2.x;
    vecOut.y = vecIn1.y + vecIn2.y;
    vecOut.z = vecIn1.z + vecIn2.z;
    return vecOut;
}

// overload - operator to provide a vector subtraction
Vector3f operator-( Vector3f vecIn1, Vector3f vecIn2)
{
    Vector3f vecOut;
    vecOut.x = vecIn1.x - vecIn2.x;
    vecOut.y = vecIn1.y - vecIn2.y;
    vecOut.z = vecIn1.z - vecIn2.z;
    return vecOut;
}

// overload * operator to provide a matrix vector product
Vector3f operator*( Mat3f matIn, Vector3f vecIn)
{
    Vector3f vecOut;
    vecOut.x = matIn.x.x*vecIn.x + matIn.x.y*vecIn.y + matIn.x.z*vecIn.z;
    vecOut.y = matIn.y.x*vecIn.x + matIn.y.y*vecIn.y + matIn.y.z*vecIn.z;
    vecOut.z = matIn.x.x*vecIn.x + matIn.z.y*vecIn.y + matIn.z.z*vecIn.z;
    return vecOut;
}

// overload * operator to provide a vector cross product
Vector3f operator*( Vector3f vecIn1, Vector3f vecIn2)
{
    Vector3f vecOut;
    vecOut.x = vecIn1.y*vecIn2.z - vecIn1.z*vecIn2.y;
    vecOut.y = vecIn1.z*vecIn2.x - vecIn1.x*vecIn2.z;
    vecOut.z = vecIn1.x*vecIn2.y - vecIn1.y*vecIn2.x;
    return vecOut;
}

// overload * operator to provide a vector scaler product
Vector3f operator*(Vector3f vecIn1, float sclIn1)
{
    Vector3f vecOut;
    vecOut.x = vecIn1.x * sclIn1;
    vecOut.y = vecIn1.y * sclIn1;
    vecOut.z = vecIn1.z * sclIn1;
    return vecOut;
}

// overload * operator to provide a vector scaler product
Vector3f operator*(float sclIn1, Vector3f vecIn1)
{
    Vector3f vecOut;
    vecOut.x = vecIn1.x * sclIn1;
    vecOut.y = vecIn1.y * sclIn1;
    vecOut.z = vecIn1.z * sclIn1;
    return vecOut;
}

void  UpdateStrapdownEquationsNED(
    Vector3f correctedDelAng, // delta angles about the xyz body axes corrected for errors (rad)
    Vector3f correctedDelVel, // delta velocities along the XYZ body axes corrected for errors (m/s)
    float accNavMag, // magnitude of navigation accel (- used to adjust GPS obs variance (m/s^2)
    Vector3f earthRateNED, // earths angular rate vector in NED (rad/s)
    Vector3f angRate, // angular rate vector in XYZ body axes measured by the IMU (rad/s)
    Vector3f accel, // acceleration vector in XYZ body axes measured by the IMU (m/s^2)
    float dt); // time lapsed since the last IMU measurement (sec)

void CovariancePrediction(
    Vector3f correctedDelAng, // delta angles about the xyz body axes corrected for errors (rad)
    Vector3f correctedDelVel, // delta velocities along the XYZ body axes corrected for errors (m/s)
    float dt, // time lapsed since the last covariance prediction step (sec)
    bool onGround, // boolean true when the flight vehicle is on the ground (not flying)
    bool useAirspeed, // boolean true if airspeed data is being used
    bool useCompass); // boolean true if magnetometer data is being used

void FuseVelPosNED(
    float innovation[6], // innovation output
    float varInnov[6], // innovation variance output
    float accNavMag, // magnitude of navigation accel (- used to adjust GPS obs variance (m/s^2)
    bool FuseGPSData, // this boolean causes the PosNE and VelNED obs to be fused
    float VelNED[3], // North, East, Down velocity obs (m/s)
    bool useVelD, // this boolean casues the D component of the VelNED vector to be used
    float PosNE[3], // North, East position obs (m)
    float StatesAtGpsTime[24], // States at the effective measurement time for PosNE and VelNED measurements
    bool FuseHgtData, // this boolean causes the HgtMea obs to be fused
    float HgtMea, //  measured height (m)
    float StatesAtHgtTime[24], // States at the effective measurement time for the HgtMea measurement
    bool useAirspeed); // this boolean indicates that airspeed measurements are also being used

void FuseMagnetometer(
    float innovation[6], // innovation output
    float varInnov[6], // innovation variance output
    bool FuseData, // boolean true when magnetometer data is to be fused
    float MagData[3], // magnetometer flux radings in X,Y,Z body axes
    float StatesAtMeasTime[24], // filter satates at the effective measurement time
    bool useCompass); // boolean true if magnetometer data is being used

void FuseAirspeed(
    float innovation, // innovation output
    float varInnov, // innovation variance output
    bool FuseData, // boolean true when airspeed data is to be fused
    float VtasMeas, // true airspeed measurement (m/s)
    float StatesAtMeasTime[24], // filter states at the effective measurement time
    bool useAirspeed); // boolean true if airspeed data is being used

void zeroRows(float covMat[24][24], uint8_t first, uint8_t last);

void zeroCols(float covMat[24][24], uint8_t first, uint8_t last);

float sq(float valIn);

void cross3D(float vecOut[3], float vecIn1[3], float vecIn2[3]);

void quatNorm(float quatOut[4], float quatIn[4]);

// store staes along with system time stamp in msces
void StoreStates(uint32_t msec);

// recall stste vector stored at closest time to the one specified by msec
void RecallStates(float statesForFusion[24], uint32_t msec);

// Global variables
#define GRAVITY_MSS 9.80665
#define deg2rad 0.017453292
#define rad2deg 57.295780
static float P[24][24]; // covariance matrix
static float states[24]; // state matrix
static float storedStates[24][50]; // state vectors stored for the last 50 time steps
uint32_t statetimeStamp[50]; // time stamp for each state vector stored

int main()
{
    uint8_t i;
    uint8_t j;

     // open IMU data file
    FILE * pImuFile;
    pImuFile = fopen ("IMU.txt","r");
    // open magnetometer data file
    FILE * pMagFile;
    pMagFile = fopen ("MAG.txt","r");
    // open GPS data file
    FILE * pGpsFile;
    pGpsFile = fopen ("GPS.txt","r");
    // open AHRS data file
    FILE * pAttFile;
    pAttFile = fopen ("ATT.txt","r");

    // IMU input data variables
    float imuIn;
    float tempImu[8];
    uint32_t IMUframe;
    float dt;
    static uint32_t IMUtime = 0;
    Vector3f angRate;
    Vector3f accel;
    int imuReadStatus;

    // GPS input data variables
    float gpsIn;
    float tempGps[14];
    uint32_t GPSframe; // col 0
    uint8_t GPSstatus; // col 1, 3 = good
    uint32_t GPStime; // col 2
    float gpsCourse; // col 11 * deg2rad
    float gpsGndSpd; // col 10
    float gpsVelD; // col 12
    float gpsLat; // col 6 * deg2rad
    float gpsLon; // col 7 * deg2rad
    float hgt; // col 8

    // Magnetometer input data variables
    float magIn;
    float tempMag[8];
    uint32_t MAGframe;
    uint32_t MAGtime;
    float MagData[3];
    float MagBias[3];

    // AHRS input data variables
    float ahrsIn;
    float ahrsMag[7];
    uint32_t AHRSframe;
    uint32_t AHRStime;
    float ahrsEul[3];

    while (imuReadStatus != EOF)
    {
        for (j=0; j<=7; j++)
        {
            imuReadStatus = fscanf (pImuFile, "%f", &imuIn);
            if (imuReadStatus != EOF) tempImu[j] = imuIn;
        }
        if (imuReadStatus != EOF)
        {
            IMUframe  = tempImu[0];
            dt        = 0.001*(tempImu[1] - IMUtime);
            IMUtime   = tempImu[1];
            angRate.x = tempImu[2];
            angRate.y = tempImu[3];
            angRate.z = tempImu[4];
            accel.x   = tempImu[5];
            accel.y   = tempImu[6];
            accel.z   = tempImu[7];
        }
    }
    fclose (pImuFile);
    fclose (pMagFile);
    fclose (pGpsFile);
    fclose (pAttFile);
}

void CovariancePrediction(
    Vector3f correctedDelAng,
    Vector3f correctedDelVel,
    float dt,
    bool onGround,
    bool useAirspeed,
    bool useCompass)
{
    // constants
    const float pi = 3.1415927;
    const float windVelSigma = dt*0.1; // 0.1 m/s/s
    const float dAngBiasSigma = dt*(0.05/3600.0*pi/180.0);
    const float dVelBiasSigma = dt*(0.01/60.0);
    const float magEarthSigma = dt*(10.0/60.0); // 10 mgauss per minute
    const float magBodySigma  = dt*(100.0/60.0); // 100 mgauss per minute - allow to adapt more rapidly due to changes in airframe
    const float daxCov = sq(dt*(50.0/60.0*pi/180.0));
    const float dayCov = sq(dt*(50.0/60.0*pi/180.0));
    const float dazCov = sq(dt*(50.0/60.0*pi/180.0));
    const float dvxCov = sq(dt*0.5);
    const float dvyCov = sq(dt*0.5);
    const float dvzCov = sq(dt*0.5);

    // integers
    uint8_t i;
    uint8_t j;

    // process noise
    float processNoise[24];
    for (i= 0; i<=9; i++) processNoise[i]  = 1.0e-9;
    for (i=10; i<=12; i++) processNoise[i] = dAngBiasSigma;
    for (i=13; i<=15; i++) processNoise[i] = dVelBiasSigma;
    for (i=16; i<=17; i++) processNoise[i] = windVelSigma;
    for (i=18; i<=20; i++) processNoise[i] = magEarthSigma;
    for (i=21; i<=23; i++) processNoise[i] = magBodySigma;
    for (i= 0; i<=23; i++) processNoise[i] = sq(processNoise[i]);

    // time varying inputs
    float dvx = correctedDelVel.x;
    float dvy = correctedDelVel.y;
    float dvz = correctedDelVel.z;
    float dax = correctedDelAng.x;
    float day = correctedDelAng.y;
    float daz = correctedDelAng.z;

    float q0 = states[0];
    float q1 = states[1];
    float q2 = states[2];
    float q3 = states[3];

    float dax_b = states[10];
    float day_b = states[11];
    float daz_b = states[12];

    float dvx_b = states[13];
    float dvy_b = states[14];
    float dvz_b = states[15];

    // arrays
    float SF[21];
    float SG[8];
    float SQ[11];
    float SPP[13];
    float nextP[24][24];

    // misc
    float temp;

    // Predicted covariance calculation
    SF[0] = dvz - dvz_b;
    SF[1] = dvy - dvy_b;
    SF[2] = dvx - dvx_b;
    SF[3] = 2*q1*SF[2] + 2*q2*SF[1] + 2*q3*SF[0];
    SF[4] = 2*q0*SF[1] - 2*q1*SF[0] + 2*q3*SF[2];
    SF[5] = 2*q0*SF[2] + 2*q2*SF[0] - 2*q3*SF[1];
    SF[6] = day*0.5 - day_b*0.5;
    SF[7] = daz*0.5 - daz_b*0.5;
    SF[8] = dax*0.5 - dax_b*0.5;
    SF[9] = dax_b*0.5 - dax*0.5;
    SF[10] = daz_b*0.5 - daz*0.5;
    SF[11] = day_b*0.5 - day*0.5;
    SF[12] = 2*q1*SF[1];
    SF[13] = 2*q0*SF[0];
    SF[14] = q1*0.5;
    SF[15] = q2*0.5;
    SF[16] = q3*0.5;
    SF[17] = sq(q3);
    SF[18] = sq(q2);
    SF[19] = sq(q1);
    SF[20] = sq(q0);

    SG[0] = q0*0.5;
    SG[1] = sq(q3);
    SG[2] = sq(q2);
    SG[3] = sq(q1);
    SG[4] = sq(q0);
    SG[5] = 2*q2*q3;
    SG[6] = 2*q1*q3;
    SG[7] = 2*q1*q2;

    SQ[0] = dvzCov*(SG[5] - 2*q0*q1)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvyCov*(SG[5] + 2*q0*q1)*(SG[1] - SG[2] + SG[3] - SG[4]) + dvxCov*(SG[6] - 2*q0*q2)*(SG[7] + 2*q0*q3);
    SQ[1] = dvzCov*(SG[6] + 2*q0*q2)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvxCov*(SG[6] - 2*q0*q2)*(SG[1] + SG[2] - SG[3] - SG[4]) + dvyCov*(SG[5] + 2*q0*q1)*(SG[7] - 2*q0*q3);
    SQ[2] = dvzCov*(SG[5] - 2*q0*q1)*(SG[6] + 2*q0*q2) - dvyCov*(SG[7] - 2*q0*q3)*(SG[1] - SG[2] + SG[3] - SG[4]) - dvxCov*(SG[7] + 2*q0*q3)*(SG[1] + SG[2] - SG[3] - SG[4]);
    SQ[3] = (dayCov*q1*SG[0])*0.5 - (dazCov*q1*SG[0])*0.5 - (daxCov*q2*q3)*0.25;
    SQ[4] = (dazCov*q2*SG[0])*0.5 - (daxCov*q2*SG[0])*0.5 - (dayCov*q1*q3)*0.25;
    SQ[5] = (daxCov*q3*SG[0])*0.5 - (dayCov*q3*SG[0])*0.5 - (dazCov*q1*q2)*0.25;
    SQ[6] = (daxCov*q1*q2)*0.25 - (dazCov*q3*SG[0])*0.5 - (dayCov*q1*q2)*0.25;
    SQ[7] = (dazCov*q1*q3)*0.25 - (daxCov*q1*q3)*0.25 - (dayCov*q2*SG[0])*0.5;
    SQ[8] = (dayCov*q2*q3)*0.25 - (daxCov*q1*SG[0])*0.5 - (dazCov*q2*q3)*0.25;
    SQ[9] = sq(SG[0]);
    SQ[10] = sq(q1);

    SPP[0] = SF[12] + SF[13] - 2*q2*SF[2];
    SPP[1] = SF[17] - SF[18] - SF[19] + SF[20];
    SPP[2] = SF[17] - SF[18] + SF[19] - SF[20];
    SPP[3] = SF[17] + SF[18] - SF[19] - SF[20];
    SPP[4] = 2*q0*SF[2] + 2*q2*SF[0] - 2*q3*SF[1];
    SPP[5] = 2*q0*SF[1] - 2*q1*SF[0] + 2*q3*SF[2];
    SPP[6] = 2*q0*q2 - 2*q1*q3;
    SPP[7] = 2*q0*q1 - 2*q2*q3;
    SPP[8] = 2*q0*q3 - 2*q1*q2;
    SPP[9] = 2*q0*q1 + 2*q2*q3;
    SPP[10] = 2*q0*q3 + 2*q1*q2;
    SPP[11] = 2*q0*q2 + 2*q1*q3;
    SPP[12] = SF[16];

    nextP[0][0] = P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12] + (daxCov*SQ[10])*0.25 + SF[9]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SF[11]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) + SF[10]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) + SF[14]*(P[0][10] + P[1][10]*SF[9] + P[2][10]*SF[11] + P[3][10]*SF[10] + P[10][10]*SF[14] + P[11][10]*SF[15] + P[12][10]*SPP[12]) + SF[15]*(P[0][11] + P[1][11]*SF[9] + P[2][11]*SF[11] + P[3][11]*SF[10] + P[10][11]*SF[14] + P[11][11]*SF[15] + P[12][11]*SPP[12]) + SPP[12]*(P[0][12] + P[1][12]*SF[9] + P[2][12]*SF[11] + P[3][12]*SF[10] + P[10][12]*SF[14] + P[11][12]*SF[15] + P[12][12]*SPP[12]) + (dayCov*sq(q2))*0.25 + (dazCov*sq(q3))*0.25;
    nextP[0][1] = P[0][1] + SQ[8] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12] + SF[8]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) + SF[7]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) + SF[11]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) - SF[15]*(P[0][12] + P[1][12]*SF[9] + P[2][12]*SF[11] + P[3][12]*SF[10] + P[10][12]*SF[14] + P[11][12]*SF[15] + P[12][12]*SPP[12]) + SPP[12]*(P[0][11] + P[1][11]*SF[9] + P[2][11]*SF[11] + P[3][11]*SF[10] + P[10][11]*SF[14] + P[11][11]*SF[15] + P[12][11]*SPP[12]) - (q0*(P[0][10] + P[1][10]*SF[9] + P[2][10]*SF[11] + P[3][10]*SF[10] + P[10][10]*SF[14] + P[11][10]*SF[15] + P[12][10]*SPP[12]))*0.5;
    nextP[0][2] = P[0][2] + SQ[7] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12] + SF[6]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) + SF[10]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SF[8]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) + SF[14]*(P[0][12] + P[1][12]*SF[9] + P[2][12]*SF[11] + P[3][12]*SF[10] + P[10][12]*SF[14] + P[11][12]*SF[15] + P[12][12]*SPP[12]) - SPP[12]*(P[0][10] + P[1][10]*SF[9] + P[2][10]*SF[11] + P[3][10]*SF[10] + P[10][10]*SF[14] + P[11][10]*SF[15] + P[12][10]*SPP[12]) - (q0*(P[0][11] + P[1][11]*SF[9] + P[2][11]*SF[11] + P[3][11]*SF[10] + P[10][11]*SF[14] + P[11][11]*SF[15] + P[12][11]*SPP[12]))*0.5;
    nextP[0][3] = P[0][3] + SQ[6] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12] + SF[7]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) + SF[6]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SF[9]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) + SF[15]*(P[0][10] + P[1][10]*SF[9] + P[2][10]*SF[11] + P[3][10]*SF[10] + P[10][10]*SF[14] + P[11][10]*SF[15] + P[12][10]*SPP[12]) - SF[14]*(P[0][11] + P[1][11]*SF[9] + P[2][11]*SF[11] + P[3][11]*SF[10] + P[10][11]*SF[14] + P[11][11]*SF[15] + P[12][11]*SPP[12]) - (q0*(P[0][12] + P[1][12]*SF[9] + P[2][12]*SF[11] + P[3][12]*SF[10] + P[10][12]*SF[14] + P[11][12]*SF[15] + P[12][12]*SPP[12]))*0.5;
    nextP[0][4] = P[0][4] + P[1][4]*SF[9] + P[2][4]*SF[11] + P[3][4]*SF[10] + P[10][4]*SF[14] + P[11][4]*SF[15] + P[12][4]*SPP[12] + SF[5]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) + SF[3]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SPP[0]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) - SPP[5]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) + SPP[3]*(P[0][13] + P[1][13]*SF[9] + P[2][13]*SF[11] + P[3][13]*SF[10] + P[10][13]*SF[14] + P[11][13]*SF[15] + P[12][13]*SPP[12]) - SPP[11]*(P[0][15] + P[1][15]*SF[9] + P[2][15]*SF[11] + P[3][15]*SF[10] + P[10][15]*SF[14] + P[11][15]*SF[15] + P[12][15]*SPP[12]) + (2*q0*q3 - 2*q1*q2)*(P[0][14] + P[1][14]*SF[9] + P[2][14]*SF[11] + P[3][14]*SF[10] + P[10][14]*SF[14] + P[11][14]*SF[15] + P[12][14]*SPP[12]);
    nextP[0][5] = P[0][5] + P[1][5]*SF[9] + P[2][5]*SF[11] + P[3][5]*SF[10] + P[10][5]*SF[14] + P[11][5]*SF[15] + P[12][5]*SPP[12] + SF[4]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) + SF[3]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) + SF[5]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) - SPP[0]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SPP[2]*(P[0][14] + P[1][14]*SF[9] + P[2][14]*SF[11] + P[3][14]*SF[10] + P[10][14]*SF[14] + P[11][14]*SF[15] + P[12][14]*SPP[12]) - SPP[10]*(P[0][13] + P[1][13]*SF[9] + P[2][13]*SF[11] + P[3][13]*SF[10] + P[10][13]*SF[14] + P[11][13]*SF[15] + P[12][13]*SPP[12]) + (2*q0*q1 - 2*q2*q3)*(P[0][15] + P[1][15]*SF[9] + P[2][15]*SF[11] + P[3][15]*SF[10] + P[10][15]*SF[14] + P[11][15]*SF[15] + P[12][15]*SPP[12]);
    nextP[0][6] = P[0][6] + P[1][6]*SF[9] + P[2][6]*SF[11] + P[3][6]*SF[10] + P[10][6]*SF[14] + P[11][6]*SF[15] + P[12][6]*SPP[12] + SF[4]*(P[0][1] + P[1][1]*SF[9] + P[2][1]*SF[11] + P[3][1]*SF[10] + P[10][1]*SF[14] + P[11][1]*SF[15] + P[12][1]*SPP[12]) + SF[3]*(P[0][3] + P[1][3]*SF[9] + P[2][3]*SF[11] + P[3][3]*SF[10] + P[10][3]*SF[14] + P[11][3]*SF[15] + P[12][3]*SPP[12]) + SPP[0]*(P[0][0] + P[1][0]*SF[9] + P[2][0]*SF[11] + P[3][0]*SF[10] + P[10][0]*SF[14] + P[11][0]*SF[15] + P[12][0]*SPP[12]) - SPP[4]*(P[0][2] + P[1][2]*SF[9] + P[2][2]*SF[11] + P[3][2]*SF[10] + P[10][2]*SF[14] + P[11][2]*SF[15] + P[12][2]*SPP[12]) + SPP[6]*(P[0][13] + P[1][13]*SF[9] + P[2][13]*SF[11] + P[3][13]*SF[10] + P[10][13]*SF[14] + P[11][13]*SF[15] + P[12][13]*SPP[12]) - SPP[1]*(P[0][15] + P[1][15]*SF[9] + P[2][15]*SF[11] + P[3][15]*SF[10] + P[10][15]*SF[14] + P[11][15]*SF[15] + P[12][15]*SPP[12]) - SPP[9]*(P[0][14] + P[1][14]*SF[9] + P[2][14]*SF[11] + P[3][14]*SF[10] + P[10][14]*SF[14] + P[11][14]*SF[15] + P[12][14]*SPP[12]);
    nextP[0][7] = P[0][7] + P[1][7]*SF[9] + P[2][7]*SF[11] + P[3][7]*SF[10] + P[10][7]*SF[14] + P[11][7]*SF[15] + P[12][7]*SPP[12] + dt*(P[0][4] + P[1][4]*SF[9] + P[2][4]*SF[11] + P[3][4]*SF[10] + P[10][4]*SF[14] + P[11][4]*SF[15] + P[12][4]*SPP[12]);
    nextP[0][8] = P[0][8] + P[1][8]*SF[9] + P[2][8]*SF[11] + P[3][8]*SF[10] + P[10][8]*SF[14] + P[11][8]*SF[15] + P[12][8]*SPP[12] + dt*(P[0][5] + P[1][5]*SF[9] + P[2][5]*SF[11] + P[3][5]*SF[10] + P[10][5]*SF[14] + P[11][5]*SF[15] + P[12][5]*SPP[12]);
    nextP[0][9] = P[0][9] + P[1][9]*SF[9] + P[2][9]*SF[11] + P[3][9]*SF[10] + P[10][9]*SF[14] + P[11][9]*SF[15] + P[12][9]*SPP[12] + dt*(P[0][6] + P[1][6]*SF[9] + P[2][6]*SF[11] + P[3][6]*SF[10] + P[10][6]*SF[14] + P[11][6]*SF[15] + P[12][6]*SPP[12]);
    nextP[0][10] = P[0][10] + P[1][10]*SF[9] + P[2][10]*SF[11] + P[3][10]*SF[10] + P[10][10]*SF[14] + P[11][10]*SF[15] + P[12][10]*SPP[12];
    nextP[0][11] = P[0][11] + P[1][11]*SF[9] + P[2][11]*SF[11] + P[3][11]*SF[10] + P[10][11]*SF[14] + P[11][11]*SF[15] + P[12][11]*SPP[12];
    nextP[0][12] = P[0][12] + P[1][12]*SF[9] + P[2][12]*SF[11] + P[3][12]*SF[10] + P[10][12]*SF[14] + P[11][12]*SF[15] + P[12][12]*SPP[12];
    nextP[0][13] = P[0][13] + P[1][13]*SF[9] + P[2][13]*SF[11] + P[3][13]*SF[10] + P[10][13]*SF[14] + P[11][13]*SF[15] + P[12][13]*SPP[12];
    nextP[0][14] = P[0][14] + P[1][14]*SF[9] + P[2][14]*SF[11] + P[3][14]*SF[10] + P[10][14]*SF[14] + P[11][14]*SF[15] + P[12][14]*SPP[12];
    nextP[0][15] = P[0][15] + P[1][15]*SF[9] + P[2][15]*SF[11] + P[3][15]*SF[10] + P[10][15]*SF[14] + P[11][15]*SF[15] + P[12][15]*SPP[12];
    nextP[0][16] = P[0][16] + P[1][16]*SF[9] + P[2][16]*SF[11] + P[3][16]*SF[10] + P[10][16]*SF[14] + P[11][16]*SF[15] + P[12][16]*SPP[12];
    nextP[0][17] = P[0][17] + P[1][17]*SF[9] + P[2][17]*SF[11] + P[3][17]*SF[10] + P[10][17]*SF[14] + P[11][17]*SF[15] + P[12][17]*SPP[12];
    nextP[0][18] = P[0][18] + P[1][18]*SF[9] + P[2][18]*SF[11] + P[3][18]*SF[10] + P[10][18]*SF[14] + P[11][18]*SF[15] + P[12][18]*SPP[12];
    nextP[0][19] = P[0][19] + P[1][19]*SF[9] + P[2][19]*SF[11] + P[3][19]*SF[10] + P[10][19]*SF[14] + P[11][19]*SF[15] + P[12][19]*SPP[12];
    nextP[0][20] = P[0][20] + P[1][20]*SF[9] + P[2][20]*SF[11] + P[3][20]*SF[10] + P[10][20]*SF[14] + P[11][20]*SF[15] + P[12][20]*SPP[12];
    nextP[0][21] = P[0][21] + P[1][21]*SF[9] + P[2][21]*SF[11] + P[3][21]*SF[10] + P[10][21]*SF[14] + P[11][21]*SF[15] + P[12][21]*SPP[12];
    nextP[0][22] = P[0][22] + P[1][22]*SF[9] + P[2][22]*SF[11] + P[3][22]*SF[10] + P[10][22]*SF[14] + P[11][22]*SF[15] + P[12][22]*SPP[12];
    nextP[0][23] = P[0][23] + P[1][23]*SF[9] + P[2][23]*SF[11] + P[3][23]*SF[10] + P[10][23]*SF[14] + P[11][23]*SF[15] + P[12][23]*SPP[12];
    nextP[1][0] = P[1][0] + SQ[8] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5 + SF[9]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SF[11]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) + SF[10]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) + SF[14]*(P[1][10] + P[0][10]*SF[8] + P[2][10]*SF[7] + P[3][10]*SF[11] - P[12][10]*SF[15] + P[11][10]*SPP[12] - (P[10][10]*q0)*0.5) + SF[15]*(P[1][11] + P[0][11]*SF[8] + P[2][11]*SF[7] + P[3][11]*SF[11] - P[12][11]*SF[15] + P[11][11]*SPP[12] - (P[10][11]*q0)*0.5) + SPP[12]*(P[1][12] + P[0][12]*SF[8] + P[2][12]*SF[7] + P[3][12]*SF[11] - P[12][12]*SF[15] + P[11][12]*SPP[12] - (P[10][12]*q0)*0.5);
    nextP[1][1] = P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] + daxCov*SQ[9] - (P[10][1]*q0)*0.5 + SF[8]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) + SF[7]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) + SF[11]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) - SF[15]*(P[1][12] + P[0][12]*SF[8] + P[2][12]*SF[7] + P[3][12]*SF[11] - P[12][12]*SF[15] + P[11][12]*SPP[12] - (P[10][12]*q0)*0.5) + SPP[12]*(P[1][11] + P[0][11]*SF[8] + P[2][11]*SF[7] + P[3][11]*SF[11] - P[12][11]*SF[15] + P[11][11]*SPP[12] - (P[10][11]*q0)*0.5) + (dayCov*sq(q3))*0.25 + (dazCov*sq(q2))*0.25 - (q0*(P[1][10] + P[0][10]*SF[8] + P[2][10]*SF[7] + P[3][10]*SF[11] - P[12][10]*SF[15] + P[11][10]*SPP[12] - (P[10][10]*q0)*0.5))*0.5;
    nextP[1][2] = P[1][2] + SQ[5] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5 + SF[6]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) + SF[10]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SF[8]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) + SF[14]*(P[1][12] + P[0][12]*SF[8] + P[2][12]*SF[7] + P[3][12]*SF[11] - P[12][12]*SF[15] + P[11][12]*SPP[12] - (P[10][12]*q0)*0.5) - SPP[12]*(P[1][10] + P[0][10]*SF[8] + P[2][10]*SF[7] + P[3][10]*SF[11] - P[12][10]*SF[15] + P[11][10]*SPP[12] - (P[10][10]*q0)*0.5) - (q0*(P[1][11] + P[0][11]*SF[8] + P[2][11]*SF[7] + P[3][11]*SF[11] - P[12][11]*SF[15] + P[11][11]*SPP[12] - (P[10][11]*q0)*0.5))*0.5;
    nextP[1][3] = P[1][3] + SQ[4] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5 + SF[7]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) + SF[6]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SF[9]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) + SF[15]*(P[1][10] + P[0][10]*SF[8] + P[2][10]*SF[7] + P[3][10]*SF[11] - P[12][10]*SF[15] + P[11][10]*SPP[12] - (P[10][10]*q0)*0.5) - SF[14]*(P[1][11] + P[0][11]*SF[8] + P[2][11]*SF[7] + P[3][11]*SF[11] - P[12][11]*SF[15] + P[11][11]*SPP[12] - (P[10][11]*q0)*0.5) - (q0*(P[1][12] + P[0][12]*SF[8] + P[2][12]*SF[7] + P[3][12]*SF[11] - P[12][12]*SF[15] + P[11][12]*SPP[12] - (P[10][12]*q0)*0.5))*0.5;
    nextP[1][4] = P[1][4] + P[0][4]*SF[8] + P[2][4]*SF[7] + P[3][4]*SF[11] - P[12][4]*SF[15] + P[11][4]*SPP[12] - (P[10][4]*q0)*0.5 + SF[5]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) + SF[3]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SPP[0]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) - SPP[5]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) + SPP[3]*(P[1][13] + P[0][13]*SF[8] + P[2][13]*SF[7] + P[3][13]*SF[11] - P[12][13]*SF[15] + P[11][13]*SPP[12] - (P[10][13]*q0)*0.5) + SPP[8]*(P[1][14] + P[0][14]*SF[8] + P[2][14]*SF[7] + P[3][14]*SF[11] - P[12][14]*SF[15] + P[11][14]*SPP[12] - (P[10][14]*q0)*0.5) - SPP[11]*(P[1][15] + P[0][15]*SF[8] + P[2][15]*SF[7] + P[3][15]*SF[11] - P[12][15]*SF[15] + P[11][15]*SPP[12] - (P[10][15]*q0)*0.5);
    nextP[1][5] = P[1][5] + P[0][5]*SF[8] + P[2][5]*SF[7] + P[3][5]*SF[11] - P[12][5]*SF[15] + P[11][5]*SPP[12] - (P[10][5]*q0)*0.5 + SF[4]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) + SF[3]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) + SF[5]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) - SPP[0]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SPP[2]*(P[1][14] + P[0][14]*SF[8] + P[2][14]*SF[7] + P[3][14]*SF[11] - P[12][14]*SF[15] + P[11][14]*SPP[12] - (P[10][14]*q0)*0.5) - SPP[10]*(P[1][13] + P[0][13]*SF[8] + P[2][13]*SF[7] + P[3][13]*SF[11] - P[12][13]*SF[15] + P[11][13]*SPP[12] - (P[10][13]*q0)*0.5) + SPP[7]*(P[1][15] + P[0][15]*SF[8] + P[2][15]*SF[7] + P[3][15]*SF[11] - P[12][15]*SF[15] + P[11][15]*SPP[12] - (P[10][15]*q0)*0.5);
    nextP[1][6] = P[1][6] + P[0][6]*SF[8] + P[2][6]*SF[7] + P[3][6]*SF[11] - P[12][6]*SF[15] + P[11][6]*SPP[12] - (P[10][6]*q0)*0.5 + SF[4]*(P[1][1] + P[0][1]*SF[8] + P[2][1]*SF[7] + P[3][1]*SF[11] - P[12][1]*SF[15] + P[11][1]*SPP[12] - (P[10][1]*q0)*0.5) + SF[3]*(P[1][3] + P[0][3]*SF[8] + P[2][3]*SF[7] + P[3][3]*SF[11] - P[12][3]*SF[15] + P[11][3]*SPP[12] - (P[10][3]*q0)*0.5) + SPP[0]*(P[1][0] + P[0][0]*SF[8] + P[2][0]*SF[7] + P[3][0]*SF[11] - P[12][0]*SF[15] + P[11][0]*SPP[12] - (P[10][0]*q0)*0.5) - SPP[4]*(P[1][2] + P[0][2]*SF[8] + P[2][2]*SF[7] + P[3][2]*SF[11] - P[12][2]*SF[15] + P[11][2]*SPP[12] - (P[10][2]*q0)*0.5) + SPP[6]*(P[1][13] + P[0][13]*SF[8] + P[2][13]*SF[7] + P[3][13]*SF[11] - P[12][13]*SF[15] + P[11][13]*SPP[12] - (P[10][13]*q0)*0.5) - SPP[1]*(P[1][15] + P[0][15]*SF[8] + P[2][15]*SF[7] + P[3][15]*SF[11] - P[12][15]*SF[15] + P[11][15]*SPP[12] - (P[10][15]*q0)*0.5) - SPP[9]*(P[1][14] + P[0][14]*SF[8] + P[2][14]*SF[7] + P[3][14]*SF[11] - P[12][14]*SF[15] + P[11][14]*SPP[12] - (P[10][14]*q0)*0.5);
    nextP[1][7] = P[1][7] + P[0][7]*SF[8] + P[2][7]*SF[7] + P[3][7]*SF[11] - P[12][7]*SF[15] + P[11][7]*SPP[12] - (P[10][7]*q0)*0.5 + dt*(P[1][4] + P[0][4]*SF[8] + P[2][4]*SF[7] + P[3][4]*SF[11] - P[12][4]*SF[15] + P[11][4]*SPP[12] - (P[10][4]*q0)*0.5);
    nextP[1][8] = P[1][8] + P[0][8]*SF[8] + P[2][8]*SF[7] + P[3][8]*SF[11] - P[12][8]*SF[15] + P[11][8]*SPP[12] - (P[10][8]*q0)*0.5 + dt*(P[1][5] + P[0][5]*SF[8] + P[2][5]*SF[7] + P[3][5]*SF[11] - P[12][5]*SF[15] + P[11][5]*SPP[12] - (P[10][5]*q0)*0.5);
    nextP[1][9] = P[1][9] + P[0][9]*SF[8] + P[2][9]*SF[7] + P[3][9]*SF[11] - P[12][9]*SF[15] + P[11][9]*SPP[12] - (P[10][9]*q0)*0.5 + dt*(P[1][6] + P[0][6]*SF[8] + P[2][6]*SF[7] + P[3][6]*SF[11] - P[12][6]*SF[15] + P[11][6]*SPP[12] - (P[10][6]*q0)*0.5);
    nextP[1][10] = P[1][10] + P[0][10]*SF[8] + P[2][10]*SF[7] + P[3][10]*SF[11] - P[12][10]*SF[15] + P[11][10]*SPP[12] - (P[10][10]*q0)*0.5;
    nextP[1][11] = P[1][11] + P[0][11]*SF[8] + P[2][11]*SF[7] + P[3][11]*SF[11] - P[12][11]*SF[15] + P[11][11]*SPP[12] - (P[10][11]*q0)*0.5;
    nextP[1][12] = P[1][12] + P[0][12]*SF[8] + P[2][12]*SF[7] + P[3][12]*SF[11] - P[12][12]*SF[15] + P[11][12]*SPP[12] - (P[10][12]*q0)*0.5;
    nextP[1][13] = P[1][13] + P[0][13]*SF[8] + P[2][13]*SF[7] + P[3][13]*SF[11] - P[12][13]*SF[15] + P[11][13]*SPP[12] - (P[10][13]*q0)*0.5;
    nextP[1][14] = P[1][14] + P[0][14]*SF[8] + P[2][14]*SF[7] + P[3][14]*SF[11] - P[12][14]*SF[15] + P[11][14]*SPP[12] - (P[10][14]*q0)*0.5;
    nextP[1][15] = P[1][15] + P[0][15]*SF[8] + P[2][15]*SF[7] + P[3][15]*SF[11] - P[12][15]*SF[15] + P[11][15]*SPP[12] - (P[10][15]*q0)*0.5;
    nextP[1][16] = P[1][16] + P[0][16]*SF[8] + P[2][16]*SF[7] + P[3][16]*SF[11] - P[12][16]*SF[15] + P[11][16]*SPP[12] - (P[10][16]*q0)*0.5;
    nextP[1][17] = P[1][17] + P[0][17]*SF[8] + P[2][17]*SF[7] + P[3][17]*SF[11] - P[12][17]*SF[15] + P[11][17]*SPP[12] - (P[10][17]*q0)*0.5;
    nextP[1][18] = P[1][18] + P[0][18]*SF[8] + P[2][18]*SF[7] + P[3][18]*SF[11] - P[12][18]*SF[15] + P[11][18]*SPP[12] - (P[10][18]*q0)*0.5;
    nextP[1][19] = P[1][19] + P[0][19]*SF[8] + P[2][19]*SF[7] + P[3][19]*SF[11] - P[12][19]*SF[15] + P[11][19]*SPP[12] - (P[10][19]*q0)*0.5;
    nextP[1][20] = P[1][20] + P[0][20]*SF[8] + P[2][20]*SF[7] + P[3][20]*SF[11] - P[12][20]*SF[15] + P[11][20]*SPP[12] - (P[10][20]*q0)*0.5;
    nextP[1][21] = P[1][21] + P[0][21]*SF[8] + P[2][21]*SF[7] + P[3][21]*SF[11] - P[12][21]*SF[15] + P[11][21]*SPP[12] - (P[10][21]*q0)*0.5;
    nextP[1][22] = P[1][22] + P[0][22]*SF[8] + P[2][22]*SF[7] + P[3][22]*SF[11] - P[12][22]*SF[15] + P[11][22]*SPP[12] - (P[10][22]*q0)*0.5;
    nextP[1][23] = P[1][23] + P[0][23]*SF[8] + P[2][23]*SF[7] + P[3][23]*SF[11] - P[12][23]*SF[15] + P[11][23]*SPP[12] - (P[10][23]*q0)*0.5;
    nextP[2][0] = P[2][0] + SQ[7] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5 + SF[9]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SF[11]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) + SF[10]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) + SF[14]*(P[2][10] + P[0][10]*SF[6] + P[1][10]*SF[10] + P[3][10]*SF[8] + P[12][10]*SF[14] - P[10][10]*SPP[12] - (P[11][10]*q0)*0.5) + SF[15]*(P[2][11] + P[0][11]*SF[6] + P[1][11]*SF[10] + P[3][11]*SF[8] + P[12][11]*SF[14] - P[10][11]*SPP[12] - (P[11][11]*q0)*0.5) + SPP[12]*(P[2][12] + P[0][12]*SF[6] + P[1][12]*SF[10] + P[3][12]*SF[8] + P[12][12]*SF[14] - P[10][12]*SPP[12] - (P[11][12]*q0)*0.5);
    nextP[2][1] = P[2][1] + SQ[5] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5 + SF[8]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) + SF[7]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) + SF[11]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) - SF[15]*(P[2][12] + P[0][12]*SF[6] + P[1][12]*SF[10] + P[3][12]*SF[8] + P[12][12]*SF[14] - P[10][12]*SPP[12] - (P[11][12]*q0)*0.5) + SPP[12]*(P[2][11] + P[0][11]*SF[6] + P[1][11]*SF[10] + P[3][11]*SF[8] + P[12][11]*SF[14] - P[10][11]*SPP[12] - (P[11][11]*q0)*0.5) - (q0*(P[2][10] + P[0][10]*SF[6] + P[1][10]*SF[10] + P[3][10]*SF[8] + P[12][10]*SF[14] - P[10][10]*SPP[12] - (P[11][10]*q0)*0.5))*0.5;
    nextP[2][2] = P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] + dayCov*SQ[9] + (dazCov*SQ[10])*0.25 - (P[11][2]*q0)*0.5 + SF[6]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) + SF[10]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SF[8]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) + SF[14]*(P[2][12] + P[0][12]*SF[6] + P[1][12]*SF[10] + P[3][12]*SF[8] + P[12][12]*SF[14] - P[10][12]*SPP[12] - (P[11][12]*q0)*0.5) - SPP[12]*(P[2][10] + P[0][10]*SF[6] + P[1][10]*SF[10] + P[3][10]*SF[8] + P[12][10]*SF[14] - P[10][10]*SPP[12] - (P[11][10]*q0)*0.5) + (daxCov*sq(q3))*0.25 - (q0*(P[2][11] + P[0][11]*SF[6] + P[1][11]*SF[10] + P[3][11]*SF[8] + P[12][11]*SF[14] - P[10][11]*SPP[12] - (P[11][11]*q0)*0.5))*0.5;
    nextP[2][3] = P[2][3] + SQ[3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5 + SF[7]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) + SF[6]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SF[9]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) + SF[15]*(P[2][10] + P[0][10]*SF[6] + P[1][10]*SF[10] + P[3][10]*SF[8] + P[12][10]*SF[14] - P[10][10]*SPP[12] - (P[11][10]*q0)*0.5) - SF[14]*(P[2][11] + P[0][11]*SF[6] + P[1][11]*SF[10] + P[3][11]*SF[8] + P[12][11]*SF[14] - P[10][11]*SPP[12] - (P[11][11]*q0)*0.5) - (q0*(P[2][12] + P[0][12]*SF[6] + P[1][12]*SF[10] + P[3][12]*SF[8] + P[12][12]*SF[14] - P[10][12]*SPP[12] - (P[11][12]*q0)*0.5))*0.5;
    nextP[2][4] = P[2][4] + P[0][4]*SF[6] + P[1][4]*SF[10] + P[3][4]*SF[8] + P[12][4]*SF[14] - P[10][4]*SPP[12] - (P[11][4]*q0)*0.5 + SF[5]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) + SF[3]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SPP[0]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) - SPP[5]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) + SPP[3]*(P[2][13] + P[0][13]*SF[6] + P[1][13]*SF[10] + P[3][13]*SF[8] + P[12][13]*SF[14] - P[10][13]*SPP[12] - (P[11][13]*q0)*0.5) + SPP[8]*(P[2][14] + P[0][14]*SF[6] + P[1][14]*SF[10] + P[3][14]*SF[8] + P[12][14]*SF[14] - P[10][14]*SPP[12] - (P[11][14]*q0)*0.5) - SPP[11]*(P[2][15] + P[0][15]*SF[6] + P[1][15]*SF[10] + P[3][15]*SF[8] + P[12][15]*SF[14] - P[10][15]*SPP[12] - (P[11][15]*q0)*0.5);
    nextP[2][5] = P[2][5] + P[0][5]*SF[6] + P[1][5]*SF[10] + P[3][5]*SF[8] + P[12][5]*SF[14] - P[10][5]*SPP[12] - (P[11][5]*q0)*0.5 + SF[4]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) + SF[3]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) + SF[5]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) - SPP[0]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SPP[2]*(P[2][14] + P[0][14]*SF[6] + P[1][14]*SF[10] + P[3][14]*SF[8] + P[12][14]*SF[14] - P[10][14]*SPP[12] - (P[11][14]*q0)*0.5) - SPP[10]*(P[2][13] + P[0][13]*SF[6] + P[1][13]*SF[10] + P[3][13]*SF[8] + P[12][13]*SF[14] - P[10][13]*SPP[12] - (P[11][13]*q0)*0.5) + SPP[7]*(P[2][15] + P[0][15]*SF[6] + P[1][15]*SF[10] + P[3][15]*SF[8] + P[12][15]*SF[14] - P[10][15]*SPP[12] - (P[11][15]*q0)*0.5);
    nextP[2][6] = P[2][6] + P[0][6]*SF[6] + P[1][6]*SF[10] + P[3][6]*SF[8] + P[12][6]*SF[14] - P[10][6]*SPP[12] - (P[11][6]*q0)*0.5 + SF[4]*(P[2][1] + P[0][1]*SF[6] + P[1][1]*SF[10] + P[3][1]*SF[8] + P[12][1]*SF[14] - P[10][1]*SPP[12] - (P[11][1]*q0)*0.5) + SF[3]*(P[2][3] + P[0][3]*SF[6] + P[1][3]*SF[10] + P[3][3]*SF[8] + P[12][3]*SF[14] - P[10][3]*SPP[12] - (P[11][3]*q0)*0.5) + SPP[0]*(P[2][0] + P[0][0]*SF[6] + P[1][0]*SF[10] + P[3][0]*SF[8] + P[12][0]*SF[14] - P[10][0]*SPP[12] - (P[11][0]*q0)*0.5) - SPP[4]*(P[2][2] + P[0][2]*SF[6] + P[1][2]*SF[10] + P[3][2]*SF[8] + P[12][2]*SF[14] - P[10][2]*SPP[12] - (P[11][2]*q0)*0.5) + SPP[6]*(P[2][13] + P[0][13]*SF[6] + P[1][13]*SF[10] + P[3][13]*SF[8] + P[12][13]*SF[14] - P[10][13]*SPP[12] - (P[11][13]*q0)*0.5) - SPP[1]*(P[2][15] + P[0][15]*SF[6] + P[1][15]*SF[10] + P[3][15]*SF[8] + P[12][15]*SF[14] - P[10][15]*SPP[12] - (P[11][15]*q0)*0.5) - SPP[9]*(P[2][14] + P[0][14]*SF[6] + P[1][14]*SF[10] + P[3][14]*SF[8] + P[12][14]*SF[14] - P[10][14]*SPP[12] - (P[11][14]*q0)*0.5);
    nextP[2][7] = P[2][7] + P[0][7]*SF[6] + P[1][7]*SF[10] + P[3][7]*SF[8] + P[12][7]*SF[14] - P[10][7]*SPP[12] - (P[11][7]*q0)*0.5 + dt*(P[2][4] + P[0][4]*SF[6] + P[1][4]*SF[10] + P[3][4]*SF[8] + P[12][4]*SF[14] - P[10][4]*SPP[12] - (P[11][4]*q0)*0.5);
    nextP[2][8] = P[2][8] + P[0][8]*SF[6] + P[1][8]*SF[10] + P[3][8]*SF[8] + P[12][8]*SF[14] - P[10][8]*SPP[12] - (P[11][8]*q0)*0.5 + dt*(P[2][5] + P[0][5]*SF[6] + P[1][5]*SF[10] + P[3][5]*SF[8] + P[12][5]*SF[14] - P[10][5]*SPP[12] - (P[11][5]*q0)*0.5);
    nextP[2][9] = P[2][9] + P[0][9]*SF[6] + P[1][9]*SF[10] + P[3][9]*SF[8] + P[12][9]*SF[14] - P[10][9]*SPP[12] - (P[11][9]*q0)*0.5 + dt*(P[2][6] + P[0][6]*SF[6] + P[1][6]*SF[10] + P[3][6]*SF[8] + P[12][6]*SF[14] - P[10][6]*SPP[12] - (P[11][6]*q0)*0.5);
    nextP[2][10] = P[2][10] + P[0][10]*SF[6] + P[1][10]*SF[10] + P[3][10]*SF[8] + P[12][10]*SF[14] - P[10][10]*SPP[12] - (P[11][10]*q0)*0.5;
    nextP[2][11] = P[2][11] + P[0][11]*SF[6] + P[1][11]*SF[10] + P[3][11]*SF[8] + P[12][11]*SF[14] - P[10][11]*SPP[12] - (P[11][11]*q0)*0.5;
    nextP[2][12] = P[2][12] + P[0][12]*SF[6] + P[1][12]*SF[10] + P[3][12]*SF[8] + P[12][12]*SF[14] - P[10][12]*SPP[12] - (P[11][12]*q0)*0.5;
    nextP[2][13] = P[2][13] + P[0][13]*SF[6] + P[1][13]*SF[10] + P[3][13]*SF[8] + P[12][13]*SF[14] - P[10][13]*SPP[12] - (P[11][13]*q0)*0.5;
    nextP[2][14] = P[2][14] + P[0][14]*SF[6] + P[1][14]*SF[10] + P[3][14]*SF[8] + P[12][14]*SF[14] - P[10][14]*SPP[12] - (P[11][14]*q0)*0.5;
    nextP[2][15] = P[2][15] + P[0][15]*SF[6] + P[1][15]*SF[10] + P[3][15]*SF[8] + P[12][15]*SF[14] - P[10][15]*SPP[12] - (P[11][15]*q0)*0.5;
    nextP[2][16] = P[2][16] + P[0][16]*SF[6] + P[1][16]*SF[10] + P[3][16]*SF[8] + P[12][16]*SF[14] - P[10][16]*SPP[12] - (P[11][16]*q0)*0.5;
    nextP[2][17] = P[2][17] + P[0][17]*SF[6] + P[1][17]*SF[10] + P[3][17]*SF[8] + P[12][17]*SF[14] - P[10][17]*SPP[12] - (P[11][17]*q0)*0.5;
    nextP[2][18] = P[2][18] + P[0][18]*SF[6] + P[1][18]*SF[10] + P[3][18]*SF[8] + P[12][18]*SF[14] - P[10][18]*SPP[12] - (P[11][18]*q0)*0.5;
    nextP[2][19] = P[2][19] + P[0][19]*SF[6] + P[1][19]*SF[10] + P[3][19]*SF[8] + P[12][19]*SF[14] - P[10][19]*SPP[12] - (P[11][19]*q0)*0.5;
    nextP[2][20] = P[2][20] + P[0][20]*SF[6] + P[1][20]*SF[10] + P[3][20]*SF[8] + P[12][20]*SF[14] - P[10][20]*SPP[12] - (P[11][20]*q0)*0.5;
    nextP[2][21] = P[2][21] + P[0][21]*SF[6] + P[1][21]*SF[10] + P[3][21]*SF[8] + P[12][21]*SF[14] - P[10][21]*SPP[12] - (P[11][21]*q0)*0.5;
    nextP[2][22] = P[2][22] + P[0][22]*SF[6] + P[1][22]*SF[10] + P[3][22]*SF[8] + P[12][22]*SF[14] - P[10][22]*SPP[12] - (P[11][22]*q0)*0.5;
    nextP[2][23] = P[2][23] + P[0][23]*SF[6] + P[1][23]*SF[10] + P[3][23]*SF[8] + P[12][23]*SF[14] - P[10][23]*SPP[12] - (P[11][23]*q0)*0.5;
    nextP[3][0] = P[3][0] + SQ[6] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5 + SF[9]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SF[11]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) + SF[10]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) + SF[14]*(P[3][10] + P[0][10]*SF[7] + P[1][10]*SF[6] + P[2][10]*SF[9] + P[10][10]*SF[15] - P[11][10]*SF[14] - (P[12][10]*q0)*0.5) + SF[15]*(P[3][11] + P[0][11]*SF[7] + P[1][11]*SF[6] + P[2][11]*SF[9] + P[10][11]*SF[15] - P[11][11]*SF[14] - (P[12][11]*q0)*0.5) + SPP[12]*(P[3][12] + P[0][12]*SF[7] + P[1][12]*SF[6] + P[2][12]*SF[9] + P[10][12]*SF[15] - P[11][12]*SF[14] - (P[12][12]*q0)*0.5);
    nextP[3][1] = P[3][1] + SQ[4] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5 + SF[8]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) + SF[7]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) + SF[11]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) - SF[15]*(P[3][12] + P[0][12]*SF[7] + P[1][12]*SF[6] + P[2][12]*SF[9] + P[10][12]*SF[15] - P[11][12]*SF[14] - (P[12][12]*q0)*0.5) + SPP[12]*(P[3][11] + P[0][11]*SF[7] + P[1][11]*SF[6] + P[2][11]*SF[9] + P[10][11]*SF[15] - P[11][11]*SF[14] - (P[12][11]*q0)*0.5) - (q0*(P[3][10] + P[0][10]*SF[7] + P[1][10]*SF[6] + P[2][10]*SF[9] + P[10][10]*SF[15] - P[11][10]*SF[14] - (P[12][10]*q0)*0.5))*0.5;
    nextP[3][2] = P[3][2] + SQ[3] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5 + SF[6]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) + SF[10]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SF[8]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) + SF[14]*(P[3][12] + P[0][12]*SF[7] + P[1][12]*SF[6] + P[2][12]*SF[9] + P[10][12]*SF[15] - P[11][12]*SF[14] - (P[12][12]*q0)*0.5) - SPP[12]*(P[3][10] + P[0][10]*SF[7] + P[1][10]*SF[6] + P[2][10]*SF[9] + P[10][10]*SF[15] - P[11][10]*SF[14] - (P[12][10]*q0)*0.5) - (q0*(P[3][11] + P[0][11]*SF[7] + P[1][11]*SF[6] + P[2][11]*SF[9] + P[10][11]*SF[15] - P[11][11]*SF[14] - (P[12][11]*q0)*0.5))*0.5;
    nextP[3][3] = P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] + (dayCov*SQ[10])*0.25 + dazCov*SQ[9] - (P[12][3]*q0)*0.5 + SF[7]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) + SF[6]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SF[9]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) + SF[15]*(P[3][10] + P[0][10]*SF[7] + P[1][10]*SF[6] + P[2][10]*SF[9] + P[10][10]*SF[15] - P[11][10]*SF[14] - (P[12][10]*q0)*0.5) - SF[14]*(P[3][11] + P[0][11]*SF[7] + P[1][11]*SF[6] + P[2][11]*SF[9] + P[10][11]*SF[15] - P[11][11]*SF[14] - (P[12][11]*q0)*0.5) + (daxCov*sq(q2))*0.25 - (q0*(P[3][12] + P[0][12]*SF[7] + P[1][12]*SF[6] + P[2][12]*SF[9] + P[10][12]*SF[15] - P[11][12]*SF[14] - (P[12][12]*q0)*0.5))*0.5;
    nextP[3][4] = P[3][4] + P[0][4]*SF[7] + P[1][4]*SF[6] + P[2][4]*SF[9] + P[10][4]*SF[15] - P[11][4]*SF[14] - (P[12][4]*q0)*0.5 + SF[5]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) + SF[3]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SPP[0]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) - SPP[5]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) + SPP[3]*(P[3][13] + P[0][13]*SF[7] + P[1][13]*SF[6] + P[2][13]*SF[9] + P[10][13]*SF[15] - P[11][13]*SF[14] - (P[12][13]*q0)*0.5) + SPP[8]*(P[3][14] + P[0][14]*SF[7] + P[1][14]*SF[6] + P[2][14]*SF[9] + P[10][14]*SF[15] - P[11][14]*SF[14] - (P[12][14]*q0)*0.5) - SPP[11]*(P[3][15] + P[0][15]*SF[7] + P[1][15]*SF[6] + P[2][15]*SF[9] + P[10][15]*SF[15] - P[11][15]*SF[14] - (P[12][15]*q0)*0.5);
    nextP[3][5] = P[3][5] + P[0][5]*SF[7] + P[1][5]*SF[6] + P[2][5]*SF[9] + P[10][5]*SF[15] - P[11][5]*SF[14] - (P[12][5]*q0)*0.5 + SF[4]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) + SF[3]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) + SF[5]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) - SPP[0]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SPP[2]*(P[3][14] + P[0][14]*SF[7] + P[1][14]*SF[6] + P[2][14]*SF[9] + P[10][14]*SF[15] - P[11][14]*SF[14] - (P[12][14]*q0)*0.5) - SPP[10]*(P[3][13] + P[0][13]*SF[7] + P[1][13]*SF[6] + P[2][13]*SF[9] + P[10][13]*SF[15] - P[11][13]*SF[14] - (P[12][13]*q0)*0.5) + SPP[7]*(P[3][15] + P[0][15]*SF[7] + P[1][15]*SF[6] + P[2][15]*SF[9] + P[10][15]*SF[15] - P[11][15]*SF[14] - (P[12][15]*q0)*0.5);
    nextP[3][6] = P[3][6] + P[0][6]*SF[7] + P[1][6]*SF[6] + P[2][6]*SF[9] + P[10][6]*SF[15] - P[11][6]*SF[14] - (P[12][6]*q0)*0.5 + SF[4]*(P[3][1] + P[0][1]*SF[7] + P[1][1]*SF[6] + P[2][1]*SF[9] + P[10][1]*SF[15] - P[11][1]*SF[14] - (P[12][1]*q0)*0.5) + SF[3]*(P[3][3] + P[0][3]*SF[7] + P[1][3]*SF[6] + P[2][3]*SF[9] + P[10][3]*SF[15] - P[11][3]*SF[14] - (P[12][3]*q0)*0.5) + SPP[0]*(P[3][0] + P[0][0]*SF[7] + P[1][0]*SF[6] + P[2][0]*SF[9] + P[10][0]*SF[15] - P[11][0]*SF[14] - (P[12][0]*q0)*0.5) - SPP[4]*(P[3][2] + P[0][2]*SF[7] + P[1][2]*SF[6] + P[2][2]*SF[9] + P[10][2]*SF[15] - P[11][2]*SF[14] - (P[12][2]*q0)*0.5) + SPP[6]*(P[3][13] + P[0][13]*SF[7] + P[1][13]*SF[6] + P[2][13]*SF[9] + P[10][13]*SF[15] - P[11][13]*SF[14] - (P[12][13]*q0)*0.5) - SPP[1]*(P[3][15] + P[0][15]*SF[7] + P[1][15]*SF[6] + P[2][15]*SF[9] + P[10][15]*SF[15] - P[11][15]*SF[14] - (P[12][15]*q0)*0.5) - SPP[9]*(P[3][14] + P[0][14]*SF[7] + P[1][14]*SF[6] + P[2][14]*SF[9] + P[10][14]*SF[15] - P[11][14]*SF[14] - (P[12][14]*q0)*0.5);
    nextP[3][7] = P[3][7] + P[0][7]*SF[7] + P[1][7]*SF[6] + P[2][7]*SF[9] + P[10][7]*SF[15] - P[11][7]*SF[14] - (P[12][7]*q0)*0.5 + dt*(P[3][4] + P[0][4]*SF[7] + P[1][4]*SF[6] + P[2][4]*SF[9] + P[10][4]*SF[15] - P[11][4]*SF[14] - (P[12][4]*q0)*0.5);
    nextP[3][8] = P[3][8] + P[0][8]*SF[7] + P[1][8]*SF[6] + P[2][8]*SF[9] + P[10][8]*SF[15] - P[11][8]*SF[14] - (P[12][8]*q0)*0.5 + dt*(P[3][5] + P[0][5]*SF[7] + P[1][5]*SF[6] + P[2][5]*SF[9] + P[10][5]*SF[15] - P[11][5]*SF[14] - (P[12][5]*q0)*0.5);
    nextP[3][9] = P[3][9] + P[0][9]*SF[7] + P[1][9]*SF[6] + P[2][9]*SF[9] + P[10][9]*SF[15] - P[11][9]*SF[14] - (P[12][9]*q0)*0.5 + dt*(P[3][6] + P[0][6]*SF[7] + P[1][6]*SF[6] + P[2][6]*SF[9] + P[10][6]*SF[15] - P[11][6]*SF[14] - (P[12][6]*q0)*0.5);
    nextP[3][10] = P[3][10] + P[0][10]*SF[7] + P[1][10]*SF[6] + P[2][10]*SF[9] + P[10][10]*SF[15] - P[11][10]*SF[14] - (P[12][10]*q0)*0.5;
    nextP[3][11] = P[3][11] + P[0][11]*SF[7] + P[1][11]*SF[6] + P[2][11]*SF[9] + P[10][11]*SF[15] - P[11][11]*SF[14] - (P[12][11]*q0)*0.5;
    nextP[3][12] = P[3][12] + P[0][12]*SF[7] + P[1][12]*SF[6] + P[2][12]*SF[9] + P[10][12]*SF[15] - P[11][12]*SF[14] - (P[12][12]*q0)*0.5;
    nextP[3][13] = P[3][13] + P[0][13]*SF[7] + P[1][13]*SF[6] + P[2][13]*SF[9] + P[10][13]*SF[15] - P[11][13]*SF[14] - (P[12][13]*q0)*0.5;
    nextP[3][14] = P[3][14] + P[0][14]*SF[7] + P[1][14]*SF[6] + P[2][14]*SF[9] + P[10][14]*SF[15] - P[11][14]*SF[14] - (P[12][14]*q0)*0.5;
    nextP[3][15] = P[3][15] + P[0][15]*SF[7] + P[1][15]*SF[6] + P[2][15]*SF[9] + P[10][15]*SF[15] - P[11][15]*SF[14] - (P[12][15]*q0)*0.5;
    nextP[3][16] = P[3][16] + P[0][16]*SF[7] + P[1][16]*SF[6] + P[2][16]*SF[9] + P[10][16]*SF[15] - P[11][16]*SF[14] - (P[12][16]*q0)*0.5;
    nextP[3][17] = P[3][17] + P[0][17]*SF[7] + P[1][17]*SF[6] + P[2][17]*SF[9] + P[10][17]*SF[15] - P[11][17]*SF[14] - (P[12][17]*q0)*0.5;
    nextP[3][18] = P[3][18] + P[0][18]*SF[7] + P[1][18]*SF[6] + P[2][18]*SF[9] + P[10][18]*SF[15] - P[11][18]*SF[14] - (P[12][18]*q0)*0.5;
    nextP[3][19] = P[3][19] + P[0][19]*SF[7] + P[1][19]*SF[6] + P[2][19]*SF[9] + P[10][19]*SF[15] - P[11][19]*SF[14] - (P[12][19]*q0)*0.5;
    nextP[3][20] = P[3][20] + P[0][20]*SF[7] + P[1][20]*SF[6] + P[2][20]*SF[9] + P[10][20]*SF[15] - P[11][20]*SF[14] - (P[12][20]*q0)*0.5;
    nextP[3][21] = P[3][21] + P[0][21]*SF[7] + P[1][21]*SF[6] + P[2][21]*SF[9] + P[10][21]*SF[15] - P[11][21]*SF[14] - (P[12][21]*q0)*0.5;
    nextP[3][22] = P[3][22] + P[0][22]*SF[7] + P[1][22]*SF[6] + P[2][22]*SF[9] + P[10][22]*SF[15] - P[11][22]*SF[14] - (P[12][22]*q0)*0.5;
    nextP[3][23] = P[3][23] + P[0][23]*SF[7] + P[1][23]*SF[6] + P[2][23]*SF[9] + P[10][23]*SF[15] - P[11][23]*SF[14] - (P[12][23]*q0)*0.5;
    nextP[4][0] = P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11] + SF[9]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SF[11]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) + SF[10]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) + SF[14]*(P[4][10] + P[0][10]*SF[5] + P[1][10]*SF[3] + P[2][10]*SPP[0] - P[3][10]*SPP[5] + P[13][10]*SPP[3] + P[14][10]*SPP[8] - P[15][10]*SPP[11]) + SF[15]*(P[4][11] + P[0][11]*SF[5] + P[1][11]*SF[3] + P[2][11]*SPP[0] - P[3][11]*SPP[5] + P[13][11]*SPP[3] + P[14][11]*SPP[8] - P[15][11]*SPP[11]) + SPP[12]*(P[4][12] + P[0][12]*SF[5] + P[1][12]*SF[3] + P[2][12]*SPP[0] - P[3][12]*SPP[5] + P[13][12]*SPP[3] + P[14][12]*SPP[8] - P[15][12]*SPP[11]);
    nextP[4][1] = P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11] + SF[8]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) + SF[7]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) + SF[11]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) - SF[15]*(P[4][12] + P[0][12]*SF[5] + P[1][12]*SF[3] + P[2][12]*SPP[0] - P[3][12]*SPP[5] + P[13][12]*SPP[3] + P[14][12]*SPP[8] - P[15][12]*SPP[11]) + SPP[12]*(P[4][11] + P[0][11]*SF[5] + P[1][11]*SF[3] + P[2][11]*SPP[0] - P[3][11]*SPP[5] + P[13][11]*SPP[3] + P[14][11]*SPP[8] - P[15][11]*SPP[11]) - (q0*(P[4][10] + P[0][10]*SF[5] + P[1][10]*SF[3] + P[2][10]*SPP[0] - P[3][10]*SPP[5] + P[13][10]*SPP[3] + P[14][10]*SPP[8] - P[15][10]*SPP[11]))*0.5;
    nextP[4][2] = P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11] + SF[6]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) + SF[10]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SF[8]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) + SF[14]*(P[4][12] + P[0][12]*SF[5] + P[1][12]*SF[3] + P[2][12]*SPP[0] - P[3][12]*SPP[5] + P[13][12]*SPP[3] + P[14][12]*SPP[8] - P[15][12]*SPP[11]) - SPP[12]*(P[4][10] + P[0][10]*SF[5] + P[1][10]*SF[3] + P[2][10]*SPP[0] - P[3][10]*SPP[5] + P[13][10]*SPP[3] + P[14][10]*SPP[8] - P[15][10]*SPP[11]) - (q0*(P[4][11] + P[0][11]*SF[5] + P[1][11]*SF[3] + P[2][11]*SPP[0] - P[3][11]*SPP[5] + P[13][11]*SPP[3] + P[14][11]*SPP[8] - P[15][11]*SPP[11]))*0.5;
    nextP[4][3] = P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11] + SF[7]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) + SF[6]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SF[9]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) + SF[15]*(P[4][10] + P[0][10]*SF[5] + P[1][10]*SF[3] + P[2][10]*SPP[0] - P[3][10]*SPP[5] + P[13][10]*SPP[3] + P[14][10]*SPP[8] - P[15][10]*SPP[11]) - SF[14]*(P[4][11] + P[0][11]*SF[5] + P[1][11]*SF[3] + P[2][11]*SPP[0] - P[3][11]*SPP[5] + P[13][11]*SPP[3] + P[14][11]*SPP[8] - P[15][11]*SPP[11]) - (q0*(P[4][12] + P[0][12]*SF[5] + P[1][12]*SF[3] + P[2][12]*SPP[0] - P[3][12]*SPP[5] + P[13][12]*SPP[3] + P[14][12]*SPP[8] - P[15][12]*SPP[11]))*0.5;
    nextP[4][4] = P[4][4] + P[0][4]*SF[5] + P[1][4]*SF[3] + P[2][4]*SPP[0] - P[3][4]*SPP[5] + P[13][4]*SPP[3] + P[14][4]*SPP[8] - P[15][4]*SPP[11] + dvyCov*sq(SG[7] - 2*q0*q3) + dvzCov*sq(SG[6] + 2*q0*q2) + SF[5]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) + SF[3]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SPP[0]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) - SPP[5]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) + SPP[3]*(P[4][13] + P[0][13]*SF[5] + P[1][13]*SF[3] + P[2][13]*SPP[0] - P[3][13]*SPP[5] + P[13][13]*SPP[3] + P[14][13]*SPP[8] - P[15][13]*SPP[11]) + SPP[8]*(P[4][14] + P[0][14]*SF[5] + P[1][14]*SF[3] + P[2][14]*SPP[0] - P[3][14]*SPP[5] + P[13][14]*SPP[3] + P[14][14]*SPP[8] - P[15][14]*SPP[11]) - SPP[11]*(P[4][15] + P[0][15]*SF[5] + P[1][15]*SF[3] + P[2][15]*SPP[0] - P[3][15]*SPP[5] + P[13][15]*SPP[3] + P[14][15]*SPP[8] - P[15][15]*SPP[11]) + dvxCov*sq(SG[1] + SG[2] - SG[3] - SG[4]);
    nextP[4][5] = P[4][5] + SQ[2] + P[0][5]*SF[5] + P[1][5]*SF[3] + P[2][5]*SPP[0] - P[3][5]*SPP[5] + P[13][5]*SPP[3] + P[14][5]*SPP[8] - P[15][5]*SPP[11] + SF[4]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) + SF[3]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) + SF[5]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) - SPP[0]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SPP[2]*(P[4][14] + P[0][14]*SF[5] + P[1][14]*SF[3] + P[2][14]*SPP[0] - P[3][14]*SPP[5] + P[13][14]*SPP[3] + P[14][14]*SPP[8] - P[15][14]*SPP[11]) - SPP[10]*(P[4][13] + P[0][13]*SF[5] + P[1][13]*SF[3] + P[2][13]*SPP[0] - P[3][13]*SPP[5] + P[13][13]*SPP[3] + P[14][13]*SPP[8] - P[15][13]*SPP[11]) + SPP[7]*(P[4][15] + P[0][15]*SF[5] + P[1][15]*SF[3] + P[2][15]*SPP[0] - P[3][15]*SPP[5] + P[13][15]*SPP[3] + P[14][15]*SPP[8] - P[15][15]*SPP[11]);
    nextP[4][6] = P[4][6] + SQ[1] + P[0][6]*SF[5] + P[1][6]*SF[3] + P[2][6]*SPP[0] - P[3][6]*SPP[5] + P[13][6]*SPP[3] + P[14][6]*SPP[8] - P[15][6]*SPP[11] + SF[4]*(P[4][1] + P[0][1]*SF[5] + P[1][1]*SF[3] + P[2][1]*SPP[0] - P[3][1]*SPP[5] + P[13][1]*SPP[3] + P[14][1]*SPP[8] - P[15][1]*SPP[11]) + SF[3]*(P[4][3] + P[0][3]*SF[5] + P[1][3]*SF[3] + P[2][3]*SPP[0] - P[3][3]*SPP[5] + P[13][3]*SPP[3] + P[14][3]*SPP[8] - P[15][3]*SPP[11]) + SPP[0]*(P[4][0] + P[0][0]*SF[5] + P[1][0]*SF[3] + P[2][0]*SPP[0] - P[3][0]*SPP[5] + P[13][0]*SPP[3] + P[14][0]*SPP[8] - P[15][0]*SPP[11]) - SPP[4]*(P[4][2] + P[0][2]*SF[5] + P[1][2]*SF[3] + P[2][2]*SPP[0] - P[3][2]*SPP[5] + P[13][2]*SPP[3] + P[14][2]*SPP[8] - P[15][2]*SPP[11]) + SPP[6]*(P[4][13] + P[0][13]*SF[5] + P[1][13]*SF[3] + P[2][13]*SPP[0] - P[3][13]*SPP[5] + P[13][13]*SPP[3] + P[14][13]*SPP[8] - P[15][13]*SPP[11]) - SPP[1]*(P[4][15] + P[0][15]*SF[5] + P[1][15]*SF[3] + P[2][15]*SPP[0] - P[3][15]*SPP[5] + P[13][15]*SPP[3] + P[14][15]*SPP[8] - P[15][15]*SPP[11]) - SPP[9]*(P[4][14] + P[0][14]*SF[5] + P[1][14]*SF[3] + P[2][14]*SPP[0] - P[3][14]*SPP[5] + P[13][14]*SPP[3] + P[14][14]*SPP[8] - P[15][14]*SPP[11]);
    nextP[4][7] = P[4][7] + P[0][7]*SF[5] + P[1][7]*SF[3] + P[2][7]*SPP[0] - P[3][7]*SPP[5] + P[13][7]*SPP[3] + P[14][7]*SPP[8] - P[15][7]*SPP[11] + dt*(P[4][4] + P[0][4]*SF[5] + P[1][4]*SF[3] + P[2][4]*SPP[0] - P[3][4]*SPP[5] + P[13][4]*SPP[3] + P[14][4]*SPP[8] - P[15][4]*SPP[11]);
    nextP[4][8] = P[4][8] + P[0][8]*SF[5] + P[1][8]*SF[3] + P[2][8]*SPP[0] - P[3][8]*SPP[5] + P[13][8]*SPP[3] + P[14][8]*SPP[8] - P[15][8]*SPP[11] + dt*(P[4][5] + P[0][5]*SF[5] + P[1][5]*SF[3] + P[2][5]*SPP[0] - P[3][5]*SPP[5] + P[13][5]*SPP[3] + P[14][5]*SPP[8] - P[15][5]*SPP[11]);
    nextP[4][9] = P[4][9] + P[0][9]*SF[5] + P[1][9]*SF[3] + P[2][9]*SPP[0] - P[3][9]*SPP[5] + P[13][9]*SPP[3] + P[14][9]*SPP[8] - P[15][9]*SPP[11] + dt*(P[4][6] + P[0][6]*SF[5] + P[1][6]*SF[3] + P[2][6]*SPP[0] - P[3][6]*SPP[5] + P[13][6]*SPP[3] + P[14][6]*SPP[8] - P[15][6]*SPP[11]);
    nextP[4][10] = P[4][10] + P[0][10]*SF[5] + P[1][10]*SF[3] + P[2][10]*SPP[0] - P[3][10]*SPP[5] + P[13][10]*SPP[3] + P[14][10]*SPP[8] - P[15][10]*SPP[11];
    nextP[4][11] = P[4][11] + P[0][11]*SF[5] + P[1][11]*SF[3] + P[2][11]*SPP[0] - P[3][11]*SPP[5] + P[13][11]*SPP[3] + P[14][11]*SPP[8] - P[15][11]*SPP[11];
    nextP[4][12] = P[4][12] + P[0][12]*SF[5] + P[1][12]*SF[3] + P[2][12]*SPP[0] - P[3][12]*SPP[5] + P[13][12]*SPP[3] + P[14][12]*SPP[8] - P[15][12]*SPP[11];
    nextP[4][13] = P[4][13] + P[0][13]*SF[5] + P[1][13]*SF[3] + P[2][13]*SPP[0] - P[3][13]*SPP[5] + P[13][13]*SPP[3] + P[14][13]*SPP[8] - P[15][13]*SPP[11];
    nextP[4][14] = P[4][14] + P[0][14]*SF[5] + P[1][14]*SF[3] + P[2][14]*SPP[0] - P[3][14]*SPP[5] + P[13][14]*SPP[3] + P[14][14]*SPP[8] - P[15][14]*SPP[11];
    nextP[4][15] = P[4][15] + P[0][15]*SF[5] + P[1][15]*SF[3] + P[2][15]*SPP[0] - P[3][15]*SPP[5] + P[13][15]*SPP[3] + P[14][15]*SPP[8] - P[15][15]*SPP[11];
    nextP[4][16] = P[4][16] + P[0][16]*SF[5] + P[1][16]*SF[3] + P[2][16]*SPP[0] - P[3][16]*SPP[5] + P[13][16]*SPP[3] + P[14][16]*SPP[8] - P[15][16]*SPP[11];
    nextP[4][17] = P[4][17] + P[0][17]*SF[5] + P[1][17]*SF[3] + P[2][17]*SPP[0] - P[3][17]*SPP[5] + P[13][17]*SPP[3] + P[14][17]*SPP[8] - P[15][17]*SPP[11];
    nextP[4][18] = P[4][18] + P[0][18]*SF[5] + P[1][18]*SF[3] + P[2][18]*SPP[0] - P[3][18]*SPP[5] + P[13][18]*SPP[3] + P[14][18]*SPP[8] - P[15][18]*SPP[11];
    nextP[4][19] = P[4][19] + P[0][19]*SF[5] + P[1][19]*SF[3] + P[2][19]*SPP[0] - P[3][19]*SPP[5] + P[13][19]*SPP[3] + P[14][19]*SPP[8] - P[15][19]*SPP[11];
    nextP[4][20] = P[4][20] + P[0][20]*SF[5] + P[1][20]*SF[3] + P[2][20]*SPP[0] - P[3][20]*SPP[5] + P[13][20]*SPP[3] + P[14][20]*SPP[8] - P[15][20]*SPP[11];
    nextP[4][21] = P[4][21] + P[0][21]*SF[5] + P[1][21]*SF[3] + P[2][21]*SPP[0] - P[3][21]*SPP[5] + P[13][21]*SPP[3] + P[14][21]*SPP[8] - P[15][21]*SPP[11];
    nextP[4][22] = P[4][22] + P[0][22]*SF[5] + P[1][22]*SF[3] + P[2][22]*SPP[0] - P[3][22]*SPP[5] + P[13][22]*SPP[3] + P[14][22]*SPP[8] - P[15][22]*SPP[11];
    nextP[4][23] = P[4][23] + P[0][23]*SF[5] + P[1][23]*SF[3] + P[2][23]*SPP[0] - P[3][23]*SPP[5] + P[13][23]*SPP[3] + P[14][23]*SPP[8] - P[15][23]*SPP[11];
    nextP[5][0] = P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7] + SF[9]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SF[11]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) + SF[10]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) + SF[14]*(P[5][10] + P[0][10]*SF[4] + P[2][10]*SF[3] + P[3][10]*SF[5] - P[1][10]*SPP[0] - P[13][10]*SPP[10] + P[14][10]*SPP[2] + P[15][10]*SPP[7]) + SF[15]*(P[5][11] + P[0][11]*SF[4] + P[2][11]*SF[3] + P[3][11]*SF[5] - P[1][11]*SPP[0] - P[13][11]*SPP[10] + P[14][11]*SPP[2] + P[15][11]*SPP[7]) + SPP[12]*(P[5][12] + P[0][12]*SF[4] + P[2][12]*SF[3] + P[3][12]*SF[5] - P[1][12]*SPP[0] - P[13][12]*SPP[10] + P[14][12]*SPP[2] + P[15][12]*SPP[7]);
    nextP[5][1] = P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7] + SF[8]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) + SF[7]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) + SF[11]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) - SF[15]*(P[5][12] + P[0][12]*SF[4] + P[2][12]*SF[3] + P[3][12]*SF[5] - P[1][12]*SPP[0] - P[13][12]*SPP[10] + P[14][12]*SPP[2] + P[15][12]*SPP[7]) + SPP[12]*(P[5][11] + P[0][11]*SF[4] + P[2][11]*SF[3] + P[3][11]*SF[5] - P[1][11]*SPP[0] - P[13][11]*SPP[10] + P[14][11]*SPP[2] + P[15][11]*SPP[7]) - (q0*(P[5][10] + P[0][10]*SF[4] + P[2][10]*SF[3] + P[3][10]*SF[5] - P[1][10]*SPP[0] - P[13][10]*SPP[10] + P[14][10]*SPP[2] + P[15][10]*SPP[7]))*0.5;
    nextP[5][2] = P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7] + SF[6]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) + SF[10]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SF[8]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) + SF[14]*(P[5][12] + P[0][12]*SF[4] + P[2][12]*SF[3] + P[3][12]*SF[5] - P[1][12]*SPP[0] - P[13][12]*SPP[10] + P[14][12]*SPP[2] + P[15][12]*SPP[7]) - SPP[12]*(P[5][10] + P[0][10]*SF[4] + P[2][10]*SF[3] + P[3][10]*SF[5] - P[1][10]*SPP[0] - P[13][10]*SPP[10] + P[14][10]*SPP[2] + P[15][10]*SPP[7]) - (q0*(P[5][11] + P[0][11]*SF[4] + P[2][11]*SF[3] + P[3][11]*SF[5] - P[1][11]*SPP[0] - P[13][11]*SPP[10] + P[14][11]*SPP[2] + P[15][11]*SPP[7]))*0.5;
    nextP[5][3] = P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7] + SF[7]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) + SF[6]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SF[9]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) + SF[15]*(P[5][10] + P[0][10]*SF[4] + P[2][10]*SF[3] + P[3][10]*SF[5] - P[1][10]*SPP[0] - P[13][10]*SPP[10] + P[14][10]*SPP[2] + P[15][10]*SPP[7]) - SF[14]*(P[5][11] + P[0][11]*SF[4] + P[2][11]*SF[3] + P[3][11]*SF[5] - P[1][11]*SPP[0] - P[13][11]*SPP[10] + P[14][11]*SPP[2] + P[15][11]*SPP[7]) - (q0*(P[5][12] + P[0][12]*SF[4] + P[2][12]*SF[3] + P[3][12]*SF[5] - P[1][12]*SPP[0] - P[13][12]*SPP[10] + P[14][12]*SPP[2] + P[15][12]*SPP[7]))*0.5;
    nextP[5][4] = P[5][4] + SQ[2] + P[0][4]*SF[4] + P[2][4]*SF[3] + P[3][4]*SF[5] - P[1][4]*SPP[0] - P[13][4]*SPP[10] + P[14][4]*SPP[2] + P[15][4]*SPP[7] + SF[5]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) + SF[3]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SPP[0]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) - SPP[5]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) + SPP[3]*(P[5][13] + P[0][13]*SF[4] + P[2][13]*SF[3] + P[3][13]*SF[5] - P[1][13]*SPP[0] - P[13][13]*SPP[10] + P[14][13]*SPP[2] + P[15][13]*SPP[7]) + SPP[8]*(P[5][14] + P[0][14]*SF[4] + P[2][14]*SF[3] + P[3][14]*SF[5] - P[1][14]*SPP[0] - P[13][14]*SPP[10] + P[14][14]*SPP[2] + P[15][14]*SPP[7]) - SPP[11]*(P[5][15] + P[0][15]*SF[4] + P[2][15]*SF[3] + P[3][15]*SF[5] - P[1][15]*SPP[0] - P[13][15]*SPP[10] + P[14][15]*SPP[2] + P[15][15]*SPP[7]);
    nextP[5][5] = P[5][5] + P[0][5]*SF[4] + P[2][5]*SF[3] + P[3][5]*SF[5] - P[1][5]*SPP[0] - P[13][5]*SPP[10] + P[14][5]*SPP[2] + P[15][5]*SPP[7] + dvxCov*sq(SG[7] + 2*q0*q3) + dvzCov*sq(SG[5] - 2*q0*q1) + SF[4]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) + SF[3]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) + SF[5]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) - SPP[0]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SPP[2]*(P[5][14] + P[0][14]*SF[4] + P[2][14]*SF[3] + P[3][14]*SF[5] - P[1][14]*SPP[0] - P[13][14]*SPP[10] + P[14][14]*SPP[2] + P[15][14]*SPP[7]) - SPP[10]*(P[5][13] + P[0][13]*SF[4] + P[2][13]*SF[3] + P[3][13]*SF[5] - P[1][13]*SPP[0] - P[13][13]*SPP[10] + P[14][13]*SPP[2] + P[15][13]*SPP[7]) + SPP[7]*(P[5][15] + P[0][15]*SF[4] + P[2][15]*SF[3] + P[3][15]*SF[5] - P[1][15]*SPP[0] - P[13][15]*SPP[10] + P[14][15]*SPP[2] + P[15][15]*SPP[7]) + dvyCov*sq(SG[1] - SG[2] + SG[3] - SG[4]);
    nextP[5][6] = P[5][6] + SQ[0] + P[0][6]*SF[4] + P[2][6]*SF[3] + P[3][6]*SF[5] - P[1][6]*SPP[0] - P[13][6]*SPP[10] + P[14][6]*SPP[2] + P[15][6]*SPP[7] + SF[4]*(P[5][1] + P[0][1]*SF[4] + P[2][1]*SF[3] + P[3][1]*SF[5] - P[1][1]*SPP[0] - P[13][1]*SPP[10] + P[14][1]*SPP[2] + P[15][1]*SPP[7]) + SF[3]*(P[5][3] + P[0][3]*SF[4] + P[2][3]*SF[3] + P[3][3]*SF[5] - P[1][3]*SPP[0] - P[13][3]*SPP[10] + P[14][3]*SPP[2] + P[15][3]*SPP[7]) + SPP[0]*(P[5][0] + P[0][0]*SF[4] + P[2][0]*SF[3] + P[3][0]*SF[5] - P[1][0]*SPP[0] - P[13][0]*SPP[10] + P[14][0]*SPP[2] + P[15][0]*SPP[7]) - SPP[4]*(P[5][2] + P[0][2]*SF[4] + P[2][2]*SF[3] + P[3][2]*SF[5] - P[1][2]*SPP[0] - P[13][2]*SPP[10] + P[14][2]*SPP[2] + P[15][2]*SPP[7]) + SPP[6]*(P[5][13] + P[0][13]*SF[4] + P[2][13]*SF[3] + P[3][13]*SF[5] - P[1][13]*SPP[0] - P[13][13]*SPP[10] + P[14][13]*SPP[2] + P[15][13]*SPP[7]) - SPP[1]*(P[5][15] + P[0][15]*SF[4] + P[2][15]*SF[3] + P[3][15]*SF[5] - P[1][15]*SPP[0] - P[13][15]*SPP[10] + P[14][15]*SPP[2] + P[15][15]*SPP[7]) - SPP[9]*(P[5][14] + P[0][14]*SF[4] + P[2][14]*SF[3] + P[3][14]*SF[5] - P[1][14]*SPP[0] - P[13][14]*SPP[10] + P[14][14]*SPP[2] + P[15][14]*SPP[7]);
    nextP[5][7] = P[5][7] + P[0][7]*SF[4] + P[2][7]*SF[3] + P[3][7]*SF[5] - P[1][7]*SPP[0] - P[13][7]*SPP[10] + P[14][7]*SPP[2] + P[15][7]*SPP[7] + dt*(P[5][4] + P[0][4]*SF[4] + P[2][4]*SF[3] + P[3][4]*SF[5] - P[1][4]*SPP[0] - P[13][4]*SPP[10] + P[14][4]*SPP[2] + P[15][4]*SPP[7]);
    nextP[5][8] = P[5][8] + P[0][8]*SF[4] + P[2][8]*SF[3] + P[3][8]*SF[5] - P[1][8]*SPP[0] - P[13][8]*SPP[10] + P[14][8]*SPP[2] + P[15][8]*SPP[7] + dt*(P[5][5] + P[0][5]*SF[4] + P[2][5]*SF[3] + P[3][5]*SF[5] - P[1][5]*SPP[0] - P[13][5]*SPP[10] + P[14][5]*SPP[2] + P[15][5]*SPP[7]);
    nextP[5][9] = P[5][9] + P[0][9]*SF[4] + P[2][9]*SF[3] + P[3][9]*SF[5] - P[1][9]*SPP[0] - P[13][9]*SPP[10] + P[14][9]*SPP[2] + P[15][9]*SPP[7] + dt*(P[5][6] + P[0][6]*SF[4] + P[2][6]*SF[3] + P[3][6]*SF[5] - P[1][6]*SPP[0] - P[13][6]*SPP[10] + P[14][6]*SPP[2] + P[15][6]*SPP[7]);
    nextP[5][10] = P[5][10] + P[0][10]*SF[4] + P[2][10]*SF[3] + P[3][10]*SF[5] - P[1][10]*SPP[0] - P[13][10]*SPP[10] + P[14][10]*SPP[2] + P[15][10]*SPP[7];
    nextP[5][11] = P[5][11] + P[0][11]*SF[4] + P[2][11]*SF[3] + P[3][11]*SF[5] - P[1][11]*SPP[0] - P[13][11]*SPP[10] + P[14][11]*SPP[2] + P[15][11]*SPP[7];
    nextP[5][12] = P[5][12] + P[0][12]*SF[4] + P[2][12]*SF[3] + P[3][12]*SF[5] - P[1][12]*SPP[0] - P[13][12]*SPP[10] + P[14][12]*SPP[2] + P[15][12]*SPP[7];
    nextP[5][13] = P[5][13] + P[0][13]*SF[4] + P[2][13]*SF[3] + P[3][13]*SF[5] - P[1][13]*SPP[0] - P[13][13]*SPP[10] + P[14][13]*SPP[2] + P[15][13]*SPP[7];
    nextP[5][14] = P[5][14] + P[0][14]*SF[4] + P[2][14]*SF[3] + P[3][14]*SF[5] - P[1][14]*SPP[0] - P[13][14]*SPP[10] + P[14][14]*SPP[2] + P[15][14]*SPP[7];
    nextP[5][15] = P[5][15] + P[0][15]*SF[4] + P[2][15]*SF[3] + P[3][15]*SF[5] - P[1][15]*SPP[0] - P[13][15]*SPP[10] + P[14][15]*SPP[2] + P[15][15]*SPP[7];
    nextP[5][16] = P[5][16] + P[0][16]*SF[4] + P[2][16]*SF[3] + P[3][16]*SF[5] - P[1][16]*SPP[0] - P[13][16]*SPP[10] + P[14][16]*SPP[2] + P[15][16]*SPP[7];
    nextP[5][17] = P[5][17] + P[0][17]*SF[4] + P[2][17]*SF[3] + P[3][17]*SF[5] - P[1][17]*SPP[0] - P[13][17]*SPP[10] + P[14][17]*SPP[2] + P[15][17]*SPP[7];
    nextP[5][18] = P[5][18] + P[0][18]*SF[4] + P[2][18]*SF[3] + P[3][18]*SF[5] - P[1][18]*SPP[0] - P[13][18]*SPP[10] + P[14][18]*SPP[2] + P[15][18]*SPP[7];
    nextP[5][19] = P[5][19] + P[0][19]*SF[4] + P[2][19]*SF[3] + P[3][19]*SF[5] - P[1][19]*SPP[0] - P[13][19]*SPP[10] + P[14][19]*SPP[2] + P[15][19]*SPP[7];
    nextP[5][20] = P[5][20] + P[0][20]*SF[4] + P[2][20]*SF[3] + P[3][20]*SF[5] - P[1][20]*SPP[0] - P[13][20]*SPP[10] + P[14][20]*SPP[2] + P[15][20]*SPP[7];
    nextP[5][21] = P[5][21] + P[0][21]*SF[4] + P[2][21]*SF[3] + P[3][21]*SF[5] - P[1][21]*SPP[0] - P[13][21]*SPP[10] + P[14][21]*SPP[2] + P[15][21]*SPP[7];
    nextP[5][22] = P[5][22] + P[0][22]*SF[4] + P[2][22]*SF[3] + P[3][22]*SF[5] - P[1][22]*SPP[0] - P[13][22]*SPP[10] + P[14][22]*SPP[2] + P[15][22]*SPP[7];
    nextP[5][23] = P[5][23] + P[0][23]*SF[4] + P[2][23]*SF[3] + P[3][23]*SF[5] - P[1][23]*SPP[0] - P[13][23]*SPP[10] + P[14][23]*SPP[2] + P[15][23]*SPP[7];
    nextP[6][0] = P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1] + SF[9]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SF[11]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) + SF[10]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) + SF[14]*(P[6][10] + P[1][10]*SF[4] + P[3][10]*SF[3] + P[0][10]*SPP[0] - P[2][10]*SPP[4] + P[13][10]*SPP[6] - P[14][10]*SPP[9] - P[15][10]*SPP[1]) + SF[15]*(P[6][11] + P[1][11]*SF[4] + P[3][11]*SF[3] + P[0][11]*SPP[0] - P[2][11]*SPP[4] + P[13][11]*SPP[6] - P[14][11]*SPP[9] - P[15][11]*SPP[1]) + SPP[12]*(P[6][12] + P[1][12]*SF[4] + P[3][12]*SF[3] + P[0][12]*SPP[0] - P[2][12]*SPP[4] + P[13][12]*SPP[6] - P[14][12]*SPP[9] - P[15][12]*SPP[1]);
    nextP[6][1] = P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1] + SF[8]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) + SF[7]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) + SF[11]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) - SF[15]*(P[6][12] + P[1][12]*SF[4] + P[3][12]*SF[3] + P[0][12]*SPP[0] - P[2][12]*SPP[4] + P[13][12]*SPP[6] - P[14][12]*SPP[9] - P[15][12]*SPP[1]) + SPP[12]*(P[6][11] + P[1][11]*SF[4] + P[3][11]*SF[3] + P[0][11]*SPP[0] - P[2][11]*SPP[4] + P[13][11]*SPP[6] - P[14][11]*SPP[9] - P[15][11]*SPP[1]) - (q0*(P[6][10] + P[1][10]*SF[4] + P[3][10]*SF[3] + P[0][10]*SPP[0] - P[2][10]*SPP[4] + P[13][10]*SPP[6] - P[14][10]*SPP[9] - P[15][10]*SPP[1]))*0.5;
    nextP[6][2] = P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1] + SF[6]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) + SF[10]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SF[8]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) + SF[14]*(P[6][12] + P[1][12]*SF[4] + P[3][12]*SF[3] + P[0][12]*SPP[0] - P[2][12]*SPP[4] + P[13][12]*SPP[6] - P[14][12]*SPP[9] - P[15][12]*SPP[1]) - SPP[12]*(P[6][10] + P[1][10]*SF[4] + P[3][10]*SF[3] + P[0][10]*SPP[0] - P[2][10]*SPP[4] + P[13][10]*SPP[6] - P[14][10]*SPP[9] - P[15][10]*SPP[1]) - (q0*(P[6][11] + P[1][11]*SF[4] + P[3][11]*SF[3] + P[0][11]*SPP[0] - P[2][11]*SPP[4] + P[13][11]*SPP[6] - P[14][11]*SPP[9] - P[15][11]*SPP[1]))*0.5;
    nextP[6][3] = P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1] + SF[7]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) + SF[6]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SF[9]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) + SF[15]*(P[6][10] + P[1][10]*SF[4] + P[3][10]*SF[3] + P[0][10]*SPP[0] - P[2][10]*SPP[4] + P[13][10]*SPP[6] - P[14][10]*SPP[9] - P[15][10]*SPP[1]) - SF[14]*(P[6][11] + P[1][11]*SF[4] + P[3][11]*SF[3] + P[0][11]*SPP[0] - P[2][11]*SPP[4] + P[13][11]*SPP[6] - P[14][11]*SPP[9] - P[15][11]*SPP[1]) - (q0*(P[6][12] + P[1][12]*SF[4] + P[3][12]*SF[3] + P[0][12]*SPP[0] - P[2][12]*SPP[4] + P[13][12]*SPP[6] - P[14][12]*SPP[9] - P[15][12]*SPP[1]))*0.5;
    nextP[6][4] = P[6][4] + SQ[1] + P[1][4]*SF[4] + P[3][4]*SF[3] + P[0][4]*SPP[0] - P[2][4]*SPP[4] + P[13][4]*SPP[6] - P[14][4]*SPP[9] - P[15][4]*SPP[1] + SF[5]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) + SF[3]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SPP[0]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) - SPP[5]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) + SPP[3]*(P[6][13] + P[1][13]*SF[4] + P[3][13]*SF[3] + P[0][13]*SPP[0] - P[2][13]*SPP[4] + P[13][13]*SPP[6] - P[14][13]*SPP[9] - P[15][13]*SPP[1]) + SPP[8]*(P[6][14] + P[1][14]*SF[4] + P[3][14]*SF[3] + P[0][14]*SPP[0] - P[2][14]*SPP[4] + P[13][14]*SPP[6] - P[14][14]*SPP[9] - P[15][14]*SPP[1]) - SPP[11]*(P[6][15] + P[1][15]*SF[4] + P[3][15]*SF[3] + P[0][15]*SPP[0] - P[2][15]*SPP[4] + P[13][15]*SPP[6] - P[14][15]*SPP[9] - P[15][15]*SPP[1]);
    nextP[6][5] = P[6][5] + SQ[0] + P[1][5]*SF[4] + P[3][5]*SF[3] + P[0][5]*SPP[0] - P[2][5]*SPP[4] + P[13][5]*SPP[6] - P[14][5]*SPP[9] - P[15][5]*SPP[1] + SF[4]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) + SF[3]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) + SF[5]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) - SPP[0]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SPP[2]*(P[6][14] + P[1][14]*SF[4] + P[3][14]*SF[3] + P[0][14]*SPP[0] - P[2][14]*SPP[4] + P[13][14]*SPP[6] - P[14][14]*SPP[9] - P[15][14]*SPP[1]) - SPP[10]*(P[6][13] + P[1][13]*SF[4] + P[3][13]*SF[3] + P[0][13]*SPP[0] - P[2][13]*SPP[4] + P[13][13]*SPP[6] - P[14][13]*SPP[9] - P[15][13]*SPP[1]) + SPP[7]*(P[6][15] + P[1][15]*SF[4] + P[3][15]*SF[3] + P[0][15]*SPP[0] - P[2][15]*SPP[4] + P[13][15]*SPP[6] - P[14][15]*SPP[9] - P[15][15]*SPP[1]);
    nextP[6][6] = P[6][6] + P[1][6]*SF[4] + P[3][6]*SF[3] + P[0][6]*SPP[0] - P[2][6]*SPP[4] + P[13][6]*SPP[6] - P[14][6]*SPP[9] - P[15][6]*SPP[1] + dvxCov*sq(SG[6] - 2*q0*q2) + dvyCov*sq(SG[5] + 2*q0*q1) + SF[4]*(P[6][1] + P[1][1]*SF[4] + P[3][1]*SF[3] + P[0][1]*SPP[0] - P[2][1]*SPP[4] + P[13][1]*SPP[6] - P[14][1]*SPP[9] - P[15][1]*SPP[1]) + SF[3]*(P[6][3] + P[1][3]*SF[4] + P[3][3]*SF[3] + P[0][3]*SPP[0] - P[2][3]*SPP[4] + P[13][3]*SPP[6] - P[14][3]*SPP[9] - P[15][3]*SPP[1]) + SPP[0]*(P[6][0] + P[1][0]*SF[4] + P[3][0]*SF[3] + P[0][0]*SPP[0] - P[2][0]*SPP[4] + P[13][0]*SPP[6] - P[14][0]*SPP[9] - P[15][0]*SPP[1]) - SPP[4]*(P[6][2] + P[1][2]*SF[4] + P[3][2]*SF[3] + P[0][2]*SPP[0] - P[2][2]*SPP[4] + P[13][2]*SPP[6] - P[14][2]*SPP[9] - P[15][2]*SPP[1]) + SPP[6]*(P[6][13] + P[1][13]*SF[4] + P[3][13]*SF[3] + P[0][13]*SPP[0] - P[2][13]*SPP[4] + P[13][13]*SPP[6] - P[14][13]*SPP[9] - P[15][13]*SPP[1]) - SPP[1]*(P[6][15] + P[1][15]*SF[4] + P[3][15]*SF[3] + P[0][15]*SPP[0] - P[2][15]*SPP[4] + P[13][15]*SPP[6] - P[14][15]*SPP[9] - P[15][15]*SPP[1]) - SPP[9]*(P[6][14] + P[1][14]*SF[4] + P[3][14]*SF[3] + P[0][14]*SPP[0] - P[2][14]*SPP[4] + P[13][14]*SPP[6] - P[14][14]*SPP[9] - P[15][14]*SPP[1]) + dvzCov*sq(SG[1] - SG[2] - SG[3] + SG[4]);
    nextP[6][7] = P[6][7] + P[1][7]*SF[4] + P[3][7]*SF[3] + P[0][7]*SPP[0] - P[2][7]*SPP[4] + P[13][7]*SPP[6] - P[14][7]*SPP[9] - P[15][7]*SPP[1] + dt*(P[6][4] + P[1][4]*SF[4] + P[3][4]*SF[3] + P[0][4]*SPP[0] - P[2][4]*SPP[4] + P[13][4]*SPP[6] - P[14][4]*SPP[9] - P[15][4]*SPP[1]);
    nextP[6][8] = P[6][8] + P[1][8]*SF[4] + P[3][8]*SF[3] + P[0][8]*SPP[0] - P[2][8]*SPP[4] + P[13][8]*SPP[6] - P[14][8]*SPP[9] - P[15][8]*SPP[1] + dt*(P[6][5] + P[1][5]*SF[4] + P[3][5]*SF[3] + P[0][5]*SPP[0] - P[2][5]*SPP[4] + P[13][5]*SPP[6] - P[14][5]*SPP[9] - P[15][5]*SPP[1]);
    nextP[6][9] = P[6][9] + P[1][9]*SF[4] + P[3][9]*SF[3] + P[0][9]*SPP[0] - P[2][9]*SPP[4] + P[13][9]*SPP[6] - P[14][9]*SPP[9] - P[15][9]*SPP[1] + dt*(P[6][6] + P[1][6]*SF[4] + P[3][6]*SF[3] + P[0][6]*SPP[0] - P[2][6]*SPP[4] + P[13][6]*SPP[6] - P[14][6]*SPP[9] - P[15][6]*SPP[1]);
    nextP[6][10] = P[6][10] + P[1][10]*SF[4] + P[3][10]*SF[3] + P[0][10]*SPP[0] - P[2][10]*SPP[4] + P[13][10]*SPP[6] - P[14][10]*SPP[9] - P[15][10]*SPP[1];
    nextP[6][11] = P[6][11] + P[1][11]*SF[4] + P[3][11]*SF[3] + P[0][11]*SPP[0] - P[2][11]*SPP[4] + P[13][11]*SPP[6] - P[14][11]*SPP[9] - P[15][11]*SPP[1];
    nextP[6][12] = P[6][12] + P[1][12]*SF[4] + P[3][12]*SF[3] + P[0][12]*SPP[0] - P[2][12]*SPP[4] + P[13][12]*SPP[6] - P[14][12]*SPP[9] - P[15][12]*SPP[1];
    nextP[6][13] = P[6][13] + P[1][13]*SF[4] + P[3][13]*SF[3] + P[0][13]*SPP[0] - P[2][13]*SPP[4] + P[13][13]*SPP[6] - P[14][13]*SPP[9] - P[15][13]*SPP[1];
    nextP[6][14] = P[6][14] + P[1][14]*SF[4] + P[3][14]*SF[3] + P[0][14]*SPP[0] - P[2][14]*SPP[4] + P[13][14]*SPP[6] - P[14][14]*SPP[9] - P[15][14]*SPP[1];
    nextP[6][15] = P[6][15] + P[1][15]*SF[4] + P[3][15]*SF[3] + P[0][15]*SPP[0] - P[2][15]*SPP[4] + P[13][15]*SPP[6] - P[14][15]*SPP[9] - P[15][15]*SPP[1];
    nextP[6][16] = P[6][16] + P[1][16]*SF[4] + P[3][16]*SF[3] + P[0][16]*SPP[0] - P[2][16]*SPP[4] + P[13][16]*SPP[6] - P[14][16]*SPP[9] - P[15][16]*SPP[1];
    nextP[6][17] = P[6][17] + P[1][17]*SF[4] + P[3][17]*SF[3] + P[0][17]*SPP[0] - P[2][17]*SPP[4] + P[13][17]*SPP[6] - P[14][17]*SPP[9] - P[15][17]*SPP[1];
    nextP[6][18] = P[6][18] + P[1][18]*SF[4] + P[3][18]*SF[3] + P[0][18]*SPP[0] - P[2][18]*SPP[4] + P[13][18]*SPP[6] - P[14][18]*SPP[9] - P[15][18]*SPP[1];
    nextP[6][19] = P[6][19] + P[1][19]*SF[4] + P[3][19]*SF[3] + P[0][19]*SPP[0] - P[2][19]*SPP[4] + P[13][19]*SPP[6] - P[14][19]*SPP[9] - P[15][19]*SPP[1];
    nextP[6][20] = P[6][20] + P[1][20]*SF[4] + P[3][20]*SF[3] + P[0][20]*SPP[0] - P[2][20]*SPP[4] + P[13][20]*SPP[6] - P[14][20]*SPP[9] - P[15][20]*SPP[1];
    nextP[6][21] = P[6][21] + P[1][21]*SF[4] + P[3][21]*SF[3] + P[0][21]*SPP[0] - P[2][21]*SPP[4] + P[13][21]*SPP[6] - P[14][21]*SPP[9] - P[15][21]*SPP[1];
    nextP[6][22] = P[6][22] + P[1][22]*SF[4] + P[3][22]*SF[3] + P[0][22]*SPP[0] - P[2][22]*SPP[4] + P[13][22]*SPP[6] - P[14][22]*SPP[9] - P[15][22]*SPP[1];
    nextP[6][23] = P[6][23] + P[1][23]*SF[4] + P[3][23]*SF[3] + P[0][23]*SPP[0] - P[2][23]*SPP[4] + P[13][23]*SPP[6] - P[14][23]*SPP[9] - P[15][23]*SPP[1];
    nextP[7][0] = P[7][0] + P[4][0]*dt + SF[9]*(P[7][1] + P[4][1]*dt) + SF[11]*(P[7][2] + P[4][2]*dt) + SF[10]*(P[7][3] + P[4][3]*dt) + SF[14]*(P[7][10] + P[4][10]*dt) + SF[15]*(P[7][11] + P[4][11]*dt) + SPP[12]*(P[7][12] + P[4][12]*dt);
    nextP[7][1] = P[7][1] + P[4][1]*dt + SF[8]*(P[7][0] + P[4][0]*dt) + SF[7]*(P[7][2] + P[4][2]*dt) + SF[11]*(P[7][3] + P[4][3]*dt) - SF[15]*(P[7][12] + P[4][12]*dt) + SPP[12]*(P[7][11] + P[4][11]*dt) - (q0*(P[7][10] + P[4][10]*dt))*0.5;
    nextP[7][2] = P[7][2] + P[4][2]*dt + SF[6]*(P[7][0] + P[4][0]*dt) + SF[10]*(P[7][1] + P[4][1]*dt) + SF[8]*(P[7][3] + P[4][3]*dt) + SF[14]*(P[7][12] + P[4][12]*dt) - SPP[12]*(P[7][10] + P[4][10]*dt) - (q0*(P[7][11] + P[4][11]*dt))*0.5;
    nextP[7][3] = P[7][3] + P[4][3]*dt + SF[7]*(P[7][0] + P[4][0]*dt) + SF[6]*(P[7][1] + P[4][1]*dt) + SF[9]*(P[7][2] + P[4][2]*dt) + SF[15]*(P[7][10] + P[4][10]*dt) - SF[14]*(P[7][11] + P[4][11]*dt) - (q0*(P[7][12] + P[4][12]*dt))*0.5;
    nextP[7][4] = P[7][4] + P[4][4]*dt + SF[3]*(P[7][1] + P[4][1]*dt) + SF[5]*(P[7][0] + P[4][0]*dt) + SPP[0]*(P[7][2] + P[4][2]*dt) - SPP[5]*(P[7][3] + P[4][3]*dt) + SPP[3]*(P[7][13] + P[4][13]*dt) + SPP[8]*(P[7][14] + P[4][14]*dt) - SPP[11]*(P[7][15] + P[4][15]*dt);
    nextP[7][5] = P[7][5] + P[4][5]*dt + SF[4]*(P[7][0] + P[4][0]*dt) + SF[3]*(P[7][2] + P[4][2]*dt) + SF[5]*(P[7][3] + P[4][3]*dt) - SPP[0]*(P[7][1] + P[4][1]*dt) + SPP[2]*(P[7][14] + P[4][14]*dt) - SPP[10]*(P[7][13] + P[4][13]*dt) + SPP[7]*(P[7][15] + P[4][15]*dt);
    nextP[7][6] = P[7][6] + P[4][6]*dt + SF[4]*(P[7][1] + P[4][1]*dt) + SF[3]*(P[7][3] + P[4][3]*dt) + SPP[0]*(P[7][0] + P[4][0]*dt) - SPP[4]*(P[7][2] + P[4][2]*dt) - SPP[1]*(P[7][15] + P[4][15]*dt) + SPP[6]*(P[7][13] + P[4][13]*dt) - SPP[9]*(P[7][14] + P[4][14]*dt);
    nextP[7][7] = P[7][7] + P[4][7]*dt + dt*(P[7][4] + P[4][4]*dt);
    nextP[7][8] = P[7][8] + P[4][8]*dt + dt*(P[7][5] + P[4][5]*dt);
    nextP[7][9] = P[7][9] + P[4][9]*dt + dt*(P[7][6] + P[4][6]*dt);
    nextP[7][10] = P[7][10] + P[4][10]*dt;
    nextP[7][11] = P[7][11] + P[4][11]*dt;
    nextP[7][12] = P[7][12] + P[4][12]*dt;
    nextP[7][13] = P[7][13] + P[4][13]*dt;
    nextP[7][14] = P[7][14] + P[4][14]*dt;
    nextP[7][15] = P[7][15] + P[4][15]*dt;
    nextP[7][16] = P[7][16] + P[4][16]*dt;
    nextP[7][17] = P[7][17] + P[4][17]*dt;
    nextP[7][18] = P[7][18] + P[4][18]*dt;
    nextP[7][19] = P[7][19] + P[4][19]*dt;
    nextP[7][20] = P[7][20] + P[4][20]*dt;
    nextP[7][21] = P[7][21] + P[4][21]*dt;
    nextP[7][22] = P[7][22] + P[4][22]*dt;
    nextP[7][23] = P[7][23] + P[4][23]*dt;
    nextP[8][0] = P[8][0] + P[5][0]*dt + SF[9]*(P[8][1] + P[5][1]*dt) + SF[11]*(P[8][2] + P[5][2]*dt) + SF[10]*(P[8][3] + P[5][3]*dt) + SF[14]*(P[8][10] + P[5][10]*dt) + SF[15]*(P[8][11] + P[5][11]*dt) + SPP[12]*(P[8][12] + P[5][12]*dt);
    nextP[8][1] = P[8][1] + P[5][1]*dt + SF[8]*(P[8][0] + P[5][0]*dt) + SF[7]*(P[8][2] + P[5][2]*dt) + SF[11]*(P[8][3] + P[5][3]*dt) - SF[15]*(P[8][12] + P[5][12]*dt) + SPP[12]*(P[8][11] + P[5][11]*dt) - (q0*(P[8][10] + P[5][10]*dt))*0.5;
    nextP[8][2] = P[8][2] + P[5][2]*dt + SF[6]*(P[8][0] + P[5][0]*dt) + SF[10]*(P[8][1] + P[5][1]*dt) + SF[8]*(P[8][3] + P[5][3]*dt) + SF[14]*(P[8][12] + P[5][12]*dt) - SPP[12]*(P[8][10] + P[5][10]*dt) - (q0*(P[8][11] + P[5][11]*dt))*0.5;
    nextP[8][3] = P[8][3] + P[5][3]*dt + SF[7]*(P[8][0] + P[5][0]*dt) + SF[6]*(P[8][1] + P[5][1]*dt) + SF[9]*(P[8][2] + P[5][2]*dt) + SF[15]*(P[8][10] + P[5][10]*dt) - SF[14]*(P[8][11] + P[5][11]*dt) - (q0*(P[8][12] + P[5][12]*dt))*0.5;
    nextP[8][4] = P[8][4] + P[5][4]*dt + SF[3]*(P[8][1] + P[5][1]*dt) + SF[5]*(P[8][0] + P[5][0]*dt) + SPP[0]*(P[8][2] + P[5][2]*dt) - SPP[5]*(P[8][3] + P[5][3]*dt) + SPP[3]*(P[8][13] + P[5][13]*dt) + SPP[8]*(P[8][14] + P[5][14]*dt) - SPP[11]*(P[8][15] + P[5][15]*dt);
    nextP[8][5] = P[8][5] + P[5][5]*dt + SF[4]*(P[8][0] + P[5][0]*dt) + SF[3]*(P[8][2] + P[5][2]*dt) + SF[5]*(P[8][3] + P[5][3]*dt) - SPP[0]*(P[8][1] + P[5][1]*dt) + SPP[2]*(P[8][14] + P[5][14]*dt) - SPP[10]*(P[8][13] + P[5][13]*dt) + SPP[7]*(P[8][15] + P[5][15]*dt);
    nextP[8][6] = P[8][6] + P[5][6]*dt + SF[4]*(P[8][1] + P[5][1]*dt) + SF[3]*(P[8][3] + P[5][3]*dt) + SPP[0]*(P[8][0] + P[5][0]*dt) - SPP[4]*(P[8][2] + P[5][2]*dt) - SPP[1]*(P[8][15] + P[5][15]*dt) + SPP[6]*(P[8][13] + P[5][13]*dt) - SPP[9]*(P[8][14] + P[5][14]*dt);
    nextP[8][7] = P[8][7] + P[5][7]*dt + dt*(P[8][4] + P[5][4]*dt);
    nextP[8][8] = P[8][8] + P[5][8]*dt + dt*(P[8][5] + P[5][5]*dt);
    nextP[8][9] = P[8][9] + P[5][9]*dt + dt*(P[8][6] + P[5][6]*dt);
    nextP[8][10] = P[8][10] + P[5][10]*dt;
    nextP[8][11] = P[8][11] + P[5][11]*dt;
    nextP[8][12] = P[8][12] + P[5][12]*dt;
    nextP[8][13] = P[8][13] + P[5][13]*dt;
    nextP[8][14] = P[8][14] + P[5][14]*dt;
    nextP[8][15] = P[8][15] + P[5][15]*dt;
    nextP[8][16] = P[8][16] + P[5][16]*dt;
    nextP[8][17] = P[8][17] + P[5][17]*dt;
    nextP[8][18] = P[8][18] + P[5][18]*dt;
    nextP[8][19] = P[8][19] + P[5][19]*dt;
    nextP[8][20] = P[8][20] + P[5][20]*dt;
    nextP[8][21] = P[8][21] + P[5][21]*dt;
    nextP[8][22] = P[8][22] + P[5][22]*dt;
    nextP[8][23] = P[8][23] + P[5][23]*dt;
    nextP[9][0] = P[9][0] + P[6][0]*dt + SF[9]*(P[9][1] + P[6][1]*dt) + SF[11]*(P[9][2] + P[6][2]*dt) + SF[10]*(P[9][3] + P[6][3]*dt) + SF[14]*(P[9][10] + P[6][10]*dt) + SF[15]*(P[9][11] + P[6][11]*dt) + SPP[12]*(P[9][12] + P[6][12]*dt);
    nextP[9][1] = P[9][1] + P[6][1]*dt + SF[8]*(P[9][0] + P[6][0]*dt) + SF[7]*(P[9][2] + P[6][2]*dt) + SF[11]*(P[9][3] + P[6][3]*dt) - SF[15]*(P[9][12] + P[6][12]*dt) + SPP[12]*(P[9][11] + P[6][11]*dt) - (q0*(P[9][10] + P[6][10]*dt))*0.5;
    nextP[9][2] = P[9][2] + P[6][2]*dt + SF[6]*(P[9][0] + P[6][0]*dt) + SF[10]*(P[9][1] + P[6][1]*dt) + SF[8]*(P[9][3] + P[6][3]*dt) + SF[14]*(P[9][12] + P[6][12]*dt) - SPP[12]*(P[9][10] + P[6][10]*dt) - (q0*(P[9][11] + P[6][11]*dt))*0.5;
    nextP[9][3] = P[9][3] + P[6][3]*dt + SF[7]*(P[9][0] + P[6][0]*dt) + SF[6]*(P[9][1] + P[6][1]*dt) + SF[9]*(P[9][2] + P[6][2]*dt) + SF[15]*(P[9][10] + P[6][10]*dt) - SF[14]*(P[9][11] + P[6][11]*dt) - (q0*(P[9][12] + P[6][12]*dt))*0.5;
    nextP[9][4] = P[9][4] + P[6][4]*dt + SF[3]*(P[9][1] + P[6][1]*dt) + SF[5]*(P[9][0] + P[6][0]*dt) + SPP[0]*(P[9][2] + P[6][2]*dt) - SPP[5]*(P[9][3] + P[6][3]*dt) + SPP[3]*(P[9][13] + P[6][13]*dt) + SPP[8]*(P[9][14] + P[6][14]*dt) - SPP[11]*(P[9][15] + P[6][15]*dt);
    nextP[9][5] = P[9][5] + P[6][5]*dt + SF[4]*(P[9][0] + P[6][0]*dt) + SF[3]*(P[9][2] + P[6][2]*dt) + SF[5]*(P[9][3] + P[6][3]*dt) - SPP[0]*(P[9][1] + P[6][1]*dt) + SPP[2]*(P[9][14] + P[6][14]*dt) - SPP[10]*(P[9][13] + P[6][13]*dt) + SPP[7]*(P[9][15] + P[6][15]*dt);
    nextP[9][6] = P[9][6] + P[6][6]*dt + SF[4]*(P[9][1] + P[6][1]*dt) + SF[3]*(P[9][3] + P[6][3]*dt) + SPP[0]*(P[9][0] + P[6][0]*dt) - SPP[4]*(P[9][2] + P[6][2]*dt) - SPP[1]*(P[9][15] + P[6][15]*dt) + SPP[6]*(P[9][13] + P[6][13]*dt) - SPP[9]*(P[9][14] + P[6][14]*dt);
    nextP[9][7] = P[9][7] + P[6][7]*dt + dt*(P[9][4] + P[6][4]*dt);
    nextP[9][8] = P[9][8] + P[6][8]*dt + dt*(P[9][5] + P[6][5]*dt);
    nextP[9][9] = P[9][9] + P[6][9]*dt + dt*(P[9][6] + P[6][6]*dt);
    nextP[9][10] = P[9][10] + P[6][10]*dt;
    nextP[9][11] = P[9][11] + P[6][11]*dt;
    nextP[9][12] = P[9][12] + P[6][12]*dt;
    nextP[9][13] = P[9][13] + P[6][13]*dt;
    nextP[9][14] = P[9][14] + P[6][14]*dt;
    nextP[9][15] = P[9][15] + P[6][15]*dt;
    nextP[9][16] = P[9][16] + P[6][16]*dt;
    nextP[9][17] = P[9][17] + P[6][17]*dt;
    nextP[9][18] = P[9][18] + P[6][18]*dt;
    nextP[9][19] = P[9][19] + P[6][19]*dt;
    nextP[9][20] = P[9][20] + P[6][20]*dt;
    nextP[9][21] = P[9][21] + P[6][21]*dt;
    nextP[9][22] = P[9][22] + P[6][22]*dt;
    nextP[9][23] = P[9][23] + P[6][23]*dt;
    nextP[10][0] = P[10][0] + P[10][1]*SF[9] + P[10][2]*SF[11] + P[10][3]*SF[10] + P[10][10]*SF[14] + P[10][11]*SF[15] + P[10][12]*SPP[12];
    nextP[10][1] = P[10][1] + P[10][0]*SF[8] + P[10][2]*SF[7] + P[10][3]*SF[11] - P[10][12]*SF[15] + P[10][11]*SPP[12] - (P[10][10]*q0)*0.5;
    nextP[10][2] = P[10][2] + P[10][0]*SF[6] + P[10][1]*SF[10] + P[10][3]*SF[8] + P[10][12]*SF[14] - P[10][10]*SPP[12] - (P[10][11]*q0)*0.5;
    nextP[10][3] = P[10][3] + P[10][0]*SF[7] + P[10][1]*SF[6] + P[10][2]*SF[9] + P[10][10]*SF[15] - P[10][11]*SF[14] - (P[10][12]*q0)*0.5;
    nextP[10][4] = P[10][4] + P[10][1]*SF[3] + P[10][0]*SF[5] + P[10][2]*SPP[0] - P[10][3]*SPP[5] + P[10][13]*SPP[3] + P[10][14]*SPP[8] - P[10][15]*SPP[11];
    nextP[10][5] = P[10][5] + P[10][0]*SF[4] + P[10][2]*SF[3] + P[10][3]*SF[5] - P[10][1]*SPP[0] + P[10][14]*SPP[2] + P[10][15]*SPP[7] - P[10][13]*SPP[10];
    nextP[10][6] = P[10][6] + P[10][1]*SF[4] + P[10][3]*SF[3] + P[10][0]*SPP[0] - P[10][2]*SPP[4] - P[10][15]*SPP[1] + P[10][13]*SPP[6] - P[10][14]*SPP[9];
    nextP[10][7] = P[10][7] + P[10][4]*dt;
    nextP[10][8] = P[10][8] + P[10][5]*dt;
    nextP[10][9] = P[10][9] + P[10][6]*dt;
    nextP[10][10] = P[10][10];
    nextP[10][11] = P[10][11];
    nextP[10][12] = P[10][12];
    nextP[10][13] = P[10][13];
    nextP[10][14] = P[10][14];
    nextP[10][15] = P[10][15];
    nextP[10][16] = P[10][16];
    nextP[10][17] = P[10][17];
    nextP[10][18] = P[10][18];
    nextP[10][19] = P[10][19];
    nextP[10][20] = P[10][20];
    nextP[10][21] = P[10][21];
    nextP[10][22] = P[10][22];
    nextP[10][23] = P[10][23];
    nextP[11][0] = P[11][0] + P[11][1]*SF[9] + P[11][2]*SF[11] + P[11][3]*SF[10] + P[11][10]*SF[14] + P[11][11]*SF[15] + P[11][12]*SPP[12];
    nextP[11][1] = P[11][1] + P[11][0]*SF[8] + P[11][2]*SF[7] + P[11][3]*SF[11] - P[11][12]*SF[15] + P[11][11]*SPP[12] - (P[11][10]*q0)*0.5;
    nextP[11][2] = P[11][2] + P[11][0]*SF[6] + P[11][1]*SF[10] + P[11][3]*SF[8] + P[11][12]*SF[14] - P[11][10]*SPP[12] - (P[11][11]*q0)*0.5;
    nextP[11][3] = P[11][3] + P[11][0]*SF[7] + P[11][1]*SF[6] + P[11][2]*SF[9] + P[11][10]*SF[15] - P[11][11]*SF[14] - (P[11][12]*q0)*0.5;
    nextP[11][4] = P[11][4] + P[11][1]*SF[3] + P[11][0]*SF[5] + P[11][2]*SPP[0] - P[11][3]*SPP[5] + P[11][13]*SPP[3] + P[11][14]*SPP[8] - P[11][15]*SPP[11];
    nextP[11][5] = P[11][5] + P[11][0]*SF[4] + P[11][2]*SF[3] + P[11][3]*SF[5] - P[11][1]*SPP[0] + P[11][14]*SPP[2] + P[11][15]*SPP[7] - P[11][13]*SPP[10];
    nextP[11][6] = P[11][6] + P[11][1]*SF[4] + P[11][3]*SF[3] + P[11][0]*SPP[0] - P[11][2]*SPP[4] - P[11][15]*SPP[1] + P[11][13]*SPP[6] - P[11][14]*SPP[9];
    nextP[11][7] = P[11][7] + P[11][4]*dt;
    nextP[11][8] = P[11][8] + P[11][5]*dt;
    nextP[11][9] = P[11][9] + P[11][6]*dt;
    nextP[11][10] = P[11][10];
    nextP[11][11] = P[11][11];
    nextP[11][12] = P[11][12];
    nextP[11][13] = P[11][13];
    nextP[11][14] = P[11][14];
    nextP[11][15] = P[11][15];
    nextP[11][16] = P[11][16];
    nextP[11][17] = P[11][17];
    nextP[11][18] = P[11][18];
    nextP[11][19] = P[11][19];
    nextP[11][20] = P[11][20];
    nextP[11][21] = P[11][21];
    nextP[11][22] = P[11][22];
    nextP[11][23] = P[11][23];
    nextP[12][0] = P[12][0] + P[12][1]*SF[9] + P[12][2]*SF[11] + P[12][3]*SF[10] + P[12][10]*SF[14] + P[12][11]*SF[15] + P[12][12]*SPP[12];
    nextP[12][1] = P[12][1] + P[12][0]*SF[8] + P[12][2]*SF[7] + P[12][3]*SF[11] - P[12][12]*SF[15] + P[12][11]*SPP[12] - (P[12][10]*q0)*0.5;
    nextP[12][2] = P[12][2] + P[12][0]*SF[6] + P[12][1]*SF[10] + P[12][3]*SF[8] + P[12][12]*SF[14] - P[12][10]*SPP[12] - (P[12][11]*q0)*0.5;
    nextP[12][3] = P[12][3] + P[12][0]*SF[7] + P[12][1]*SF[6] + P[12][2]*SF[9] + P[12][10]*SF[15] - P[12][11]*SF[14] - (P[12][12]*q0)*0.5;
    nextP[12][4] = P[12][4] + P[12][1]*SF[3] + P[12][0]*SF[5] + P[12][2]*SPP[0] - P[12][3]*SPP[5] + P[12][13]*SPP[3] + P[12][14]*SPP[8] - P[12][15]*SPP[11];
    nextP[12][5] = P[12][5] + P[12][0]*SF[4] + P[12][2]*SF[3] + P[12][3]*SF[5] - P[12][1]*SPP[0] + P[12][14]*SPP[2] + P[12][15]*SPP[7] - P[12][13]*SPP[10];
    nextP[12][6] = P[12][6] + P[12][1]*SF[4] + P[12][3]*SF[3] + P[12][0]*SPP[0] - P[12][2]*SPP[4] - P[12][15]*SPP[1] + P[12][13]*SPP[6] - P[12][14]*SPP[9];
    nextP[12][7] = P[12][7] + P[12][4]*dt;
    nextP[12][8] = P[12][8] + P[12][5]*dt;
    nextP[12][9] = P[12][9] + P[12][6]*dt;
    nextP[12][10] = P[12][10];
    nextP[12][11] = P[12][11];
    nextP[12][12] = P[12][12];
    nextP[12][13] = P[12][13];
    nextP[12][14] = P[12][14];
    nextP[12][15] = P[12][15];
    nextP[12][16] = P[12][16];
    nextP[12][17] = P[12][17];
    nextP[12][18] = P[12][18];
    nextP[12][19] = P[12][19];
    nextP[12][20] = P[12][20];
    nextP[12][21] = P[12][21];
    nextP[12][22] = P[12][22];
    nextP[12][23] = P[12][23];
    nextP[13][0] = P[13][0] + P[13][1]*SF[9] + P[13][2]*SF[11] + P[13][3]*SF[10] + P[13][10]*SF[14] + P[13][11]*SF[15] + P[13][12]*SPP[12];
    nextP[13][1] = P[13][1] + P[13][0]*SF[8] + P[13][2]*SF[7] + P[13][3]*SF[11] - P[13][12]*SF[15] + P[13][11]*SPP[12] - (P[13][10]*q0)*0.5;
    nextP[13][2] = P[13][2] + P[13][0]*SF[6] + P[13][1]*SF[10] + P[13][3]*SF[8] + P[13][12]*SF[14] - P[13][10]*SPP[12] - (P[13][11]*q0)*0.5;
    nextP[13][3] = P[13][3] + P[13][0]*SF[7] + P[13][1]*SF[6] + P[13][2]*SF[9] + P[13][10]*SF[15] - P[13][11]*SF[14] - (P[13][12]*q0)*0.5;
    nextP[13][4] = P[13][4] + P[13][1]*SF[3] + P[13][0]*SF[5] + P[13][2]*SPP[0] - P[13][3]*SPP[5] + P[13][13]*SPP[3] + P[13][14]*SPP[8] - P[13][15]*SPP[11];
    nextP[13][5] = P[13][5] + P[13][0]*SF[4] + P[13][2]*SF[3] + P[13][3]*SF[5] - P[13][1]*SPP[0] + P[13][14]*SPP[2] + P[13][15]*SPP[7] - P[13][13]*SPP[10];
    nextP[13][6] = P[13][6] + P[13][1]*SF[4] + P[13][3]*SF[3] + P[13][0]*SPP[0] - P[13][2]*SPP[4] - P[13][15]*SPP[1] + P[13][13]*SPP[6] - P[13][14]*SPP[9];
    nextP[13][7] = P[13][7] + P[13][4]*dt;
    nextP[13][8] = P[13][8] + P[13][5]*dt;
    nextP[13][9] = P[13][9] + P[13][6]*dt;
    nextP[13][10] = P[13][10];
    nextP[13][11] = P[13][11];
    nextP[13][12] = P[13][12];
    nextP[13][13] = P[13][13];
    nextP[13][14] = P[13][14];
    nextP[13][15] = P[13][15];
    nextP[13][16] = P[13][16];
    nextP[13][17] = P[13][17];
    nextP[13][18] = P[13][18];
    nextP[13][19] = P[13][19];
    nextP[13][20] = P[13][20];
    nextP[13][21] = P[13][21];
    nextP[13][22] = P[13][22];
    nextP[13][23] = P[13][23];
    nextP[14][0] = P[14][0] + P[14][1]*SF[9] + P[14][2]*SF[11] + P[14][3]*SF[10] + P[14][10]*SF[14] + P[14][11]*SF[15] + P[14][12]*SPP[12];
    nextP[14][1] = P[14][1] + P[14][0]*SF[8] + P[14][2]*SF[7] + P[14][3]*SF[11] - P[14][12]*SF[15] + P[14][11]*SPP[12] - (P[14][10]*q0)*0.5;
    nextP[14][2] = P[14][2] + P[14][0]*SF[6] + P[14][1]*SF[10] + P[14][3]*SF[8] + P[14][12]*SF[14] - P[14][10]*SPP[12] - (P[14][11]*q0)*0.5;
    nextP[14][3] = P[14][3] + P[14][0]*SF[7] + P[14][1]*SF[6] + P[14][2]*SF[9] + P[14][10]*SF[15] - P[14][11]*SF[14] - (P[14][12]*q0)*0.5;
    nextP[14][4] = P[14][4] + P[14][1]*SF[3] + P[14][0]*SF[5] + P[14][2]*SPP[0] - P[14][3]*SPP[5] + P[14][13]*SPP[3] + P[14][14]*SPP[8] - P[14][15]*SPP[11];
    nextP[14][5] = P[14][5] + P[14][0]*SF[4] + P[14][2]*SF[3] + P[14][3]*SF[5] - P[14][1]*SPP[0] + P[14][14]*SPP[2] + P[14][15]*SPP[7] - P[14][13]*SPP[10];
    nextP[14][6] = P[14][6] + P[14][1]*SF[4] + P[14][3]*SF[3] + P[14][0]*SPP[0] - P[14][2]*SPP[4] - P[14][15]*SPP[1] + P[14][13]*SPP[6] - P[14][14]*SPP[9];
    nextP[14][7] = P[14][7] + P[14][4]*dt;
    nextP[14][8] = P[14][8] + P[14][5]*dt;
    nextP[14][9] = P[14][9] + P[14][6]*dt;
    nextP[14][10] = P[14][10];
    nextP[14][11] = P[14][11];
    nextP[14][12] = P[14][12];
    nextP[14][13] = P[14][13];
    nextP[14][14] = P[14][14];
    nextP[14][15] = P[14][15];
    nextP[14][16] = P[14][16];
    nextP[14][17] = P[14][17];
    nextP[14][18] = P[14][18];
    nextP[14][19] = P[14][19];
    nextP[14][20] = P[14][20];
    nextP[14][21] = P[14][21];
    nextP[14][22] = P[14][22];
    nextP[14][23] = P[14][23];
    nextP[15][0] = P[15][0] + P[15][1]*SF[9] + P[15][2]*SF[11] + P[15][3]*SF[10] + P[15][10]*SF[14] + P[15][11]*SF[15] + P[15][12]*SPP[12];
    nextP[15][1] = P[15][1] + P[15][0]*SF[8] + P[15][2]*SF[7] + P[15][3]*SF[11] - P[15][12]*SF[15] + P[15][11]*SPP[12] - (P[15][10]*q0)*0.5;
    nextP[15][2] = P[15][2] + P[15][0]*SF[6] + P[15][1]*SF[10] + P[15][3]*SF[8] + P[15][12]*SF[14] - P[15][10]*SPP[12] - (P[15][11]*q0)*0.5;
    nextP[15][3] = P[15][3] + P[15][0]*SF[7] + P[15][1]*SF[6] + P[15][2]*SF[9] + P[15][10]*SF[15] - P[15][11]*SF[14] - (P[15][12]*q0)*0.5;
    nextP[15][4] = P[15][4] + P[15][1]*SF[3] + P[15][0]*SF[5] + P[15][2]*SPP[0] - P[15][3]*SPP[5] + P[15][13]*SPP[3] + P[15][14]*SPP[8] - P[15][15]*SPP[11];
    nextP[15][5] = P[15][5] + P[15][0]*SF[4] + P[15][2]*SF[3] + P[15][3]*SF[5] - P[15][1]*SPP[0] + P[15][14]*SPP[2] + P[15][15]*SPP[7] - P[15][13]*SPP[10];
    nextP[15][6] = P[15][6] + P[15][1]*SF[4] + P[15][3]*SF[3] + P[15][0]*SPP[0] - P[15][2]*SPP[4] - P[15][15]*SPP[1] + P[15][13]*SPP[6] - P[15][14]*SPP[9];
    nextP[15][7] = P[15][7] + P[15][4]*dt;
    nextP[15][8] = P[15][8] + P[15][5]*dt;
    nextP[15][9] = P[15][9] + P[15][6]*dt;
    nextP[15][10] = P[15][10];
    nextP[15][11] = P[15][11];
    nextP[15][12] = P[15][12];
    nextP[15][13] = P[15][13];
    nextP[15][14] = P[15][14];
    nextP[15][15] = P[15][15];
    nextP[15][16] = P[15][16];
    nextP[15][17] = P[15][17];
    nextP[15][18] = P[15][18];
    nextP[15][19] = P[15][19];
    nextP[15][20] = P[15][20];
    nextP[15][21] = P[15][21];
    nextP[15][22] = P[15][22];
    nextP[15][23] = P[15][23];
    nextP[16][0] = P[16][0] + P[16][1]*SF[9] + P[16][2]*SF[11] + P[16][3]*SF[10] + P[16][10]*SF[14] + P[16][11]*SF[15] + P[16][12]*SPP[12];
    nextP[16][1] = P[16][1] + P[16][0]*SF[8] + P[16][2]*SF[7] + P[16][3]*SF[11] - P[16][12]*SF[15] + P[16][11]*SPP[12] - (P[16][10]*q0)*0.5;
    nextP[16][2] = P[16][2] + P[16][0]*SF[6] + P[16][1]*SF[10] + P[16][3]*SF[8] + P[16][12]*SF[14] - P[16][10]*SPP[12] - (P[16][11]*q0)*0.5;
    nextP[16][3] = P[16][3] + P[16][0]*SF[7] + P[16][1]*SF[6] + P[16][2]*SF[9] + P[16][10]*SF[15] - P[16][11]*SF[14] - (P[16][12]*q0)*0.5;
    nextP[16][4] = P[16][4] + P[16][1]*SF[3] + P[16][0]*SF[5] + P[16][2]*SPP[0] - P[16][3]*SPP[5] + P[16][13]*SPP[3] + P[16][14]*SPP[8] - P[16][15]*SPP[11];
    nextP[16][5] = P[16][5] + P[16][0]*SF[4] + P[16][2]*SF[3] + P[16][3]*SF[5] - P[16][1]*SPP[0] + P[16][14]*SPP[2] + P[16][15]*SPP[7] - P[16][13]*SPP[10];
    nextP[16][6] = P[16][6] + P[16][1]*SF[4] + P[16][3]*SF[3] + P[16][0]*SPP[0] - P[16][2]*SPP[4] - P[16][15]*SPP[1] + P[16][13]*SPP[6] - P[16][14]*SPP[9];
    nextP[16][7] = P[16][7] + P[16][4]*dt;
    nextP[16][8] = P[16][8] + P[16][5]*dt;
    nextP[16][9] = P[16][9] + P[16][6]*dt;
    nextP[16][10] = P[16][10];
    nextP[16][11] = P[16][11];
    nextP[16][12] = P[16][12];
    nextP[16][13] = P[16][13];
    nextP[16][14] = P[16][14];
    nextP[16][15] = P[16][15];
    nextP[16][16] = P[16][16];
    nextP[16][17] = P[16][17];
    nextP[16][18] = P[16][18];
    nextP[16][19] = P[16][19];
    nextP[16][20] = P[16][20];
    nextP[16][21] = P[16][21];
    nextP[16][22] = P[16][22];
    nextP[16][23] = P[16][23];
    nextP[17][0] = P[17][0] + P[17][1]*SF[9] + P[17][2]*SF[11] + P[17][3]*SF[10] + P[17][10]*SF[14] + P[17][11]*SF[15] + P[17][12]*SPP[12];
    nextP[17][1] = P[17][1] + P[17][0]*SF[8] + P[17][2]*SF[7] + P[17][3]*SF[11] - P[17][12]*SF[15] + P[17][11]*SPP[12] - (P[17][10]*q0)*0.5;
    nextP[17][2] = P[17][2] + P[17][0]*SF[6] + P[17][1]*SF[10] + P[17][3]*SF[8] + P[17][12]*SF[14] - P[17][10]*SPP[12] - (P[17][11]*q0)*0.5;
    nextP[17][3] = P[17][3] + P[17][0]*SF[7] + P[17][1]*SF[6] + P[17][2]*SF[9] + P[17][10]*SF[15] - P[17][11]*SF[14] - (P[17][12]*q0)*0.5;
    nextP[17][4] = P[17][4] + P[17][1]*SF[3] + P[17][0]*SF[5] + P[17][2]*SPP[0] - P[17][3]*SPP[5] + P[17][13]*SPP[3] + P[17][14]*SPP[8] - P[17][15]*SPP[11];
    nextP[17][5] = P[17][5] + P[17][0]*SF[4] + P[17][2]*SF[3] + P[17][3]*SF[5] - P[17][1]*SPP[0] + P[17][14]*SPP[2] + P[17][15]*SPP[7] - P[17][13]*SPP[10];
    nextP[17][6] = P[17][6] + P[17][1]*SF[4] + P[17][3]*SF[3] + P[17][0]*SPP[0] - P[17][2]*SPP[4] - P[17][15]*SPP[1] + P[17][13]*SPP[6] - P[17][14]*SPP[9];
    nextP[17][7] = P[17][7] + P[17][4]*dt;
    nextP[17][8] = P[17][8] + P[17][5]*dt;
    nextP[17][9] = P[17][9] + P[17][6]*dt;
    nextP[17][10] = P[17][10];
    nextP[17][11] = P[17][11];
    nextP[17][12] = P[17][12];
    nextP[17][13] = P[17][13];
    nextP[17][14] = P[17][14];
    nextP[17][15] = P[17][15];
    nextP[17][16] = P[17][16];
    nextP[17][17] = P[17][17];
    nextP[17][18] = P[17][18];
    nextP[17][19] = P[17][19];
    nextP[17][20] = P[17][20];
    nextP[17][21] = P[17][21];
    nextP[17][22] = P[17][22];
    nextP[17][23] = P[17][23];
    nextP[18][0] = P[18][0] + P[18][1]*SF[9] + P[18][2]*SF[11] + P[18][3]*SF[10] + P[18][10]*SF[14] + P[18][11]*SF[15] + P[18][12]*SPP[12];
    nextP[18][1] = P[18][1] + P[18][0]*SF[8] + P[18][2]*SF[7] + P[18][3]*SF[11] - P[18][12]*SF[15] + P[18][11]*SPP[12] - (P[18][10]*q0)*0.5;
    nextP[18][2] = P[18][2] + P[18][0]*SF[6] + P[18][1]*SF[10] + P[18][3]*SF[8] + P[18][12]*SF[14] - P[18][10]*SPP[12] - (P[18][11]*q0)*0.5;
    nextP[18][3] = P[18][3] + P[18][0]*SF[7] + P[18][1]*SF[6] + P[18][2]*SF[9] + P[18][10]*SF[15] - P[18][11]*SF[14] - (P[18][12]*q0)*0.5;
    nextP[18][4] = P[18][4] + P[18][1]*SF[3] + P[18][0]*SF[5] + P[18][2]*SPP[0] - P[18][3]*SPP[5] + P[18][13]*SPP[3] + P[18][14]*SPP[8] - P[18][15]*SPP[11];
    nextP[18][5] = P[18][5] + P[18][0]*SF[4] + P[18][2]*SF[3] + P[18][3]*SF[5] - P[18][1]*SPP[0] + P[18][14]*SPP[2] + P[18][15]*SPP[7] - P[18][13]*SPP[10];
    nextP[18][6] = P[18][6] + P[18][1]*SF[4] + P[18][3]*SF[3] + P[18][0]*SPP[0] - P[18][2]*SPP[4] - P[18][15]*SPP[1] + P[18][13]*SPP[6] - P[18][14]*SPP[9];
    nextP[18][7] = P[18][7] + P[18][4]*dt;
    nextP[18][8] = P[18][8] + P[18][5]*dt;
    nextP[18][9] = P[18][9] + P[18][6]*dt;
    nextP[18][10] = P[18][10];
    nextP[18][11] = P[18][11];
    nextP[18][12] = P[18][12];
    nextP[18][13] = P[18][13];
    nextP[18][14] = P[18][14];
    nextP[18][15] = P[18][15];
    nextP[18][16] = P[18][16];
    nextP[18][17] = P[18][17];
    nextP[18][18] = P[18][18];
    nextP[18][19] = P[18][19];
    nextP[18][20] = P[18][20];
    nextP[18][21] = P[18][21];
    nextP[18][22] = P[18][22];
    nextP[18][23] = P[18][23];
    nextP[19][0] = P[19][0] + P[19][1]*SF[9] + P[19][2]*SF[11] + P[19][3]*SF[10] + P[19][10]*SF[14] + P[19][11]*SF[15] + P[19][12]*SPP[12];
    nextP[19][1] = P[19][1] + P[19][0]*SF[8] + P[19][2]*SF[7] + P[19][3]*SF[11] - P[19][12]*SF[15] + P[19][11]*SPP[12] - (P[19][10]*q0)*0.5;
    nextP[19][2] = P[19][2] + P[19][0]*SF[6] + P[19][1]*SF[10] + P[19][3]*SF[8] + P[19][12]*SF[14] - P[19][10]*SPP[12] - (P[19][11]*q0)*0.5;
    nextP[19][3] = P[19][3] + P[19][0]*SF[7] + P[19][1]*SF[6] + P[19][2]*SF[9] + P[19][10]*SF[15] - P[19][11]*SF[14] - (P[19][12]*q0)*0.5;
    nextP[19][4] = P[19][4] + P[19][1]*SF[3] + P[19][0]*SF[5] + P[19][2]*SPP[0] - P[19][3]*SPP[5] + P[19][13]*SPP[3] + P[19][14]*SPP[8] - P[19][15]*SPP[11];
    nextP[19][5] = P[19][5] + P[19][0]*SF[4] + P[19][2]*SF[3] + P[19][3]*SF[5] - P[19][1]*SPP[0] + P[19][14]*SPP[2] + P[19][15]*SPP[7] - P[19][13]*SPP[10];
    nextP[19][6] = P[19][6] + P[19][1]*SF[4] + P[19][3]*SF[3] + P[19][0]*SPP[0] - P[19][2]*SPP[4] - P[19][15]*SPP[1] + P[19][13]*SPP[6] - P[19][14]*SPP[9];
    nextP[19][7] = P[19][7] + P[19][4]*dt;
    nextP[19][8] = P[19][8] + P[19][5]*dt;
    nextP[19][9] = P[19][9] + P[19][6]*dt;
    nextP[19][10] = P[19][10];
    nextP[19][11] = P[19][11];
    nextP[19][12] = P[19][12];
    nextP[19][13] = P[19][13];
    nextP[19][14] = P[19][14];
    nextP[19][15] = P[19][15];
    nextP[19][16] = P[19][16];
    nextP[19][17] = P[19][17];
    nextP[19][18] = P[19][18];
    nextP[19][19] = P[19][19];
    nextP[19][20] = P[19][20];
    nextP[19][21] = P[19][21];
    nextP[19][22] = P[19][22];
    nextP[19][23] = P[19][23];
    nextP[20][0] = P[20][0] + P[20][1]*SF[9] + P[20][2]*SF[11] + P[20][3]*SF[10] + P[20][10]*SF[14] + P[20][11]*SF[15] + P[20][12]*SPP[12];
    nextP[20][1] = P[20][1] + P[20][0]*SF[8] + P[20][2]*SF[7] + P[20][3]*SF[11] - P[20][12]*SF[15] + P[20][11]*SPP[12] - (P[20][10]*q0)*0.5;
    nextP[20][2] = P[20][2] + P[20][0]*SF[6] + P[20][1]*SF[10] + P[20][3]*SF[8] + P[20][12]*SF[14] - P[20][10]*SPP[12] - (P[20][11]*q0)*0.5;
    nextP[20][3] = P[20][3] + P[20][0]*SF[7] + P[20][1]*SF[6] + P[20][2]*SF[9] + P[20][10]*SF[15] - P[20][11]*SF[14] - (P[20][12]*q0)*0.5;
    nextP[20][4] = P[20][4] + P[20][1]*SF[3] + P[20][0]*SF[5] + P[20][2]*SPP[0] - P[20][3]*SPP[5] + P[20][13]*SPP[3] + P[20][14]*SPP[8] - P[20][15]*SPP[11];
    nextP[20][5] = P[20][5] + P[20][0]*SF[4] + P[20][2]*SF[3] + P[20][3]*SF[5] - P[20][1]*SPP[0] + P[20][14]*SPP[2] + P[20][15]*SPP[7] - P[20][13]*SPP[10];
    nextP[20][6] = P[20][6] + P[20][1]*SF[4] + P[20][3]*SF[3] + P[20][0]*SPP[0] - P[20][2]*SPP[4] - P[20][15]*SPP[1] + P[20][13]*SPP[6] - P[20][14]*SPP[9];
    nextP[20][7] = P[20][7] + P[20][4]*dt;
    nextP[20][8] = P[20][8] + P[20][5]*dt;
    nextP[20][9] = P[20][9] + P[20][6]*dt;
    nextP[20][10] = P[20][10];
    nextP[20][11] = P[20][11];
    nextP[20][12] = P[20][12];
    nextP[20][13] = P[20][13];
    nextP[20][14] = P[20][14];
    nextP[20][15] = P[20][15];
    nextP[20][16] = P[20][16];
    nextP[20][17] = P[20][17];
    nextP[20][18] = P[20][18];
    nextP[20][19] = P[20][19];
    nextP[20][20] = P[20][20];
    nextP[20][21] = P[20][21];
    nextP[20][22] = P[20][22];
    nextP[20][23] = P[20][23];
    nextP[21][0] = P[21][0] + P[21][1]*SF[9] + P[21][2]*SF[11] + P[21][3]*SF[10] + P[21][10]*SF[14] + P[21][11]*SF[15] + P[21][12]*SPP[12];
    nextP[21][1] = P[21][1] + P[21][0]*SF[8] + P[21][2]*SF[7] + P[21][3]*SF[11] - P[21][12]*SF[15] + P[21][11]*SPP[12] - (P[21][10]*q0)*0.5;
    nextP[21][2] = P[21][2] + P[21][0]*SF[6] + P[21][1]*SF[10] + P[21][3]*SF[8] + P[21][12]*SF[14] - P[21][10]*SPP[12] - (P[21][11]*q0)*0.5;
    nextP[21][3] = P[21][3] + P[21][0]*SF[7] + P[21][1]*SF[6] + P[21][2]*SF[9] + P[21][10]*SF[15] - P[21][11]*SF[14] - (P[21][12]*q0)*0.5;
    nextP[21][4] = P[21][4] + P[21][1]*SF[3] + P[21][0]*SF[5] + P[21][2]*SPP[0] - P[21][3]*SPP[5] + P[21][13]*SPP[3] + P[21][14]*SPP[8] - P[21][15]*SPP[11];
    nextP[21][5] = P[21][5] + P[21][0]*SF[4] + P[21][2]*SF[3] + P[21][3]*SF[5] - P[21][1]*SPP[0] + P[21][14]*SPP[2] + P[21][15]*SPP[7] - P[21][13]*SPP[10];
    nextP[21][6] = P[21][6] + P[21][1]*SF[4] + P[21][3]*SF[3] + P[21][0]*SPP[0] - P[21][2]*SPP[4] - P[21][15]*SPP[1] + P[21][13]*SPP[6] - P[21][14]*SPP[9];
    nextP[21][7] = P[21][7] + P[21][4]*dt;
    nextP[21][8] = P[21][8] + P[21][5]*dt;
    nextP[21][9] = P[21][9] + P[21][6]*dt;
    nextP[21][10] = P[21][10];
    nextP[21][11] = P[21][11];
    nextP[21][12] = P[21][12];
    nextP[21][13] = P[21][13];
    nextP[21][14] = P[21][14];
    nextP[21][15] = P[21][15];
    nextP[21][16] = P[21][16];
    nextP[21][17] = P[21][17];
    nextP[21][18] = P[21][18];
    nextP[21][19] = P[21][19];
    nextP[21][20] = P[21][20];
    nextP[21][21] = P[21][21];
    nextP[21][22] = P[21][22];
    nextP[21][23] = P[21][23];
    nextP[22][0] = P[22][0] + P[22][1]*SF[9] + P[22][2]*SF[11] + P[22][3]*SF[10] + P[22][10]*SF[14] + P[22][11]*SF[15] + P[22][12]*SPP[12];
    nextP[22][1] = P[22][1] + P[22][0]*SF[8] + P[22][2]*SF[7] + P[22][3]*SF[11] - P[22][12]*SF[15] + P[22][11]*SPP[12] - (P[22][10]*q0)*0.5;
    nextP[22][2] = P[22][2] + P[22][0]*SF[6] + P[22][1]*SF[10] + P[22][3]*SF[8] + P[22][12]*SF[14] - P[22][10]*SPP[12] - (P[22][11]*q0)*0.5;
    nextP[22][3] = P[22][3] + P[22][0]*SF[7] + P[22][1]*SF[6] + P[22][2]*SF[9] + P[22][10]*SF[15] - P[22][11]*SF[14] - (P[22][12]*q0)*0.5;
    nextP[22][4] = P[22][4] + P[22][1]*SF[3] + P[22][0]*SF[5] + P[22][2]*SPP[0] - P[22][3]*SPP[5] + P[22][13]*SPP[3] + P[22][14]*SPP[8] - P[22][15]*SPP[11];
    nextP[22][5] = P[22][5] + P[22][0]*SF[4] + P[22][2]*SF[3] + P[22][3]*SF[5] - P[22][1]*SPP[0] + P[22][14]*SPP[2] + P[22][15]*SPP[7] - P[22][13]*SPP[10];
    nextP[22][6] = P[22][6] + P[22][1]*SF[4] + P[22][3]*SF[3] + P[22][0]*SPP[0] - P[22][2]*SPP[4] - P[22][15]*SPP[1] + P[22][13]*SPP[6] - P[22][14]*SPP[9];
    nextP[22][7] = P[22][7] + P[22][4]*dt;
    nextP[22][8] = P[22][8] + P[22][5]*dt;
    nextP[22][9] = P[22][9] + P[22][6]*dt;
    nextP[22][10] = P[22][10];
    nextP[22][11] = P[22][11];
    nextP[22][12] = P[22][12];
    nextP[22][13] = P[22][13];
    nextP[22][14] = P[22][14];
    nextP[22][15] = P[22][15];
    nextP[22][16] = P[22][16];
    nextP[22][17] = P[22][17];
    nextP[22][18] = P[22][18];
    nextP[22][19] = P[22][19];
    nextP[22][20] = P[22][20];
    nextP[22][21] = P[22][21];
    nextP[22][22] = P[22][22];
    nextP[22][23] = P[22][23];
    nextP[23][0] = P[23][0] + P[23][1]*SF[9] + P[23][2]*SF[11] + P[23][3]*SF[10] + P[23][10]*SF[14] + P[23][11]*SF[15] + P[23][12]*SPP[12];
    nextP[23][1] = P[23][1] + P[23][0]*SF[8] + P[23][2]*SF[7] + P[23][3]*SF[11] - P[23][12]*SF[15] + P[23][11]*SPP[12] - (P[23][10]*q0)*0.5;
    nextP[23][2] = P[23][2] + P[23][0]*SF[6] + P[23][1]*SF[10] + P[23][3]*SF[8] + P[23][12]*SF[14] - P[23][10]*SPP[12] - (P[23][11]*q0)*0.5;
    nextP[23][3] = P[23][3] + P[23][0]*SF[7] + P[23][1]*SF[6] + P[23][2]*SF[9] + P[23][10]*SF[15] - P[23][11]*SF[14] - (P[23][12]*q0)*0.5;
    nextP[23][4] = P[23][4] + P[23][1]*SF[3] + P[23][0]*SF[5] + P[23][2]*SPP[0] - P[23][3]*SPP[5] + P[23][13]*SPP[3] + P[23][14]*SPP[8] - P[23][15]*SPP[11];
    nextP[23][5] = P[23][5] + P[23][0]*SF[4] + P[23][2]*SF[3] + P[23][3]*SF[5] - P[23][1]*SPP[0] + P[23][14]*SPP[2] + P[23][15]*SPP[7] - P[23][13]*SPP[10];
    nextP[23][6] = P[23][6] + P[23][1]*SF[4] + P[23][3]*SF[3] + P[23][0]*SPP[0] - P[23][2]*SPP[4] - P[23][15]*SPP[1] + P[23][13]*SPP[6] - P[23][14]*SPP[9];
    nextP[23][7] = P[23][7] + P[23][4]*dt;
    nextP[23][8] = P[23][8] + P[23][5]*dt;
    nextP[23][9] = P[23][9] + P[23][6]*dt;
    nextP[23][10] = P[23][10];
    nextP[23][11] = P[23][11];
    nextP[23][12] = P[23][12];
    nextP[23][13] = P[23][13];
    nextP[23][14] = P[23][14];
    nextP[23][15] = P[23][15];
    nextP[23][16] = P[23][16];
    nextP[23][17] = P[23][17];
    nextP[23][18] = P[23][18];
    nextP[23][19] = P[23][19];
    nextP[23][20] = P[23][20];
    nextP[23][21] = P[23][21];
    nextP[23][22] = P[23][22];
    nextP[23][23] = P[23][23];

    for (i=0; i<= 23; i++)
    {
        nextP[i][i] = nextP[i][i] + processNoise[i];
    }

    // If on ground, inhibit accelerometer bias, wind  and magnetic field state updates by
    // setting the corresponding covariance terms to zero
    // This is a quick hack - for efficiency should really not calculate terms above when this condition is true
    if (onGround)
    {
        zeroRows(nextP,13,15);
        zeroCols(nextP,13,15);
        zeroRows(nextP,16,17);
        zeroCols(nextP,16,17);
        zeroRows(nextP,18,23);
        zeroCols(nextP,18,23);
    }
    // If on ground, inhibit accelerometer bias updates by setting the coresponding
    // covariance terms to zero. To be efficient, the corresponding terms should
    // also not be calculated above
    if (onGround)
    {
        zeroRows(nextP,13,15);
        zeroCols(nextP,13,15);
    }

    // If on ground or no magnetometer fitted, inhibit magnetometer bias updates by
    // setting the coresponding covariance terms to zero. To be efficient, the
    // corresponding terms should also not be calculated above
    if (onGround || !useCompass)
    {
        zeroRows(nextP,18,23);
        zeroCols(nextP,18,23);
    }

    // If on ground or not using airspeed sensing, inhibit wind velocity
    // covariance growth. To be efficient, the corresponding terms should also
    // not be calculated above
    if (onGround || ~useAirspeed)
    {
        for (i=16; i<=17; i++)
        {
            for (j=0; j<=23; j++)
            {
                nextP[i][j] = P[i][j];
                nextP[j][i] = P[j][i];
            }
        }
    }

    // If the total position variance exceds 1E6 (1000m), then stop covariance
    // growth by setting the predicted to the previous values
    // This prevent an ill conditioned matrix from occurring for long periods
    // without GPS
    if ((P[7][7] + P[8][8]) > 1E6)
    {
        for (i=7; i<=8; i++)
        {
            for (j=0; j<=23; j++)
            {
                nextP[i][j] = P[i][j];
                nextP[j][i] = P[j][i];
            }
        }
    }

    // Force symmetry on the covariance matrix to prevent ill-conditioning
    // of the matrix which would cause the filter to blow-up
    for (i=0; i<=23; i++) P[i][i] = nextP[i][i];
    for (i=1; i<=23; i++)
    {
        for (j=i-1; j<=22; j++)
        {
            temp = 0.5*(nextP[i][j] + nextP[j][i]);
            P[i][j] = temp;
            P[j][i] = temp;
        }
    }

}

void  UpdateStrapdownEquationsNED(
    Vector3f correctedDelAng,
    Vector3f correctedDelVel,
    float accNavMag,
    Vector3f earthRateNED,
    Vector3f angRate,
    Vector3f accel,
    float dt)
{
    static Vector3f prevAngRate;
    static Vector3f prevAccel;
    static Vector3f  prevDelAng;
    static Vector3f  prevDelVel;
    Vector3f dAng;
    Vector3f dVel;
    Vector3f delVelNav;
    float q00 = sq(states[0]);
    float q11;
    float q22;
    float q33;
    float q01;
    float q02;
    float q03;
    float q12;
    float q13;
    float q23;
    static Mat3f  Tbn;
    static Mat3f  Tnb;
    float rotationMag;
    float rotScaler;
    float qUpdated[4];
    float quatMag;
    float quatMagInv;
    float deltaQuat[4];
    float lastVelocity[3];
    const Vector3f gravityNED = {0.0,0.0,GRAVITY_MSS};

// Convert IMU data to delta angles and velocities using trapezoidal
// integration
    dAng = (angRate + prevAngRate)*dt*0.5;
    dVel = (accel + prevAccel)*dt*0.5;
    prevAngRate = angRate;
    prevAccel   = accel;

// Remove sensor bias errors
    dAng.x = dAng.x - states[10];
    dAng.y = dAng.y - states[11];
    dAng.z = dAng.z - states[12];
    dVel.x = dVel.x - states[13];
    dVel.y = dVel.y - states[14];
    dVel.z = dVel.z - states[15];

// Apply rotational and skulling corrections
// * and + operators have been overloaded

    correctedDelVel = dVel + 0.5*(prevDelAng + dAng)*(prevDelVel + dVel) + 0.1666667*((prevDelAng + dAng)*((prevDelAng + dAng)*(prevDelVel + dVel))) + 0.08333333*((prevDelAng*dVel) + (prevDelVel*dAng));

// Apply corrections for earths rotation rate and coning errors
// * and + operators have been overloaded
    correctedDelAng   = dAng - Tnb*earthRateNED*dt + 0.08333333*prevDelAng*dAng;

// Save current measurements
    prevDelAng = dAng;
    prevDelVel = dVel;

// Convert the rotation vector to its equivalent quaternion
    rotationMag = sqrt(sq(correctedDelAng.x) + sq(correctedDelAng.y) + sq(correctedDelAng.z));
    if (rotationMag < 1e-12)
    {
        deltaQuat[0] = 1.0;
        deltaQuat[1] = 0.0;
        deltaQuat[2] = 0.0;
        deltaQuat[3] = 0.0;
    }
    else
    {
        deltaQuat[0] = cos(0.5*rotationMag);
        rotScaler = (sin(0.5*rotationMag))/rotationMag;
        deltaQuat[1] = correctedDelAng.x*rotScaler;
        deltaQuat[2] = correctedDelAng.y*rotScaler;
        deltaQuat[3] = correctedDelAng.z*rotScaler;
    }

// Update the quaternions by rotating from the previous attitude through
// the delta angle rotation quaternion
    qUpdated[0] = states[0]*deltaQuat[0] - states[1]*deltaQuat[1] - states[2]*deltaQuat[2] - states[3]*deltaQuat[3];
    qUpdated[1] = states[0]*deltaQuat[1] + states[1]*deltaQuat[0] - states[2]*deltaQuat[3] + states[3]*deltaQuat[2];
    qUpdated[2] = states[0]*deltaQuat[2] + states[1]*deltaQuat[3] + states[2]*deltaQuat[0] - states[3]*deltaQuat[1];
    qUpdated[3] = states[0]*deltaQuat[3] - states[1]*deltaQuat[2] + states[2]*deltaQuat[1] + states[3]*deltaQuat[0];

// Normalise the quaternions and update the quaternion states
    quatMag = sqrt(sq(qUpdated[0]) + sq(qUpdated[1]) + sq(qUpdated[2]) + sq(qUpdated[3]));
    if (quatMag > 1e-16)
    {
        quatMagInv = 1.0/quatMag;
        states[0] = quatMagInv*qUpdated[0];
        states[1] = quatMagInv*qUpdated[1];
        states[2] = quatMagInv*qUpdated[2];
        states[3] = quatMagInv*qUpdated[3];
    }

// Calculate the body to nav cosine matrix
    q00 = sq(states[0]);
    q11 = sq(states[1]);
    q22 = sq(states[2]);
    q33 = sq(states[3]);
    q01 =  states[0]*states[1];
    q02 =  states[0]*states[2];
    q03 =  states[0]*states[3];
    q12 =  states[1]*states[2];
    q13 =  states[1]*states[3];
    q23 =  states[2]*states[3];

    Tnb.x.x = q00 + q11 - q22 - q33;
    Tnb.y.y = q00 - q11 + q22 - q33;
    Tnb.z.z = q00 - q11 - q22 + q33;
    Tnb.x.y = 2.0*(q12 + q03);
    Tnb.x.z = 2.0*(q13 - q02);
    Tnb.y.x = 2.0*(q12 - q03);
    Tnb.y.z = 2.0*(q23 + q01);
    Tnb.z.x = 2.0*(q13 + q02);
    Tnb.z.y = 2.0*(q23 - q01);

    Tbn.x.x = Tnb.x.x;
    Tbn.y.y = Tnb.y.y;
    Tbn.z.z = Tnb.z.z;
    Tbn.x.y = Tnb.y.x;
    Tbn.x.z = Tnb.z.x;
    Tbn.y.x = Tnb.x.y;
    Tbn.y.z = Tnb.z.y;
    Tbn.z.x = Tnb.x.z;
    Tbn.z.y = Tnb.y.z;

// transform body delta velocities to delta velocities in the nav frame
// * and + operators have been overloaded
    delVelNav = Tbn*correctedDelVel + gravityNED*dt;

// calculate the magnitude of the nav acceleration (required for GPS
// variance estimation)
    accNavMag = sqrt(delVelNav.x*delVelNav.x + delVelNav.y*delVelNav.y + delVelNav.z*delVelNav.z)/dt;

// If calculating position save previous velocity
    lastVelocity[0] = states[4];
    lastVelocity[1] = states[5];
    lastVelocity[2] = states[6];

// Sum delta velocities to get velocity after removing gravitational skew and correcting for transport rate
    states[4] = states[4] + delVelNav.x;
    states[5] = states[5] + delVelNav.y;
    states[6] = states[6] + delVelNav.z;

// If calculating postions, do a trapezoidal integration for position
    states[7] = states[7] + 0.5*(states[4] + lastVelocity[0])*dt;
    states[8] = states[8] + 0.5*(states[5] + lastVelocity[1])*dt;
    states[9] = states[9] + 0.5*(states[6] + lastVelocity[2])*dt;
}

void FuseVelPosNED(
    float innovation[6], // innovation output
    float varInnov[6], // innovation variance output
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

// declare variables used to check measurement errors
    float velInnov[3] = {0.0,0.0,0.0};
    float posInnov[2] = {0.0,0.0};
    float hgtInnov = 0.0;

// declare indices used to access arrays
    uint8_t stateIndex;
    uint8_t obsIndex;
    uint8_t i;
    uint8_t j;
    uint8_t iMax;

// declare variables used by state and covariance update calculations
    float velErr;
    float posErr;
    float R_OBS[6];
    float observation[6];
    float KHP[24][24];
    float SK;
    float K[24];
    float quatMag;
    uint8_t startIndex;
    uint8_t endIndex;

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
            if (useVelD) iMax = 2;
            else iMax = 1;
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
}

void FuseMagnetometer(
    float innovation[3], // innovation output
    float varInnov[3], // innovation variance output
    bool FuseData, // boolean true when magnetometer data is to be fused
    float MagData[3], // magnetometer flux radings in X,Y,Z body axes
    float StatesAtMeasTime[24], // filter satates at the effective measurement time
    bool useCompass) // boolean true if magnetometer data is being used
{

    static float q0 = 1.0;
    static float q1 = 0.0;
    static float q2 = 0.0;
    static float q3 = 0.0;
    static float magN = 0.0;
    static float magE = 0.0;
    static float magD = 0.0;
    static float magXbias = 0.0;
    static float magYbias = 0.0;
    static float magZbias = 0.0;
    static uint8_t obsIndex = 0;
    static float SH_MAG[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    static float Tnb[3][3] =
    {
        {1.0,0.0,0.0} ,
        {0.0,1.0,0.0} ,
        {0.0,0.0,1.0}
    };
    static float MagPred[3] = {0.0,0.0,0.0};
    uint8_t i;
    uint8_t j;
    uint8_t k;
    const float R_MAG = 2500.0;
    float H_MAG[24];
    float SK_MX[6];
    float SK_MY[5];
    float SK_MZ[6];
    float K_MAG[24];
    float KH[24][24];
    float KHP[24][24];

// Perform sequential fusion of Magnetometer measurements.
// This assumes that the errors in the different components are
// uncorrelated which is not true, however in the absence of covariance
// data fit is the only assumption we can make
// so we might as well take advantage of the computational efficiencies
// associated with sequential fusion
    if (useCompass && (FuseData || obsIndex == 1 || obsIndex == 2))
    {

        // Sequential fusion of XYZ components to spread processing load across
        // three prediction time steps.

        // Calculate observation jacobians and Kalman gains
        if (FuseData)
        {
            // Copy required states to local variable names
            q0       = StatesAtMeasTime[0];
            q1       = StatesAtMeasTime[1];
            q2       = StatesAtMeasTime[2];
            q3       = StatesAtMeasTime[3];
            magN     = StatesAtMeasTime[18];
            magE     = StatesAtMeasTime[19];
            magD     = StatesAtMeasTime[20];
            magXbias = StatesAtMeasTime[21];
            magYbias = StatesAtMeasTime[22];
            magZbias = StatesAtMeasTime[23];
            // rotate predicted earth components into body axes and calculate
            // predicted measurments
            Tnb[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
            Tnb[0][1] = 2.0*(q1*q2 + q0*q3);
            Tnb[0][2] = 2.0*(q1*q3-q0*q2);
            Tnb[1][0] = 2.0*(q1*q2 - q0*q3);
            Tnb[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
            Tnb[1][2] = 2.0*(q2*q3 + q0*q1);
            Tnb[2][0] = 2.0*(q1*q3 + q0*q2);
            Tnb[2][1] = 2.0*(q2*q3 - q0*q1);
            Tnb[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
            MagPred[0] = Tnb[0][0]*magN + Tnb[0][1]*magE  + Tnb[0][2]*magD + magXbias;
            MagPred[1] = Tnb[1][0]*magN + Tnb[1][1]*magE  + Tnb[1][2]*magD + magYbias;
            MagPred[2] = Tnb[2][0]*magN + Tnb[2][1]*magE  + Tnb[0][2]*magD + magZbias;
            // Calculate observation jacobians
            SH_MAG[0] = 2*magD*q3 + 2*magE*q2 + 2*magN*q1;
            SH_MAG[1] = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
            SH_MAG[2] = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
            SH_MAG[3] = sq(q3);
            SH_MAG[4] = sq(q2);
            SH_MAG[5] = sq(q1);
            SH_MAG[6] = sq(q0);
            SH_MAG[7] = 2*magN*q0;
            SH_MAG[8] = 2*magE*q3;
            H_MAG[0] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            H_MAG[1] = SH_MAG[0];
            H_MAG[2] = 2*magE*q1 - 2*magD*q0 - 2*magN*q2;
            H_MAG[3] = SH_MAG[2];
            H_MAG[18] = SH_MAG[5] - SH_MAG[4] - SH_MAG[3] + SH_MAG[6];
            H_MAG[19] = 2*q0*q3 + 2*q1*q2;
            H_MAG[20] = 2*q1*q3 - 2*q0*q2;
            H_MAG[21] = 1;
            // Calculate Kalman gain
            SK_MX[0] = 1/(P[21][21] + R_MAG + P[1][21]*SH_MAG[0] + P[3][21]*SH_MAG[2] - P[18][21]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) - (2*magD*q0 - 2*magE*q1 + 2*magN*q2)*(P[21][2] + P[1][2]*SH_MAG[0] + P[3][2]*SH_MAG[2] - P[18][2]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][2]*(2*q0*q3 + 2*q1*q2) - P[20][2]*(2*q0*q2 - 2*q1*q3) - P[2][2]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[21][0] + P[1][0]*SH_MAG[0] + P[3][0]*SH_MAG[2] - P[18][0]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][0]*(2*q0*q3 + 2*q1*q2) - P[20][0]*(2*q0*q2 - 2*q1*q3) - P[2][0]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[21][1] + P[1][1]*SH_MAG[0] + P[3][1]*SH_MAG[2] - P[18][1]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][1]*(2*q0*q3 + 2*q1*q2) - P[20][1]*(2*q0*q2 - 2*q1*q3) - P[2][1]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[2]*(P[21][3] + P[1][3]*SH_MAG[0] + P[3][3]*SH_MAG[2] - P[18][3]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][3]*(2*q0*q3 + 2*q1*q2) - P[20][3]*(2*q0*q2 - 2*q1*q3) - P[2][3]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6])*(P[21][18] + P[1][18]*SH_MAG[0] + P[3][18]*SH_MAG[2] - P[18][18]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][18]*(2*q0*q3 + 2*q1*q2) - P[20][18]*(2*q0*q2 - 2*q1*q3) - P[2][18]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][18]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[19][21]*(2*q0*q3 + 2*q1*q2) - P[20][21]*(2*q0*q2 - 2*q1*q3) - P[2][21]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + (2*q0*q3 + 2*q1*q2)*(P[21][19] + P[1][19]*SH_MAG[0] + P[3][19]*SH_MAG[2] - P[18][19]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][19]*(2*q0*q3 + 2*q1*q2) - P[20][19]*(2*q0*q2 - 2*q1*q3) - P[2][19]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][19]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (2*q0*q2 - 2*q1*q3)*(P[21][20] + P[1][20]*SH_MAG[0] + P[3][20]*SH_MAG[2] - P[18][20]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[19][20]*(2*q0*q3 + 2*q1*q2) - P[20][20]*(2*q0*q2 - 2*q1*q3) - P[2][20]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][20]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[0][21]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MX[1] = SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6];
            SK_MX[2] = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
            SK_MX[3] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MX[4] = 2*q0*q2 - 2*q1*q3;
            SK_MX[5] = 2*q0*q3 + 2*q1*q2;
            K_MAG[0] = SK_MX[0]*(P[0][21] + P[0][1]*SH_MAG[0] + P[0][3]*SH_MAG[2] + P[0][0]*SK_MX[3] - P[0][2]*SK_MX[2] - P[0][18]*SK_MX[1] + P[0][19]*SK_MX[5] - P[0][20]*SK_MX[4]);
            K_MAG[1] = SK_MX[0]*(P[1][21] + P[1][1]*SH_MAG[0] + P[1][3]*SH_MAG[2] + P[1][0]*SK_MX[3] - P[1][2]*SK_MX[2] - P[1][18]*SK_MX[1] + P[1][19]*SK_MX[5] - P[1][20]*SK_MX[4]);
            K_MAG[2] = SK_MX[0]*(P[2][21] + P[2][1]*SH_MAG[0] + P[2][3]*SH_MAG[2] + P[2][0]*SK_MX[3] - P[2][2]*SK_MX[2] - P[2][18]*SK_MX[1] + P[2][19]*SK_MX[5] - P[2][20]*SK_MX[4]);
            K_MAG[3] = SK_MX[0]*(P[3][21] + P[3][1]*SH_MAG[0] + P[3][3]*SH_MAG[2] + P[3][0]*SK_MX[3] - P[3][2]*SK_MX[2] - P[3][18]*SK_MX[1] + P[3][19]*SK_MX[5] - P[3][20]*SK_MX[4]);
            K_MAG[4] = SK_MX[0]*(P[4][21] + P[4][1]*SH_MAG[0] + P[4][3]*SH_MAG[2] + P[4][0]*SK_MX[3] - P[4][2]*SK_MX[2] - P[4][18]*SK_MX[1] + P[4][19]*SK_MX[5] - P[4][20]*SK_MX[4]);
            K_MAG[5] = SK_MX[0]*(P[5][21] + P[5][1]*SH_MAG[0] + P[5][3]*SH_MAG[2] + P[5][0]*SK_MX[3] - P[5][2]*SK_MX[2] - P[5][18]*SK_MX[1] + P[5][19]*SK_MX[5] - P[5][20]*SK_MX[4]);
            K_MAG[6] = SK_MX[0]*(P[6][21] + P[6][1]*SH_MAG[0] + P[6][3]*SH_MAG[2] + P[6][0]*SK_MX[3] - P[6][2]*SK_MX[2] - P[6][18]*SK_MX[1] + P[6][19]*SK_MX[5] - P[6][20]*SK_MX[4]);
            K_MAG[7] = SK_MX[0]*(P[7][21] + P[7][1]*SH_MAG[0] + P[7][3]*SH_MAG[2] + P[7][0]*SK_MX[3] - P[7][2]*SK_MX[2] - P[7][18]*SK_MX[1] + P[7][19]*SK_MX[5] - P[7][20]*SK_MX[4]);
            K_MAG[8] = SK_MX[0]*(P[8][21] + P[8][1]*SH_MAG[0] + P[8][3]*SH_MAG[2] + P[8][0]*SK_MX[3] - P[8][2]*SK_MX[2] - P[8][18]*SK_MX[1] + P[8][19]*SK_MX[5] - P[8][20]*SK_MX[4]);
            K_MAG[9] = SK_MX[0]*(P[9][21] + P[9][1]*SH_MAG[0] + P[9][3]*SH_MAG[2] + P[9][0]*SK_MX[3] - P[9][2]*SK_MX[2] - P[9][18]*SK_MX[1] + P[9][19]*SK_MX[5] - P[9][20]*SK_MX[4]);
            K_MAG[10] = SK_MX[0]*(P[10][21] + P[10][1]*SH_MAG[0] + P[10][3]*SH_MAG[2] + P[10][0]*SK_MX[3] - P[10][2]*SK_MX[2] - P[10][18]*SK_MX[1] + P[10][19]*SK_MX[5] - P[10][20]*SK_MX[4]);
            K_MAG[11] = SK_MX[0]*(P[11][21] + P[11][1]*SH_MAG[0] + P[11][3]*SH_MAG[2] + P[11][0]*SK_MX[3] - P[11][2]*SK_MX[2] - P[11][18]*SK_MX[1] + P[11][19]*SK_MX[5] - P[11][20]*SK_MX[4]);
            K_MAG[12] = SK_MX[0]*(P[12][21] + P[12][1]*SH_MAG[0] + P[12][3]*SH_MAG[2] + P[12][0]*SK_MX[3] - P[12][2]*SK_MX[2] - P[12][18]*SK_MX[1] + P[12][19]*SK_MX[5] - P[12][20]*SK_MX[4]);
            K_MAG[13] = SK_MX[0]*(P[13][21] + P[13][1]*SH_MAG[0] + P[13][3]*SH_MAG[2] + P[13][0]*SK_MX[3] - P[13][2]*SK_MX[2] - P[13][18]*SK_MX[1] + P[13][19]*SK_MX[5] - P[13][20]*SK_MX[4]);
            K_MAG[14] = SK_MX[0]*(P[14][21] + P[14][1]*SH_MAG[0] + P[14][3]*SH_MAG[2] + P[14][0]*SK_MX[3] - P[14][2]*SK_MX[2] - P[14][18]*SK_MX[1] + P[14][19]*SK_MX[5] - P[14][20]*SK_MX[4]);
            K_MAG[15] = SK_MX[0]*(P[15][21] + P[15][1]*SH_MAG[0] + P[15][3]*SH_MAG[2] + P[15][0]*SK_MX[3] - P[15][2]*SK_MX[2] - P[15][18]*SK_MX[1] + P[15][19]*SK_MX[5] - P[15][20]*SK_MX[4]);
            K_MAG[16] = SK_MX[0]*(P[16][21] + P[16][1]*SH_MAG[0] + P[16][3]*SH_MAG[2] + P[16][0]*SK_MX[3] - P[16][2]*SK_MX[2] - P[16][18]*SK_MX[1] + P[16][19]*SK_MX[5] - P[16][20]*SK_MX[4]);
            K_MAG[17] = SK_MX[0]*(P[17][21] + P[17][1]*SH_MAG[0] + P[17][3]*SH_MAG[2] + P[17][0]*SK_MX[3] - P[17][2]*SK_MX[2] - P[17][18]*SK_MX[1] + P[17][19]*SK_MX[5] - P[17][20]*SK_MX[4]);
            K_MAG[18] = SK_MX[0]*(P[18][21] + P[18][1]*SH_MAG[0] + P[18][3]*SH_MAG[2] + P[18][0]*SK_MX[3] - P[18][2]*SK_MX[2] - P[18][18]*SK_MX[1] + P[18][19]*SK_MX[5] - P[18][20]*SK_MX[4]);
            K_MAG[19] = SK_MX[0]*(P[19][21] + P[19][1]*SH_MAG[0] + P[19][3]*SH_MAG[2] + P[19][0]*SK_MX[3] - P[19][2]*SK_MX[2] - P[19][18]*SK_MX[1] + P[19][19]*SK_MX[5] - P[19][20]*SK_MX[4]);
            K_MAG[20] = SK_MX[0]*(P[20][21] + P[20][1]*SH_MAG[0] + P[20][3]*SH_MAG[2] + P[20][0]*SK_MX[3] - P[20][2]*SK_MX[2] - P[20][18]*SK_MX[1] + P[20][19]*SK_MX[5] - P[20][20]*SK_MX[4]);
            K_MAG[21] = SK_MX[0]*(P[21][21] + P[21][1]*SH_MAG[0] + P[21][3]*SH_MAG[2] + P[21][0]*SK_MX[3] - P[21][2]*SK_MX[2] - P[21][18]*SK_MX[1] + P[21][19]*SK_MX[5] - P[21][20]*SK_MX[4]);
            K_MAG[22] = SK_MX[0]*(P[22][21] + P[22][1]*SH_MAG[0] + P[22][3]*SH_MAG[2] + P[22][0]*SK_MX[3] - P[22][2]*SK_MX[2] - P[22][18]*SK_MX[1] + P[22][19]*SK_MX[5] - P[22][20]*SK_MX[4]);
            K_MAG[23] = SK_MX[0]*(P[23][21] + P[23][1]*SH_MAG[0] + P[23][3]*SH_MAG[2] + P[23][0]*SK_MX[3] - P[23][2]*SK_MX[2] - P[23][18]*SK_MX[1] + P[23][19]*SK_MX[5] - P[23][20]*SK_MX[4]);
            varInnov[0] = 1.0/SK_MX[0];
            // reset the observation index to 0 (we start by fusing the X
            // measurement)
            obsIndex = 0;
        }
        else if (obsIndex == 1) // we are now fusing the Y measurement
        {
            // Calculate observation jacobians
            H_MAG[0] = SH_MAG[2];
            H_MAG[1] = SH_MAG[1];
            H_MAG[2] = SH_MAG[0];
            H_MAG[3] = 2*magD*q2 - SH_MAG[8] - SH_MAG[7];
            H_MAG[18] = 2*q1*q2 - 2*q0*q3;
            H_MAG[19] = SH_MAG[4] - SH_MAG[3] - SH_MAG[5] + SH_MAG[6];
            H_MAG[20] = 2*q0*q1 + 2*q2*q3;
            H_MAG[22] = 1;
            // Calculate Kalman gain
            SK_MY[0] = 1/(P[22][22] + R_MAG + P[0][22]*SH_MAG[2] + P[1][22]*SH_MAG[1] + P[2][22]*SH_MAG[0] - P[19][22]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - (2*q0*q3 - 2*q1*q2)*(P[22][18] + P[0][18]*SH_MAG[2] + P[1][18]*SH_MAG[1] + P[2][18]*SH_MAG[0] - P[19][18]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][18]*(2*q0*q3 - 2*q1*q2) + P[20][18]*(2*q0*q1 + 2*q2*q3) - P[3][18]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (2*q0*q1 + 2*q2*q3)*(P[22][20] + P[0][20]*SH_MAG[2] + P[1][20]*SH_MAG[1] + P[2][20]*SH_MAG[0] - P[19][20]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][20]*(2*q0*q3 - 2*q1*q2) + P[20][20]*(2*q0*q1 + 2*q2*q3) - P[3][20]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[22][3] + P[0][3]*SH_MAG[2] + P[1][3]*SH_MAG[1] + P[2][3]*SH_MAG[0] - P[19][3]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][3]*(2*q0*q3 - 2*q1*q2) + P[20][3]*(2*q0*q1 + 2*q2*q3) - P[3][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - P[18][22]*(2*q0*q3 - 2*q1*q2) + P[20][22]*(2*q0*q1 + 2*q2*q3) + SH_MAG[2]*(P[22][0] + P[0][0]*SH_MAG[2] + P[1][0]*SH_MAG[1] + P[2][0]*SH_MAG[0] - P[19][0]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][0]*(2*q0*q3 - 2*q1*q2) + P[20][0]*(2*q0*q1 + 2*q2*q3) - P[3][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[1]*(P[22][1] + P[0][1]*SH_MAG[2] + P[1][1]*SH_MAG[1] + P[2][1]*SH_MAG[0] - P[19][1]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][1]*(2*q0*q3 - 2*q1*q2) + P[20][1]*(2*q0*q1 + 2*q2*q3) - P[3][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[22][2] + P[0][2]*SH_MAG[2] + P[1][2]*SH_MAG[1] + P[2][2]*SH_MAG[0] - P[19][2]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][2]*(2*q0*q3 - 2*q1*q2) + P[20][2]*(2*q0*q1 + 2*q2*q3) - P[3][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6])*(P[22][19] + P[0][19]*SH_MAG[2] + P[1][19]*SH_MAG[1] + P[2][19]*SH_MAG[0] - P[19][19]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[18][19]*(2*q0*q3 - 2*q1*q2) + P[20][19]*(2*q0*q1 + 2*q2*q3) - P[3][19]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - P[3][22]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MY[1] = SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6];
            SK_MY[2] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MY[3] = 2*q0*q3 - 2*q1*q2;
            SK_MY[4] = 2*q0*q1 + 2*q2*q3;
            K_MAG[0] = SK_MY[0]*(P[0][22] + P[0][0]*SH_MAG[2] + P[0][1]*SH_MAG[1] + P[0][2]*SH_MAG[0] - P[0][3]*SK_MY[2] - P[0][19]*SK_MY[1] - P[0][18]*SK_MY[3] + P[0][20]*SK_MY[4]);
            K_MAG[1] = SK_MY[0]*(P[1][22] + P[1][0]*SH_MAG[2] + P[1][1]*SH_MAG[1] + P[1][2]*SH_MAG[0] - P[1][3]*SK_MY[2] - P[1][19]*SK_MY[1] - P[1][18]*SK_MY[3] + P[1][20]*SK_MY[4]);
            K_MAG[2] = SK_MY[0]*(P[2][22] + P[2][0]*SH_MAG[2] + P[2][1]*SH_MAG[1] + P[2][2]*SH_MAG[0] - P[2][3]*SK_MY[2] - P[2][19]*SK_MY[1] - P[2][18]*SK_MY[3] + P[2][20]*SK_MY[4]);
            K_MAG[3] = SK_MY[0]*(P[3][22] + P[3][0]*SH_MAG[2] + P[3][1]*SH_MAG[1] + P[3][2]*SH_MAG[0] - P[3][3]*SK_MY[2] - P[3][19]*SK_MY[1] - P[3][18]*SK_MY[3] + P[3][20]*SK_MY[4]);
            K_MAG[4] = SK_MY[0]*(P[4][22] + P[4][0]*SH_MAG[2] + P[4][1]*SH_MAG[1] + P[4][2]*SH_MAG[0] - P[4][3]*SK_MY[2] - P[4][19]*SK_MY[1] - P[4][18]*SK_MY[3] + P[4][20]*SK_MY[4]);
            K_MAG[5] = SK_MY[0]*(P[5][22] + P[5][0]*SH_MAG[2] + P[5][1]*SH_MAG[1] + P[5][2]*SH_MAG[0] - P[5][3]*SK_MY[2] - P[5][19]*SK_MY[1] - P[5][18]*SK_MY[3] + P[5][20]*SK_MY[4]);
            K_MAG[6] = SK_MY[0]*(P[6][22] + P[6][0]*SH_MAG[2] + P[6][1]*SH_MAG[1] + P[6][2]*SH_MAG[0] - P[6][3]*SK_MY[2] - P[6][19]*SK_MY[1] - P[6][18]*SK_MY[3] + P[6][20]*SK_MY[4]);
            K_MAG[7] = SK_MY[0]*(P[7][22] + P[7][0]*SH_MAG[2] + P[7][1]*SH_MAG[1] + P[7][2]*SH_MAG[0] - P[7][3]*SK_MY[2] - P[7][19]*SK_MY[1] - P[7][18]*SK_MY[3] + P[7][20]*SK_MY[4]);
            K_MAG[8] = SK_MY[0]*(P[8][22] + P[8][0]*SH_MAG[2] + P[8][1]*SH_MAG[1] + P[8][2]*SH_MAG[0] - P[8][3]*SK_MY[2] - P[8][19]*SK_MY[1] - P[8][18]*SK_MY[3] + P[8][20]*SK_MY[4]);
            K_MAG[9] = SK_MY[0]*(P[9][22] + P[9][0]*SH_MAG[2] + P[9][1]*SH_MAG[1] + P[9][2]*SH_MAG[0] - P[9][3]*SK_MY[2] - P[9][19]*SK_MY[1] - P[9][18]*SK_MY[3] + P[9][20]*SK_MY[4]);
            K_MAG[10] = SK_MY[0]*(P[10][22] + P[10][0]*SH_MAG[2] + P[10][1]*SH_MAG[1] + P[10][2]*SH_MAG[0] - P[10][3]*SK_MY[2] - P[10][19]*SK_MY[1] - P[10][18]*SK_MY[3] + P[10][20]*SK_MY[4]);
            K_MAG[11] = SK_MY[0]*(P[11][22] + P[11][0]*SH_MAG[2] + P[11][1]*SH_MAG[1] + P[11][2]*SH_MAG[0] - P[11][3]*SK_MY[2] - P[11][19]*SK_MY[1] - P[11][18]*SK_MY[3] + P[11][20]*SK_MY[4]);
            K_MAG[12] = SK_MY[0]*(P[12][22] + P[12][0]*SH_MAG[2] + P[12][1]*SH_MAG[1] + P[12][2]*SH_MAG[0] - P[12][3]*SK_MY[2] - P[12][19]*SK_MY[1] - P[12][18]*SK_MY[3] + P[12][20]*SK_MY[4]);
            K_MAG[13] = SK_MY[0]*(P[13][22] + P[13][0]*SH_MAG[2] + P[13][1]*SH_MAG[1] + P[13][2]*SH_MAG[0] - P[13][3]*SK_MY[2] - P[13][19]*SK_MY[1] - P[13][18]*SK_MY[3] + P[13][20]*SK_MY[4]);
            K_MAG[14] = SK_MY[0]*(P[14][22] + P[14][0]*SH_MAG[2] + P[14][1]*SH_MAG[1] + P[14][2]*SH_MAG[0] - P[14][3]*SK_MY[2] - P[14][19]*SK_MY[1] - P[14][18]*SK_MY[3] + P[14][20]*SK_MY[4]);
            K_MAG[15] = SK_MY[0]*(P[15][22] + P[15][0]*SH_MAG[2] + P[15][1]*SH_MAG[1] + P[15][2]*SH_MAG[0] - P[15][3]*SK_MY[2] - P[15][19]*SK_MY[1] - P[15][18]*SK_MY[3] + P[15][20]*SK_MY[4]);
            K_MAG[16] = SK_MY[0]*(P[16][22] + P[16][0]*SH_MAG[2] + P[16][1]*SH_MAG[1] + P[16][2]*SH_MAG[0] - P[16][3]*SK_MY[2] - P[16][19]*SK_MY[1] - P[16][18]*SK_MY[3] + P[16][20]*SK_MY[4]);
            K_MAG[17] = SK_MY[0]*(P[17][22] + P[17][0]*SH_MAG[2] + P[17][1]*SH_MAG[1] + P[17][2]*SH_MAG[0] - P[17][3]*SK_MY[2] - P[17][19]*SK_MY[1] - P[17][18]*SK_MY[3] + P[17][20]*SK_MY[4]);
            K_MAG[18] = SK_MY[0]*(P[18][22] + P[18][0]*SH_MAG[2] + P[18][1]*SH_MAG[1] + P[18][2]*SH_MAG[0] - P[18][3]*SK_MY[2] - P[18][19]*SK_MY[1] - P[18][18]*SK_MY[3] + P[18][20]*SK_MY[4]);
            K_MAG[19] = SK_MY[0]*(P[19][22] + P[19][0]*SH_MAG[2] + P[19][1]*SH_MAG[1] + P[19][2]*SH_MAG[0] - P[19][3]*SK_MY[2] - P[19][19]*SK_MY[1] - P[19][18]*SK_MY[3] + P[19][20]*SK_MY[4]);
            K_MAG[20] = SK_MY[0]*(P[20][22] + P[20][0]*SH_MAG[2] + P[20][1]*SH_MAG[1] + P[20][2]*SH_MAG[0] - P[20][3]*SK_MY[2] - P[20][19]*SK_MY[1] - P[20][18]*SK_MY[3] + P[20][20]*SK_MY[4]);
            K_MAG[21] = SK_MY[0]*(P[21][22] + P[21][0]*SH_MAG[2] + P[21][1]*SH_MAG[1] + P[21][2]*SH_MAG[0] - P[21][3]*SK_MY[2] - P[21][19]*SK_MY[1] - P[21][18]*SK_MY[3] + P[21][20]*SK_MY[4]);
            K_MAG[22] = SK_MY[0]*(P[22][22] + P[22][0]*SH_MAG[2] + P[22][1]*SH_MAG[1] + P[22][2]*SH_MAG[0] - P[22][3]*SK_MY[2] - P[22][19]*SK_MY[1] - P[22][18]*SK_MY[3] + P[22][20]*SK_MY[4]);
            K_MAG[23] = SK_MY[0]*(P[23][22] + P[23][0]*SH_MAG[2] + P[23][1]*SH_MAG[1] + P[23][2]*SH_MAG[0] - P[23][3]*SK_MY[2] - P[23][19]*SK_MY[1] - P[23][18]*SK_MY[3] + P[23][20]*SK_MY[4]);
            varInnov[1] = 1.0/SK_MY[0];
        }
        else if (obsIndex == 2) // we are now fusing the Z measurement
        {
            // Calculate observation jacobians
            H_MAG[0] = SH_MAG[1];
            H_MAG[1] = 2*magN*q3 - 2*magE*q0 - 2*magD*q1;
            H_MAG[2] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            H_MAG[3] = SH_MAG[0];
            H_MAG[18] = 2*q0*q2 + 2*q1*q3;
            H_MAG[19] = 2*q2*q3 - 2*q0*q1;
            H_MAG[20] = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
            H_MAG[23] = 1;
            // Calculate Kalman gain
            SK_MZ[0] = 1/(P[23][23] + R_MAG + P[0][23]*SH_MAG[1] + P[3][23]*SH_MAG[0] + P[20][23]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) - (2*magD*q1 + 2*magE*q0 - 2*magN*q3)*(P[23][1] + P[0][1]*SH_MAG[1] + P[3][1]*SH_MAG[0] + P[20][1]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][1]*(2*q0*q2 + 2*q1*q3) - P[19][1]*(2*q0*q1 - 2*q2*q3) - P[1][1]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[23][2] + P[0][2]*SH_MAG[1] + P[3][2]*SH_MAG[0] + P[20][2]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][2]*(2*q0*q2 + 2*q1*q3) - P[19][2]*(2*q0*q1 - 2*q2*q3) - P[1][2]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[1]*(P[23][0] + P[0][0]*SH_MAG[1] + P[3][0]*SH_MAG[0] + P[20][0]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][0]*(2*q0*q2 + 2*q1*q3) - P[19][0]*(2*q0*q1 - 2*q2*q3) - P[1][0]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[23][3] + P[0][3]*SH_MAG[1] + P[3][3]*SH_MAG[0] + P[20][3]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][3]*(2*q0*q2 + 2*q1*q3) - P[19][3]*(2*q0*q1 - 2*q2*q3) - P[1][3]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6])*(P[23][20] + P[0][20]*SH_MAG[1] + P[3][20]*SH_MAG[0] + P[20][20]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][20]*(2*q0*q2 + 2*q1*q3) - P[19][20]*(2*q0*q1 - 2*q2*q3) - P[1][20]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][20]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[18][23]*(2*q0*q2 + 2*q1*q3) - P[19][23]*(2*q0*q1 - 2*q2*q3) - P[1][23]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + (2*q0*q2 + 2*q1*q3)*(P[23][18] + P[0][18]*SH_MAG[1] + P[3][18]*SH_MAG[0] + P[20][18]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][18]*(2*q0*q2 + 2*q1*q3) - P[19][18]*(2*q0*q1 - 2*q2*q3) - P[1][18]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][18]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (2*q0*q1 - 2*q2*q3)*(P[23][19] + P[0][19]*SH_MAG[1] + P[3][19]*SH_MAG[0] + P[20][19]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[18][19]*(2*q0*q2 + 2*q1*q3) - P[19][19]*(2*q0*q1 - 2*q2*q3) - P[1][19]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][19]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[2][23]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MZ[1] = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
            SK_MZ[2] = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
            SK_MZ[3] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MZ[4] = 2*q0*q1 - 2*q2*q3;
            SK_MZ[5] = 2*q0*q2 + 2*q1*q3;
            K_MAG[0] = SK_MZ[0]*(P[0][23] + P[0][0]*SH_MAG[1] + P[0][3]*SH_MAG[0] - P[0][1]*SK_MZ[2] + P[0][2]*SK_MZ[3] + P[0][20]*SK_MZ[1] + P[0][18]*SK_MZ[5] - P[0][19]*SK_MZ[4]);
            K_MAG[1] = SK_MZ[0]*(P[1][23] + P[1][0]*SH_MAG[1] + P[1][3]*SH_MAG[0] - P[1][1]*SK_MZ[2] + P[1][2]*SK_MZ[3] + P[1][20]*SK_MZ[1] + P[1][18]*SK_MZ[5] - P[1][19]*SK_MZ[4]);
            K_MAG[2] = SK_MZ[0]*(P[2][23] + P[2][0]*SH_MAG[1] + P[2][3]*SH_MAG[0] - P[2][1]*SK_MZ[2] + P[2][2]*SK_MZ[3] + P[2][20]*SK_MZ[1] + P[2][18]*SK_MZ[5] - P[2][19]*SK_MZ[4]);
            K_MAG[3] = SK_MZ[0]*(P[3][23] + P[3][0]*SH_MAG[1] + P[3][3]*SH_MAG[0] - P[3][1]*SK_MZ[2] + P[3][2]*SK_MZ[3] + P[3][20]*SK_MZ[1] + P[3][18]*SK_MZ[5] - P[3][19]*SK_MZ[4]);
            K_MAG[4] = SK_MZ[0]*(P[4][23] + P[4][0]*SH_MAG[1] + P[4][3]*SH_MAG[0] - P[4][1]*SK_MZ[2] + P[4][2]*SK_MZ[3] + P[4][20]*SK_MZ[1] + P[4][18]*SK_MZ[5] - P[4][19]*SK_MZ[4]);
            K_MAG[5] = SK_MZ[0]*(P[5][23] + P[5][0]*SH_MAG[1] + P[5][3]*SH_MAG[0] - P[5][1]*SK_MZ[2] + P[5][2]*SK_MZ[3] + P[5][20]*SK_MZ[1] + P[5][18]*SK_MZ[5] - P[5][19]*SK_MZ[4]);
            K_MAG[6] = SK_MZ[0]*(P[6][23] + P[6][0]*SH_MAG[1] + P[6][3]*SH_MAG[0] - P[6][1]*SK_MZ[2] + P[6][2]*SK_MZ[3] + P[6][20]*SK_MZ[1] + P[6][18]*SK_MZ[5] - P[6][19]*SK_MZ[4]);
            K_MAG[7] = SK_MZ[0]*(P[7][23] + P[7][0]*SH_MAG[1] + P[7][3]*SH_MAG[0] - P[7][1]*SK_MZ[2] + P[7][2]*SK_MZ[3] + P[7][20]*SK_MZ[1] + P[7][18]*SK_MZ[5] - P[7][19]*SK_MZ[4]);
            K_MAG[8] = SK_MZ[0]*(P[8][23] + P[8][0]*SH_MAG[1] + P[8][3]*SH_MAG[0] - P[8][1]*SK_MZ[2] + P[8][2]*SK_MZ[3] + P[8][20]*SK_MZ[1] + P[8][18]*SK_MZ[5] - P[8][19]*SK_MZ[4]);
            K_MAG[9] = SK_MZ[0]*(P[9][23] + P[9][0]*SH_MAG[1] + P[9][3]*SH_MAG[0] - P[9][1]*SK_MZ[2] + P[9][2]*SK_MZ[3] + P[9][20]*SK_MZ[1] + P[9][18]*SK_MZ[5] - P[9][19]*SK_MZ[4]);
            K_MAG[10] = SK_MZ[0]*(P[10][23] + P[10][0]*SH_MAG[1] + P[10][3]*SH_MAG[0] - P[10][1]*SK_MZ[2] + P[10][2]*SK_MZ[3] + P[10][20]*SK_MZ[1] + P[10][18]*SK_MZ[5] - P[10][19]*SK_MZ[4]);
            K_MAG[11] = SK_MZ[0]*(P[11][23] + P[11][0]*SH_MAG[1] + P[11][3]*SH_MAG[0] - P[11][1]*SK_MZ[2] + P[11][2]*SK_MZ[3] + P[11][20]*SK_MZ[1] + P[11][18]*SK_MZ[5] - P[11][19]*SK_MZ[4]);
            K_MAG[12] = SK_MZ[0]*(P[12][23] + P[12][0]*SH_MAG[1] + P[12][3]*SH_MAG[0] - P[12][1]*SK_MZ[2] + P[12][2]*SK_MZ[3] + P[12][20]*SK_MZ[1] + P[12][18]*SK_MZ[5] - P[12][19]*SK_MZ[4]);
            K_MAG[13] = SK_MZ[0]*(P[13][23] + P[13][0]*SH_MAG[1] + P[13][3]*SH_MAG[0] - P[13][1]*SK_MZ[2] + P[13][2]*SK_MZ[3] + P[13][20]*SK_MZ[1] + P[13][18]*SK_MZ[5] - P[13][19]*SK_MZ[4]);
            K_MAG[14] = SK_MZ[0]*(P[14][23] + P[14][0]*SH_MAG[1] + P[14][3]*SH_MAG[0] - P[14][1]*SK_MZ[2] + P[14][2]*SK_MZ[3] + P[14][20]*SK_MZ[1] + P[14][18]*SK_MZ[5] - P[14][19]*SK_MZ[4]);
            K_MAG[15] = SK_MZ[0]*(P[15][23] + P[15][0]*SH_MAG[1] + P[15][3]*SH_MAG[0] - P[15][1]*SK_MZ[2] + P[15][2]*SK_MZ[3] + P[15][20]*SK_MZ[1] + P[15][18]*SK_MZ[5] - P[15][19]*SK_MZ[4]);
            K_MAG[16] = SK_MZ[0]*(P[16][23] + P[16][0]*SH_MAG[1] + P[16][3]*SH_MAG[0] - P[16][1]*SK_MZ[2] + P[16][2]*SK_MZ[3] + P[16][20]*SK_MZ[1] + P[16][18]*SK_MZ[5] - P[16][19]*SK_MZ[4]);
            K_MAG[17] = SK_MZ[0]*(P[17][23] + P[17][0]*SH_MAG[1] + P[17][3]*SH_MAG[0] - P[17][1]*SK_MZ[2] + P[17][2]*SK_MZ[3] + P[17][20]*SK_MZ[1] + P[17][18]*SK_MZ[5] - P[17][19]*SK_MZ[4]);
            K_MAG[18] = SK_MZ[0]*(P[18][23] + P[18][0]*SH_MAG[1] + P[18][3]*SH_MAG[0] - P[18][1]*SK_MZ[2] + P[18][2]*SK_MZ[3] + P[18][20]*SK_MZ[1] + P[18][18]*SK_MZ[5] - P[18][19]*SK_MZ[4]);
            K_MAG[19] = SK_MZ[0]*(P[19][23] + P[19][0]*SH_MAG[1] + P[19][3]*SH_MAG[0] - P[19][1]*SK_MZ[2] + P[19][2]*SK_MZ[3] + P[19][20]*SK_MZ[1] + P[19][18]*SK_MZ[5] - P[19][19]*SK_MZ[4]);
            K_MAG[20] = SK_MZ[0]*(P[20][23] + P[20][0]*SH_MAG[1] + P[20][3]*SH_MAG[0] - P[20][1]*SK_MZ[2] + P[20][2]*SK_MZ[3] + P[20][20]*SK_MZ[1] + P[20][18]*SK_MZ[5] - P[20][19]*SK_MZ[4]);
            K_MAG[21] = SK_MZ[0]*(P[21][23] + P[21][0]*SH_MAG[1] + P[21][3]*SH_MAG[0] - P[21][1]*SK_MZ[2] + P[21][2]*SK_MZ[3] + P[21][20]*SK_MZ[1] + P[21][18]*SK_MZ[5] - P[21][19]*SK_MZ[4]);
            K_MAG[22] = SK_MZ[0]*(P[22][23] + P[22][0]*SH_MAG[1] + P[22][3]*SH_MAG[0] - P[22][1]*SK_MZ[2] + P[22][2]*SK_MZ[3] + P[22][20]*SK_MZ[1] + P[22][18]*SK_MZ[5] - P[22][19]*SK_MZ[4]);
            K_MAG[23] = SK_MZ[0]*(P[23][23] + P[23][0]*SH_MAG[1] + P[23][3]*SH_MAG[0] - P[23][1]*SK_MZ[2] + P[23][2]*SK_MZ[3] + P[23][20]*SK_MZ[1] + P[23][18]*SK_MZ[5] - P[23][19]*SK_MZ[4]);
            varInnov[2] = 1.0/SK_MZ[0];
        }
        // Calculate the measurement innovation
        innovation[obsIndex] = MagPred[obsIndex] - MagData[obsIndex];
        // Check the innovation for consistency and don't fuse if > 5Sigma
        if ((innovation[obsIndex]*innovation[obsIndex]/varInnov[obsIndex]) < 25.0)
        {
            // correct the state vector
            for (j= 0; j<=23; j++)
            {
                states[j] = states[j] - K_MAG[j] * innovation[obsIndex];
            }
            // normalise the quaternion states
            float quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
            if (quatMag > 1e-12)
            {
                for (j= 0; j<=3; j++)
                {
                    float quatMagInv = 1.0/quatMag;
                    states[j] = states[j] * quatMagInv;
                }
            }
            // correct the covariance P = (I - K*H)*P
            // take advantage of the empty columns in KH to reduce the
            // number of operations
            for (i = 0; i<=23; i++)
            {
                for (j = 0; j<=3; j++)
                {
                    KH[i][j] = K_MAG[i] * H_MAG[j];
                }
                for (j = 18; j<=23; j++)
                {
                    KH[i][j] = K_MAG[i] * H_MAG[j];
                }
            }
            for (i = 0; i<=23; i++)
            {
                for (j = 0; j<=23; j++)
                {
                    for (k = 0; k<=3; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][j] * P[k][j];
                    }
                    for (k = 18; k<=23; k++)
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
            obsIndex = obsIndex + 1;
        }
    }
}

void FuseAirspeed(
    float innovation, // innovation output
    float varInnov, // innovation variance output
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
    uint8_t i;
    uint8_t j;
    uint8_t k;
    const float R_TAS = 2.0;
    float SH_TAS[3];
    float SK_TAS[2];
    float H_TAS[24];
    float K_TAS[24];
    float KH[24][24];
    float KHP[24][24];
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
    }
}

void zeroRows(float covMat[24][24], uint8_t first, uint8_t last)
{
    uint8_t row;
    uint8_t col;
    for (row=first; row<=last; row++)
    {
        for (col=0; col<=23; col++)
        {
            covMat[row][col] = 0.0;
        }
    }
}

void zeroCols(float covMat[24][24], uint8_t first, uint8_t last)
{
    uint8_t row;
    uint8_t col;
    for (col=first; col<=last; col++)
    {
        for (row=1; row<=24; row++)
        {
            covMat[row][col] = 0.0;
        }
    }
}

float sq(float valIn)
{
    return valIn*valIn;
}

// Store states in a history array along with time stamp
void StoreStates(uint32_t msec)
{
    static uint8_t storeIndex = 0;
    uint8_t i;
    if (storeIndex > 49) storeIndex = 0;
    for (i=0; i<=23; i++) storedStates[i][storeIndex] = states[i];
    statetimeStamp[storeIndex] = msec;
    storeIndex = storeIndex + 1;
}

// Output the state vector stored at the time that best matches that specified by msec
void RecallStates(float statesForFusion[24], uint32_t msec)
{
    long int timeDelta;
    long int bestTimeDelta = 200;
    uint8_t storeIndex;
    uint8_t bestStoreIndex = 0;
    uint8_t i;
    for (storeIndex=0; storeIndex<=49; storeIndex++)
    {
        timeDelta = msec - statetimeStamp[storeIndex];
        if (timeDelta < 0) timeDelta = -timeDelta;
        if (timeDelta < bestTimeDelta)
        {
            bestStoreIndex = storeIndex;
            bestTimeDelta = timeDelta;
        }
    }
    if (bestTimeDelta <= 200) // only output stored state if < 200 msec retrieval error
    {
        for (i=0; i<=23; i++) statesForFusion[i] = storedStates[i][bestStoreIndex];
    }
    else // otherwise output current state
    {
        for (i=0; i<=23; i++) statesForFusion[i] = states[i];
    }
}
