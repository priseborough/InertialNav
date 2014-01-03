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

    float length(void) const;
    Vector3f zero(void) const;
};

float Vector3f::length(void) const
{
    return sqrt(x*x + y*y + z*z);
}

Vector3f Vector3f::zero(void) const
{
    Vector3f ret = *this;
    ret.x = 0.0;
    ret.y = 0.0;
    ret.z = 0.0;
    return ret;
}

class Mat3f
{
private:
public:
    Vector3f x;
    Vector3f y;
    Vector3f z;

    Mat3f transpose(void) const;
};

static void swap_var(float &d1, float &d2)
{
    float tmp = d1;
    d1 = d2;
    d2 = tmp;
}

Mat3f Mat3f::transpose(void) const
{
    Mat3f ret = *this;
    swap_var(ret.x.y, ret.y.x);
    swap_var(ret.x.z, ret.z.x);
    swap_var(ret.y.z, ret.z.y);
    return ret;
}

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

// overload % operator to provide a vector cross product
Vector3f operator%( Vector3f vecIn1, Vector3f vecIn2)
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

void  UpdateStrapdownEquationsNED();

void CovariancePrediction();

void FuseVelposNED();

void FuseMagnetometer();

void FuseAirspeed();

void zeroRows(float covMat[21][21], uint8_t first, uint8_t last);

void zeroCols(float covMat[21][21], uint8_t first, uint8_t last);

float sq(float valIn);

void quatNorm(float quatOut[4], float quatIn[4]);

uint32_t millis();

// store staes along with system time stamp in msces
void StoreStates();

// recall stste vector stored at closest time to the one specified by msec
void RecallStates(float statesForFusion[21], uint32_t msec);

void quat2Tnb(Mat3f &Tnb, float quat[4]);

void quat2Tbn(Mat3f &Tbn, float quat[4]);

void calcEarthRateNED(Vector3f &omega, float latitude);

void eul2quat(float quat[4], float eul[3]);

void quat2eul(float eul[3],float quat[4]);

void calcvelNED(float velNED[3], float gpsCourse, float gpsGndSpd, float gpsVelD);

void calcposNED(float posNED[3], float lat, float lon, float hgt, float latRef, float lonRef, float hgtRef);

void calcLLH(float posNED[3], float lat, float lon, float hgt, float latRef, float lonRef, float hgtRef);

void OnGroundCheck();

void CovarianceInit();

void readIMUData();

void readGpsData();

void readMagData();

void readAirSpdData();

void readAhrsData();

void readTimingData();

void InitialiseFilter();

void WriteFilterOutput();

void CloseFiles();

// Global variables
#define GRAVITY_MSS 9.80665f
#define deg2rad 0.017453292f
#define rad2deg 57.295780f
#define pi 3.141592657f
#define earthRate 0.000072921f
#define earthRadius 6378145.0f
#define earthRadiusInv  1.5678540e-7f
float KH[21][21]; //  intermediate result used for covariance updates
float KHP[21][21]; // intermediate result used for covariance updates
static float P[21][21]; // covariance matrix
float Kfusion[21]; // Kalman gains
static float states[21]; // state matrix
static float storedStates[21][50]; // state vectors stored for the last 50 time steps
uint32_t statetimeStamp[50]; // time stamp for each state vector stored
Vector3f correctedDelAng; // delta angles about the xyz body axes corrected for errors (rad)
Vector3f summedDelAng; // summed delta angles about the xyz body axes corrected for errors (rad)
Vector3f summedDelVel; // summed delta velocities along the XYZ body axes corrected for errors (m/s)
float accNavMag; // magnitude of navigation accel (- used to adjust GPS obs variance (m/s^2)
Vector3f earthRateNED; // earths angular rate vector in NED (rad/s)
Vector3f angRate; // angular rate vector in XYZ body axes measured by the IMU (rad/s)
Vector3f accel; // acceleration vector in XYZ body axes measured by the IMU (m/s^2)
Vector3f dVelIMU;
Vector3f dAngIMU;
float dtIMU; // time lapsed since the last IMU measurement or covariance update (sec)
float dt; // time lapsed since last covariance prediction
bool onGround = true; // boolean true when the flight vehicle is on the ground (not flying)
bool useAirspeed = true; // boolean true if airspeed data is being used
bool useCompass = true; // boolean true if magnetometer data is being used
uint8_t fusionModeGPS = 0; // 0 = GPS outputs 3D velocity, 1 = GPS outputs 2D velocity, 2 = GPS outputs no velocity
float innovVelPos[6]; // innovation output
float varInnovVelPos[6]; // innovation variance output
bool fuseVelData = false; // this boolean causes the posNE and velNED obs to be fused
bool fusePosData = false; // this boolean causes the posNE and velNED obs to be fused
bool fuseHgtData = false; // this boolean causes the hgtMea obs to be fused
float velNED[3]; // North, East, Down velocity obs (m/s)
float posNE[2]; // North, East position obs (m)
float hgtMea; //  measured height (m)
float posNED[3]; // North, East Down position (m)
float statesAtVelTime[21]; // States at the effective measurement time for posNE and velNED measurements
float statesAtPosTime[21]; // States at the effective measurement time for posNE and velNED measurements
float statesAtHgtTime[21]; // States at the effective measurement time for the hgtMea measurement
float innovMag[3]; // innovation output
float varInnovMag[3]; // innovation variance output
bool fuseMagData = false; // boolean true when magnetometer data is to be fused
Vector3f magData; // magnetometer flux radings in X,Y,Z body axes
float statesAtMagMeasTime[21]; // filter satates at the effective measurement time
float innovVtas; // innovation output
float varInnovVtas; // innovation variance output
bool fuseVtasData = false; // boolean true when airspeed data is to be fused
float VtasMeas; // true airspeed measurement (m/s)
float statesAtVtasMeasTime[21]; // filter states at the effective measurement time
float latRef; // WGS-84 latitude of reference point (rad)
float lonRef; // WGS-84 longitude of reference point (rad)
float hgtRef; // WGS-84 height of reference point (m)
Vector3f magBias; // states representing magnetometer bias vector in XYZ body axes
float eulerEst[3]; // Euler angles calculated from filter states
float eulerDif[3]; // difference between Euler angle estimated by EKF and the AHRS solution
uint8_t covSkipCount = 0; // Number of state prediction frames (IMU daya updates to skip before doing the covariance prediction
const float covTimeStepMax = 0.07f; // maximum time allowed between covariance predictions
const float covDelAngMax = 0.05f; // maximum delta angle between covariance predictions
float EAS2TAS = 1.0f; // ratio f true to equivalent airspeed
bool endOfData = false; //boolean set to true when all files have returned data

// Estimated time delays (msec)
uint32_t msecVelDelay = 230;
uint32_t msecPosDelay = 210;
uint32_t msecHgtDelay = 350;
uint32_t msecMagDelay = 30;
uint32_t msecTasDelay = 210;

// IMU input data variables
float imuIn;
float tempImu[8];
float IMUtimestamp;
static uint32_t IMUmsec = 0;

// GPS input data variables
float gpsCourse;
float gpsGndSpd;
float gpsVelD;
float gpsLat;
float gpsLon;
float gpsHgt;
bool newDataGps;
uint8_t GPSstatus;

// Magnetometer input data variables
float magIn;
float tempMag[8];
float tempMagPrev[8];
float MAGtimestamp = 0;
uint32_t MAGmsec = 0;
uint32_t lastMAGmsec = 0;
bool newDataMag = false;

// AHRS input data variables
float ahrsIn;
float tempAhrs[7];
float tempAhrsPrev[7];
float AHRStimestamp = 0;
uint32_t AHRSmsec = 0;
uint32_t lastAHRStime = 0;
float ahrsEul[3];
float ahrsErrorRP;
float ahrsErrorYaw;

// ADS input data variables
float adsIn;
float tempAds[10];
float tempAdsPrev[10];
float ADStimestamp = 0;
uint32_t ADSmsec = 0;
uint32_t lastADSmsec = 0;
float Veas;
bool newAdsData = false;

// input data timing
uint32_t msecAlignTime;
uint32_t msecStartTime;
uint32_t msecEndTime;

// Data file identifiers
FILE * pImuFile;
FILE * pMagFile;
FILE * pGpsFile;
FILE * pAhrsFile;
FILE * pAdsFile;
FILE * pStateOutFile;
FILE * pEulOutFile;
FILE * pCovOutFile;
FILE * pRefPosVelOutFile;
FILE * pVelPosFuseFile;
FILE * pMagFuseFile;
FILE * pTasFuseFile;
FILE * pTimeFile;

bool statesInitialised = false;

int main()
{

    // open data files
    pImuFile = fopen ("IMU.txt","r");
    pMagFile = fopen ("MAG.txt","r");
    pGpsFile = fopen ("GPS.txt","r");
    pAhrsFile = fopen ("ATT.txt","r");
    pAdsFile = fopen ("NTUN.txt","r");
    pTimeFile = fopen ("timing.txt","r");
    pStateOutFile = fopen ("StateDataOut.txt","w");
    pEulOutFile = fopen ("EulDataOut.txt","w");
    pCovOutFile = fopen ("CovDataOut.txt","w");
    pRefPosVelOutFile = fopen ("RefVelPosDataOut.txt","w");
    pVelPosFuseFile = fopen ("VelPosFuse.txt","w");
    pMagFuseFile = fopen ("MagFuse.txt","w");
    pTasFuseFile = fopen ("TasFuse.txt","w");

    // read test data from files for first timestamp
    readIMUData();
    readGpsData();
    readMagData();
    readAirSpdData();
    readAhrsData();
    readTimingData();

    while (!endOfData)
    {
        if ((IMUmsec >= msecStartTime) && (IMUmsec <= msecEndTime))
        {
            // Initialise states, covariance and other data
            if ((IMUmsec > msecAlignTime) && !statesInitialised && (GPSstatus == 3))
            {
                InitialiseFilter();
            }

            // If valid IMU data and states initialised, predict states and covariances
            if (statesInitialised)
            {
                // Run the strapdown INS equations every IMU update
                UpdateStrapdownEquationsNED();
                // debug code - could be tunred into a filter mnitoring/watchdog function
                float tempQuat[4];
                for (uint8_t j=0; j<=3; j++) tempQuat[j] = states[j];
                quat2eul(eulerEst, tempQuat);
                for (uint8_t j=0; j<=2; j++) eulerDif[j] = eulerEst[j] - ahrsEul[j];
                if (eulerDif[2] > pi) eulerDif[2] -= 2*pi;
                if (eulerDif[2] < -pi) eulerDif[2] += 2*pi;
                // store the predicted states for subsequent use by measurement fusion
                StoreStates();
                // Check if on ground - status is used by covariance prediction
                OnGroundCheck();
                // sum delta angles and time used by covariance prediction
                summedDelAng = summedDelAng + correctedDelAng;
                summedDelVel = summedDelVel + dVelIMU;
                dt += dtIMU;
                // perform a covariance prediction if the total delta angle has exceeded the limit
                // or the time limit will be exceeded at the next IMU update
                if ((dt >= (covTimeStepMax - dtIMU)) || (summedDelAng.length() > covDelAngMax))
                {
                    CovariancePrediction();
                    summedDelAng = summedDelAng.zero();
                    summedDelVel = summedDelVel.zero();
                    dt = 0.0f;
                }
            }

            // Fuse GPS Measurements
            if (newDataGps && statesInitialised)
            {
                // Convert GPS measurements to Pos NE, hgt and Vel NED
                calcvelNED(velNED, gpsCourse, gpsGndSpd, gpsVelD);
                calcposNED(posNED, gpsLat, gpsLon, gpsHgt, latRef, lonRef, hgtRef);
                posNE[0] = posNED[0];
                posNE[1] = posNED[1];
                hgtMea =  -posNED[2];
                // set fusion flags
                fuseVelData = true;
                fusePosData = true;
                fuseHgtData = true;
                // recall states stored at time of measurement after adjusting for delays
                RecallStates(statesAtVelTime, (IMUmsec - msecVelDelay));
                RecallStates(statesAtPosTime, (IMUmsec - msecPosDelay));
                RecallStates(statesAtHgtTime, (IMUmsec - msecHgtDelay));
                // run the fusion step
                FuseVelposNED();
            }
            else
            {
                fuseVelData = false;
                fusePosData = false;
                fuseHgtData = false;
            }

            // Fuse Magnetometer Measurements
            if (newDataMag && statesInitialised)
            {
                fuseMagData = true;
                RecallStates(statesAtMagMeasTime, (IMUmsec - msecMagDelay)); // Assume 50 msec avg delay for magnetometer data
            }
            else
            {
                fuseMagData = false;
            }
            if (statesInitialised) FuseMagnetometer();

            // Fuse Airspeed Measurements
            if (newAdsData && statesInitialised && VtasMeas > 8.0f)
            {
                fuseVtasData = true;
                RecallStates(statesAtVtasMeasTime, (IMUmsec - msecTasDelay)); // assume 100 msec avg delay for airspeed data
                FuseAirspeed();
            }
            else
            {
                fuseVtasData = false;
            }

            // debug output
            //printf("Euler Angle Difference = %3.1f , %3.1f , %3.1f deg\n", rad2deg*eulerDif[0],rad2deg*eulerDif[1],rad2deg*eulerDif[2]);
            WriteFilterOutput();
        }
        // read test data from files for next timestamp
        readIMUData();
        readGpsData();
        readMagData();
        readAirSpdData();
        readAhrsData();
        if (IMUmsec > msecEndTime)
        {
            CloseFiles();
            endOfData = true;
        }
    }
}

void  UpdateStrapdownEquationsNED()
{
    static Vector3f  prevDelAng;
    Vector3f delVelNav;
    float q00;
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
    static float lastVelocity[3];
    const Vector3f gravityNED = {0.0,0.0,GRAVITY_MSS};

// Remove sensor bias errors
    correctedDelAng.x = dAngIMU.x - states[10];
    correctedDelAng.y = dAngIMU.y - states[11];
    correctedDelAng.z = dAngIMU.z - states[12];
    dVelIMU.x = dVelIMU.x;
    dVelIMU.y = dVelIMU.y;
    dVelIMU.z = dVelIMU.z;

// Save current measurements
    prevDelAng = correctedDelAng;

// Apply corrections for earths rotation rate and coning errors
// * and + operators have been overloaded
    correctedDelAng   = correctedDelAng - Tnb*earthRateNED*dtIMU + 8.333333333333333e-2f*(prevDelAng % correctedDelAng);

// Convert the rotation vector to its equivalent quaternion
    rotationMag = correctedDelAng.length();
    if (rotationMag < 1e-12f)
    {
        deltaQuat[0] = 1.0;
        deltaQuat[1] = 0.0;
        deltaQuat[2] = 0.0;
        deltaQuat[3] = 0.0;
    }
    else
    {
        deltaQuat[0] = cos(0.5f*rotationMag);
        rotScaler = (sin(0.5f*rotationMag))/rotationMag;
        deltaQuat[1] = correctedDelAng.x*rotScaler;
        deltaQuat[2] = correctedDelAng.y*rotScaler;
        deltaQuat[3] = correctedDelAng.z*rotScaler;
    }

// Update the quaternions by rotating from the previous attitude through
// the delta angle rotation quaternion
    qUpdated[0] = states[0]*deltaQuat[0] - states[1]*deltaQuat[1] - states[2]*deltaQuat[2] - states[3]*deltaQuat[3];
    qUpdated[1] = states[0]*deltaQuat[1] + states[1]*deltaQuat[0] + states[2]*deltaQuat[3] - states[3]*deltaQuat[2];
    qUpdated[2] = states[0]*deltaQuat[2] + states[2]*deltaQuat[0] + states[3]*deltaQuat[1] - states[1]*deltaQuat[3];
    qUpdated[3] = states[0]*deltaQuat[3] + states[3]*deltaQuat[0] + states[1]*deltaQuat[2] - states[2]*deltaQuat[1];

// Normalise the quaternions and update the quaternion states
    quatMag = sqrt(sq(qUpdated[0]) + sq(qUpdated[1]) + sq(qUpdated[2]) + sq(qUpdated[3]));
    if (quatMag > 1e-16f)
    {
        quatMagInv = 1.0f/quatMag;
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

    Tbn.x.x = q00 + q11 - q22 - q33;
    Tbn.y.y = q00 - q11 + q22 - q33;
    Tbn.z.z = q00 - q11 - q22 + q33;
    Tbn.x.y = 2*(q12 - q03);
    Tbn.x.z = 2*(q13 + q02);
    Tbn.y.x = 2*(q12 + q03);
    Tbn.y.z = 2*(q23 - q01);
    Tbn.z.x = 2*(q13 - q02);
    Tbn.z.y = 2*(q23 + q01);

    Tnb = Tbn.transpose();

// transform body delta velocities to delta velocities in the nav frame
// * and + operators have been overloaded
    //delVelNav = Tbn*dVelIMU + gravityNED*dtIMU;
    delVelNav.x = Tbn.x.x*dVelIMU.x + Tbn.x.y*dVelIMU.y + Tbn.x.z*dVelIMU.z + gravityNED.x*dtIMU;
    delVelNav.y = Tbn.y.x*dVelIMU.x + Tbn.y.y*dVelIMU.y + Tbn.y.z*dVelIMU.z + gravityNED.y*dtIMU;
    delVelNav.z = Tbn.z.x*dVelIMU.x + Tbn.z.y*dVelIMU.y + Tbn.z.z*dVelIMU.z + gravityNED.z*dtIMU;

// calculate the magnitude of the nav acceleration (required for GPS
// variance estimation)
    accNavMag = delVelNav.length()/dtIMU;

// If calculating position save previous velocity
    lastVelocity[0] = states[4];
    lastVelocity[1] = states[5];
    lastVelocity[2] = states[6];

// Sum delta velocities to get velocity
    states[4] = states[4] + delVelNav.x;
    states[5] = states[5] + delVelNav.y;
    states[6] = states[6] + delVelNav.z;

// If calculating postions, do a trapezoidal integration for position
    states[7] = states[7] + 0.5f*(states[4] + lastVelocity[0])*dtIMU;
    states[8] = states[8] + 0.5f*(states[5] + lastVelocity[1])*dtIMU;
    states[9] = states[9] + 0.5f*(states[6] + lastVelocity[2])*dtIMU;

}

void CovariancePrediction()
{
    // scalars
    float windVelSigma;
    float dAngBiasSigma;
    float dVelBiasSigma;
    float magEarthSigma;
    float magBodySigma;
    float daxCov;
    float dayCov;
    float dazCov;
    float dvxCov;
    float dvyCov;
    float dvzCov;
    float dvx;
    float dvy;
    float dvz;
    float dax;
    float day;
    float daz;
    float q0;
    float q1;
    float q2;
    float q3;
    float dax_b;
    float day_b;
    float daz_b;

    // arrays
    float processNoise[21];
    float SF[14];
    float SG[8];
    float SQ[11];
    float SPP[13];
    float nextP[21][21];

    // calculate covariance prediction process noise
    const float yawVarScale = 10.0f;
    windVelSigma = dt*0.1f;
    dAngBiasSigma = dt*5.0e-7f;
    dVelBiasSigma = dt*1.0e-4f;
    magEarthSigma = dt*1.5e-4f;
    magBodySigma  = dt*1.5e-4f;
    for (uint8_t i= 0; i<=9; i++) processNoise[i]  = 1.0e-9f;
    for (uint8_t i=10; i<=12; i++) processNoise[i] = dAngBiasSigma;
    if (onGround) processNoise[12] = dAngBiasSigma * yawVarScale;
    for (uint8_t i=13; i<=14; i++) processNoise[i] = windVelSigma;
    for (uint8_t i=15; i<=17; i++) processNoise[i] = magEarthSigma;
    for (uint8_t i=18; i<=20; i++) processNoise[i] = magBodySigma;
    for (uint8_t i= 0; i<=20; i++) processNoise[i] = sq(processNoise[i]);

    // set variables used to calculate covariance growth
    dvx = summedDelVel.x;
    dvy = summedDelVel.y;
    dvz = summedDelVel.z;
    dax = summedDelAng.x;
    day = summedDelAng.y;
    daz = summedDelAng.z;
    q0 = states[0];
    q1 = states[1];
    q2 = states[2];
    q3 = states[3];
    dax_b = states[10];
    day_b = states[11];
    daz_b = states[12];
    daxCov = sq(dt*1.4544411e-2f);
    dayCov = sq(dt*1.4544411e-2f);
    dazCov = sq(dt*1.4544411e-2f);
    if (onGround) dazCov = dazCov * sq(yawVarScale);
    dvxCov = sq(dt*0.5f);
    dvyCov = sq(dt*0.5f);
    dvzCov = sq(dt*0.5f);

    // Predicted covariance calculation
    SF[0] = 2*dvx*q1 + 2*dvy*q2 + 2*dvz*q3;
    SF[1] = 2*dvx*q3 + 2*dvy*q0 - 2*dvz*q1;
    SF[2] = 2*dvx*q0 - 2*dvy*q3 + 2*dvz*q2;
    SF[3] = day/2 - day_b/2;
    SF[4] = daz/2 - daz_b/2;
    SF[5] = dax/2 - dax_b/2;
    SF[6] = dax_b/2 - dax/2;
    SF[7] = daz_b/2 - daz/2;
    SF[8] = day_b/2 - day/2;
    SF[9] = q1/2;
    SF[10] = q2/2;
    SF[11] = q3/2;
    SF[12] = 2*dvz*q0;
    SF[13] = 2*dvy*q1;

    SG[0] = q0/2;
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
    SQ[3] = (dayCov*q1*SG[0])/2 - (dazCov*q1*SG[0])/2 - (daxCov*q2*q3)/4;
    SQ[4] = (dazCov*q2*SG[0])/2 - (daxCov*q2*SG[0])/2 - (dayCov*q1*q3)/4;
    SQ[5] = (daxCov*q3*SG[0])/2 - (dayCov*q3*SG[0])/2 - (dazCov*q1*q2)/4;
    SQ[6] = (daxCov*q1*q2)/4 - (dazCov*q3*SG[0])/2 - (dayCov*q1*q2)/4;
    SQ[7] = (dazCov*q1*q3)/4 - (daxCov*q1*q3)/4 - (dayCov*q2*SG[0])/2;
    SQ[8] = (dayCov*q2*q3)/4 - (daxCov*q1*SG[0])/2 - (dazCov*q2*q3)/4;
    SQ[9] = sq(SG[0]);
    SQ[10] = sq(q1);

    SPP[0] = SF[12] + SF[13] - 2*dvx*q2;
    SPP[1] = 2*dvx*q0 - 2*dvy*q3 + 2*dvz*q2;
    SPP[2] = 2*dvx*q3 + 2*dvy*q0 - 2*dvz*q1;
    SPP[3] = SF[11];
    SPP[4] = SF[10];
    SPP[5] = SF[9];
    SPP[6] = SF[7];
    SPP[7] = SF[8];

    nextP[0][0] = P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3] + (daxCov*SQ[10])/4 + SF[6]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]) + SPP[7]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]) + SPP[6]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]) + SPP[5]*(P[0][10] + P[1][10]*SF[6] + P[2][10]*SPP[7] + P[3][10]*SPP[6] + P[10][10]*SPP[5] + P[11][10]*SPP[4] + P[12][10]*SPP[3]) + SPP[4]*(P[0][11] + P[1][11]*SF[6] + P[2][11]*SPP[7] + P[3][11]*SPP[6] + P[10][11]*SPP[5] + P[11][11]*SPP[4] + P[12][11]*SPP[3]) + SPP[3]*(P[0][12] + P[1][12]*SF[6] + P[2][12]*SPP[7] + P[3][12]*SPP[6] + P[10][12]*SPP[5] + P[11][12]*SPP[4] + P[12][12]*SPP[3]) + (dayCov*sq(q2))/4 + (dazCov*sq(q3))/4;
    nextP[0][1] = P[0][1] + SQ[8] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3] + SF[5]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) + SF[4]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]) + SPP[7]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]) + SPP[3]*(P[0][11] + P[1][11]*SF[6] + P[2][11]*SPP[7] + P[3][11]*SPP[6] + P[10][11]*SPP[5] + P[11][11]*SPP[4] + P[12][11]*SPP[3]) - SPP[4]*(P[0][12] + P[1][12]*SF[6] + P[2][12]*SPP[7] + P[3][12]*SPP[6] + P[10][12]*SPP[5] + P[11][12]*SPP[4] + P[12][12]*SPP[3]) - (q0*(P[0][10] + P[1][10]*SF[6] + P[2][10]*SPP[7] + P[3][10]*SPP[6] + P[10][10]*SPP[5] + P[11][10]*SPP[4] + P[12][10]*SPP[3]))/2;
    nextP[0][2] = P[0][2] + SQ[7] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3] + SF[3]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) + SF[5]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]) + SPP[6]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]) - SPP[3]*(P[0][10] + P[1][10]*SF[6] + P[2][10]*SPP[7] + P[3][10]*SPP[6] + P[10][10]*SPP[5] + P[11][10]*SPP[4] + P[12][10]*SPP[3]) + SPP[5]*(P[0][12] + P[1][12]*SF[6] + P[2][12]*SPP[7] + P[3][12]*SPP[6] + P[10][12]*SPP[5] + P[11][12]*SPP[4] + P[12][12]*SPP[3]) - (q0*(P[0][11] + P[1][11]*SF[6] + P[2][11]*SPP[7] + P[3][11]*SPP[6] + P[10][11]*SPP[5] + P[11][11]*SPP[4] + P[12][11]*SPP[3]))/2;
    nextP[0][3] = P[0][3] + SQ[6] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3] + SF[4]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) + SF[3]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]) + SF[6]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]) + SPP[4]*(P[0][10] + P[1][10]*SF[6] + P[2][10]*SPP[7] + P[3][10]*SPP[6] + P[10][10]*SPP[5] + P[11][10]*SPP[4] + P[12][10]*SPP[3]) - SPP[5]*(P[0][11] + P[1][11]*SF[6] + P[2][11]*SPP[7] + P[3][11]*SPP[6] + P[10][11]*SPP[5] + P[11][11]*SPP[4] + P[12][11]*SPP[3]) - (q0*(P[0][12] + P[1][12]*SF[6] + P[2][12]*SPP[7] + P[3][12]*SPP[6] + P[10][12]*SPP[5] + P[11][12]*SPP[4] + P[12][12]*SPP[3]))/2;
    nextP[0][4] = P[0][4] + P[1][4]*SF[6] + P[2][4]*SPP[7] + P[3][4]*SPP[6] + P[10][4]*SPP[5] + P[11][4]*SPP[4] + P[12][4]*SPP[3] + SF[2]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) + SF[0]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]) + SPP[0]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]) - SPP[2]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]);
    nextP[0][5] = P[0][5] + P[1][5]*SF[6] + P[2][5]*SPP[7] + P[3][5]*SPP[6] + P[10][5]*SPP[5] + P[11][5]*SPP[4] + P[12][5]*SPP[3] + SF[1]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) + SF[0]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]) + SF[2]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]) - SPP[0]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]);
    nextP[0][6] = P[0][6] + P[1][6]*SF[6] + P[2][6]*SPP[7] + P[3][6]*SPP[6] + P[10][6]*SPP[5] + P[11][6]*SPP[4] + P[12][6]*SPP[3] + SF[1]*(P[0][1] + P[1][1]*SF[6] + P[2][1]*SPP[7] + P[3][1]*SPP[6] + P[10][1]*SPP[5] + P[11][1]*SPP[4] + P[12][1]*SPP[3]) + SF[0]*(P[0][3] + P[1][3]*SF[6] + P[2][3]*SPP[7] + P[3][3]*SPP[6] + P[10][3]*SPP[5] + P[11][3]*SPP[4] + P[12][3]*SPP[3]) + SPP[0]*(P[0][0] + P[1][0]*SF[6] + P[2][0]*SPP[7] + P[3][0]*SPP[6] + P[10][0]*SPP[5] + P[11][0]*SPP[4] + P[12][0]*SPP[3]) - SPP[1]*(P[0][2] + P[1][2]*SF[6] + P[2][2]*SPP[7] + P[3][2]*SPP[6] + P[10][2]*SPP[5] + P[11][2]*SPP[4] + P[12][2]*SPP[3]);
    nextP[0][7] = P[0][7] + P[1][7]*SF[6] + P[2][7]*SPP[7] + P[3][7]*SPP[6] + P[10][7]*SPP[5] + P[11][7]*SPP[4] + P[12][7]*SPP[3] + dt*(P[0][4] + P[1][4]*SF[6] + P[2][4]*SPP[7] + P[3][4]*SPP[6] + P[10][4]*SPP[5] + P[11][4]*SPP[4] + P[12][4]*SPP[3]);
    nextP[0][8] = P[0][8] + P[1][8]*SF[6] + P[2][8]*SPP[7] + P[3][8]*SPP[6] + P[10][8]*SPP[5] + P[11][8]*SPP[4] + P[12][8]*SPP[3] + dt*(P[0][5] + P[1][5]*SF[6] + P[2][5]*SPP[7] + P[3][5]*SPP[6] + P[10][5]*SPP[5] + P[11][5]*SPP[4] + P[12][5]*SPP[3]);
    nextP[0][9] = P[0][9] + P[1][9]*SF[6] + P[2][9]*SPP[7] + P[3][9]*SPP[6] + P[10][9]*SPP[5] + P[11][9]*SPP[4] + P[12][9]*SPP[3] + dt*(P[0][6] + P[1][6]*SF[6] + P[2][6]*SPP[7] + P[3][6]*SPP[6] + P[10][6]*SPP[5] + P[11][6]*SPP[4] + P[12][6]*SPP[3]);
    nextP[0][10] = P[0][10] + P[1][10]*SF[6] + P[2][10]*SPP[7] + P[3][10]*SPP[6] + P[10][10]*SPP[5] + P[11][10]*SPP[4] + P[12][10]*SPP[3];
    nextP[0][11] = P[0][11] + P[1][11]*SF[6] + P[2][11]*SPP[7] + P[3][11]*SPP[6] + P[10][11]*SPP[5] + P[11][11]*SPP[4] + P[12][11]*SPP[3];
    nextP[0][12] = P[0][12] + P[1][12]*SF[6] + P[2][12]*SPP[7] + P[3][12]*SPP[6] + P[10][12]*SPP[5] + P[11][12]*SPP[4] + P[12][12]*SPP[3];
    nextP[0][13] = P[0][13] + P[1][13]*SF[6] + P[2][13]*SPP[7] + P[3][13]*SPP[6] + P[10][13]*SPP[5] + P[11][13]*SPP[4] + P[12][13]*SPP[3];
    nextP[0][14] = P[0][14] + P[1][14]*SF[6] + P[2][14]*SPP[7] + P[3][14]*SPP[6] + P[10][14]*SPP[5] + P[11][14]*SPP[4] + P[12][14]*SPP[3];
    nextP[0][15] = P[0][15] + P[1][15]*SF[6] + P[2][15]*SPP[7] + P[3][15]*SPP[6] + P[10][15]*SPP[5] + P[11][15]*SPP[4] + P[12][15]*SPP[3];
    nextP[0][16] = P[0][16] + P[1][16]*SF[6] + P[2][16]*SPP[7] + P[3][16]*SPP[6] + P[10][16]*SPP[5] + P[11][16]*SPP[4] + P[12][16]*SPP[3];
    nextP[0][17] = P[0][17] + P[1][17]*SF[6] + P[2][17]*SPP[7] + P[3][17]*SPP[6] + P[10][17]*SPP[5] + P[11][17]*SPP[4] + P[12][17]*SPP[3];
    nextP[0][18] = P[0][18] + P[1][18]*SF[6] + P[2][18]*SPP[7] + P[3][18]*SPP[6] + P[10][18]*SPP[5] + P[11][18]*SPP[4] + P[12][18]*SPP[3];
    nextP[0][19] = P[0][19] + P[1][19]*SF[6] + P[2][19]*SPP[7] + P[3][19]*SPP[6] + P[10][19]*SPP[5] + P[11][19]*SPP[4] + P[12][19]*SPP[3];
    nextP[0][20] = P[0][20] + P[1][20]*SF[6] + P[2][20]*SPP[7] + P[3][20]*SPP[6] + P[10][20]*SPP[5] + P[11][20]*SPP[4] + P[12][20]*SPP[3];
    nextP[1][0] = P[1][0] + SQ[8] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2 + SF[6]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2) + SPP[7]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2) + SPP[6]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2) + SPP[5]*(P[1][10] + P[0][10]*SF[5] + P[2][10]*SF[4] + P[3][10]*SPP[7] + P[11][10]*SPP[3] - P[12][10]*SPP[4] - (P[10][10]*q0)/2) + SPP[4]*(P[1][11] + P[0][11]*SF[5] + P[2][11]*SF[4] + P[3][11]*SPP[7] + P[11][11]*SPP[3] - P[12][11]*SPP[4] - (P[10][11]*q0)/2) + SPP[3]*(P[1][12] + P[0][12]*SF[5] + P[2][12]*SF[4] + P[3][12]*SPP[7] + P[11][12]*SPP[3] - P[12][12]*SPP[4] - (P[10][12]*q0)/2);
    nextP[1][1] = P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] + daxCov*SQ[9] - (P[10][1]*q0)/2 + SF[5]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) + SF[4]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2) + SPP[7]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2) + SPP[3]*(P[1][11] + P[0][11]*SF[5] + P[2][11]*SF[4] + P[3][11]*SPP[7] + P[11][11]*SPP[3] - P[12][11]*SPP[4] - (P[10][11]*q0)/2) - SPP[4]*(P[1][12] + P[0][12]*SF[5] + P[2][12]*SF[4] + P[3][12]*SPP[7] + P[11][12]*SPP[3] - P[12][12]*SPP[4] - (P[10][12]*q0)/2) + (dayCov*sq(q3))/4 + (dazCov*sq(q2))/4 - (q0*(P[1][10] + P[0][10]*SF[5] + P[2][10]*SF[4] + P[3][10]*SPP[7] + P[11][10]*SPP[3] - P[12][10]*SPP[4] - (P[10][10]*q0)/2))/2;
    nextP[1][2] = P[1][2] + SQ[5] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2 + SF[3]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) + SF[5]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2) + SPP[6]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2) - SPP[3]*(P[1][10] + P[0][10]*SF[5] + P[2][10]*SF[4] + P[3][10]*SPP[7] + P[11][10]*SPP[3] - P[12][10]*SPP[4] - (P[10][10]*q0)/2) + SPP[5]*(P[1][12] + P[0][12]*SF[5] + P[2][12]*SF[4] + P[3][12]*SPP[7] + P[11][12]*SPP[3] - P[12][12]*SPP[4] - (P[10][12]*q0)/2) - (q0*(P[1][11] + P[0][11]*SF[5] + P[2][11]*SF[4] + P[3][11]*SPP[7] + P[11][11]*SPP[3] - P[12][11]*SPP[4] - (P[10][11]*q0)/2))/2;
    nextP[1][3] = P[1][3] + SQ[4] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2 + SF[4]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) + SF[3]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2) + SF[6]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2) + SPP[4]*(P[1][10] + P[0][10]*SF[5] + P[2][10]*SF[4] + P[3][10]*SPP[7] + P[11][10]*SPP[3] - P[12][10]*SPP[4] - (P[10][10]*q0)/2) - SPP[5]*(P[1][11] + P[0][11]*SF[5] + P[2][11]*SF[4] + P[3][11]*SPP[7] + P[11][11]*SPP[3] - P[12][11]*SPP[4] - (P[10][11]*q0)/2) - (q0*(P[1][12] + P[0][12]*SF[5] + P[2][12]*SF[4] + P[3][12]*SPP[7] + P[11][12]*SPP[3] - P[12][12]*SPP[4] - (P[10][12]*q0)/2))/2;
    nextP[1][4] = P[1][4] + P[0][4]*SF[5] + P[2][4]*SF[4] + P[3][4]*SPP[7] + P[11][4]*SPP[3] - P[12][4]*SPP[4] - (P[10][4]*q0)/2 + SF[2]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) + SF[0]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2) + SPP[0]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2) - SPP[2]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2);
    nextP[1][5] = P[1][5] + P[0][5]*SF[5] + P[2][5]*SF[4] + P[3][5]*SPP[7] + P[11][5]*SPP[3] - P[12][5]*SPP[4] - (P[10][5]*q0)/2 + SF[1]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) + SF[0]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2) + SF[2]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2) - SPP[0]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2);
    nextP[1][6] = P[1][6] + P[0][6]*SF[5] + P[2][6]*SF[4] + P[3][6]*SPP[7] + P[11][6]*SPP[3] - P[12][6]*SPP[4] - (P[10][6]*q0)/2 + SF[1]*(P[1][1] + P[0][1]*SF[5] + P[2][1]*SF[4] + P[3][1]*SPP[7] + P[11][1]*SPP[3] - P[12][1]*SPP[4] - (P[10][1]*q0)/2) + SF[0]*(P[1][3] + P[0][3]*SF[5] + P[2][3]*SF[4] + P[3][3]*SPP[7] + P[11][3]*SPP[3] - P[12][3]*SPP[4] - (P[10][3]*q0)/2) + SPP[0]*(P[1][0] + P[0][0]*SF[5] + P[2][0]*SF[4] + P[3][0]*SPP[7] + P[11][0]*SPP[3] - P[12][0]*SPP[4] - (P[10][0]*q0)/2) - SPP[1]*(P[1][2] + P[0][2]*SF[5] + P[2][2]*SF[4] + P[3][2]*SPP[7] + P[11][2]*SPP[3] - P[12][2]*SPP[4] - (P[10][2]*q0)/2);
    nextP[1][7] = P[1][7] + P[0][7]*SF[5] + P[2][7]*SF[4] + P[3][7]*SPP[7] + P[11][7]*SPP[3] - P[12][7]*SPP[4] - (P[10][7]*q0)/2 + dt*(P[1][4] + P[0][4]*SF[5] + P[2][4]*SF[4] + P[3][4]*SPP[7] + P[11][4]*SPP[3] - P[12][4]*SPP[4] - (P[10][4]*q0)/2);
    nextP[1][8] = P[1][8] + P[0][8]*SF[5] + P[2][8]*SF[4] + P[3][8]*SPP[7] + P[11][8]*SPP[3] - P[12][8]*SPP[4] - (P[10][8]*q0)/2 + dt*(P[1][5] + P[0][5]*SF[5] + P[2][5]*SF[4] + P[3][5]*SPP[7] + P[11][5]*SPP[3] - P[12][5]*SPP[4] - (P[10][5]*q0)/2);
    nextP[1][9] = P[1][9] + P[0][9]*SF[5] + P[2][9]*SF[4] + P[3][9]*SPP[7] + P[11][9]*SPP[3] - P[12][9]*SPP[4] - (P[10][9]*q0)/2 + dt*(P[1][6] + P[0][6]*SF[5] + P[2][6]*SF[4] + P[3][6]*SPP[7] + P[11][6]*SPP[3] - P[12][6]*SPP[4] - (P[10][6]*q0)/2);
    nextP[1][10] = P[1][10] + P[0][10]*SF[5] + P[2][10]*SF[4] + P[3][10]*SPP[7] + P[11][10]*SPP[3] - P[12][10]*SPP[4] - (P[10][10]*q0)/2;
    nextP[1][11] = P[1][11] + P[0][11]*SF[5] + P[2][11]*SF[4] + P[3][11]*SPP[7] + P[11][11]*SPP[3] - P[12][11]*SPP[4] - (P[10][11]*q0)/2;
    nextP[1][12] = P[1][12] + P[0][12]*SF[5] + P[2][12]*SF[4] + P[3][12]*SPP[7] + P[11][12]*SPP[3] - P[12][12]*SPP[4] - (P[10][12]*q0)/2;
    nextP[1][13] = P[1][13] + P[0][13]*SF[5] + P[2][13]*SF[4] + P[3][13]*SPP[7] + P[11][13]*SPP[3] - P[12][13]*SPP[4] - (P[10][13]*q0)/2;
    nextP[1][14] = P[1][14] + P[0][14]*SF[5] + P[2][14]*SF[4] + P[3][14]*SPP[7] + P[11][14]*SPP[3] - P[12][14]*SPP[4] - (P[10][14]*q0)/2;
    nextP[1][15] = P[1][15] + P[0][15]*SF[5] + P[2][15]*SF[4] + P[3][15]*SPP[7] + P[11][15]*SPP[3] - P[12][15]*SPP[4] - (P[10][15]*q0)/2;
    nextP[1][16] = P[1][16] + P[0][16]*SF[5] + P[2][16]*SF[4] + P[3][16]*SPP[7] + P[11][16]*SPP[3] - P[12][16]*SPP[4] - (P[10][16]*q0)/2;
    nextP[1][17] = P[1][17] + P[0][17]*SF[5] + P[2][17]*SF[4] + P[3][17]*SPP[7] + P[11][17]*SPP[3] - P[12][17]*SPP[4] - (P[10][17]*q0)/2;
    nextP[1][18] = P[1][18] + P[0][18]*SF[5] + P[2][18]*SF[4] + P[3][18]*SPP[7] + P[11][18]*SPP[3] - P[12][18]*SPP[4] - (P[10][18]*q0)/2;
    nextP[1][19] = P[1][19] + P[0][19]*SF[5] + P[2][19]*SF[4] + P[3][19]*SPP[7] + P[11][19]*SPP[3] - P[12][19]*SPP[4] - (P[10][19]*q0)/2;
    nextP[1][20] = P[1][20] + P[0][20]*SF[5] + P[2][20]*SF[4] + P[3][20]*SPP[7] + P[11][20]*SPP[3] - P[12][20]*SPP[4] - (P[10][20]*q0)/2;
    nextP[2][0] = P[2][0] + SQ[7] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2 + SF[6]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2) + SPP[7]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2) + SPP[6]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2) + SPP[5]*(P[2][10] + P[0][10]*SF[3] + P[3][10]*SF[5] + P[1][10]*SPP[6] - P[10][10]*SPP[3] + P[12][10]*SPP[5] - (P[11][10]*q0)/2) + SPP[4]*(P[2][11] + P[0][11]*SF[3] + P[3][11]*SF[5] + P[1][11]*SPP[6] - P[10][11]*SPP[3] + P[12][11]*SPP[5] - (P[11][11]*q0)/2) + SPP[3]*(P[2][12] + P[0][12]*SF[3] + P[3][12]*SF[5] + P[1][12]*SPP[6] - P[10][12]*SPP[3] + P[12][12]*SPP[5] - (P[11][12]*q0)/2);
    nextP[2][1] = P[2][1] + SQ[5] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2 + SF[5]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) + SF[4]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2) + SPP[7]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2) + SPP[3]*(P[2][11] + P[0][11]*SF[3] + P[3][11]*SF[5] + P[1][11]*SPP[6] - P[10][11]*SPP[3] + P[12][11]*SPP[5] - (P[11][11]*q0)/2) - SPP[4]*(P[2][12] + P[0][12]*SF[3] + P[3][12]*SF[5] + P[1][12]*SPP[6] - P[10][12]*SPP[3] + P[12][12]*SPP[5] - (P[11][12]*q0)/2) - (q0*(P[2][10] + P[0][10]*SF[3] + P[3][10]*SF[5] + P[1][10]*SPP[6] - P[10][10]*SPP[3] + P[12][10]*SPP[5] - (P[11][10]*q0)/2))/2;
    nextP[2][2] = P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] + dayCov*SQ[9] + (dazCov*SQ[10])/4 - (P[11][2]*q0)/2 + SF[3]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) + SF[5]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2) + SPP[6]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2) - SPP[3]*(P[2][10] + P[0][10]*SF[3] + P[3][10]*SF[5] + P[1][10]*SPP[6] - P[10][10]*SPP[3] + P[12][10]*SPP[5] - (P[11][10]*q0)/2) + SPP[5]*(P[2][12] + P[0][12]*SF[3] + P[3][12]*SF[5] + P[1][12]*SPP[6] - P[10][12]*SPP[3] + P[12][12]*SPP[5] - (P[11][12]*q0)/2) + (daxCov*sq(q3))/4 - (q0*(P[2][11] + P[0][11]*SF[3] + P[3][11]*SF[5] + P[1][11]*SPP[6] - P[10][11]*SPP[3] + P[12][11]*SPP[5] - (P[11][11]*q0)/2))/2;
    nextP[2][3] = P[2][3] + SQ[3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2 + SF[4]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) + SF[3]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2) + SF[6]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2) + SPP[4]*(P[2][10] + P[0][10]*SF[3] + P[3][10]*SF[5] + P[1][10]*SPP[6] - P[10][10]*SPP[3] + P[12][10]*SPP[5] - (P[11][10]*q0)/2) - SPP[5]*(P[2][11] + P[0][11]*SF[3] + P[3][11]*SF[5] + P[1][11]*SPP[6] - P[10][11]*SPP[3] + P[12][11]*SPP[5] - (P[11][11]*q0)/2) - (q0*(P[2][12] + P[0][12]*SF[3] + P[3][12]*SF[5] + P[1][12]*SPP[6] - P[10][12]*SPP[3] + P[12][12]*SPP[5] - (P[11][12]*q0)/2))/2;
    nextP[2][4] = P[2][4] + P[0][4]*SF[3] + P[3][4]*SF[5] + P[1][4]*SPP[6] - P[10][4]*SPP[3] + P[12][4]*SPP[5] - (P[11][4]*q0)/2 + SF[2]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) + SF[0]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2) + SPP[0]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2) - SPP[2]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2);
    nextP[2][5] = P[2][5] + P[0][5]*SF[3] + P[3][5]*SF[5] + P[1][5]*SPP[6] - P[10][5]*SPP[3] + P[12][5]*SPP[5] - (P[11][5]*q0)/2 + SF[1]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) + SF[0]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2) + SF[2]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2) - SPP[0]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2);
    nextP[2][6] = P[2][6] + P[0][6]*SF[3] + P[3][6]*SF[5] + P[1][6]*SPP[6] - P[10][6]*SPP[3] + P[12][6]*SPP[5] - (P[11][6]*q0)/2 + SF[1]*(P[2][1] + P[0][1]*SF[3] + P[3][1]*SF[5] + P[1][1]*SPP[6] - P[10][1]*SPP[3] + P[12][1]*SPP[5] - (P[11][1]*q0)/2) + SF[0]*(P[2][3] + P[0][3]*SF[3] + P[3][3]*SF[5] + P[1][3]*SPP[6] - P[10][3]*SPP[3] + P[12][3]*SPP[5] - (P[11][3]*q0)/2) + SPP[0]*(P[2][0] + P[0][0]*SF[3] + P[3][0]*SF[5] + P[1][0]*SPP[6] - P[10][0]*SPP[3] + P[12][0]*SPP[5] - (P[11][0]*q0)/2) - SPP[1]*(P[2][2] + P[0][2]*SF[3] + P[3][2]*SF[5] + P[1][2]*SPP[6] - P[10][2]*SPP[3] + P[12][2]*SPP[5] - (P[11][2]*q0)/2);
    nextP[2][7] = P[2][7] + P[0][7]*SF[3] + P[3][7]*SF[5] + P[1][7]*SPP[6] - P[10][7]*SPP[3] + P[12][7]*SPP[5] - (P[11][7]*q0)/2 + dt*(P[2][4] + P[0][4]*SF[3] + P[3][4]*SF[5] + P[1][4]*SPP[6] - P[10][4]*SPP[3] + P[12][4]*SPP[5] - (P[11][4]*q0)/2);
    nextP[2][8] = P[2][8] + P[0][8]*SF[3] + P[3][8]*SF[5] + P[1][8]*SPP[6] - P[10][8]*SPP[3] + P[12][8]*SPP[5] - (P[11][8]*q0)/2 + dt*(P[2][5] + P[0][5]*SF[3] + P[3][5]*SF[5] + P[1][5]*SPP[6] - P[10][5]*SPP[3] + P[12][5]*SPP[5] - (P[11][5]*q0)/2);
    nextP[2][9] = P[2][9] + P[0][9]*SF[3] + P[3][9]*SF[5] + P[1][9]*SPP[6] - P[10][9]*SPP[3] + P[12][9]*SPP[5] - (P[11][9]*q0)/2 + dt*(P[2][6] + P[0][6]*SF[3] + P[3][6]*SF[5] + P[1][6]*SPP[6] - P[10][6]*SPP[3] + P[12][6]*SPP[5] - (P[11][6]*q0)/2);
    nextP[2][10] = P[2][10] + P[0][10]*SF[3] + P[3][10]*SF[5] + P[1][10]*SPP[6] - P[10][10]*SPP[3] + P[12][10]*SPP[5] - (P[11][10]*q0)/2;
    nextP[2][11] = P[2][11] + P[0][11]*SF[3] + P[3][11]*SF[5] + P[1][11]*SPP[6] - P[10][11]*SPP[3] + P[12][11]*SPP[5] - (P[11][11]*q0)/2;
    nextP[2][12] = P[2][12] + P[0][12]*SF[3] + P[3][12]*SF[5] + P[1][12]*SPP[6] - P[10][12]*SPP[3] + P[12][12]*SPP[5] - (P[11][12]*q0)/2;
    nextP[2][13] = P[2][13] + P[0][13]*SF[3] + P[3][13]*SF[5] + P[1][13]*SPP[6] - P[10][13]*SPP[3] + P[12][13]*SPP[5] - (P[11][13]*q0)/2;
    nextP[2][14] = P[2][14] + P[0][14]*SF[3] + P[3][14]*SF[5] + P[1][14]*SPP[6] - P[10][14]*SPP[3] + P[12][14]*SPP[5] - (P[11][14]*q0)/2;
    nextP[2][15] = P[2][15] + P[0][15]*SF[3] + P[3][15]*SF[5] + P[1][15]*SPP[6] - P[10][15]*SPP[3] + P[12][15]*SPP[5] - (P[11][15]*q0)/2;
    nextP[2][16] = P[2][16] + P[0][16]*SF[3] + P[3][16]*SF[5] + P[1][16]*SPP[6] - P[10][16]*SPP[3] + P[12][16]*SPP[5] - (P[11][16]*q0)/2;
    nextP[2][17] = P[2][17] + P[0][17]*SF[3] + P[3][17]*SF[5] + P[1][17]*SPP[6] - P[10][17]*SPP[3] + P[12][17]*SPP[5] - (P[11][17]*q0)/2;
    nextP[2][18] = P[2][18] + P[0][18]*SF[3] + P[3][18]*SF[5] + P[1][18]*SPP[6] - P[10][18]*SPP[3] + P[12][18]*SPP[5] - (P[11][18]*q0)/2;
    nextP[2][19] = P[2][19] + P[0][19]*SF[3] + P[3][19]*SF[5] + P[1][19]*SPP[6] - P[10][19]*SPP[3] + P[12][19]*SPP[5] - (P[11][19]*q0)/2;
    nextP[2][20] = P[2][20] + P[0][20]*SF[3] + P[3][20]*SF[5] + P[1][20]*SPP[6] - P[10][20]*SPP[3] + P[12][20]*SPP[5] - (P[11][20]*q0)/2;
    nextP[3][0] = P[3][0] + SQ[6] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2 + SF[6]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2) + SPP[7]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2) + SPP[6]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2) + SPP[5]*(P[3][10] + P[0][10]*SF[4] + P[1][10]*SF[3] + P[2][10]*SF[6] + P[10][10]*SPP[4] - P[11][10]*SPP[5] - (P[12][10]*q0)/2) + SPP[4]*(P[3][11] + P[0][11]*SF[4] + P[1][11]*SF[3] + P[2][11]*SF[6] + P[10][11]*SPP[4] - P[11][11]*SPP[5] - (P[12][11]*q0)/2) + SPP[3]*(P[3][12] + P[0][12]*SF[4] + P[1][12]*SF[3] + P[2][12]*SF[6] + P[10][12]*SPP[4] - P[11][12]*SPP[5] - (P[12][12]*q0)/2);
    nextP[3][1] = P[3][1] + SQ[4] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2 + SF[5]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) + SF[4]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2) + SPP[7]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2) + SPP[3]*(P[3][11] + P[0][11]*SF[4] + P[1][11]*SF[3] + P[2][11]*SF[6] + P[10][11]*SPP[4] - P[11][11]*SPP[5] - (P[12][11]*q0)/2) - SPP[4]*(P[3][12] + P[0][12]*SF[4] + P[1][12]*SF[3] + P[2][12]*SF[6] + P[10][12]*SPP[4] - P[11][12]*SPP[5] - (P[12][12]*q0)/2) - (q0*(P[3][10] + P[0][10]*SF[4] + P[1][10]*SF[3] + P[2][10]*SF[6] + P[10][10]*SPP[4] - P[11][10]*SPP[5] - (P[12][10]*q0)/2))/2;
    nextP[3][2] = P[3][2] + SQ[3] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2 + SF[3]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) + SF[5]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2) + SPP[6]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2) - SPP[3]*(P[3][10] + P[0][10]*SF[4] + P[1][10]*SF[3] + P[2][10]*SF[6] + P[10][10]*SPP[4] - P[11][10]*SPP[5] - (P[12][10]*q0)/2) + SPP[5]*(P[3][12] + P[0][12]*SF[4] + P[1][12]*SF[3] + P[2][12]*SF[6] + P[10][12]*SPP[4] - P[11][12]*SPP[5] - (P[12][12]*q0)/2) - (q0*(P[3][11] + P[0][11]*SF[4] + P[1][11]*SF[3] + P[2][11]*SF[6] + P[10][11]*SPP[4] - P[11][11]*SPP[5] - (P[12][11]*q0)/2))/2;
    nextP[3][3] = P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] + (dayCov*SQ[10])/4 + dazCov*SQ[9] - (P[12][3]*q0)/2 + SF[4]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) + SF[3]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2) + SF[6]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2) + SPP[4]*(P[3][10] + P[0][10]*SF[4] + P[1][10]*SF[3] + P[2][10]*SF[6] + P[10][10]*SPP[4] - P[11][10]*SPP[5] - (P[12][10]*q0)/2) - SPP[5]*(P[3][11] + P[0][11]*SF[4] + P[1][11]*SF[3] + P[2][11]*SF[6] + P[10][11]*SPP[4] - P[11][11]*SPP[5] - (P[12][11]*q0)/2) + (daxCov*sq(q2))/4 - (q0*(P[3][12] + P[0][12]*SF[4] + P[1][12]*SF[3] + P[2][12]*SF[6] + P[10][12]*SPP[4] - P[11][12]*SPP[5] - (P[12][12]*q0)/2))/2;
    nextP[3][4] = P[3][4] + P[0][4]*SF[4] + P[1][4]*SF[3] + P[2][4]*SF[6] + P[10][4]*SPP[4] - P[11][4]*SPP[5] - (P[12][4]*q0)/2 + SF[2]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) + SF[0]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2) + SPP[0]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2) - SPP[2]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2);
    nextP[3][5] = P[3][5] + P[0][5]*SF[4] + P[1][5]*SF[3] + P[2][5]*SF[6] + P[10][5]*SPP[4] - P[11][5]*SPP[5] - (P[12][5]*q0)/2 + SF[1]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) + SF[0]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2) + SF[2]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2) - SPP[0]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2);
    nextP[3][6] = P[3][6] + P[0][6]*SF[4] + P[1][6]*SF[3] + P[2][6]*SF[6] + P[10][6]*SPP[4] - P[11][6]*SPP[5] - (P[12][6]*q0)/2 + SF[1]*(P[3][1] + P[0][1]*SF[4] + P[1][1]*SF[3] + P[2][1]*SF[6] + P[10][1]*SPP[4] - P[11][1]*SPP[5] - (P[12][1]*q0)/2) + SF[0]*(P[3][3] + P[0][3]*SF[4] + P[1][3]*SF[3] + P[2][3]*SF[6] + P[10][3]*SPP[4] - P[11][3]*SPP[5] - (P[12][3]*q0)/2) + SPP[0]*(P[3][0] + P[0][0]*SF[4] + P[1][0]*SF[3] + P[2][0]*SF[6] + P[10][0]*SPP[4] - P[11][0]*SPP[5] - (P[12][0]*q0)/2) - SPP[1]*(P[3][2] + P[0][2]*SF[4] + P[1][2]*SF[3] + P[2][2]*SF[6] + P[10][2]*SPP[4] - P[11][2]*SPP[5] - (P[12][2]*q0)/2);
    nextP[3][7] = P[3][7] + P[0][7]*SF[4] + P[1][7]*SF[3] + P[2][7]*SF[6] + P[10][7]*SPP[4] - P[11][7]*SPP[5] - (P[12][7]*q0)/2 + dt*(P[3][4] + P[0][4]*SF[4] + P[1][4]*SF[3] + P[2][4]*SF[6] + P[10][4]*SPP[4] - P[11][4]*SPP[5] - (P[12][4]*q0)/2);
    nextP[3][8] = P[3][8] + P[0][8]*SF[4] + P[1][8]*SF[3] + P[2][8]*SF[6] + P[10][8]*SPP[4] - P[11][8]*SPP[5] - (P[12][8]*q0)/2 + dt*(P[3][5] + P[0][5]*SF[4] + P[1][5]*SF[3] + P[2][5]*SF[6] + P[10][5]*SPP[4] - P[11][5]*SPP[5] - (P[12][5]*q0)/2);
    nextP[3][9] = P[3][9] + P[0][9]*SF[4] + P[1][9]*SF[3] + P[2][9]*SF[6] + P[10][9]*SPP[4] - P[11][9]*SPP[5] - (P[12][9]*q0)/2 + dt*(P[3][6] + P[0][6]*SF[4] + P[1][6]*SF[3] + P[2][6]*SF[6] + P[10][6]*SPP[4] - P[11][6]*SPP[5] - (P[12][6]*q0)/2);
    nextP[3][10] = P[3][10] + P[0][10]*SF[4] + P[1][10]*SF[3] + P[2][10]*SF[6] + P[10][10]*SPP[4] - P[11][10]*SPP[5] - (P[12][10]*q0)/2;
    nextP[3][11] = P[3][11] + P[0][11]*SF[4] + P[1][11]*SF[3] + P[2][11]*SF[6] + P[10][11]*SPP[4] - P[11][11]*SPP[5] - (P[12][11]*q0)/2;
    nextP[3][12] = P[3][12] + P[0][12]*SF[4] + P[1][12]*SF[3] + P[2][12]*SF[6] + P[10][12]*SPP[4] - P[11][12]*SPP[5] - (P[12][12]*q0)/2;
    nextP[3][13] = P[3][13] + P[0][13]*SF[4] + P[1][13]*SF[3] + P[2][13]*SF[6] + P[10][13]*SPP[4] - P[11][13]*SPP[5] - (P[12][13]*q0)/2;
    nextP[3][14] = P[3][14] + P[0][14]*SF[4] + P[1][14]*SF[3] + P[2][14]*SF[6] + P[10][14]*SPP[4] - P[11][14]*SPP[5] - (P[12][14]*q0)/2;
    nextP[3][15] = P[3][15] + P[0][15]*SF[4] + P[1][15]*SF[3] + P[2][15]*SF[6] + P[10][15]*SPP[4] - P[11][15]*SPP[5] - (P[12][15]*q0)/2;
    nextP[3][16] = P[3][16] + P[0][16]*SF[4] + P[1][16]*SF[3] + P[2][16]*SF[6] + P[10][16]*SPP[4] - P[11][16]*SPP[5] - (P[12][16]*q0)/2;
    nextP[3][17] = P[3][17] + P[0][17]*SF[4] + P[1][17]*SF[3] + P[2][17]*SF[6] + P[10][17]*SPP[4] - P[11][17]*SPP[5] - (P[12][17]*q0)/2;
    nextP[3][18] = P[3][18] + P[0][18]*SF[4] + P[1][18]*SF[3] + P[2][18]*SF[6] + P[10][18]*SPP[4] - P[11][18]*SPP[5] - (P[12][18]*q0)/2;
    nextP[3][19] = P[3][19] + P[0][19]*SF[4] + P[1][19]*SF[3] + P[2][19]*SF[6] + P[10][19]*SPP[4] - P[11][19]*SPP[5] - (P[12][19]*q0)/2;
    nextP[3][20] = P[3][20] + P[0][20]*SF[4] + P[1][20]*SF[3] + P[2][20]*SF[6] + P[10][20]*SPP[4] - P[11][20]*SPP[5] - (P[12][20]*q0)/2;
    nextP[4][0] = P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2] + SF[6]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]) + SPP[7]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]) + SPP[6]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) + SPP[5]*(P[4][10] + P[0][10]*SF[2] + P[1][10]*SF[0] + P[2][10]*SPP[0] - P[3][10]*SPP[2]) + SPP[4]*(P[4][11] + P[0][11]*SF[2] + P[1][11]*SF[0] + P[2][11]*SPP[0] - P[3][11]*SPP[2]) + SPP[3]*(P[4][12] + P[0][12]*SF[2] + P[1][12]*SF[0] + P[2][12]*SPP[0] - P[3][12]*SPP[2]);
    nextP[4][1] = P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2] + SF[5]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) + SF[4]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]) + SPP[7]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) + SPP[3]*(P[4][11] + P[0][11]*SF[2] + P[1][11]*SF[0] + P[2][11]*SPP[0] - P[3][11]*SPP[2]) - SPP[4]*(P[4][12] + P[0][12]*SF[2] + P[1][12]*SF[0] + P[2][12]*SPP[0] - P[3][12]*SPP[2]) - (q0*(P[4][10] + P[0][10]*SF[2] + P[1][10]*SF[0] + P[2][10]*SPP[0] - P[3][10]*SPP[2]))/2;
    nextP[4][2] = P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2] + SF[3]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) + SF[5]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) + SPP[6]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]) - SPP[3]*(P[4][10] + P[0][10]*SF[2] + P[1][10]*SF[0] + P[2][10]*SPP[0] - P[3][10]*SPP[2]) + SPP[5]*(P[4][12] + P[0][12]*SF[2] + P[1][12]*SF[0] + P[2][12]*SPP[0] - P[3][12]*SPP[2]) - (q0*(P[4][11] + P[0][11]*SF[2] + P[1][11]*SF[0] + P[2][11]*SPP[0] - P[3][11]*SPP[2]))/2;
    nextP[4][3] = P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2] + SF[4]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) + SF[3]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]) + SF[6]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]) + SPP[4]*(P[4][10] + P[0][10]*SF[2] + P[1][10]*SF[0] + P[2][10]*SPP[0] - P[3][10]*SPP[2]) - SPP[5]*(P[4][11] + P[0][11]*SF[2] + P[1][11]*SF[0] + P[2][11]*SPP[0] - P[3][11]*SPP[2]) - (q0*(P[4][12] + P[0][12]*SF[2] + P[1][12]*SF[0] + P[2][12]*SPP[0] - P[3][12]*SPP[2]))/2;
    nextP[4][4] = P[4][4] + P[0][4]*SF[2] + P[1][4]*SF[0] + P[2][4]*SPP[0] - P[3][4]*SPP[2] + dvyCov*sq(SG[7] - 2*q0*q3) + dvzCov*sq(SG[6] + 2*q0*q2) + SF[2]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) + SF[0]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]) + SPP[0]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]) - SPP[2]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) + dvxCov*sq(SG[1] + SG[2] - SG[3] - SG[4]);
    nextP[4][5] = P[4][5] + SQ[2] + P[0][5]*SF[2] + P[1][5]*SF[0] + P[2][5]*SPP[0] - P[3][5]*SPP[2] + SF[1]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) + SF[0]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]) + SF[2]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) - SPP[0]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]);
    nextP[4][6] = P[4][6] + SQ[1] + P[0][6]*SF[2] + P[1][6]*SF[0] + P[2][6]*SPP[0] - P[3][6]*SPP[2] + SF[1]*(P[4][1] + P[0][1]*SF[2] + P[1][1]*SF[0] + P[2][1]*SPP[0] - P[3][1]*SPP[2]) + SF[0]*(P[4][3] + P[0][3]*SF[2] + P[1][3]*SF[0] + P[2][3]*SPP[0] - P[3][3]*SPP[2]) + SPP[0]*(P[4][0] + P[0][0]*SF[2] + P[1][0]*SF[0] + P[2][0]*SPP[0] - P[3][0]*SPP[2]) - SPP[1]*(P[4][2] + P[0][2]*SF[2] + P[1][2]*SF[0] + P[2][2]*SPP[0] - P[3][2]*SPP[2]);
    nextP[4][7] = P[4][7] + P[0][7]*SF[2] + P[1][7]*SF[0] + P[2][7]*SPP[0] - P[3][7]*SPP[2] + dt*(P[4][4] + P[0][4]*SF[2] + P[1][4]*SF[0] + P[2][4]*SPP[0] - P[3][4]*SPP[2]);
    nextP[4][8] = P[4][8] + P[0][8]*SF[2] + P[1][8]*SF[0] + P[2][8]*SPP[0] - P[3][8]*SPP[2] + dt*(P[4][5] + P[0][5]*SF[2] + P[1][5]*SF[0] + P[2][5]*SPP[0] - P[3][5]*SPP[2]);
    nextP[4][9] = P[4][9] + P[0][9]*SF[2] + P[1][9]*SF[0] + P[2][9]*SPP[0] - P[3][9]*SPP[2] + dt*(P[4][6] + P[0][6]*SF[2] + P[1][6]*SF[0] + P[2][6]*SPP[0] - P[3][6]*SPP[2]);
    nextP[4][10] = P[4][10] + P[0][10]*SF[2] + P[1][10]*SF[0] + P[2][10]*SPP[0] - P[3][10]*SPP[2];
    nextP[4][11] = P[4][11] + P[0][11]*SF[2] + P[1][11]*SF[0] + P[2][11]*SPP[0] - P[3][11]*SPP[2];
    nextP[4][12] = P[4][12] + P[0][12]*SF[2] + P[1][12]*SF[0] + P[2][12]*SPP[0] - P[3][12]*SPP[2];
    nextP[4][13] = P[4][13] + P[0][13]*SF[2] + P[1][13]*SF[0] + P[2][13]*SPP[0] - P[3][13]*SPP[2];
    nextP[4][14] = P[4][14] + P[0][14]*SF[2] + P[1][14]*SF[0] + P[2][14]*SPP[0] - P[3][14]*SPP[2];
    nextP[4][15] = P[4][15] + P[0][15]*SF[2] + P[1][15]*SF[0] + P[2][15]*SPP[0] - P[3][15]*SPP[2];
    nextP[4][16] = P[4][16] + P[0][16]*SF[2] + P[1][16]*SF[0] + P[2][16]*SPP[0] - P[3][16]*SPP[2];
    nextP[4][17] = P[4][17] + P[0][17]*SF[2] + P[1][17]*SF[0] + P[2][17]*SPP[0] - P[3][17]*SPP[2];
    nextP[4][18] = P[4][18] + P[0][18]*SF[2] + P[1][18]*SF[0] + P[2][18]*SPP[0] - P[3][18]*SPP[2];
    nextP[4][19] = P[4][19] + P[0][19]*SF[2] + P[1][19]*SF[0] + P[2][19]*SPP[0] - P[3][19]*SPP[2];
    nextP[4][20] = P[4][20] + P[0][20]*SF[2] + P[1][20]*SF[0] + P[2][20]*SPP[0] - P[3][20]*SPP[2];
    nextP[5][0] = P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0] + SF[6]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) + SPP[7]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]) + SPP[6]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]) + SPP[5]*(P[5][10] + P[0][10]*SF[1] + P[2][10]*SF[0] + P[3][10]*SF[2] - P[1][10]*SPP[0]) + SPP[4]*(P[5][11] + P[0][11]*SF[1] + P[2][11]*SF[0] + P[3][11]*SF[2] - P[1][11]*SPP[0]) + SPP[3]*(P[5][12] + P[0][12]*SF[1] + P[2][12]*SF[0] + P[3][12]*SF[2] - P[1][12]*SPP[0]);
    nextP[5][1] = P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0] + SF[5]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) + SF[4]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]) + SPP[7]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]) + SPP[3]*(P[5][11] + P[0][11]*SF[1] + P[2][11]*SF[0] + P[3][11]*SF[2] - P[1][11]*SPP[0]) - SPP[4]*(P[5][12] + P[0][12]*SF[1] + P[2][12]*SF[0] + P[3][12]*SF[2] - P[1][12]*SPP[0]) - (q0*(P[5][10] + P[0][10]*SF[1] + P[2][10]*SF[0] + P[3][10]*SF[2] - P[1][10]*SPP[0]))/2;
    nextP[5][2] = P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0] + SF[3]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) + SF[5]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]) + SPP[6]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) - SPP[3]*(P[5][10] + P[0][10]*SF[1] + P[2][10]*SF[0] + P[3][10]*SF[2] - P[1][10]*SPP[0]) + SPP[5]*(P[5][12] + P[0][12]*SF[1] + P[2][12]*SF[0] + P[3][12]*SF[2] - P[1][12]*SPP[0]) - (q0*(P[5][11] + P[0][11]*SF[1] + P[2][11]*SF[0] + P[3][11]*SF[2] - P[1][11]*SPP[0]))/2;
    nextP[5][3] = P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0] + SF[4]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) + SF[3]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) + SF[6]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]) + SPP[4]*(P[5][10] + P[0][10]*SF[1] + P[2][10]*SF[0] + P[3][10]*SF[2] - P[1][10]*SPP[0]) - SPP[5]*(P[5][11] + P[0][11]*SF[1] + P[2][11]*SF[0] + P[3][11]*SF[2] - P[1][11]*SPP[0]) - (q0*(P[5][12] + P[0][12]*SF[1] + P[2][12]*SF[0] + P[3][12]*SF[2] - P[1][12]*SPP[0]))/2;
    nextP[5][4] = P[5][4] + SQ[2] + P[0][4]*SF[1] + P[2][4]*SF[0] + P[3][4]*SF[2] - P[1][4]*SPP[0] + SF[2]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) + SF[0]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) + SPP[0]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]) - SPP[2]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]);
    nextP[5][5] = P[5][5] + P[0][5]*SF[1] + P[2][5]*SF[0] + P[3][5]*SF[2] - P[1][5]*SPP[0] + dvxCov*sq(SG[7] + 2*q0*q3) + dvzCov*sq(SG[5] - 2*q0*q1) + SF[1]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) + SF[0]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]) + SF[2]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]) - SPP[0]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) + dvyCov*sq(SG[1] - SG[2] + SG[3] - SG[4]);
    nextP[5][6] = P[5][6] + SQ[0] + P[0][6]*SF[1] + P[2][6]*SF[0] + P[3][6]*SF[2] - P[1][6]*SPP[0] + SF[1]*(P[5][1] + P[0][1]*SF[1] + P[2][1]*SF[0] + P[3][1]*SF[2] - P[1][1]*SPP[0]) + SF[0]*(P[5][3] + P[0][3]*SF[1] + P[2][3]*SF[0] + P[3][3]*SF[2] - P[1][3]*SPP[0]) + SPP[0]*(P[5][0] + P[0][0]*SF[1] + P[2][0]*SF[0] + P[3][0]*SF[2] - P[1][0]*SPP[0]) - SPP[1]*(P[5][2] + P[0][2]*SF[1] + P[2][2]*SF[0] + P[3][2]*SF[2] - P[1][2]*SPP[0]);
    nextP[5][7] = P[5][7] + P[0][7]*SF[1] + P[2][7]*SF[0] + P[3][7]*SF[2] - P[1][7]*SPP[0] + dt*(P[5][4] + P[0][4]*SF[1] + P[2][4]*SF[0] + P[3][4]*SF[2] - P[1][4]*SPP[0]);
    nextP[5][8] = P[5][8] + P[0][8]*SF[1] + P[2][8]*SF[0] + P[3][8]*SF[2] - P[1][8]*SPP[0] + dt*(P[5][5] + P[0][5]*SF[1] + P[2][5]*SF[0] + P[3][5]*SF[2] - P[1][5]*SPP[0]);
    nextP[5][9] = P[5][9] + P[0][9]*SF[1] + P[2][9]*SF[0] + P[3][9]*SF[2] - P[1][9]*SPP[0] + dt*(P[5][6] + P[0][6]*SF[1] + P[2][6]*SF[0] + P[3][6]*SF[2] - P[1][6]*SPP[0]);
    nextP[5][10] = P[5][10] + P[0][10]*SF[1] + P[2][10]*SF[0] + P[3][10]*SF[2] - P[1][10]*SPP[0];
    nextP[5][11] = P[5][11] + P[0][11]*SF[1] + P[2][11]*SF[0] + P[3][11]*SF[2] - P[1][11]*SPP[0];
    nextP[5][12] = P[5][12] + P[0][12]*SF[1] + P[2][12]*SF[0] + P[3][12]*SF[2] - P[1][12]*SPP[0];
    nextP[5][13] = P[5][13] + P[0][13]*SF[1] + P[2][13]*SF[0] + P[3][13]*SF[2] - P[1][13]*SPP[0];
    nextP[5][14] = P[5][14] + P[0][14]*SF[1] + P[2][14]*SF[0] + P[3][14]*SF[2] - P[1][14]*SPP[0];
    nextP[5][15] = P[5][15] + P[0][15]*SF[1] + P[2][15]*SF[0] + P[3][15]*SF[2] - P[1][15]*SPP[0];
    nextP[5][16] = P[5][16] + P[0][16]*SF[1] + P[2][16]*SF[0] + P[3][16]*SF[2] - P[1][16]*SPP[0];
    nextP[5][17] = P[5][17] + P[0][17]*SF[1] + P[2][17]*SF[0] + P[3][17]*SF[2] - P[1][17]*SPP[0];
    nextP[5][18] = P[5][18] + P[0][18]*SF[1] + P[2][18]*SF[0] + P[3][18]*SF[2] - P[1][18]*SPP[0];
    nextP[5][19] = P[5][19] + P[0][19]*SF[1] + P[2][19]*SF[0] + P[3][19]*SF[2] - P[1][19]*SPP[0];
    nextP[5][20] = P[5][20] + P[0][20]*SF[1] + P[2][20]*SF[0] + P[3][20]*SF[2] - P[1][20]*SPP[0];
    nextP[6][0] = P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1] + SF[6]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]) + SPP[7]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) + SPP[6]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]) + SPP[5]*(P[6][10] + P[1][10]*SF[1] + P[3][10]*SF[0] + P[0][10]*SPP[0] - P[2][10]*SPP[1]) + SPP[4]*(P[6][11] + P[1][11]*SF[1] + P[3][11]*SF[0] + P[0][11]*SPP[0] - P[2][11]*SPP[1]) + SPP[3]*(P[6][12] + P[1][12]*SF[1] + P[3][12]*SF[0] + P[0][12]*SPP[0] - P[2][12]*SPP[1]);
    nextP[6][1] = P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1] + SF[5]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) + SF[4]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) + SPP[7]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]) + SPP[3]*(P[6][11] + P[1][11]*SF[1] + P[3][11]*SF[0] + P[0][11]*SPP[0] - P[2][11]*SPP[1]) - SPP[4]*(P[6][12] + P[1][12]*SF[1] + P[3][12]*SF[0] + P[0][12]*SPP[0] - P[2][12]*SPP[1]) - (q0*(P[6][10] + P[1][10]*SF[1] + P[3][10]*SF[0] + P[0][10]*SPP[0] - P[2][10]*SPP[1]))/2;
    nextP[6][2] = P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1] + SF[3]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) + SF[5]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]) + SPP[6]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]) - SPP[3]*(P[6][10] + P[1][10]*SF[1] + P[3][10]*SF[0] + P[0][10]*SPP[0] - P[2][10]*SPP[1]) + SPP[5]*(P[6][12] + P[1][12]*SF[1] + P[3][12]*SF[0] + P[0][12]*SPP[0] - P[2][12]*SPP[1]) - (q0*(P[6][11] + P[1][11]*SF[1] + P[3][11]*SF[0] + P[0][11]*SPP[0] - P[2][11]*SPP[1]))/2;
    nextP[6][3] = P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1] + SF[4]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) + SF[3]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]) + SF[6]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) + SPP[4]*(P[6][10] + P[1][10]*SF[1] + P[3][10]*SF[0] + P[0][10]*SPP[0] - P[2][10]*SPP[1]) - SPP[5]*(P[6][11] + P[1][11]*SF[1] + P[3][11]*SF[0] + P[0][11]*SPP[0] - P[2][11]*SPP[1]) - (q0*(P[6][12] + P[1][12]*SF[1] + P[3][12]*SF[0] + P[0][12]*SPP[0] - P[2][12]*SPP[1]))/2;
    nextP[6][4] = P[6][4] + SQ[1] + P[1][4]*SF[1] + P[3][4]*SF[0] + P[0][4]*SPP[0] - P[2][4]*SPP[1] + SF[2]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) + SF[0]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]) + SPP[0]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) - SPP[2]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]);
    nextP[6][5] = P[6][5] + SQ[0] + P[1][5]*SF[1] + P[3][5]*SF[0] + P[0][5]*SPP[0] - P[2][5]*SPP[1] + SF[1]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) + SF[0]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) + SF[2]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]) - SPP[0]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]);
    nextP[6][6] = P[6][6] + P[1][6]*SF[1] + P[3][6]*SF[0] + P[0][6]*SPP[0] - P[2][6]*SPP[1] + dvxCov*sq(SG[6] - 2*q0*q2) + dvyCov*sq(SG[5] + 2*q0*q1) + SF[1]*(P[6][1] + P[1][1]*SF[1] + P[3][1]*SF[0] + P[0][1]*SPP[0] - P[2][1]*SPP[1]) + SF[0]*(P[6][3] + P[1][3]*SF[1] + P[3][3]*SF[0] + P[0][3]*SPP[0] - P[2][3]*SPP[1]) + SPP[0]*(P[6][0] + P[1][0]*SF[1] + P[3][0]*SF[0] + P[0][0]*SPP[0] - P[2][0]*SPP[1]) - SPP[1]*(P[6][2] + P[1][2]*SF[1] + P[3][2]*SF[0] + P[0][2]*SPP[0] - P[2][2]*SPP[1]) + dvzCov*sq(SG[1] - SG[2] - SG[3] + SG[4]);
    nextP[6][7] = P[6][7] + P[1][7]*SF[1] + P[3][7]*SF[0] + P[0][7]*SPP[0] - P[2][7]*SPP[1] + dt*(P[6][4] + P[1][4]*SF[1] + P[3][4]*SF[0] + P[0][4]*SPP[0] - P[2][4]*SPP[1]);
    nextP[6][8] = P[6][8] + P[1][8]*SF[1] + P[3][8]*SF[0] + P[0][8]*SPP[0] - P[2][8]*SPP[1] + dt*(P[6][5] + P[1][5]*SF[1] + P[3][5]*SF[0] + P[0][5]*SPP[0] - P[2][5]*SPP[1]);
    nextP[6][9] = P[6][9] + P[1][9]*SF[1] + P[3][9]*SF[0] + P[0][9]*SPP[0] - P[2][9]*SPP[1] + dt*(P[6][6] + P[1][6]*SF[1] + P[3][6]*SF[0] + P[0][6]*SPP[0] - P[2][6]*SPP[1]);
    nextP[6][10] = P[6][10] + P[1][10]*SF[1] + P[3][10]*SF[0] + P[0][10]*SPP[0] - P[2][10]*SPP[1];
    nextP[6][11] = P[6][11] + P[1][11]*SF[1] + P[3][11]*SF[0] + P[0][11]*SPP[0] - P[2][11]*SPP[1];
    nextP[6][12] = P[6][12] + P[1][12]*SF[1] + P[3][12]*SF[0] + P[0][12]*SPP[0] - P[2][12]*SPP[1];
    nextP[6][13] = P[6][13] + P[1][13]*SF[1] + P[3][13]*SF[0] + P[0][13]*SPP[0] - P[2][13]*SPP[1];
    nextP[6][14] = P[6][14] + P[1][14]*SF[1] + P[3][14]*SF[0] + P[0][14]*SPP[0] - P[2][14]*SPP[1];
    nextP[6][15] = P[6][15] + P[1][15]*SF[1] + P[3][15]*SF[0] + P[0][15]*SPP[0] - P[2][15]*SPP[1];
    nextP[6][16] = P[6][16] + P[1][16]*SF[1] + P[3][16]*SF[0] + P[0][16]*SPP[0] - P[2][16]*SPP[1];
    nextP[6][17] = P[6][17] + P[1][17]*SF[1] + P[3][17]*SF[0] + P[0][17]*SPP[0] - P[2][17]*SPP[1];
    nextP[6][18] = P[6][18] + P[1][18]*SF[1] + P[3][18]*SF[0] + P[0][18]*SPP[0] - P[2][18]*SPP[1];
    nextP[6][19] = P[6][19] + P[1][19]*SF[1] + P[3][19]*SF[0] + P[0][19]*SPP[0] - P[2][19]*SPP[1];
    nextP[6][20] = P[6][20] + P[1][20]*SF[1] + P[3][20]*SF[0] + P[0][20]*SPP[0] - P[2][20]*SPP[1];
    nextP[7][0] = P[7][0] + P[4][0]*dt + SF[6]*(P[7][1] + P[4][1]*dt) + SPP[7]*(P[7][2] + P[4][2]*dt) + SPP[6]*(P[7][3] + P[4][3]*dt) + SPP[5]*(P[7][10] + P[4][10]*dt) + SPP[4]*(P[7][11] + P[4][11]*dt) + SPP[3]*(P[7][12] + P[4][12]*dt);
    nextP[7][1] = P[7][1] + P[4][1]*dt + SF[5]*(P[7][0] + P[4][0]*dt) + SF[4]*(P[7][2] + P[4][2]*dt) + SPP[7]*(P[7][3] + P[4][3]*dt) + SPP[3]*(P[7][11] + P[4][11]*dt) - SPP[4]*(P[7][12] + P[4][12]*dt) - (q0*(P[7][10] + P[4][10]*dt))/2;
    nextP[7][2] = P[7][2] + P[4][2]*dt + SF[3]*(P[7][0] + P[4][0]*dt) + SF[5]*(P[7][3] + P[4][3]*dt) + SPP[6]*(P[7][1] + P[4][1]*dt) - SPP[3]*(P[7][10] + P[4][10]*dt) + SPP[5]*(P[7][12] + P[4][12]*dt) - (q0*(P[7][11] + P[4][11]*dt))/2;
    nextP[7][3] = P[7][3] + P[4][3]*dt + SF[4]*(P[7][0] + P[4][0]*dt) + SF[3]*(P[7][1] + P[4][1]*dt) + SF[6]*(P[7][2] + P[4][2]*dt) + SPP[4]*(P[7][10] + P[4][10]*dt) - SPP[5]*(P[7][11] + P[4][11]*dt) - (q0*(P[7][12] + P[4][12]*dt))/2;
    nextP[7][4] = P[7][4] + P[4][4]*dt + SF[0]*(P[7][1] + P[4][1]*dt) + SF[2]*(P[7][0] + P[4][0]*dt) + SPP[0]*(P[7][2] + P[4][2]*dt) - SPP[2]*(P[7][3] + P[4][3]*dt);
    nextP[7][5] = P[7][5] + P[4][5]*dt + SF[1]*(P[7][0] + P[4][0]*dt) + SF[0]*(P[7][2] + P[4][2]*dt) + SF[2]*(P[7][3] + P[4][3]*dt) - SPP[0]*(P[7][1] + P[4][1]*dt);
    nextP[7][6] = P[7][6] + P[4][6]*dt + SF[1]*(P[7][1] + P[4][1]*dt) + SF[0]*(P[7][3] + P[4][3]*dt) + SPP[0]*(P[7][0] + P[4][0]*dt) - SPP[1]*(P[7][2] + P[4][2]*dt);
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
    nextP[8][0] = P[8][0] + P[5][0]*dt + SF[6]*(P[8][1] + P[5][1]*dt) + SPP[7]*(P[8][2] + P[5][2]*dt) + SPP[6]*(P[8][3] + P[5][3]*dt) + SPP[5]*(P[8][10] + P[5][10]*dt) + SPP[4]*(P[8][11] + P[5][11]*dt) + SPP[3]*(P[8][12] + P[5][12]*dt);
    nextP[8][1] = P[8][1] + P[5][1]*dt + SF[5]*(P[8][0] + P[5][0]*dt) + SF[4]*(P[8][2] + P[5][2]*dt) + SPP[7]*(P[8][3] + P[5][3]*dt) + SPP[3]*(P[8][11] + P[5][11]*dt) - SPP[4]*(P[8][12] + P[5][12]*dt) - (q0*(P[8][10] + P[5][10]*dt))/2;
    nextP[8][2] = P[8][2] + P[5][2]*dt + SF[3]*(P[8][0] + P[5][0]*dt) + SF[5]*(P[8][3] + P[5][3]*dt) + SPP[6]*(P[8][1] + P[5][1]*dt) - SPP[3]*(P[8][10] + P[5][10]*dt) + SPP[5]*(P[8][12] + P[5][12]*dt) - (q0*(P[8][11] + P[5][11]*dt))/2;
    nextP[8][3] = P[8][3] + P[5][3]*dt + SF[4]*(P[8][0] + P[5][0]*dt) + SF[3]*(P[8][1] + P[5][1]*dt) + SF[6]*(P[8][2] + P[5][2]*dt) + SPP[4]*(P[8][10] + P[5][10]*dt) - SPP[5]*(P[8][11] + P[5][11]*dt) - (q0*(P[8][12] + P[5][12]*dt))/2;
    nextP[8][4] = P[8][4] + P[5][4]*dt + SF[0]*(P[8][1] + P[5][1]*dt) + SF[2]*(P[8][0] + P[5][0]*dt) + SPP[0]*(P[8][2] + P[5][2]*dt) - SPP[2]*(P[8][3] + P[5][3]*dt);
    nextP[8][5] = P[8][5] + P[5][5]*dt + SF[1]*(P[8][0] + P[5][0]*dt) + SF[0]*(P[8][2] + P[5][2]*dt) + SF[2]*(P[8][3] + P[5][3]*dt) - SPP[0]*(P[8][1] + P[5][1]*dt);
    nextP[8][6] = P[8][6] + P[5][6]*dt + SF[1]*(P[8][1] + P[5][1]*dt) + SF[0]*(P[8][3] + P[5][3]*dt) + SPP[0]*(P[8][0] + P[5][0]*dt) - SPP[1]*(P[8][2] + P[5][2]*dt);
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
    nextP[9][0] = P[9][0] + P[6][0]*dt + SF[6]*(P[9][1] + P[6][1]*dt) + SPP[7]*(P[9][2] + P[6][2]*dt) + SPP[6]*(P[9][3] + P[6][3]*dt) + SPP[5]*(P[9][10] + P[6][10]*dt) + SPP[4]*(P[9][11] + P[6][11]*dt) + SPP[3]*(P[9][12] + P[6][12]*dt);
    nextP[9][1] = P[9][1] + P[6][1]*dt + SF[5]*(P[9][0] + P[6][0]*dt) + SF[4]*(P[9][2] + P[6][2]*dt) + SPP[7]*(P[9][3] + P[6][3]*dt) + SPP[3]*(P[9][11] + P[6][11]*dt) - SPP[4]*(P[9][12] + P[6][12]*dt) - (q0*(P[9][10] + P[6][10]*dt))/2;
    nextP[9][2] = P[9][2] + P[6][2]*dt + SF[3]*(P[9][0] + P[6][0]*dt) + SF[5]*(P[9][3] + P[6][3]*dt) + SPP[6]*(P[9][1] + P[6][1]*dt) - SPP[3]*(P[9][10] + P[6][10]*dt) + SPP[5]*(P[9][12] + P[6][12]*dt) - (q0*(P[9][11] + P[6][11]*dt))/2;
    nextP[9][3] = P[9][3] + P[6][3]*dt + SF[4]*(P[9][0] + P[6][0]*dt) + SF[3]*(P[9][1] + P[6][1]*dt) + SF[6]*(P[9][2] + P[6][2]*dt) + SPP[4]*(P[9][10] + P[6][10]*dt) - SPP[5]*(P[9][11] + P[6][11]*dt) - (q0*(P[9][12] + P[6][12]*dt))/2;
    nextP[9][4] = P[9][4] + P[6][4]*dt + SF[0]*(P[9][1] + P[6][1]*dt) + SF[2]*(P[9][0] + P[6][0]*dt) + SPP[0]*(P[9][2] + P[6][2]*dt) - SPP[2]*(P[9][3] + P[6][3]*dt);
    nextP[9][5] = P[9][5] + P[6][5]*dt + SF[1]*(P[9][0] + P[6][0]*dt) + SF[0]*(P[9][2] + P[6][2]*dt) + SF[2]*(P[9][3] + P[6][3]*dt) - SPP[0]*(P[9][1] + P[6][1]*dt);
    nextP[9][6] = P[9][6] + P[6][6]*dt + SF[1]*(P[9][1] + P[6][1]*dt) + SF[0]*(P[9][3] + P[6][3]*dt) + SPP[0]*(P[9][0] + P[6][0]*dt) - SPP[1]*(P[9][2] + P[6][2]*dt);
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
    nextP[10][0] = P[10][0] + P[10][1]*SF[6] + P[10][2]*SPP[7] + P[10][3]*SPP[6] + P[10][10]*SPP[5] + P[10][11]*SPP[4] + P[10][12]*SPP[3];
    nextP[10][1] = P[10][1] + P[10][0]*SF[5] + P[10][2]*SF[4] + P[10][3]*SPP[7] + P[10][11]*SPP[3] - P[10][12]*SPP[4] - (P[10][10]*q0)/2;
    nextP[10][2] = P[10][2] + P[10][0]*SF[3] + P[10][3]*SF[5] + P[10][1]*SPP[6] - P[10][10]*SPP[3] + P[10][12]*SPP[5] - (P[10][11]*q0)/2;
    nextP[10][3] = P[10][3] + P[10][0]*SF[4] + P[10][1]*SF[3] + P[10][2]*SF[6] + P[10][10]*SPP[4] - P[10][11]*SPP[5] - (P[10][12]*q0)/2;
    nextP[10][4] = P[10][4] + P[10][1]*SF[0] + P[10][0]*SF[2] + P[10][2]*SPP[0] - P[10][3]*SPP[2];
    nextP[10][5] = P[10][5] + P[10][0]*SF[1] + P[10][2]*SF[0] + P[10][3]*SF[2] - P[10][1]*SPP[0];
    nextP[10][6] = P[10][6] + P[10][1]*SF[1] + P[10][3]*SF[0] + P[10][0]*SPP[0] - P[10][2]*SPP[1];
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
    nextP[11][0] = P[11][0] + P[11][1]*SF[6] + P[11][2]*SPP[7] + P[11][3]*SPP[6] + P[11][10]*SPP[5] + P[11][11]*SPP[4] + P[11][12]*SPP[3];
    nextP[11][1] = P[11][1] + P[11][0]*SF[5] + P[11][2]*SF[4] + P[11][3]*SPP[7] + P[11][11]*SPP[3] - P[11][12]*SPP[4] - (P[11][10]*q0)/2;
    nextP[11][2] = P[11][2] + P[11][0]*SF[3] + P[11][3]*SF[5] + P[11][1]*SPP[6] - P[11][10]*SPP[3] + P[11][12]*SPP[5] - (P[11][11]*q0)/2;
    nextP[11][3] = P[11][3] + P[11][0]*SF[4] + P[11][1]*SF[3] + P[11][2]*SF[6] + P[11][10]*SPP[4] - P[11][11]*SPP[5] - (P[11][12]*q0)/2;
    nextP[11][4] = P[11][4] + P[11][1]*SF[0] + P[11][0]*SF[2] + P[11][2]*SPP[0] - P[11][3]*SPP[2];
    nextP[11][5] = P[11][5] + P[11][0]*SF[1] + P[11][2]*SF[0] + P[11][3]*SF[2] - P[11][1]*SPP[0];
    nextP[11][6] = P[11][6] + P[11][1]*SF[1] + P[11][3]*SF[0] + P[11][0]*SPP[0] - P[11][2]*SPP[1];
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
    nextP[12][0] = P[12][0] + P[12][1]*SF[6] + P[12][2]*SPP[7] + P[12][3]*SPP[6] + P[12][10]*SPP[5] + P[12][11]*SPP[4] + P[12][12]*SPP[3];
    nextP[12][1] = P[12][1] + P[12][0]*SF[5] + P[12][2]*SF[4] + P[12][3]*SPP[7] + P[12][11]*SPP[3] - P[12][12]*SPP[4] - (P[12][10]*q0)/2;
    nextP[12][2] = P[12][2] + P[12][0]*SF[3] + P[12][3]*SF[5] + P[12][1]*SPP[6] - P[12][10]*SPP[3] + P[12][12]*SPP[5] - (P[12][11]*q0)/2;
    nextP[12][3] = P[12][3] + P[12][0]*SF[4] + P[12][1]*SF[3] + P[12][2]*SF[6] + P[12][10]*SPP[4] - P[12][11]*SPP[5] - (P[12][12]*q0)/2;
    nextP[12][4] = P[12][4] + P[12][1]*SF[0] + P[12][0]*SF[2] + P[12][2]*SPP[0] - P[12][3]*SPP[2];
    nextP[12][5] = P[12][5] + P[12][0]*SF[1] + P[12][2]*SF[0] + P[12][3]*SF[2] - P[12][1]*SPP[0];
    nextP[12][6] = P[12][6] + P[12][1]*SF[1] + P[12][3]*SF[0] + P[12][0]*SPP[0] - P[12][2]*SPP[1];
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
    nextP[13][0] = P[13][0] + P[13][1]*SF[6] + P[13][2]*SPP[7] + P[13][3]*SPP[6] + P[13][10]*SPP[5] + P[13][11]*SPP[4] + P[13][12]*SPP[3];
    nextP[13][1] = P[13][1] + P[13][0]*SF[5] + P[13][2]*SF[4] + P[13][3]*SPP[7] + P[13][11]*SPP[3] - P[13][12]*SPP[4] - (P[13][10]*q0)/2;
    nextP[13][2] = P[13][2] + P[13][0]*SF[3] + P[13][3]*SF[5] + P[13][1]*SPP[6] - P[13][10]*SPP[3] + P[13][12]*SPP[5] - (P[13][11]*q0)/2;
    nextP[13][3] = P[13][3] + P[13][0]*SF[4] + P[13][1]*SF[3] + P[13][2]*SF[6] + P[13][10]*SPP[4] - P[13][11]*SPP[5] - (P[13][12]*q0)/2;
    nextP[13][4] = P[13][4] + P[13][1]*SF[0] + P[13][0]*SF[2] + P[13][2]*SPP[0] - P[13][3]*SPP[2];
    nextP[13][5] = P[13][5] + P[13][0]*SF[1] + P[13][2]*SF[0] + P[13][3]*SF[2] - P[13][1]*SPP[0];
    nextP[13][6] = P[13][6] + P[13][1]*SF[1] + P[13][3]*SF[0] + P[13][0]*SPP[0] - P[13][2]*SPP[1];
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
    nextP[14][0] = P[14][0] + P[14][1]*SF[6] + P[14][2]*SPP[7] + P[14][3]*SPP[6] + P[14][10]*SPP[5] + P[14][11]*SPP[4] + P[14][12]*SPP[3];
    nextP[14][1] = P[14][1] + P[14][0]*SF[5] + P[14][2]*SF[4] + P[14][3]*SPP[7] + P[14][11]*SPP[3] - P[14][12]*SPP[4] - (P[14][10]*q0)/2;
    nextP[14][2] = P[14][2] + P[14][0]*SF[3] + P[14][3]*SF[5] + P[14][1]*SPP[6] - P[14][10]*SPP[3] + P[14][12]*SPP[5] - (P[14][11]*q0)/2;
    nextP[14][3] = P[14][3] + P[14][0]*SF[4] + P[14][1]*SF[3] + P[14][2]*SF[6] + P[14][10]*SPP[4] - P[14][11]*SPP[5] - (P[14][12]*q0)/2;
    nextP[14][4] = P[14][4] + P[14][1]*SF[0] + P[14][0]*SF[2] + P[14][2]*SPP[0] - P[14][3]*SPP[2];
    nextP[14][5] = P[14][5] + P[14][0]*SF[1] + P[14][2]*SF[0] + P[14][3]*SF[2] - P[14][1]*SPP[0];
    nextP[14][6] = P[14][6] + P[14][1]*SF[1] + P[14][3]*SF[0] + P[14][0]*SPP[0] - P[14][2]*SPP[1];
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
    nextP[15][0] = P[15][0] + P[15][1]*SF[6] + P[15][2]*SPP[7] + P[15][3]*SPP[6] + P[15][10]*SPP[5] + P[15][11]*SPP[4] + P[15][12]*SPP[3];
    nextP[15][1] = P[15][1] + P[15][0]*SF[5] + P[15][2]*SF[4] + P[15][3]*SPP[7] + P[15][11]*SPP[3] - P[15][12]*SPP[4] - (P[15][10]*q0)/2;
    nextP[15][2] = P[15][2] + P[15][0]*SF[3] + P[15][3]*SF[5] + P[15][1]*SPP[6] - P[15][10]*SPP[3] + P[15][12]*SPP[5] - (P[15][11]*q0)/2;
    nextP[15][3] = P[15][3] + P[15][0]*SF[4] + P[15][1]*SF[3] + P[15][2]*SF[6] + P[15][10]*SPP[4] - P[15][11]*SPP[5] - (P[15][12]*q0)/2;
    nextP[15][4] = P[15][4] + P[15][1]*SF[0] + P[15][0]*SF[2] + P[15][2]*SPP[0] - P[15][3]*SPP[2];
    nextP[15][5] = P[15][5] + P[15][0]*SF[1] + P[15][2]*SF[0] + P[15][3]*SF[2] - P[15][1]*SPP[0];
    nextP[15][6] = P[15][6] + P[15][1]*SF[1] + P[15][3]*SF[0] + P[15][0]*SPP[0] - P[15][2]*SPP[1];
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
    nextP[16][0] = P[16][0] + P[16][1]*SF[6] + P[16][2]*SPP[7] + P[16][3]*SPP[6] + P[16][10]*SPP[5] + P[16][11]*SPP[4] + P[16][12]*SPP[3];
    nextP[16][1] = P[16][1] + P[16][0]*SF[5] + P[16][2]*SF[4] + P[16][3]*SPP[7] + P[16][11]*SPP[3] - P[16][12]*SPP[4] - (P[16][10]*q0)/2;
    nextP[16][2] = P[16][2] + P[16][0]*SF[3] + P[16][3]*SF[5] + P[16][1]*SPP[6] - P[16][10]*SPP[3] + P[16][12]*SPP[5] - (P[16][11]*q0)/2;
    nextP[16][3] = P[16][3] + P[16][0]*SF[4] + P[16][1]*SF[3] + P[16][2]*SF[6] + P[16][10]*SPP[4] - P[16][11]*SPP[5] - (P[16][12]*q0)/2;
    nextP[16][4] = P[16][4] + P[16][1]*SF[0] + P[16][0]*SF[2] + P[16][2]*SPP[0] - P[16][3]*SPP[2];
    nextP[16][5] = P[16][5] + P[16][0]*SF[1] + P[16][2]*SF[0] + P[16][3]*SF[2] - P[16][1]*SPP[0];
    nextP[16][6] = P[16][6] + P[16][1]*SF[1] + P[16][3]*SF[0] + P[16][0]*SPP[0] - P[16][2]*SPP[1];
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
    nextP[17][0] = P[17][0] + P[17][1]*SF[6] + P[17][2]*SPP[7] + P[17][3]*SPP[6] + P[17][10]*SPP[5] + P[17][11]*SPP[4] + P[17][12]*SPP[3];
    nextP[17][1] = P[17][1] + P[17][0]*SF[5] + P[17][2]*SF[4] + P[17][3]*SPP[7] + P[17][11]*SPP[3] - P[17][12]*SPP[4] - (P[17][10]*q0)/2;
    nextP[17][2] = P[17][2] + P[17][0]*SF[3] + P[17][3]*SF[5] + P[17][1]*SPP[6] - P[17][10]*SPP[3] + P[17][12]*SPP[5] - (P[17][11]*q0)/2;
    nextP[17][3] = P[17][3] + P[17][0]*SF[4] + P[17][1]*SF[3] + P[17][2]*SF[6] + P[17][10]*SPP[4] - P[17][11]*SPP[5] - (P[17][12]*q0)/2;
    nextP[17][4] = P[17][4] + P[17][1]*SF[0] + P[17][0]*SF[2] + P[17][2]*SPP[0] - P[17][3]*SPP[2];
    nextP[17][5] = P[17][5] + P[17][0]*SF[1] + P[17][2]*SF[0] + P[17][3]*SF[2] - P[17][1]*SPP[0];
    nextP[17][6] = P[17][6] + P[17][1]*SF[1] + P[17][3]*SF[0] + P[17][0]*SPP[0] - P[17][2]*SPP[1];
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
    nextP[18][0] = P[18][0] + P[18][1]*SF[6] + P[18][2]*SPP[7] + P[18][3]*SPP[6] + P[18][10]*SPP[5] + P[18][11]*SPP[4] + P[18][12]*SPP[3];
    nextP[18][1] = P[18][1] + P[18][0]*SF[5] + P[18][2]*SF[4] + P[18][3]*SPP[7] + P[18][11]*SPP[3] - P[18][12]*SPP[4] - (P[18][10]*q0)/2;
    nextP[18][2] = P[18][2] + P[18][0]*SF[3] + P[18][3]*SF[5] + P[18][1]*SPP[6] - P[18][10]*SPP[3] + P[18][12]*SPP[5] - (P[18][11]*q0)/2;
    nextP[18][3] = P[18][3] + P[18][0]*SF[4] + P[18][1]*SF[3] + P[18][2]*SF[6] + P[18][10]*SPP[4] - P[18][11]*SPP[5] - (P[18][12]*q0)/2;
    nextP[18][4] = P[18][4] + P[18][1]*SF[0] + P[18][0]*SF[2] + P[18][2]*SPP[0] - P[18][3]*SPP[2];
    nextP[18][5] = P[18][5] + P[18][0]*SF[1] + P[18][2]*SF[0] + P[18][3]*SF[2] - P[18][1]*SPP[0];
    nextP[18][6] = P[18][6] + P[18][1]*SF[1] + P[18][3]*SF[0] + P[18][0]*SPP[0] - P[18][2]*SPP[1];
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
    nextP[19][0] = P[19][0] + P[19][1]*SF[6] + P[19][2]*SPP[7] + P[19][3]*SPP[6] + P[19][10]*SPP[5] + P[19][11]*SPP[4] + P[19][12]*SPP[3];
    nextP[19][1] = P[19][1] + P[19][0]*SF[5] + P[19][2]*SF[4] + P[19][3]*SPP[7] + P[19][11]*SPP[3] - P[19][12]*SPP[4] - (P[19][10]*q0)/2;
    nextP[19][2] = P[19][2] + P[19][0]*SF[3] + P[19][3]*SF[5] + P[19][1]*SPP[6] - P[19][10]*SPP[3] + P[19][12]*SPP[5] - (P[19][11]*q0)/2;
    nextP[19][3] = P[19][3] + P[19][0]*SF[4] + P[19][1]*SF[3] + P[19][2]*SF[6] + P[19][10]*SPP[4] - P[19][11]*SPP[5] - (P[19][12]*q0)/2;
    nextP[19][4] = P[19][4] + P[19][1]*SF[0] + P[19][0]*SF[2] + P[19][2]*SPP[0] - P[19][3]*SPP[2];
    nextP[19][5] = P[19][5] + P[19][0]*SF[1] + P[19][2]*SF[0] + P[19][3]*SF[2] - P[19][1]*SPP[0];
    nextP[19][6] = P[19][6] + P[19][1]*SF[1] + P[19][3]*SF[0] + P[19][0]*SPP[0] - P[19][2]*SPP[1];
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
    nextP[20][0] = P[20][0] + P[20][1]*SF[6] + P[20][2]*SPP[7] + P[20][3]*SPP[6] + P[20][10]*SPP[5] + P[20][11]*SPP[4] + P[20][12]*SPP[3];
    nextP[20][1] = P[20][1] + P[20][0]*SF[5] + P[20][2]*SF[4] + P[20][3]*SPP[7] + P[20][11]*SPP[3] - P[20][12]*SPP[4] - (P[20][10]*q0)/2;
    nextP[20][2] = P[20][2] + P[20][0]*SF[3] + P[20][3]*SF[5] + P[20][1]*SPP[6] - P[20][10]*SPP[3] + P[20][12]*SPP[5] - (P[20][11]*q0)/2;
    nextP[20][3] = P[20][3] + P[20][0]*SF[4] + P[20][1]*SF[3] + P[20][2]*SF[6] + P[20][10]*SPP[4] - P[20][11]*SPP[5] - (P[20][12]*q0)/2;
    nextP[20][4] = P[20][4] + P[20][1]*SF[0] + P[20][0]*SF[2] + P[20][2]*SPP[0] - P[20][3]*SPP[2];
    nextP[20][5] = P[20][5] + P[20][0]*SF[1] + P[20][2]*SF[0] + P[20][3]*SF[2] - P[20][1]*SPP[0];
    nextP[20][6] = P[20][6] + P[20][1]*SF[1] + P[20][3]*SF[0] + P[20][0]*SPP[0] - P[20][2]*SPP[1];
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

    for (uint8_t i=0; i<= 20; i++)
    {
        nextP[i][i] = nextP[i][i] + processNoise[i];
    }

    // If on ground or no magnetometer fitted, inhibit magnetometer bias updates by
    // setting the coresponding covariance terms to zero.
    if (onGround || !useCompass)
    {
        zeroRows(nextP,15,20);
        zeroCols(nextP,15,20);
    }

    // If on ground or not using airspeed sensing, inhibit wind velocity
    // covariance growth.
    if (onGround || !useAirspeed)
    {
        zeroRows(nextP,13,14);
        zeroCols(nextP,13,14);
    }

    // If the total position variance exceds 1E6 (1000m), then stop covariance
    // growth by setting the predicted to the previous values
    // This prevent an ill conditioned matrix from occurring for long periods
    // without GPS
    if ((P[7][7] + P[8][8]) > 1E6f)
    {
        for (uint8_t i=7; i<=8; i++)
        {
            for (uint8_t j=0; j<=20; j++)
            {
                nextP[i][j] = P[i][j];
                nextP[j][i] = P[j][i];
            }
        }
    }

    // Force symmetry on the covariance matrix to prevent ill-conditioning
    // of the matrix which would cause the filter to blow-up
    for (uint8_t i=0; i<=20; i++) P[i][i] = nextP[i][i];
    for (uint8_t i=1; i<=20; i++)
    {
        for (uint8_t j=0; j<=i-1; j++)
        {
            P[i][j] = 0.5f*(nextP[i][j] + nextP[j][i]);
            P[j][i] = P[i][j];
        }
    }
}

void FuseVelposNED()
{

// declare variables used by fault isolation logic
    uint32_t gpsRetryTime = 30000; // time in msec before GPS fusion will be retried following innovation consistency failure
    uint32_t gpsRetryTimeNoTAS = 5000; // retry time if no TAS measurement available
    uint32_t hgtRetryTime = 5000; // height measurement retry time
    uint32_t horizRetryTime;
    static uint32_t velFailTime;
    static uint32_t posFailTime;
    static uint32_t hgtFailTime;
    bool velHealth;
    bool posHealth;
    bool hgtHealth;
    bool velTimeout;
    bool posTimeout;
    bool hgtTimeout;

// declare variables used to check measurement errors
    float velInnov[3] = {0.0,0.0,0.0};
    float posInnov[2] = {0.0,0.0};
    float hgtInnov = 0.0;

// declare variables used to control access to arrays
    bool fuseData[6] = {false,false,false,false,false,false};
    uint8_t stateIndex;
    uint8_t obsIndex;
    uint8_t indexLimit;

// declare variables used by state and covariance update calculations
    float velErr;
    float posErr;
    float R_OBS[6];
    float observation[6];
    float SK;
    float quatMag;

// Perform sequential fusion of GPS measurements. This assumes that the
// errors in the different velocity and position components are
// uncorrelated which is not true, however in the absence of covariance
// data from the GPS receiver it is the only assumption we can make
// so we might as well take advantage of the computational efficiencies
// associated with sequential fusion
    if (fuseVelData || fusePosData || fuseHgtData)
    {
        // set the GPS data timeout depending on whether airspeed data is present
        if (useAirspeed) horizRetryTime = gpsRetryTime;
        else horizRetryTime = gpsRetryTimeNoTAS;

        // Form the observation vector
        for (uint8_t i=0; i<=2; i++) observation[i] = velNED[i];
        for (uint8_t i=3; i<=4; i++) observation[i] = posNE[i-3];
        observation[5] = -(hgtMea);

        // Estimate the GPS Velocity, GPS horiz position and height measurement variances.
        velErr = 0.2*accNavMag; // additional error in GPS velocities caused by manoeuvring
        posErr = 0.2*accNavMag; // additional error in GPS position caused by manoeuvring
        R_OBS[0] = 0.04f + sq(velErr);
        R_OBS[1] = R_OBS[0];
        R_OBS[2] = 0.08f + sq(velErr);
        R_OBS[3] = R_OBS[2];
        R_OBS[4] = 4.0f  + sq(posErr);
        R_OBS[5] = 4.0f;

        // Set innovation variances to zero default
        for (uint8_t i = 0; i<=5; i++)
        {
            varInnovVelPos[i] = 0.0f;
        }
        // calculate innovations and check GPS data validity using an innovation consistency check
        if (fuseVelData)
        {
            // test velocity measurements
            uint8_t imax = 2;
            if (fusionModeGPS == 1) imax = 1;
            for (uint8_t i = 0; i<=imax; i++)
            {
                velInnov[i] = statesAtVelTime[i+4] - velNED[i];
                stateIndex = 4 + i;
                varInnovVelPos[i] = P[stateIndex][stateIndex] + R_OBS[i];
            }
            // apply a 5-sigma threshold
            velHealth = (sq(velInnov[0]) + sq(velInnov[1]) + sq(velInnov[2])) < 25.0*(varInnovVelPos[0] + varInnovVelPos[1] + varInnovVelPos[2]);
            velTimeout = (millis() - velFailTime) > horizRetryTime;
            if (velHealth || velTimeout)
            {
                velHealth = true;
                velFailTime = millis();
            }
            else
            {
                velHealth = false;
            }
        }
        if (fusePosData)
        {
            // test horizontal position measurements
            posInnov[0] = statesAtPosTime[7] - posNE[0];
            posInnov[1] = statesAtPosTime[8] - posNE[1];
            varInnovVelPos[3] = P[7][7] + R_OBS[3];
            varInnovVelPos[4] = P[8][8] + R_OBS[4];
            // apply a 10-sigma threshold
            posHealth = (sq(posInnov[0]) + sq(posInnov[1])) < 100.0f*(varInnovVelPos[3] + varInnovVelPos[4]);
            posTimeout = (millis() - posFailTime) > horizRetryTime;
            if (posHealth || posTimeout)
            {
                posHealth = true;
                posFailTime = millis();
            }
            else
            {
                posHealth = false;
            }
        }
        // test height measurements
        if (fuseHgtData)
        {
            hgtInnov = statesAtHgtTime[9] + hgtMea;
            varInnovVelPos[5] = P[9][9] + R_OBS[5];
            // apply a 10-sigma threshold
            hgtHealth = sq(hgtInnov) < 100.0f*varInnovVelPos[5];
            hgtTimeout = (millis() - hgtFailTime) > hgtRetryTime;
            if (hgtHealth || hgtTimeout)
            {
                hgtHealth = true;
                hgtFailTime = millis();
            }
            else
            {
                hgtHealth = false;
            }
        }
        // Set range for sequential fusion of velocity and position measurements depending
        // on which data is available and its health
        if (fuseVelData && fusionModeGPS == 0 && velHealth)
        {
            fuseData[0] = true;
            fuseData[1] = true;
            fuseData[2] = true;
        }
        if (fuseVelData && fusionModeGPS == 1 && velHealth)
        {
            fuseData[0] = true;
            fuseData[1] = true;
        }
        if (fusePosData && fusionModeGPS <= 2 && posHealth)
        {
            fuseData[3] = true;
            fuseData[4] = true;
        }
        if (fuseHgtData && hgtHealth)
        {
            fuseData[5] = true;
        }
        // Limit range of states modified when on ground
        if(!onGround)
        {
            indexLimit = 20;
        }
        else
        {
            indexLimit = 12;
        }
        // Fuse measurements sequentially
        for (obsIndex=0; obsIndex<=5; obsIndex++)
        {
            if (fuseData[obsIndex])
            {
                stateIndex = 4 + obsIndex;
                // Calculate the measurement innovation, using states from a
                // different time coordinate if fusing height data
                if (obsIndex >= 0 && obsIndex <= 2)
                {
                    innovVelPos[obsIndex] = statesAtVelTime[stateIndex] - observation[obsIndex];
                }
                else if (obsIndex == 3 || obsIndex == 4)
                {
                    innovVelPos[obsIndex] = statesAtPosTime[stateIndex] - observation[obsIndex];
                }
                else if (obsIndex == 5)
                {
                    innovVelPos[obsIndex] = statesAtHgtTime[stateIndex] - observation[obsIndex];
                }
                // Calculate the Kalman Gain
                // Calculate innovation variances - also used for data logging
                varInnovVelPos[obsIndex] = P[stateIndex][stateIndex] + R_OBS[obsIndex];
                SK = 1.0/varInnovVelPos[obsIndex];
                for (uint8_t i= 0; i<=indexLimit; i++)
                {
                    Kfusion[i] = P[i][stateIndex]*SK;
                }
                // Calculate state corrections and re-normalise the quaternions
                for (uint8_t i = 0; i<=indexLimit; i++)
                {
                    states[i] = states[i] - Kfusion[i] * innovVelPos[obsIndex];
                }
                quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
                if (quatMag > 1e-12f) // divide by  0 protection
                {
                    for (uint8_t i = 0; i<=3; i++)
                    {
                        states[i] = states[i] / quatMag;
                    }
                }
                // Update the covariance - take advantage of direct observation of a
                // single state at index = stateIndex to reduce computations
                // Optimised implementation of standard equation P = (I - K*H)*P;
                for (uint8_t i= 0; i<=indexLimit; i++)
                {
                    for (uint8_t j= 0; j<=indexLimit; j++)
                    {
                        KHP[i][j] = Kfusion[i] * P[stateIndex][j];
                    }
                }
                for (uint8_t i= 0; i<=indexLimit; i++)
                {
                    for (uint8_t j= 0; j<=indexLimit; j++)
                    {
                        P[i][j] = P[i][j] - KHP[i][j];
                    }
                }
            }
        }
    }
}

void FuseMagnetometer()
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
    static uint8_t obsIndex;
    uint8_t indexLimit;
    float DCM[3][3] =
    {
        {1.0,0.0,0.0} ,
        {0.0,1.0,0.0} ,
        {0.0,0.0,1.0}
    };
    static float MagPred[3] = {0.0,0.0,0.0};
    static float R_MAG;
    static float SH_MAG[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    float H_MAG[21];
    float SK_MX[6];
    float SK_MY[5];
    float SK_MZ[6];

// Perform sequential fusion of Magnetometer measurements.
// This assumes that the errors in the different components are
// uncorrelated which is not true, however in the absence of covariance
// data fit is the only assumption we can make
// so we might as well take advantage of the computational efficiencies
// associated with sequential fusion
    if (useCompass && (fuseMagData || obsIndex == 1 || obsIndex == 2))
    {
        // Limit range of states modified when on ground
        if(!onGround)
        {
            indexLimit = 20;
        }
        else
        {
            indexLimit = 12;
        }

        // Sequential fusion of XYZ components to spread processing load across
        // three prediction time steps.

        // Calculate observation jacobians and Kalman gains
        if (fuseMagData)
        {
            // Copy required states to local variable names
            q0       = statesAtMagMeasTime[0];
            q1       = statesAtMagMeasTime[1];
            q2       = statesAtMagMeasTime[2];
            q3       = statesAtMagMeasTime[3];
            magN     = statesAtMagMeasTime[15];
            magE     = statesAtMagMeasTime[16];
            magD     = statesAtMagMeasTime[17];
            magXbias = statesAtMagMeasTime[18];
            magYbias = statesAtMagMeasTime[19];
            magZbias = statesAtMagMeasTime[20];

            // rotate predicted earth components into body axes and calculate
            // predicted measurments
            DCM[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
            DCM[0][1] = 2*(q1*q2 + q0*q3);
            DCM[0][2] = 2*(q1*q3-q0*q2);
            DCM[1][0] = 2*(q1*q2 - q0*q3);
            DCM[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
            DCM[1][2] = 2*(q2*q3 + q0*q1);
            DCM[2][0] = 2*(q1*q3 + q0*q2);
            DCM[2][1] = 2*(q2*q3 - q0*q1);
            DCM[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
            MagPred[0] = DCM[0][0]*magN + DCM[0][1]*magE  + DCM[0][2]*magD + magXbias;
            MagPred[1] = DCM[1][0]*magN + DCM[1][1]*magE  + DCM[1][2]*magD + magYbias;
            MagPred[2] = DCM[2][0]*magN + DCM[2][1]*magE  + DCM[2][2]*magD + magZbias;

            // scale magnetometer observation error with total angular rate
            R_MAG = 0.0025f + sq(0.05f*dAngIMU.length()/dtIMU);

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
            for (uint8_t i=0; i<=20; i++) H_MAG[i] = 0;
            H_MAG[0] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            H_MAG[1] = SH_MAG[0];
            H_MAG[2] = 2*magE*q1 - 2*magD*q0 - 2*magN*q2;
            H_MAG[3] = SH_MAG[2];
            H_MAG[15] = SH_MAG[5] - SH_MAG[4] - SH_MAG[3] + SH_MAG[6];
            H_MAG[16] = 2*q0*q3 + 2*q1*q2;
            H_MAG[17] = 2*q1*q3 - 2*q0*q2;
            H_MAG[18] = 1;

            // Calculate Kalman gain
            SK_MX[0] = 1/(P[18][18] + R_MAG + P[1][18]*SH_MAG[0] + P[3][18]*SH_MAG[2] - P[15][18]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) - (2*magD*q0 - 2*magE*q1 + 2*magN*q2)*(P[18][2] + P[1][2]*SH_MAG[0] + P[3][2]*SH_MAG[2] - P[15][2]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][2]*(2*q0*q3 + 2*q1*q2) - P[17][2]*(2*q0*q2 - 2*q1*q3) - P[2][2]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[18][0] + P[1][0]*SH_MAG[0] + P[3][0]*SH_MAG[2] - P[15][0]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][0]*(2*q0*q3 + 2*q1*q2) - P[17][0]*(2*q0*q2 - 2*q1*q3) - P[2][0]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[18][1] + P[1][1]*SH_MAG[0] + P[3][1]*SH_MAG[2] - P[15][1]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][1]*(2*q0*q3 + 2*q1*q2) - P[17][1]*(2*q0*q2 - 2*q1*q3) - P[2][1]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[2]*(P[18][3] + P[1][3]*SH_MAG[0] + P[3][3]*SH_MAG[2] - P[15][3]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][3]*(2*q0*q3 + 2*q1*q2) - P[17][3]*(2*q0*q2 - 2*q1*q3) - P[2][3]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6])*(P[18][15] + P[1][15]*SH_MAG[0] + P[3][15]*SH_MAG[2] - P[15][15]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][15]*(2*q0*q3 + 2*q1*q2) - P[17][15]*(2*q0*q2 - 2*q1*q3) - P[2][15]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][15]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[16][18]*(2*q0*q3 + 2*q1*q2) - P[17][18]*(2*q0*q2 - 2*q1*q3) - P[2][18]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + (2*q0*q3 + 2*q1*q2)*(P[18][16] + P[1][16]*SH_MAG[0] + P[3][16]*SH_MAG[2] - P[15][16]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][16]*(2*q0*q3 + 2*q1*q2) - P[17][16]*(2*q0*q2 - 2*q1*q3) - P[2][16]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][16]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (2*q0*q2 - 2*q1*q3)*(P[18][17] + P[1][17]*SH_MAG[0] + P[3][17]*SH_MAG[2] - P[15][17]*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P[16][17]*(2*q0*q3 + 2*q1*q2) - P[17][17]*(2*q0*q2 - 2*q1*q3) - P[2][17]*(2*magD*q0 - 2*magE*q1 + 2*magN*q2) + P[0][17]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[0][18]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MX[1] = SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6];
            SK_MX[2] = 2*magD*q0 - 2*magE*q1 + 2*magN*q2;
            SK_MX[3] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MX[4] = 2*q0*q2 - 2*q1*q3;
            SK_MX[5] = 2*q0*q3 + 2*q1*q2;
            Kfusion[0] = SK_MX[0]*(P[0][18] + P[0][1]*SH_MAG[0] + P[0][3]*SH_MAG[2] + P[0][0]*SK_MX[3] - P[0][2]*SK_MX[2] - P[0][15]*SK_MX[1] + P[0][16]*SK_MX[5] - P[0][17]*SK_MX[4]);
            Kfusion[1] = SK_MX[0]*(P[1][18] + P[1][1]*SH_MAG[0] + P[1][3]*SH_MAG[2] + P[1][0]*SK_MX[3] - P[1][2]*SK_MX[2] - P[1][15]*SK_MX[1] + P[1][16]*SK_MX[5] - P[1][17]*SK_MX[4]);
            Kfusion[2] = SK_MX[0]*(P[2][18] + P[2][1]*SH_MAG[0] + P[2][3]*SH_MAG[2] + P[2][0]*SK_MX[3] - P[2][2]*SK_MX[2] - P[2][15]*SK_MX[1] + P[2][16]*SK_MX[5] - P[2][17]*SK_MX[4]);
            Kfusion[3] = SK_MX[0]*(P[3][18] + P[3][1]*SH_MAG[0] + P[3][3]*SH_MAG[2] + P[3][0]*SK_MX[3] - P[3][2]*SK_MX[2] - P[3][15]*SK_MX[1] + P[3][16]*SK_MX[5] - P[3][17]*SK_MX[4]);
            Kfusion[4] = SK_MX[0]*(P[4][18] + P[4][1]*SH_MAG[0] + P[4][3]*SH_MAG[2] + P[4][0]*SK_MX[3] - P[4][2]*SK_MX[2] - P[4][15]*SK_MX[1] + P[4][16]*SK_MX[5] - P[4][17]*SK_MX[4]);
            Kfusion[5] = SK_MX[0]*(P[5][18] + P[5][1]*SH_MAG[0] + P[5][3]*SH_MAG[2] + P[5][0]*SK_MX[3] - P[5][2]*SK_MX[2] - P[5][15]*SK_MX[1] + P[5][16]*SK_MX[5] - P[5][17]*SK_MX[4]);
            Kfusion[6] = SK_MX[0]*(P[6][18] + P[6][1]*SH_MAG[0] + P[6][3]*SH_MAG[2] + P[6][0]*SK_MX[3] - P[6][2]*SK_MX[2] - P[6][15]*SK_MX[1] + P[6][16]*SK_MX[5] - P[6][17]*SK_MX[4]);
            Kfusion[7] = SK_MX[0]*(P[7][18] + P[7][1]*SH_MAG[0] + P[7][3]*SH_MAG[2] + P[7][0]*SK_MX[3] - P[7][2]*SK_MX[2] - P[7][15]*SK_MX[1] + P[7][16]*SK_MX[5] - P[7][17]*SK_MX[4]);
            Kfusion[8] = SK_MX[0]*(P[8][18] + P[8][1]*SH_MAG[0] + P[8][3]*SH_MAG[2] + P[8][0]*SK_MX[3] - P[8][2]*SK_MX[2] - P[8][15]*SK_MX[1] + P[8][16]*SK_MX[5] - P[8][17]*SK_MX[4]);
            Kfusion[9] = SK_MX[0]*(P[9][18] + P[9][1]*SH_MAG[0] + P[9][3]*SH_MAG[2] + P[9][0]*SK_MX[3] - P[9][2]*SK_MX[2] - P[9][15]*SK_MX[1] + P[9][16]*SK_MX[5] - P[9][17]*SK_MX[4]);
            Kfusion[10] = SK_MX[0]*(P[10][18] + P[10][1]*SH_MAG[0] + P[10][3]*SH_MAG[2] + P[10][0]*SK_MX[3] - P[10][2]*SK_MX[2] - P[10][15]*SK_MX[1] + P[10][16]*SK_MX[5] - P[10][17]*SK_MX[4]);
            Kfusion[11] = SK_MX[0]*(P[11][18] + P[11][1]*SH_MAG[0] + P[11][3]*SH_MAG[2] + P[11][0]*SK_MX[3] - P[11][2]*SK_MX[2] - P[11][15]*SK_MX[1] + P[11][16]*SK_MX[5] - P[11][17]*SK_MX[4]);
            Kfusion[12] = SK_MX[0]*(P[12][18] + P[12][1]*SH_MAG[0] + P[12][3]*SH_MAG[2] + P[12][0]*SK_MX[3] - P[12][2]*SK_MX[2] - P[12][15]*SK_MX[1] + P[12][16]*SK_MX[5] - P[12][17]*SK_MX[4]);
            Kfusion[13] = SK_MX[0]*(P[13][18] + P[13][1]*SH_MAG[0] + P[13][3]*SH_MAG[2] + P[13][0]*SK_MX[3] - P[13][2]*SK_MX[2] - P[13][15]*SK_MX[1] + P[13][16]*SK_MX[5] - P[13][17]*SK_MX[4]);
            Kfusion[14] = SK_MX[0]*(P[14][18] + P[14][1]*SH_MAG[0] + P[14][3]*SH_MAG[2] + P[14][0]*SK_MX[3] - P[14][2]*SK_MX[2] - P[14][15]*SK_MX[1] + P[14][16]*SK_MX[5] - P[14][17]*SK_MX[4]);
            Kfusion[15] = SK_MX[0]*(P[15][18] + P[15][1]*SH_MAG[0] + P[15][3]*SH_MAG[2] + P[15][0]*SK_MX[3] - P[15][2]*SK_MX[2] - P[15][15]*SK_MX[1] + P[15][16]*SK_MX[5] - P[15][17]*SK_MX[4]);
            Kfusion[16] = SK_MX[0]*(P[16][18] + P[16][1]*SH_MAG[0] + P[16][3]*SH_MAG[2] + P[16][0]*SK_MX[3] - P[16][2]*SK_MX[2] - P[16][15]*SK_MX[1] + P[16][16]*SK_MX[5] - P[16][17]*SK_MX[4]);
            Kfusion[17] = SK_MX[0]*(P[17][18] + P[17][1]*SH_MAG[0] + P[17][3]*SH_MAG[2] + P[17][0]*SK_MX[3] - P[17][2]*SK_MX[2] - P[17][15]*SK_MX[1] + P[17][16]*SK_MX[5] - P[17][17]*SK_MX[4]);
            Kfusion[18] = SK_MX[0]*(P[18][18] + P[18][1]*SH_MAG[0] + P[18][3]*SH_MAG[2] + P[18][0]*SK_MX[3] - P[18][2]*SK_MX[2] - P[18][15]*SK_MX[1] + P[18][16]*SK_MX[5] - P[18][17]*SK_MX[4]);
            Kfusion[19] = SK_MX[0]*(P[19][18] + P[19][1]*SH_MAG[0] + P[19][3]*SH_MAG[2] + P[19][0]*SK_MX[3] - P[19][2]*SK_MX[2] - P[19][15]*SK_MX[1] + P[19][16]*SK_MX[5] - P[19][17]*SK_MX[4]);
            Kfusion[20] = SK_MX[0]*(P[20][18] + P[20][1]*SH_MAG[0] + P[20][3]*SH_MAG[2] + P[20][0]*SK_MX[3] - P[20][2]*SK_MX[2] - P[20][15]*SK_MX[1] + P[20][16]*SK_MX[5] - P[20][17]*SK_MX[4]);
            varInnovMag[0] = 1.0f/SK_MX[0];
            innovMag[0] = MagPred[0] - magData.x;

            // reset the observation index to 0 (we start by fusing the X
            // measurement)
            obsIndex = 0;
        }
        else if (obsIndex == 1) // we are now fusing the Y measurement
        {
            // Calculate observation jacobians
            for (uint8_t i=0; i<=20; i++) H_MAG[i] = 0;
            H_MAG[0] = SH_MAG[2];
            H_MAG[1] = SH_MAG[1];
            H_MAG[2] = SH_MAG[0];
            H_MAG[3] = 2*magD*q2 - SH_MAG[8] - SH_MAG[7];
            H_MAG[15] = 2*q1*q2 - 2*q0*q3;
            H_MAG[16] = SH_MAG[4] - SH_MAG[3] - SH_MAG[5] + SH_MAG[6];
            H_MAG[17] = 2*q0*q1 + 2*q2*q3;
            H_MAG[19] = 1;

            // Calculate Kalman gain
            SK_MY[0] = 1/(P[19][19] + R_MAG + P[0][19]*SH_MAG[2] + P[1][19]*SH_MAG[1] + P[2][19]*SH_MAG[0] - P[16][19]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - (2*q0*q3 - 2*q1*q2)*(P[19][15] + P[0][15]*SH_MAG[2] + P[1][15]*SH_MAG[1] + P[2][15]*SH_MAG[0] - P[16][15]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][15]*(2*q0*q3 - 2*q1*q2) + P[17][15]*(2*q0*q1 + 2*q2*q3) - P[3][15]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (2*q0*q1 + 2*q2*q3)*(P[19][17] + P[0][17]*SH_MAG[2] + P[1][17]*SH_MAG[1] + P[2][17]*SH_MAG[0] - P[16][17]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][17]*(2*q0*q3 - 2*q1*q2) + P[17][17]*(2*q0*q1 + 2*q2*q3) - P[3][17]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[19][3] + P[0][3]*SH_MAG[2] + P[1][3]*SH_MAG[1] + P[2][3]*SH_MAG[0] - P[16][3]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][3]*(2*q0*q3 - 2*q1*q2) + P[17][3]*(2*q0*q1 + 2*q2*q3) - P[3][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - P[15][19]*(2*q0*q3 - 2*q1*q2) + P[17][19]*(2*q0*q1 + 2*q2*q3) + SH_MAG[2]*(P[19][0] + P[0][0]*SH_MAG[2] + P[1][0]*SH_MAG[1] + P[2][0]*SH_MAG[0] - P[16][0]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][0]*(2*q0*q3 - 2*q1*q2) + P[17][0]*(2*q0*q1 + 2*q2*q3) - P[3][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[1]*(P[19][1] + P[0][1]*SH_MAG[2] + P[1][1]*SH_MAG[1] + P[2][1]*SH_MAG[0] - P[16][1]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][1]*(2*q0*q3 - 2*q1*q2) + P[17][1]*(2*q0*q1 + 2*q2*q3) - P[3][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[19][2] + P[0][2]*SH_MAG[2] + P[1][2]*SH_MAG[1] + P[2][2]*SH_MAG[0] - P[16][2]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][2]*(2*q0*q3 - 2*q1*q2) + P[17][2]*(2*q0*q1 + 2*q2*q3) - P[3][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6])*(P[19][16] + P[0][16]*SH_MAG[2] + P[1][16]*SH_MAG[1] + P[2][16]*SH_MAG[0] - P[16][16]*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P[15][16]*(2*q0*q3 - 2*q1*q2) + P[17][16]*(2*q0*q1 + 2*q2*q3) - P[3][16]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - P[3][19]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MY[1] = SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6];
            SK_MY[2] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MY[3] = 2*q0*q3 - 2*q1*q2;
            SK_MY[4] = 2*q0*q1 + 2*q2*q3;
            Kfusion[0] = SK_MY[0]*(P[0][19] + P[0][0]*SH_MAG[2] + P[0][1]*SH_MAG[1] + P[0][2]*SH_MAG[0] - P[0][3]*SK_MY[2] - P[0][16]*SK_MY[1] - P[0][15]*SK_MY[3] + P[0][17]*SK_MY[4]);
            Kfusion[1] = SK_MY[0]*(P[1][19] + P[1][0]*SH_MAG[2] + P[1][1]*SH_MAG[1] + P[1][2]*SH_MAG[0] - P[1][3]*SK_MY[2] - P[1][16]*SK_MY[1] - P[1][15]*SK_MY[3] + P[1][17]*SK_MY[4]);
            Kfusion[2] = SK_MY[0]*(P[2][19] + P[2][0]*SH_MAG[2] + P[2][1]*SH_MAG[1] + P[2][2]*SH_MAG[0] - P[2][3]*SK_MY[2] - P[2][16]*SK_MY[1] - P[2][15]*SK_MY[3] + P[2][17]*SK_MY[4]);
            Kfusion[3] = SK_MY[0]*(P[3][19] + P[3][0]*SH_MAG[2] + P[3][1]*SH_MAG[1] + P[3][2]*SH_MAG[0] - P[3][3]*SK_MY[2] - P[3][16]*SK_MY[1] - P[3][15]*SK_MY[3] + P[3][17]*SK_MY[4]);
            Kfusion[4] = SK_MY[0]*(P[4][19] + P[4][0]*SH_MAG[2] + P[4][1]*SH_MAG[1] + P[4][2]*SH_MAG[0] - P[4][3]*SK_MY[2] - P[4][16]*SK_MY[1] - P[4][15]*SK_MY[3] + P[4][17]*SK_MY[4]);
            Kfusion[5] = SK_MY[0]*(P[5][19] + P[5][0]*SH_MAG[2] + P[5][1]*SH_MAG[1] + P[5][2]*SH_MAG[0] - P[5][3]*SK_MY[2] - P[5][16]*SK_MY[1] - P[5][15]*SK_MY[3] + P[5][17]*SK_MY[4]);
            Kfusion[6] = SK_MY[0]*(P[6][19] + P[6][0]*SH_MAG[2] + P[6][1]*SH_MAG[1] + P[6][2]*SH_MAG[0] - P[6][3]*SK_MY[2] - P[6][16]*SK_MY[1] - P[6][15]*SK_MY[3] + P[6][17]*SK_MY[4]);
            Kfusion[7] = SK_MY[0]*(P[7][19] + P[7][0]*SH_MAG[2] + P[7][1]*SH_MAG[1] + P[7][2]*SH_MAG[0] - P[7][3]*SK_MY[2] - P[7][16]*SK_MY[1] - P[7][15]*SK_MY[3] + P[7][17]*SK_MY[4]);
            Kfusion[8] = SK_MY[0]*(P[8][19] + P[8][0]*SH_MAG[2] + P[8][1]*SH_MAG[1] + P[8][2]*SH_MAG[0] - P[8][3]*SK_MY[2] - P[8][16]*SK_MY[1] - P[8][15]*SK_MY[3] + P[8][17]*SK_MY[4]);
            Kfusion[9] = SK_MY[0]*(P[9][19] + P[9][0]*SH_MAG[2] + P[9][1]*SH_MAG[1] + P[9][2]*SH_MAG[0] - P[9][3]*SK_MY[2] - P[9][16]*SK_MY[1] - P[9][15]*SK_MY[3] + P[9][17]*SK_MY[4]);
            Kfusion[10] = SK_MY[0]*(P[10][19] + P[10][0]*SH_MAG[2] + P[10][1]*SH_MAG[1] + P[10][2]*SH_MAG[0] - P[10][3]*SK_MY[2] - P[10][16]*SK_MY[1] - P[10][15]*SK_MY[3] + P[10][17]*SK_MY[4]);
            Kfusion[11] = SK_MY[0]*(P[11][19] + P[11][0]*SH_MAG[2] + P[11][1]*SH_MAG[1] + P[11][2]*SH_MAG[0] - P[11][3]*SK_MY[2] - P[11][16]*SK_MY[1] - P[11][15]*SK_MY[3] + P[11][17]*SK_MY[4]);
            Kfusion[12] = SK_MY[0]*(P[12][19] + P[12][0]*SH_MAG[2] + P[12][1]*SH_MAG[1] + P[12][2]*SH_MAG[0] - P[12][3]*SK_MY[2] - P[12][16]*SK_MY[1] - P[12][15]*SK_MY[3] + P[12][17]*SK_MY[4]);
            Kfusion[13] = SK_MY[0]*(P[13][19] + P[13][0]*SH_MAG[2] + P[13][1]*SH_MAG[1] + P[13][2]*SH_MAG[0] - P[13][3]*SK_MY[2] - P[13][16]*SK_MY[1] - P[13][15]*SK_MY[3] + P[13][17]*SK_MY[4]);
            Kfusion[14] = SK_MY[0]*(P[14][19] + P[14][0]*SH_MAG[2] + P[14][1]*SH_MAG[1] + P[14][2]*SH_MAG[0] - P[14][3]*SK_MY[2] - P[14][16]*SK_MY[1] - P[14][15]*SK_MY[3] + P[14][17]*SK_MY[4]);
            Kfusion[15] = SK_MY[0]*(P[15][19] + P[15][0]*SH_MAG[2] + P[15][1]*SH_MAG[1] + P[15][2]*SH_MAG[0] - P[15][3]*SK_MY[2] - P[15][16]*SK_MY[1] - P[15][15]*SK_MY[3] + P[15][17]*SK_MY[4]);
            Kfusion[16] = SK_MY[0]*(P[16][19] + P[16][0]*SH_MAG[2] + P[16][1]*SH_MAG[1] + P[16][2]*SH_MAG[0] - P[16][3]*SK_MY[2] - P[16][16]*SK_MY[1] - P[16][15]*SK_MY[3] + P[16][17]*SK_MY[4]);
            Kfusion[17] = SK_MY[0]*(P[17][19] + P[17][0]*SH_MAG[2] + P[17][1]*SH_MAG[1] + P[17][2]*SH_MAG[0] - P[17][3]*SK_MY[2] - P[17][16]*SK_MY[1] - P[17][15]*SK_MY[3] + P[17][17]*SK_MY[4]);
            Kfusion[18] = SK_MY[0]*(P[18][19] + P[18][0]*SH_MAG[2] + P[18][1]*SH_MAG[1] + P[18][2]*SH_MAG[0] - P[18][3]*SK_MY[2] - P[18][16]*SK_MY[1] - P[18][15]*SK_MY[3] + P[18][17]*SK_MY[4]);
            Kfusion[19] = SK_MY[0]*(P[19][19] + P[19][0]*SH_MAG[2] + P[19][1]*SH_MAG[1] + P[19][2]*SH_MAG[0] - P[19][3]*SK_MY[2] - P[19][16]*SK_MY[1] - P[19][15]*SK_MY[3] + P[19][17]*SK_MY[4]);
            Kfusion[20] = SK_MY[0]*(P[20][19] + P[20][0]*SH_MAG[2] + P[20][1]*SH_MAG[1] + P[20][2]*SH_MAG[0] - P[20][3]*SK_MY[2] - P[20][16]*SK_MY[1] - P[20][15]*SK_MY[3] + P[20][17]*SK_MY[4]);
            varInnovMag[1] = 1.0f/SK_MY[0];
            innovMag[1] = MagPred[1] - magData.y;
        }
        else if (obsIndex == 2) // we are now fusing the Z measurement
        {
            // Calculate observation jacobians
            for (uint8_t i=0; i<=20; i++) H_MAG[i] = 0;
            H_MAG[0] = SH_MAG[1];
            H_MAG[1] = 2*magN*q3 - 2*magE*q0 - 2*magD*q1;
            H_MAG[2] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            H_MAG[3] = SH_MAG[0];
            H_MAG[15] = 2*q0*q2 + 2*q1*q3;
            H_MAG[16] = 2*q2*q3 - 2*q0*q1;
            H_MAG[17] = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
            H_MAG[20] = 1;

            // Calculate Kalman gain
            SK_MZ[0] = 1/(P[20][20] + R_MAG + P[0][20]*SH_MAG[1] + P[3][20]*SH_MAG[0] + P[17][20]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) - (2*magD*q1 + 2*magE*q0 - 2*magN*q3)*(P[20][1] + P[0][1]*SH_MAG[1] + P[3][1]*SH_MAG[0] + P[17][1]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][1]*(2*q0*q2 + 2*q1*q3) - P[16][1]*(2*q0*q1 - 2*q2*q3) - P[1][1]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][1]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2*magD*q2)*(P[20][2] + P[0][2]*SH_MAG[1] + P[3][2]*SH_MAG[0] + P[17][2]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][2]*(2*q0*q2 + 2*q1*q3) - P[16][2]*(2*q0*q1 - 2*q2*q3) - P[1][2]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][2]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[1]*(P[20][0] + P[0][0]*SH_MAG[1] + P[3][0]*SH_MAG[0] + P[17][0]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][0]*(2*q0*q2 + 2*q1*q3) - P[16][0]*(2*q0*q1 - 2*q2*q3) - P[1][0]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][0]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + SH_MAG[0]*(P[20][3] + P[0][3]*SH_MAG[1] + P[3][3]*SH_MAG[0] + P[17][3]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][3]*(2*q0*q2 + 2*q1*q3) - P[16][3]*(2*q0*q1 - 2*q2*q3) - P[1][3]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][3]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + (SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6])*(P[20][17] + P[0][17]*SH_MAG[1] + P[3][17]*SH_MAG[0] + P[17][17]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][17]*(2*q0*q2 + 2*q1*q3) - P[16][17]*(2*q0*q1 - 2*q2*q3) - P[1][17]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][17]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[15][20]*(2*q0*q2 + 2*q1*q3) - P[16][20]*(2*q0*q1 - 2*q2*q3) - P[1][20]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + (2*q0*q2 + 2*q1*q3)*(P[20][15] + P[0][15]*SH_MAG[1] + P[3][15]*SH_MAG[0] + P[17][15]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][15]*(2*q0*q2 + 2*q1*q3) - P[16][15]*(2*q0*q1 - 2*q2*q3) - P[1][15]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][15]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) - (2*q0*q1 - 2*q2*q3)*(P[20][16] + P[0][16]*SH_MAG[1] + P[3][16]*SH_MAG[0] + P[17][16]*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P[15][16]*(2*q0*q2 + 2*q1*q3) - P[16][16]*(2*q0*q1 - 2*q2*q3) - P[1][16]*(2*magD*q1 + 2*magE*q0 - 2*magN*q3) + P[2][16]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2)) + P[2][20]*(SH_MAG[7] + SH_MAG[8] - 2*magD*q2));
            SK_MZ[1] = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
            SK_MZ[2] = 2*magD*q1 + 2*magE*q0 - 2*magN*q3;
            SK_MZ[3] = SH_MAG[7] + SH_MAG[8] - 2*magD*q2;
            SK_MZ[4] = 2*q0*q1 - 2*q2*q3;
            SK_MZ[5] = 2*q0*q2 + 2*q1*q3;
            Kfusion[0] = SK_MZ[0]*(P[0][20] + P[0][0]*SH_MAG[1] + P[0][3]*SH_MAG[0] - P[0][1]*SK_MZ[2] + P[0][2]*SK_MZ[3] + P[0][17]*SK_MZ[1] + P[0][15]*SK_MZ[5] - P[0][16]*SK_MZ[4]);
            Kfusion[1] = SK_MZ[0]*(P[1][20] + P[1][0]*SH_MAG[1] + P[1][3]*SH_MAG[0] - P[1][1]*SK_MZ[2] + P[1][2]*SK_MZ[3] + P[1][17]*SK_MZ[1] + P[1][15]*SK_MZ[5] - P[1][16]*SK_MZ[4]);
            Kfusion[2] = SK_MZ[0]*(P[2][20] + P[2][0]*SH_MAG[1] + P[2][3]*SH_MAG[0] - P[2][1]*SK_MZ[2] + P[2][2]*SK_MZ[3] + P[2][17]*SK_MZ[1] + P[2][15]*SK_MZ[5] - P[2][16]*SK_MZ[4]);
            Kfusion[3] = SK_MZ[0]*(P[3][20] + P[3][0]*SH_MAG[1] + P[3][3]*SH_MAG[0] - P[3][1]*SK_MZ[2] + P[3][2]*SK_MZ[3] + P[3][17]*SK_MZ[1] + P[3][15]*SK_MZ[5] - P[3][16]*SK_MZ[4]);
            Kfusion[4] = SK_MZ[0]*(P[4][20] + P[4][0]*SH_MAG[1] + P[4][3]*SH_MAG[0] - P[4][1]*SK_MZ[2] + P[4][2]*SK_MZ[3] + P[4][17]*SK_MZ[1] + P[4][15]*SK_MZ[5] - P[4][16]*SK_MZ[4]);
            Kfusion[5] = SK_MZ[0]*(P[5][20] + P[5][0]*SH_MAG[1] + P[5][3]*SH_MAG[0] - P[5][1]*SK_MZ[2] + P[5][2]*SK_MZ[3] + P[5][17]*SK_MZ[1] + P[5][15]*SK_MZ[5] - P[5][16]*SK_MZ[4]);
            Kfusion[6] = SK_MZ[0]*(P[6][20] + P[6][0]*SH_MAG[1] + P[6][3]*SH_MAG[0] - P[6][1]*SK_MZ[2] + P[6][2]*SK_MZ[3] + P[6][17]*SK_MZ[1] + P[6][15]*SK_MZ[5] - P[6][16]*SK_MZ[4]);
            Kfusion[7] = SK_MZ[0]*(P[7][20] + P[7][0]*SH_MAG[1] + P[7][3]*SH_MAG[0] - P[7][1]*SK_MZ[2] + P[7][2]*SK_MZ[3] + P[7][17]*SK_MZ[1] + P[7][15]*SK_MZ[5] - P[7][16]*SK_MZ[4]);
            Kfusion[8] = SK_MZ[0]*(P[8][20] + P[8][0]*SH_MAG[1] + P[8][3]*SH_MAG[0] - P[8][1]*SK_MZ[2] + P[8][2]*SK_MZ[3] + P[8][17]*SK_MZ[1] + P[8][15]*SK_MZ[5] - P[8][16]*SK_MZ[4]);
            Kfusion[9] = SK_MZ[0]*(P[9][20] + P[9][0]*SH_MAG[1] + P[9][3]*SH_MAG[0] - P[9][1]*SK_MZ[2] + P[9][2]*SK_MZ[3] + P[9][17]*SK_MZ[1] + P[9][15]*SK_MZ[5] - P[9][16]*SK_MZ[4]);
            Kfusion[10] = SK_MZ[0]*(P[10][20] + P[10][0]*SH_MAG[1] + P[10][3]*SH_MAG[0] - P[10][1]*SK_MZ[2] + P[10][2]*SK_MZ[3] + P[10][17]*SK_MZ[1] + P[10][15]*SK_MZ[5] - P[10][16]*SK_MZ[4]);
            Kfusion[11] = SK_MZ[0]*(P[11][20] + P[11][0]*SH_MAG[1] + P[11][3]*SH_MAG[0] - P[11][1]*SK_MZ[2] + P[11][2]*SK_MZ[3] + P[11][17]*SK_MZ[1] + P[11][15]*SK_MZ[5] - P[11][16]*SK_MZ[4]);
            Kfusion[12] = SK_MZ[0]*(P[12][20] + P[12][0]*SH_MAG[1] + P[12][3]*SH_MAG[0] - P[12][1]*SK_MZ[2] + P[12][2]*SK_MZ[3] + P[12][17]*SK_MZ[1] + P[12][15]*SK_MZ[5] - P[12][16]*SK_MZ[4]);
            Kfusion[13] = SK_MZ[0]*(P[13][20] + P[13][0]*SH_MAG[1] + P[13][3]*SH_MAG[0] - P[13][1]*SK_MZ[2] + P[13][2]*SK_MZ[3] + P[13][17]*SK_MZ[1] + P[13][15]*SK_MZ[5] - P[13][16]*SK_MZ[4]);
            Kfusion[14] = SK_MZ[0]*(P[14][20] + P[14][0]*SH_MAG[1] + P[14][3]*SH_MAG[0] - P[14][1]*SK_MZ[2] + P[14][2]*SK_MZ[3] + P[14][17]*SK_MZ[1] + P[14][15]*SK_MZ[5] - P[14][16]*SK_MZ[4]);
            Kfusion[15] = SK_MZ[0]*(P[15][20] + P[15][0]*SH_MAG[1] + P[15][3]*SH_MAG[0] - P[15][1]*SK_MZ[2] + P[15][2]*SK_MZ[3] + P[15][17]*SK_MZ[1] + P[15][15]*SK_MZ[5] - P[15][16]*SK_MZ[4]);
            Kfusion[16] = SK_MZ[0]*(P[16][20] + P[16][0]*SH_MAG[1] + P[16][3]*SH_MAG[0] - P[16][1]*SK_MZ[2] + P[16][2]*SK_MZ[3] + P[16][17]*SK_MZ[1] + P[16][15]*SK_MZ[5] - P[16][16]*SK_MZ[4]);
            Kfusion[17] = SK_MZ[0]*(P[17][20] + P[17][0]*SH_MAG[1] + P[17][3]*SH_MAG[0] - P[17][1]*SK_MZ[2] + P[17][2]*SK_MZ[3] + P[17][17]*SK_MZ[1] + P[17][15]*SK_MZ[5] - P[17][16]*SK_MZ[4]);
            Kfusion[18] = SK_MZ[0]*(P[18][20] + P[18][0]*SH_MAG[1] + P[18][3]*SH_MAG[0] - P[18][1]*SK_MZ[2] + P[18][2]*SK_MZ[3] + P[18][17]*SK_MZ[1] + P[18][15]*SK_MZ[5] - P[18][16]*SK_MZ[4]);
            Kfusion[19] = SK_MZ[0]*(P[19][20] + P[19][0]*SH_MAG[1] + P[19][3]*SH_MAG[0] - P[19][1]*SK_MZ[2] + P[19][2]*SK_MZ[3] + P[19][17]*SK_MZ[1] + P[19][15]*SK_MZ[5] - P[19][16]*SK_MZ[4]);
            Kfusion[20] = SK_MZ[0]*(P[20][20] + P[20][0]*SH_MAG[1] + P[20][3]*SH_MAG[0] - P[20][1]*SK_MZ[2] + P[20][2]*SK_MZ[3] + P[20][17]*SK_MZ[1] + P[20][15]*SK_MZ[5] - P[20][16]*SK_MZ[4]);
            varInnovMag[2] = 1.0f/SK_MZ[0];
            innovMag[2] = MagPred[2] - magData.z;

        }

        // Check the innovation for consistency and don't fuse if > 5Sigma
        if ((innovMag[obsIndex]*innovMag[obsIndex]/varInnovMag[obsIndex]) < 25.0)
        {
            // correct the state vector
            for (uint8_t j= 0; j<=indexLimit; j++)
            {
                states[j] = states[j] - Kfusion[j] * innovMag[obsIndex];
            }
            // normalise the quaternion states
            float quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
            if (quatMag > 1e-12)
            {
                for (uint8_t j= 0; j<=3; j++)
                {
                    float quatMagInv = 1.0f/quatMag;
                    states[j] = states[j] * quatMagInv;
                }
            }
            // correct the covariance P = (I - K*H)*P
            // take advantage of the empty columns in KH to reduce the
            // number of operations
            for (uint8_t i = 0; i<=indexLimit; i++)
            {
                for (uint8_t j = 0; j<=3; j++)
                {
                    KH[i][j] = Kfusion[i] * H_MAG[j];
                }
                for (uint8_t j = 4; j<=17; j++) KH[i][j] = 0.0f;
                if (!onGround)
                {
                    for (uint8_t j = 15; j<=20; j++)
                    {
                        KH[i][j] = Kfusion[i] * H_MAG[j];
                    }
                }
                else
                {
                    for (uint8_t j = 15; j<=20; j++)
                    {
                        KH[i][j] = 0.0f;
                    }
                }
            }
            for (uint8_t i = 0; i<=indexLimit; i++)
            {
                for (uint8_t j = 0; j<=indexLimit; j++)
                {
                    KHP[i][j] = 0.0f;
                    for (uint8_t k = 0; k<=3; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][k] * P[k][j];
                    }
                    if (!onGround)
                    {
                        for (uint8_t k = 15; k<=20; k++)
                        {
                            KHP[i][j] = KHP[i][j] + KH[i][k] * P[k][j];
                        }
                    }
                }
            }
        }
        for (uint8_t i = 0; i<=indexLimit; i++)
        {
            for (uint8_t j = 0; j<=indexLimit; j++)
            {
                P[i][j] = P[i][j] - KHP[i][j];
            }
        }
    }
    obsIndex = obsIndex + 1;
}

void FuseAirspeed()
{
    float vn;
    float ve;
    float vd;
    float vwn;
    float vwe;
    const float R_TAS = 2.0f;
    float SH_TAS[3];
    float SK_TAS;
    float H_TAS[21];
    float Kfusion[21];
    float VtasPred;
    float quatMag;

    // Copy required states to local variable names
    vn = statesAtVtasMeasTime[4];
    ve = statesAtVtasMeasTime[5];
    vd = statesAtVtasMeasTime[6];
    vwn = statesAtVtasMeasTime[13];
    vwe = statesAtVtasMeasTime[14];

    // Need to check that it is flying before fusing airspeed data
    // Calculate the predicted airspeed
    VtasPred = sqrt((ve - vwe)*(ve - vwe) + (vn - vwn)*(vn - vwn) + vd*vd);
    // Perform fusion of True Airspeed measurement
    if (useAirspeed && fuseVtasData && (VtasPred > 1.0) && (VtasMeas > 8.0))
    {
        // Calculate observation jacobians
        SH_TAS[0] = 1/(sqrt(sq(ve - vwe) + sq(vn - vwn) + sq(vd)));
        SH_TAS[1] = (SH_TAS[0]*(2*ve - 2*vwe))/2;
        SH_TAS[2] = (SH_TAS[0]*(2*vn - 2*vwn))/2;
        for (uint8_t i=0; i<=20; i++) H_TAS[i] = 0.0f;
        H_TAS[4] = SH_TAS[2];
        H_TAS[5] = SH_TAS[1];
        H_TAS[6] = vd*SH_TAS[0];
        H_TAS[13] = -SH_TAS[2];
        H_TAS[14] = -SH_TAS[1];

        // Calculate Kalman gains
        SK_TAS = 1.0f/(R_TAS + SH_TAS[2]*(P[4][4]*SH_TAS[2] + P[5][4]*SH_TAS[1] - P[13][4]*SH_TAS[2] - P[14][4]*SH_TAS[1] + P[6][4]*vd*SH_TAS[0]) + SH_TAS[1]*(P[4][5]*SH_TAS[2] + P[5][5]*SH_TAS[1] - P[13][5]*SH_TAS[2] - P[14][5]*SH_TAS[1] + P[6][5]*vd*SH_TAS[0]) - SH_TAS[2]*(P[4][13]*SH_TAS[2] + P[5][13]*SH_TAS[1] - P[13][13]*SH_TAS[2] - P[14][13]*SH_TAS[1] + P[6][13]*vd*SH_TAS[0]) - SH_TAS[1]*(P[4][14]*SH_TAS[2] + P[5][14]*SH_TAS[1] - P[13][14]*SH_TAS[2] - P[14][14]*SH_TAS[1] + P[6][14]*vd*SH_TAS[0]) + vd*SH_TAS[0]*(P[4][6]*SH_TAS[2] + P[5][6]*SH_TAS[1] - P[13][6]*SH_TAS[2] - P[14][6]*SH_TAS[1] + P[6][6]*vd*SH_TAS[0]));
        Kfusion[0] = SK_TAS*(P[0][4]*SH_TAS[2] - P[0][13]*SH_TAS[2] + P[0][5]*SH_TAS[1] - P[0][14]*SH_TAS[1] + P[0][6]*vd*SH_TAS[0]);
        Kfusion[1] = SK_TAS*(P[1][4]*SH_TAS[2] - P[1][13]*SH_TAS[2] + P[1][5]*SH_TAS[1] - P[1][14]*SH_TAS[1] + P[1][6]*vd*SH_TAS[0]);
        Kfusion[2] = SK_TAS*(P[2][4]*SH_TAS[2] - P[2][13]*SH_TAS[2] + P[2][5]*SH_TAS[1] - P[2][14]*SH_TAS[1] + P[2][6]*vd*SH_TAS[0]);
        Kfusion[3] = SK_TAS*(P[3][4]*SH_TAS[2] - P[3][13]*SH_TAS[2] + P[3][5]*SH_TAS[1] - P[3][14]*SH_TAS[1] + P[3][6]*vd*SH_TAS[0]);
        Kfusion[4] = SK_TAS*(P[4][4]*SH_TAS[2] - P[4][13]*SH_TAS[2] + P[4][5]*SH_TAS[1] - P[4][14]*SH_TAS[1] + P[4][6]*vd*SH_TAS[0]);
        Kfusion[5] = SK_TAS*(P[5][4]*SH_TAS[2] - P[5][13]*SH_TAS[2] + P[5][5]*SH_TAS[1] - P[5][14]*SH_TAS[1] + P[5][6]*vd*SH_TAS[0]);
        Kfusion[6] = SK_TAS*(P[6][4]*SH_TAS[2] - P[6][13]*SH_TAS[2] + P[6][5]*SH_TAS[1] - P[6][14]*SH_TAS[1] + P[6][6]*vd*SH_TAS[0]);
        Kfusion[7] = SK_TAS*(P[7][4]*SH_TAS[2] - P[7][13]*SH_TAS[2] + P[7][5]*SH_TAS[1] - P[7][14]*SH_TAS[1] + P[7][6]*vd*SH_TAS[0]);
        Kfusion[8] = SK_TAS*(P[8][4]*SH_TAS[2] - P[8][13]*SH_TAS[2] + P[8][5]*SH_TAS[1] - P[8][14]*SH_TAS[1] + P[8][6]*vd*SH_TAS[0]);
        Kfusion[9] = SK_TAS*(P[9][4]*SH_TAS[2] - P[9][13]*SH_TAS[2] + P[9][5]*SH_TAS[1] - P[9][14]*SH_TAS[1] + P[9][6]*vd*SH_TAS[0]);
        Kfusion[10] = SK_TAS*(P[10][4]*SH_TAS[2] - P[10][13]*SH_TAS[2] + P[10][5]*SH_TAS[1] - P[10][14]*SH_TAS[1] + P[10][6]*vd*SH_TAS[0]);
        Kfusion[11] = SK_TAS*(P[11][4]*SH_TAS[2] - P[11][13]*SH_TAS[2] + P[11][5]*SH_TAS[1] - P[11][14]*SH_TAS[1] + P[11][6]*vd*SH_TAS[0]);
        Kfusion[12] = SK_TAS*(P[12][4]*SH_TAS[2] - P[12][13]*SH_TAS[2] + P[12][5]*SH_TAS[1] - P[12][14]*SH_TAS[1] + P[12][6]*vd*SH_TAS[0]);
        Kfusion[13] = SK_TAS*(P[13][4]*SH_TAS[2] - P[13][13]*SH_TAS[2] + P[13][5]*SH_TAS[1] - P[13][14]*SH_TAS[1] + P[13][6]*vd*SH_TAS[0]);
        Kfusion[14] = SK_TAS*(P[14][4]*SH_TAS[2] - P[14][13]*SH_TAS[2] + P[14][5]*SH_TAS[1] - P[14][14]*SH_TAS[1] + P[14][6]*vd*SH_TAS[0]);
        Kfusion[15] = SK_TAS*(P[15][4]*SH_TAS[2] - P[15][13]*SH_TAS[2] + P[15][5]*SH_TAS[1] - P[15][14]*SH_TAS[1] + P[15][6]*vd*SH_TAS[0]);
        Kfusion[16] = SK_TAS*(P[16][4]*SH_TAS[2] - P[16][13]*SH_TAS[2] + P[16][5]*SH_TAS[1] - P[16][14]*SH_TAS[1] + P[16][6]*vd*SH_TAS[0]);
        Kfusion[17] = SK_TAS*(P[17][4]*SH_TAS[2] - P[17][13]*SH_TAS[2] + P[17][5]*SH_TAS[1] - P[17][14]*SH_TAS[1] + P[17][6]*vd*SH_TAS[0]);
        Kfusion[18] = SK_TAS*(P[18][4]*SH_TAS[2] - P[18][13]*SH_TAS[2] + P[18][5]*SH_TAS[1] - P[18][14]*SH_TAS[1] + P[18][6]*vd*SH_TAS[0]);
        Kfusion[19] = SK_TAS*(P[19][4]*SH_TAS[2] - P[19][13]*SH_TAS[2] + P[19][5]*SH_TAS[1] - P[19][14]*SH_TAS[1] + P[19][6]*vd*SH_TAS[0]);
        Kfusion[20] = SK_TAS*(P[20][4]*SH_TAS[2] - P[20][13]*SH_TAS[2] + P[20][5]*SH_TAS[1] - P[20][14]*SH_TAS[1] + P[20][6]*vd*SH_TAS[0]);
        varInnovVtas = 1.0f/SK_TAS;

        // Calculate the measurement innovation
        innovVtas = VtasPred - VtasMeas;
        // Check the innovation for consistency and don't fuse if > 5Sigma
        if ((innovVtas*innovVtas*SK_TAS) < 25.0)
        {
            // correct the state vector
            for (uint8_t j=0; j<=20; j++)
            {
                states[j] = states[j] - Kfusion[j] * innovVtas;
            }
            // normalise the quaternion states
            quatMag = sqrt(states[0]*states[0] + states[1]*states[1] + states[2]*states[2] + states[3]*states[3]);
            if (quatMag > 1e-12f)
            {
                for (uint8_t j= 0; j<=3; j++)
                {
                    float quatMagInv = 1.0f/quatMag;
                    states[j] = states[j] * quatMagInv;
                }
            }
            // correct the covariance P = (I - K*H)*P
            // take advantage of the empty columns in H to reduce the
            // number of operations
            for (uint8_t i = 0; i<=20; i++)
            {
                for (uint8_t j = 0; j<=3; j++) KH[i][j] = 0.0;
                for (uint8_t j = 4; j<=6; j++)
                {
                    KH[i][j] = Kfusion[i] * H_TAS[j];
                }
                for (uint8_t j = 7; j<=12; j++) KH[i][j] = 0.0;
                for (uint8_t j = 13; j<=14; j++)
                {
                    KH[i][j] = Kfusion[i] * H_TAS[j];
                }
                for (uint8_t j = 15; j<=20; j++) KH[i][j] = 0.0;
            }
            for (uint8_t i = 0; i<=20; i++)
            {
                for (uint8_t j = 0; j<=20; j++)
                {
                    KHP[i][j] = 0.0;
                    for (uint8_t k = 4; k<=6; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][k] * P[k][j];
                    }
                    for (uint8_t k = 13; k<=14; k++)
                    {
                        KHP[i][j] = KHP[i][j] + KH[i][k] * P[k][j];
                    }
                }
            }
            for (uint8_t i = 0; i<=20; i++)
            {
                for (uint8_t j = 0; j<=20; j++)
                {
                    P[i][j] = P[i][j] - KHP[i][j];
                }
            }
        }
    }
}

uint32_t millis()
{
    return IMUmsec;
}

void zeroRows(float covMat[21][21], uint8_t first, uint8_t last)
{
    uint8_t row;
    uint8_t col;
    for (row=first; row<=last; row++)
    {
        for (col=0; col<=20; col++)
        {
            covMat[row][col] = 0.0;
        }
    }
}

void zeroCols(float covMat[21][21], uint8_t first, uint8_t last)
{
    uint8_t row;
    uint8_t col;
    for (col=first; col<=last; col++)
    {
        for (row=0; row<=20; row++)
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
void StoreStates()
{
    static uint8_t storeIndex = 0;
    if (storeIndex > 49) storeIndex = 0;
    for (uint8_t i=0; i<=20; i++) storedStates[i][storeIndex] = states[i];
    statetimeStamp[storeIndex] = millis();
    storeIndex = storeIndex + 1;
}

// Output the state vector stored at the time that best matches that specified by msec
void RecallStates(float statesForFusion[21], uint32_t msec)
{
    long int timeDelta;
    long int bestTimeDelta = 200;
    uint8_t storeIndex;
    uint8_t bestStoreIndex = 0;
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
    if (bestTimeDelta < 200) // only output stored state if < 200 msec retrieval error
    {
        for (uint8_t i=0; i<=20; i++) statesForFusion[i] = storedStates[i][bestStoreIndex];
    }
    else // otherwise output current state
    {
        for (uint8_t i=0; i<=20; i++) statesForFusion[i] = states[i];
    }
}

void quat2Tnb(Mat3f &Tnb, float quat[4])
{
    // Calculate the nav to body cosine matrix
    float q00 = sq(quat[0]);
    float q11 = sq(quat[1]);
    float q22 = sq(quat[2]);
    float q33 = sq(quat[3]);
    float q01 =  quat[0]*quat[1];
    float q02 =  quat[0]*quat[2];
    float q03 =  quat[0]*quat[3];
    float q12 =  quat[1]*quat[2];
    float q13 =  quat[1]*quat[3];
    float q23 =  quat[2]*quat[3];

    Tnb.x.x = q00 + q11 - q22 - q33;
    Tnb.y.y = q00 - q11 + q22 - q33;
    Tnb.z.z = q00 - q11 - q22 + q33;
    Tnb.y.x = 2*(q12 - q03);
    Tnb.z.x = 2*(q13 + q02);
    Tnb.x.y = 2*(q12 + q03);
    Tnb.z.y = 2*(q23 - q01);
    Tnb.x.z = 2*(q13 - q02);
    Tnb.y.z = 2*(q23 + q01);
}

void quat2Tbn(Mat3f &Tbn, float quat[4])
{
    // Calculate the body to nav cosine matrix
    float q00 = sq(quat[0]);
    float q11 = sq(quat[1]);
    float q22 = sq(quat[2]);
    float q33 = sq(quat[3]);
    float q01 =  quat[0]*quat[1];
    float q02 =  quat[0]*quat[2];
    float q03 =  quat[0]*quat[3];
    float q12 =  quat[1]*quat[2];
    float q13 =  quat[1]*quat[3];
    float q23 =  quat[2]*quat[3];

    Tbn.x.x = q00 + q11 - q22 - q33;
    Tbn.y.y = q00 - q11 + q22 - q33;
    Tbn.z.z = q00 - q11 - q22 + q33;
    Tbn.x.y = 2*(q12 - q03);
    Tbn.x.z = 2*(q13 + q02);
    Tbn.y.x = 2*(q12 + q03);
    Tbn.y.z = 2*(q23 - q01);
    Tbn.z.x = 2*(q13 - q02);
    Tbn.z.y = 2*(q23 + q01);
}

void eul2quat(float quat[4], float eul[3])
{
    float u1 = cos(0.5f*eul[0]);
    float u2 = cos(0.5f*eul[1]);
    float u3 = cos(0.5f*eul[2]);
    float u4 = sin(0.5f*eul[0]);
    float u5 = sin(0.5f*eul[1]);
    float u6 = sin(0.5f*eul[2]);
    quat[0] = u1*u2*u3+u4*u5*u6;
    quat[1] = u4*u2*u3-u1*u5*u6;
    quat[2] = u1*u5*u3+u4*u2*u6;
    quat[3] = u1*u2*u6-u4*u5*u3;
}

void quat2eul(float y[3], float u[4])
{
    y[0] = atan2((2.0*(u[2]*u[3]+u[0]*u[1])) , (u[0]*u[0]-u[1]*u[1]-u[2]*u[2]+u[3]*u[3]));
    y[1] = -asin(2.0*(u[1]*u[3]-u[0]*u[2]));
    y[2] = atan2((2.0*(u[1]*u[2]+u[0]*u[3])) , (u[0]*u[0]+u[1]*u[1]-u[2]*u[2]-u[3]*u[3]));
}

void calcvelNED(float velNED[3], float gpsCourse, float gpsGndSpd, float gpsVelD)
{
    velNED[0] = gpsGndSpd*cos(gpsCourse);
    velNED[1] = gpsGndSpd*sin(gpsCourse);
    velNED[2] = gpsVelD;
}

void calcposNED(float posNED[3], float lat, float lon, float hgt, float latRef, float lonRef, float hgtRef)
{
    posNED[0] = earthRadius * (lat - latRef);
    posNED[1] = earthRadius * cos(latRef) * (lon - lonRef);
    posNED[2] = -(hgt - hgtRef);
}

void calcLLH(float posNED[3], float lat, float lon, float hgt, float latRef, float lonRef, float hgtRef)
{
    lat = latRef + posNED[0] * earthRadiusInv;
    lon = lonRef + posNED[1] * earthRadiusInv / cos(latRef);
    hgt = hgtRef - posNED[2];
}

void OnGroundCheck()
{
    onGround = (((sq(velNED[0]) + sq(velNED[1]) + sq(velNED[2])) < 4.0f) && (VtasMeas < 8.0f));
}

void CovarianceInit()
{
    // Calculate the initial covariance matrix P
    P[1][1]   = 0.25f*sq(1.0f*deg2rad);
    P[2][2]   = 0.25f*sq(1.0f*deg2rad);
    P[3][3]   = 0.25f*sq(10.0f*deg2rad);
    P[4][4]   = sq(0.7);
    P[5][5]   = P[4][4];
    P[6][6]   = sq(0.7);
    P[7][7]   = sq(15.0);
    P[8][8]   = P[7][7];
    P[9][9]   = sq(5.0);
    P[10][10] = sq(0.1*deg2rad*dtIMU);
    P[11][11] = P[10][10];
    P[12][12] = P[10][10];
    P[13][13] = sq(8.0f);
    P[14][4]  = P[13][13];
    P[15][15] = sq(0.02f);
    P[16][16] = P[15][15];
    P[17][17] = P[15][15];
    P[18][18] = sq(0.02f);
    P[19][19] = P[18][18];
    P[20][20] = P[18][18];
}

void readIMUData()
{
    static Vector3f lastAngRate;
    static Vector3f lastAccel;
    for (uint8_t j=0; j<=7; j++)
    {
        if (fscanf(pImuFile, "%f", &imuIn) != EOF) tempImu[j] = imuIn;
        else endOfData = true;
    }
    if (!endOfData)
    {
        IMUtimestamp  = tempImu[0];
        dtIMU     = 0.001f*(tempImu[1] - IMUmsec);
        IMUmsec   = tempImu[1];
        angRate.x = tempImu[2];
        angRate.y = tempImu[3];
        angRate.z = tempImu[4];
        accel.x   = tempImu[5];
        accel.y   = tempImu[6];
        accel.z   = tempImu[7];
        dAngIMU = 0.5f*(angRate + lastAngRate)*dtIMU;
        lastAngRate = angRate;
        dVelIMU = 0.5f*(accel + lastAccel)*dtIMU;
        lastAccel = accel;
    }
}

void readGpsData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update

    float gpsIn;

    static uint32_t GPSmsec = 0;
    static uint32_t lastGPSmsec = 0;
    static float tempGps[14];
    static float tempGpsPrev[14];
    static float GPStimestamp = 0;

    while (GPStimestamp <= IMUtimestamp)
    {
        for (uint8_t j=0; j<=13; j++)
        {
            tempGpsPrev[j] = tempGps[j];
            if (fscanf(pGpsFile, "%f", &gpsIn) != EOF) tempGps[j] = gpsIn;
            else endOfData = true;
        }
        if (!endOfData && (tempGps[1] == 3))
        {
            GPStimestamp  = tempGps[0];
            GPSmsec = tempGpsPrev[2];
            GPSstatus = tempGpsPrev[1];
            gpsCourse = deg2rad*tempGpsPrev[11];
            gpsGndSpd = tempGpsPrev[10];
            gpsVelD = tempGpsPrev[12];
            gpsLat = deg2rad*tempGpsPrev[6];
            gpsLon = deg2rad*tempGpsPrev[7] - pi;
            gpsHgt = tempGpsPrev[8];
        }
    }
    if (GPSmsec > lastGPSmsec)
    {
        lastGPSmsec = GPSmsec;
        newDataGps = true;
    }
    else
    {
        newDataGps = false;
    }
}

void readMagData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (MAGtimestamp <= IMUtimestamp)
    {
        for (uint8_t j=0; j<=7; j++)
        {
            tempMagPrev[j] = tempMag[j];
            if (fscanf(pMagFile, "%f", &magIn) != EOF) tempMag[j] = magIn;
            else endOfData = true;
        }
        if (!endOfData)
        {
            MAGtimestamp  = tempMag[0];
            MAGmsec = tempMagPrev[1];
            magData.x =  0.001f*(tempMagPrev[2] - tempMagPrev[5]);
            magBias.x = -0.001f*tempMagPrev[5];
            magData.y =  0.001f*(tempMagPrev[3] - tempMagPrev[6]);
            magBias.y = -0.001f*tempMagPrev[6];
            magData.z =  0.001f*(tempMagPrev[4] - tempMagPrev[7]);
            magBias.z = -0.001f*tempMagPrev[7];
        }
    }
    if (MAGmsec > lastMAGmsec)
    {
        lastMAGmsec = MAGmsec;
        newDataMag = true;
    }
    else
    {
        newDataMag = false;
    }
}

void readAirSpdData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (ADStimestamp <= IMUtimestamp)
    {
        for (uint8_t j=0; j<=9; j++)
        {
            tempAdsPrev[j] = tempAds[j];
            if (fscanf(pAdsFile, "%f", &adsIn) != EOF) tempAds[j] = adsIn;
            else endOfData = true;
        }
        if (!endOfData)
        {
            ADStimestamp  = tempAds[0];
            ADSmsec = tempAdsPrev[1];
            VtasMeas = EAS2TAS*tempAdsPrev[7];
        }
    }
    if (ADSmsec > lastADSmsec)
    {
        lastADSmsec = ADSmsec;
        newAdsData = true;
    }
    else
    {
        newAdsData = false;
    }
}

void readAhrsData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (AHRStimestamp <= IMUtimestamp)
    {
        for (uint8_t j=0; j<=6; j++)
        {
            tempAhrsPrev[j] = tempAhrs[j];
            if (fscanf(pAhrsFile, "%f", &ahrsIn) != EOF) tempAhrs[j] = ahrsIn;
            else endOfData = true;
        }
        if (!endOfData != EOF)
        {
            AHRStimestamp  = tempAhrs[0];
            AHRSmsec = tempAhrsPrev[1];
            for (uint8_t j=0; j<=2; j++)
            {
                ahrsEul[j] = deg2rad*tempAhrsPrev[j+2];
            }
            ahrsErrorRP = tempAhrs[5];
            ahrsErrorYaw = tempAhrs[6];
        }
    }
}

void InitialiseFilter()
{
    // Calculate initial filter quaternion states from ahrs solution
    float initQuat[4];
    eul2quat(initQuat, ahrsEul);



    // Calculate initial Tbn matrix and rotate Mag measurements into NED
    // to set initial NED magnetic field states
    Mat3f DCM;
    quat2Tbn(DCM, initQuat);
    Vector3f initMagNED;
    Vector3f initMagXYZ;
    initMagXYZ   = magData - magBias;
    initMagNED.x = DCM.x.x*initMagXYZ.x + DCM.x.y*initMagXYZ.y + DCM.x.z*initMagXYZ.z;
    initMagNED.y = DCM.y.x*initMagXYZ.x + DCM.y.y*initMagXYZ.y + DCM.y.z*initMagXYZ.z;
    initMagNED.z = DCM.z.x*initMagXYZ.x + DCM.z.y*initMagXYZ.y + DCM.z.z*initMagXYZ.z;

    // calculate initial velocities
    float initvelNED[3];
    calcvelNED(initvelNED, gpsCourse, gpsGndSpd, gpsVelD);

    //store initial lat,long and height
    latRef = gpsLat;
    lonRef = gpsLon;
    hgtRef = gpsHgt;

    // write to state vector
    for (uint8_t j=0; j<=3; j++) states[j] = initQuat[j]; // quaternions
    for (uint8_t j=0; j<=2; j++) states[j+4] = initvelNED[j]; // velocities
    for (uint8_t j=0; j<=7; j++) states[j+7] = 0.0f; // positiions, dAngBias, windVel
    states[15] = initMagNED.x; // Magnetic Field North
    states[16] = initMagNED.y; // Magnetic Field East
    states[17] = initMagNED.z; // Magnetic Field Down
    states[18] = magBias.x; // Magnetic Field Bias X
    states[19] = magBias.y; // Magnetic Field Bias Y
    states[20] = magBias.z; // Magnetic Field Bias Z

    statesInitialised = true;

    // initialise the covariance matrix
    CovarianceInit();

    //Define Earth rotation vector in the NED navigation frame
    calcEarthRateNED(earthRateNED, latRef);

    //Initialise summed variables used by covariance prediction
    summedDelAng.x = 0.0f;
    summedDelAng.y = 0.0f;
    summedDelAng.z = 0.0f;
    summedDelVel.x = 0.0f;
    summedDelVel.y = 0.0f;
    summedDelVel.z = 0.0f;
    dt = 0.0f;
}

void WriteFilterOutput()
{
    // filter states
    fprintf(pStateOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=20; i++)
    {
        fprintf(pStateOutFile," %e", states[i]);
    }
    fprintf(pStateOutFile,"\n");
    // Euler angles from filter states, AHRS euler angles and AHRS error RP and error Yaw
    fprintf(pEulOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=2; i++)
    {
        fprintf(pEulOutFile," %e %e", eulerEst[i], ahrsEul[i]);
    }
    fprintf(pEulOutFile," %e %e", ahrsErrorRP, ahrsErrorYaw);
    fprintf(pEulOutFile,"\n");
    // covariance matrix diagonals
    fprintf(pCovOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=20; i++)
    {
        fprintf(pCovOutFile," %e", P[i][i]);
    }
    fprintf(pCovOutFile,"\n");
    // velocity, position and height observations used by the filter
    fprintf(pRefPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pRefPosVelOutFile," %e %e %e %e %e %e", velNED[0], velNED[1], velNED[2], posNE[0], posNE[1], hgtMea);
    fprintf(pRefPosVelOutFile,"\n");
    // velocity and position innovations and innovation variances
    fprintf(pVelPosFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=5; i++)
    {
        fprintf(pVelPosFuseFile," %e %e", innovVelPos[i], varInnovVelPos[i]);
    }
    fprintf(pVelPosFuseFile,"\n");
    // magnetometer innovations and innovation variances
    fprintf(pMagFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=2; i++)
    {
        fprintf(pMagFuseFile," %e %e", innovMag[i], varInnovMag[i]);
    }
    fprintf(pMagFuseFile,"\n");
    // airspeed innovation and innovation variance
    fprintf(pTasFuseFile," %e", float(IMUmsec*0.001f));
    fprintf(pTasFuseFile," %e %e", innovVtas, varInnovVtas);
    fprintf(pTasFuseFile,"\n");

}

void calcEarthRateNED(Vector3f &omega, float latitude)
{
    //Define Earth rotation vector in the NED navigation frame
    omega.x  = earthRate*cos(latitude);
    omega.y  = 0.0;
    omega.z  = -earthRate*sin(latitude);
}

void readTimingData()
{
    float timeDataIn;
    float timeArray[9];
    for (uint8_t j=0; j<=8; j++)
    {
        if (fscanf(pTimeFile, "%f", &timeDataIn) != EOF)
        {
            timeArray[j] = timeDataIn;
        }
    }
    msecAlignTime = uint32_t(1000*timeArray[0]);
    msecStartTime = uint32_t(1000*timeArray[1]);
    msecEndTime   = uint32_t(1000*timeArray[2]);
    msecVelDelay  = uint32_t(timeArray[3]);
    msecPosDelay  = uint32_t(timeArray[4]);
    msecHgtDelay  = uint32_t(timeArray[5]);
    msecMagDelay  = uint32_t(timeArray[6]);
    msecTasDelay  = uint32_t(timeArray[7]);
    EAS2TAS       = timeArray[8];
}

void CloseFiles()
{
    fclose (pImuFile);
    fclose (pMagFile);
    fclose (pGpsFile);
    fclose (pAhrsFile);
    fclose (pAdsFile);
    fclose (pStateOutFile);
    fclose (pEulOutFile);
    fclose (pCovOutFile);
    fclose (pRefPosVelOutFile);
    fclose (pVelPosFuseFile);
    fclose (pMagFuseFile);
    fclose (pTasFuseFile);
    fclose (pTimeFile);
}
