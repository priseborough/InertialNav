#include "estimator_23states.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readIMUData();

void readGpsData();

void readMagData();

void readAirData();

void readRngData();

void readOptFlowData();

void readAhrsData();

void readTimingData();

void readOnboardData();

void WriteFilterOutput();

void CloseFiles();

float ConstrainFloat(float val, float min, float max);

bool endOfData = false; //boolean set to true when all files have returned data

// Estimated time delays (msec)
uint32_t msecVelDelay = 230;
uint32_t msecPosDelay = 210;
uint32_t msecHgtDelay = 350;
uint32_t msecRngDelay = 100;
uint32_t msecMagDelay = 30;
uint32_t msecTasDelay = 210;
uint32_t msecOptFlowDelay = 230;

// IMU input data variables
float imuIn;
float tempImu[8];
float IMUtimestamp;
static uint32_t IMUmsec = 0;

// Magnetometer input data variables
float magIn;
float tempMag[8];
float tempMagPrev[8];
float posNED[3];
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
float eulerEst[3]; // Euler angles calculated from filter states
float eulerDif[3]; // difference between Euler angle estimated by EKF and the AHRS solution
float gpsRaw[7];

// ADS input data variables
float adsIn;
float tempAds[10];
float tempAdsPrev[10];
float ADStimestamp = 0;
uint32_t ADSmsec = 0;
uint32_t lastADSmsec = 0;
float Veas;
bool newAdsData = false;
bool newDataGps = false;
bool newRngData = false;
bool newOptFlowData = false;

float onboardTimestamp = 0;
uint32_t onboardMsec = 0;
uint32_t lastOnboardMsec = 0;
bool newOnboardData = false;

float onboardPosNED[3];
float onboardVelNED[3];
float onboardLat;
float onboardLon;
float onboardHgt;

// input data timing
uint64_t msecAlignTime;
uint64_t msecStartTime;
uint64_t msecEndTime;

float gpsGndSpd;

AttPosEKF                   *_ekf;

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
FILE * pRngFuseFile;
FILE * pOptFlowFuseFile;
FILE * pTimeFile;
FILE * pGpsRawOUTFile;
FILE * pGpsRawINFile;
FILE * validationOutFile;
FILE * pOnboardPosVelOutFile;
FILE * pOnboardFile;

FILE * open_with_exit(const char* filename, const char* flags)
{
    FILE *f = fopen(filename, flags);

    if (!f) {
        printf("FAILED TO OPEN FILE: %s\n", filename);
        exit(1);
    }

    return f;
}

int printstates() {
    printf("States:\n");
    unsigned i = 0;
    printf("Quaternion:\n");
    for (; i<4; i++)
    {
        printf(" %e", _ekf->states[i]);
    }
    printf("\n");
    for (; i<4+6; i++)
    {
        printf(" %e", _ekf->states[i]);
    }
    printf("\n");
    for (; i<4+6+6; i++)
    {
        printf(" %e", _ekf->states[i]);
    }
    printf("\n");
    for (; i<n_states; i++)
    {
        printf(" %e", _ekf->states[i]);
    }
    printf("\n");

    return 0;
}

int main(int argc, char *argv[])
{
    // Instantiate EKF
    _ekf = new AttPosEKF();

    // open data files
    pImuFile = open_with_exit ("IMU.txt","r");
    pMagFile = open_with_exit ("MAG.txt","r");
    pGpsFile = open_with_exit ("GPS.txt","r");
    pAhrsFile = open_with_exit ("ATT.txt","r");
    pAdsFile = open_with_exit ("NTUN.txt","r");
    pTimeFile = open_with_exit ("timing.txt","r");
    pStateOutFile = open_with_exit ("StateDataOut.txt","w");
    pEulOutFile = open_with_exit ("EulDataOut.txt","w");
    pCovOutFile = open_with_exit ("CovDataOut.txt","w");
    pRefPosVelOutFile = open_with_exit ("RefVelPosDataOut.txt","w");
    pVelPosFuseFile = open_with_exit ("VelPosFuse.txt","w");
    pMagFuseFile = open_with_exit ("MagFuse.txt","w");
    pTasFuseFile = open_with_exit ("TasFuse.txt","w");
    pRngFuseFile = open_with_exit ("RngFuse.txt","w");
    pOptFlowFuseFile = open_with_exit ("OptFlowFuse.txt","w");
    pGpsRawINFile = fopen ("GPSraw.txt","r");
    pGpsRawOUTFile = open_with_exit ("GPSrawOut.txt","w");
    validationOutFile = fopen("ValidationOut.txt", "w");
    pOnboardFile = fopen ("GPOSonboard.txt","r");
    pOnboardPosVelOutFile = open_with_exit ("OnboardVelPosDataOut.txt","w");

    printf("Filter start\n");

    memset(gpsRaw, 0, sizeof(gpsRaw));

    readTimingData();

    printf("First data instances loaded\n");

    float dt = 0.0f; // time lapsed since last covariance prediction

    bool resetTests = false;

    // Test resets
    if (argc > 1 && (strcmp(argv[1], "reset") == 0)) {
        resetTests = true;
    }

    bool timeoutTested = false;
    bool nanTested = false;

    while (true) {

        // read test data from files for next timestamp
        unsigned nreads = 1;

        // Decide wether to perform any reset tests

        if (resetTests) {

            // Trigger a NaN reset after 25% of the log
            if (!nanTested && (IMUmsec > (msecEndTime - msecStartTime) / 4)) {
                _ekf->states[0] = 0.0f / 0.0f;
                _ekf->states[9] = 0.0f / 0.0f;
                nanTested = true;
                printf("WARNING: TRIGGERING NAN STATE ON PURPOSE!\n");
            }

            // Trigger a timeout at half the log
            if (!timeoutTested && (IMUmsec > (msecEndTime - msecStartTime) / 2)) {
                nreads = 20;
                timeoutTested = true;
                printf("WARNING: TRIGGERING SENSOR TIMEOUT ON PURPOSE!\n");
            }
        }

        // Supporting multiple reads at once to support timeout simulation.
        // The default is however to only read one dataset at a time

        // We need to re-do the dtIMU calculation here so that
        // dtIMU correctly skips if we skip readings.
        uint64_t IMUmsecPrev = IMUmsec;
        for (unsigned i = 0; i < nreads; i++) {
            readIMUData();
            readGpsData();
            readOptFlowData();
            readMagData();
            readAirData();
            readRngData();
            readAhrsData();
            readOnboardData();
        }

        // Apply dtIMU here after 1 or more reads, simulating skipped sensor
        // readings if desired.
        _ekf->dtIMU     = 0.001f*(IMUmsec - IMUmsecPrev);

        if (IMUmsec > msecEndTime || endOfData)
        {
            printf("Reached end at %8.4f seconds (end of logfile reached: %s)", IMUmsec/1000.0f, (endOfData) ? "YES" : "NO");
            CloseFiles();
            break;
        }

        if ((IMUmsec >= msecStartTime) && (IMUmsec <= msecEndTime))
        {
            // Initialise states, covariance and other data
            if ((IMUmsec > msecAlignTime) && !_ekf->statesInitialised && (_ekf->GPSstatus == 3))
            {
                if (pGpsRawINFile > 0)
                {
                    _ekf->velNED[0] = gpsRaw[4];
                    _ekf->velNED[1] = gpsRaw[5];
                    _ekf->velNED[2] = gpsRaw[6];
                }
                else
                {
                    _ekf->calcvelNED(_ekf->velNED, _ekf->gpsCourse, gpsGndSpd, _ekf->gpsVelD);
                }
                _ekf->InitialiseFilter(_ekf->velNED, _ekf->gpsLat, _ekf->gpsLon, _ekf->gpsHgt, 0.0f);
            }

            // If valid IMU data and states initialised, predict states and covariances
            if (_ekf->statesInitialised)
            {
                // Run the strapdown INS equations every IMU update
                _ekf->UpdateStrapdownEquationsNED();
                #if 1
                // debug code - could be turned into a filter monitoring/watchdog function
                float tempQuat[4];
                for (uint8_t j=0; j<4; j++) tempQuat[j] = _ekf->states[j];
                _ekf->quat2eul(eulerEst, tempQuat);
                for (uint8_t j=0; j<=2; j++) eulerDif[j] = eulerEst[j] - ahrsEul[j];
                if (eulerDif[2] > pi) eulerDif[2] -= 2*pi;
                if (eulerDif[2] < -pi) eulerDif[2] += 2*pi;
                #endif
                // store the predicted states for subsequent use by measurement fusion
                _ekf->StoreStates(IMUmsec);
                // Check if on ground - status is used by covariance prediction
                _ekf->OnGroundCheck();
                // sum delta angles and time used by covariance prediction
                _ekf->summedDelAng = _ekf->summedDelAng + _ekf->correctedDelAng;
                _ekf->summedDelVel = _ekf->summedDelVel + _ekf->dVelIMU;
                dt += _ekf->dtIMU;
                // perform a covariance prediction if the total delta angle has exceeded the limit
                // or the time limit will be exceeded at the next IMU update
                if ((dt >= (_ekf->covTimeStepMax - _ekf->dtIMU)) || (_ekf->summedDelAng.length() > _ekf->covDelAngMax))
                {
                    _ekf->CovariancePrediction(dt);
                    _ekf->summedDelAng.zero();
                    _ekf->summedDelVel.zero();
                    dt = 0.0f;
                }
            }

            // Fuse GPS Measurements
            if (newDataGps && _ekf->statesInitialised)
            {
                // Convert GPS measurements to Pos NE, hgt and Vel NED
                if (pGpsRawINFile > 0)
                {
                    _ekf->velNED[0] = gpsRaw[4];
                    _ekf->velNED[1] = gpsRaw[5];
                    _ekf->velNED[2] = gpsRaw[6];
                }
                else
                {
                    _ekf->calcvelNED(_ekf->velNED, _ekf->gpsCourse, gpsGndSpd, _ekf->gpsVelD);
                }
                _ekf->calcposNED(posNED, _ekf->gpsLat, _ekf->gpsLon, _ekf->gpsHgt, _ekf->latRef, _ekf->lonRef, _ekf->hgtRef);

                if (pOnboardFile > 0) {
                    _ekf->calcposNED(onboardPosNED, onboardLat, onboardLon, onboardHgt, _ekf->latRef, _ekf->lonRef, _ekf->hgtRef);

                }

                _ekf->posNE[0] = posNED[0];
                _ekf->posNE[1] = posNED[1];
                 // set fusion flags
                _ekf->fuseVelData = true;
                _ekf->fusePosData = true;
                // recall states stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtVelTime, (IMUmsec - msecVelDelay));
                _ekf->RecallStates(_ekf->statesAtPosTime, (IMUmsec - msecPosDelay));
                // run the fusion step
                _ekf->FuseVelposNED();
            }
            else
            {
                _ekf->fuseVelData = false;
                _ekf->fusePosData = false;
            }

            // Fuse Optical Flow Measurements
            if (newOptFlowData && _ekf->statesInitialised)
            {
                // recall states stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtOptFlowTime, (IMUmsec - msecOptFlowDelay));
                _ekf->fuseOptFlowData = true;
                _ekf->FuseOptFlow();
                _ekf->FuseOptFlow();
            }
            else
            {
                _ekf->fuseOptFlowData = false;
            }

            if (newAdsData && _ekf->statesInitialised)
            {
                // Could use a blend of GPS and baro alt data if desired
                _ekf->hgtMea = 1.0f*_ekf->baroHgt + 0.0f*_ekf->gpsHgt - _ekf->hgtRef - _ekf->baroHgtOffset;
                _ekf->fuseHgtData = true;
                // recall states stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtHgtTime, (IMUmsec - msecHgtDelay));
                // run the fusion step
                _ekf->FuseVelposNED();
            }
            else
            {
                _ekf->fuseHgtData = false;
            }

            // Fuse RangeFinder Measurements
            if (newRngData && _ekf->statesInitialised)
            {
                // recall states stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtRngTime, (IMUmsec - msecRngDelay));
                _ekf->fuseRngData = true;
                _ekf->FuseRangeFinder();
            }
            else
            {
                _ekf->fuseRngData = false;
            }

            // Fuse Magnetometer Measurements
            if (newDataMag && _ekf->statesInitialised)
            {
                _ekf->fuseMagData = true;
                _ekf->RecallStates(_ekf->statesAtMagMeasTime, (IMUmsec - msecMagDelay)); // Assume 50 msec avg delay for magnetometer data
                _ekf->magstate.obsIndex = 0;
                _ekf->FuseMagnetometer();
                _ekf->FuseMagnetometer();
                _ekf->FuseMagnetometer();
            }
            else
            {
                _ekf->fuseMagData = false;
            }

            // Fuse Airspeed Measurements
            if (newAdsData && _ekf->statesInitialised && _ekf->VtasMeas > 8.0f)
            {
                _ekf->fuseVtasData = true;
                _ekf->RecallStates(_ekf->statesAtVtasMeasTime, (IMUmsec - msecTasDelay)); // assume 100 msec avg delay for airspeed data
                _ekf->FuseAirspeed();
            }
            else
            {
                _ekf->fuseVtasData = false;
            }

            struct ekf_status_report ekf_report;

            /*
             *    CHECK IF THE INPUT DATA IS SANE
             */
            int check = _ekf->CheckAndBound(&ekf_report);

            switch (check) {
                case 0:
                    /* all ok */
                    break;
                case 1:
                {
                    printf("NaN in states, resetting\n");
                    printf("fail states: ");
                    for (unsigned i = 0; i < ekf_report.n_states; i++) {
                        printf("%f ",ekf_report.states[i]);
                    }
                    printf("\n");

                    printf("states after reset: ");
                    for (unsigned i = 0; i < ekf_report.n_states; i++) {
                        printf("%f ",_ekf->states[i]);
                    }
                    printf("\n");
                    break;
                }
                case 2:
                {
                    printf("stale IMU data, resetting\n");
                    break;
                }
                case 3:
                {
                    printf("switching to dynamic state\n");
                    break;
                }
                case 4:
                {
                    printf("excessive gyro offsets\n");
                    break;
                }
                case 5:
                {
                    printf("GPS velocity diversion\n");
                    break;
                }
                case 6:
                {
                    printf("Excessive covariances\n");
                    break;
                }


                default:
                {
                    printf("unknown reset condition\n");
                }
            }

            if (check) {
                printf("RESET OCCURED AT %d milliseconds\n", (int)IMUmsec);

                if (!ekf_report.velHealth || !ekf_report.posHealth || !ekf_report.hgtHealth || ekf_report.gyroOffsetsExcessive) {
                printf("health: VEL:%s POS:%s HGT:%s OFFS:%s\n",
                    ((ekf_report.velHealth) ? "OK" : "ERR"),
                    ((ekf_report.posHealth) ? "OK" : "ERR"),
                    ((ekf_report.hgtHealth) ? "OK" : "ERR"),
                    ((!ekf_report.gyroOffsetsExcessive) ? "OK" : "ERR"));
                }

                if (ekf_report.velTimeout || ekf_report.posTimeout || ekf_report.hgtTimeout || ekf_report.imuTimeout) {
                    printf("timeout: %s%s%s%s\n",
                        ((ekf_report.velTimeout) ? "VEL " : ""),
                        ((ekf_report.posTimeout) ? "POS " : ""),
                        ((ekf_report.hgtTimeout) ? "HGT " : ""),
                        ((ekf_report.imuTimeout) ? "IMU " : ""));
                }
            }

            // debug output
            //printf("Euler Angle Difference = %3.1f , %3.1f , %3.1f deg\n", rad2deg*eulerDif[0],rad2deg*eulerDif[1],rad2deg*eulerDif[2]);
            WriteFilterOutput();

            // State vector:
            // 0-3: quaternions (q0, q1, q2, q3)
            // 4-6: Velocity - m/sec (North, East, Down)
            // 7-9: Position - m (North, East, Down)
            // 10-12: Delta Angle bias - rad (X,Y,Z)
            // 13: Delta Velocity Z bias -m/s
            // 14-15: Wind Vector  - m/sec (North,East)
            // 16-18: Earth Magnetic Field Vector - milligauss (North, East, Down)
            // 19-21: Body Magnetic Field Vector - milligauss (X,Y,Z)
            // 22: Terrain Vertical Offset - m

            // printf("\n");
            // printf("dtIMU: %8.6f, dt: %8.6f, imuMsec: %u\n", _ekf->dtIMU, dt, IMUmsec);
            // printf("posNED: %8.4f, %8.4f, %8.4f, velNED: %8.4f, %8.4f, %8.4f\n", (double)_ekf->posNED[0], (double)_ekf->posNED[1], (double)_ekf->posNED[2],
            //     (double)_ekf->velNED[0], (double)_ekf->velNED[1], (double)_ekf->velNED[2]);
            // printf("vTAS: %8.4f baro alt: %8.4f\n", _ekf->VtasMeas, _ekf->hgtMea);
            // printf("mag: %8.4f, %8.4f, %8.4f\n", (double)_ekf->magData.x, (double)_ekf->magData.y, (double)_ekf->magData.z);
            // printf("states (quat)        [1-4]: %8.4f, %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[0], (double)_ekf->states[1], (double)_ekf->states[2], (double)_ekf->states[3]);
            // printf("states (vel m/s)     [5-7]: %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[4], (double)_ekf->states[5], (double)_ekf->states[6]);
            // printf("states (pos m)      [8-10]: %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[7], (double)_ekf->states[8], (double)_ekf->states[9]);
            // printf("states (delta ang) [11-13]: %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[10], (double)_ekf->states[11], (double)_ekf->states[12]);
            // printf("states (delta vel) [14]: %8.4ff\n", (double)_ekf->states[13]);
            // printf("states (wind)      [15-16]: %8.4f, %8.4f\n", (double)_ekf->states[14], (double)_ekf->states[15]);
            // printf("states (earth mag) [17-19]: %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[16], (double)_ekf->states[17], (double)_ekf->states[18]);
            // printf("states (body mag)  [20-22]: %8.4f, %8.4f, %8.4f\n", (double)_ekf->states[19], (double)_ekf->states[20], (double)_ekf->states[21]);
            // printf("states (terain offset) [23]: %8.4ff\n", (double)_ekf->states[22]);
            // printf("states: %s %s %s %s %s %s %s %s %s\n",
            //     (_ekf->statesInitialised) ? "INITIALIZED" : "NON_INIT",
            //     (_ekf->onGround) ? "ON_GROUND" : "AIRBORNE",
            //     (_ekf->fuseVelData) ? "FUSE_VEL" : "INH_VEL",
            //     (_ekf->fusePosData) ? "FUSE_POS" : "INH_POS",
            //     (_ekf->fuseHgtData) ? "FUSE_HGT" : "INH_HGT",
            //     (_ekf->fuseMagData) ? "FUSE_MAG" : "INH_MAG",
            //     (_ekf->fuseVtasData) ? "FUSE_VTAS" : "INH_VTAS",
            //     (_ekf->useAirspeed) ? "USE_AIRSPD" : "IGN_AIRSPD",
            //     (_ekf->useCompass) ? "USE_COMPASS" : "IGN_COMPASS");
        }
    }

    printf("\n\nSuccess: Finished processing complete dataset. Text files written.\n");
}

uint32_t millis()
{
    return IMUmsec;
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
        _ekf->dtIMU     = 0.001f*(tempImu[1] - IMUmsec);
        IMUmsec   = tempImu[1];
        _ekf->angRate.x = tempImu[2];
        _ekf->angRate.y = tempImu[3];
        _ekf->angRate.z = tempImu[4];
        _ekf->accel.x   = tempImu[5];
        _ekf->accel.y   = tempImu[6];
        _ekf->accel.z   = tempImu[7];
        _ekf->dAngIMU = 0.5f*(_ekf->angRate + lastAngRate)*_ekf->dtIMU;
        lastAngRate = _ekf->angRate;
        _ekf->dVelIMU = 0.5f*(_ekf->accel + lastAccel)*_ekf->dtIMU;
        lastAccel = _ekf->accel;
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

    while (GPStimestamp <= IMUtimestamp && !endOfData)
    {

        // Load APM GPS file format
        for (unsigned j = 0; j < 14; j++)
        {
            tempGpsPrev[j] = tempGps[j];
            if (fscanf(pGpsFile, "%f", &gpsIn) != EOF) {
                tempGps[j] = gpsIn;
            }
            else
            {
                endOfData = true;
                break;
            }
        }

        if (pGpsRawINFile > 0) {
            // Load RAW GPS file format in addition
            for (unsigned j = 0; j < sizeof(gpsRaw) / sizeof(gpsRaw[0]); j++)
            {
                if (fscanf(pGpsRawINFile, "%f", &gpsIn) != EOF) {
                    gpsRaw[j] = gpsIn;
                }
                else
                {
                    endOfData = true;
                    break;
                }
            }
        }

        if (!endOfData && (tempGps[1] == 3))
        {
            GPStimestamp  = tempGps[0];
            GPSmsec = tempGpsPrev[2];
            _ekf->GPSstatus = tempGpsPrev[1];
            _ekf->gpsCourse = deg2rad*tempGpsPrev[11];
            gpsGndSpd = tempGpsPrev[10];
            _ekf->gpsVelD = tempGpsPrev[12];
            _ekf->gpsLat = deg2rad*tempGpsPrev[6];
            _ekf->gpsLon = deg2rad*tempGpsPrev[7] - pi;
            _ekf->gpsHgt = tempGpsPrev[8];
        } else if (endOfData) {
            break;
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

void readOptFlowData()
{
    // currently synthesize optical flow measurements from GPS velocities and estimated angles
    if (newDataGps) {
        float q0 = 0.0f;
        float q1 = 0.0f;
        float q2 = 0.0f;
        float q3 = 1.0f;
        Vector3f relVelSensor;
        // Transformation matrix from nav to body axes
        Mat3f Tnb;
        // Transformation matrix from body to sensor axes
        // assume camera is aligned with Z body axis plus a misalignment
        // defined by 3 small angles about X, Y and Z body axis
        Mat3f Tbs;
        // Transformation matrix from navigation to sensor axes
        Mat3f Tns;
        // Copy required states to local variable names
        q0       = _ekf->statesAtVelTime[0];
        q1       = _ekf->statesAtVelTime[1];
        q2       = _ekf->statesAtVelTime[2];
        q3       = _ekf->statesAtVelTime[3];

        // Define rotation from body to sensor axes
        Tbs.x.y =  _ekf->a3;
        Tbs.y.x = -_ekf->a3;
        Tbs.x.z = -_ekf->a2;
        Tbs.z.x =  _ekf->a2;
        Tbs.y.z =  _ekf->a1;
        Tbs.z.y = -_ekf->a1;

        // calculate rotation from NED to body axes
        float q00 = q0*q0;
        float q11 = q1*q1;
        float q22 = q2*q2;
        float q33 = q3*q3;
        float q01 = q0 * q1;
        float q02 = q0 * q2;
        float q03 = q0 * q3;
        float q12 = q1 * q2;
        float q13 = q1 * q3;
        float q23 = q2 * q3;
        Tnb.x.x = q00 + q11 - q22 - q33;
        Tnb.y.y = q00 - q11 + q22 - q33;
        Tnb.z.z = q00 - q11 - q22 + q33;
        Tnb.y.x = 2*(q12 - q03);
        Tnb.z.x = 2*(q13 + q02);
        Tnb.x.y = 2*(q12 + q03);
        Tnb.z.y = 2*(q23 - q01);
        Tnb.x.z = 2*(q13 - q02);
        Tnb.y.z = 2*(q23 + q01);

        // calculate transformation from NED to sensor axes
        Tns = Tbs*Tnb;

        // calculate range from ground plain to centre of sensor fov assuming flat earth
        float range = ConstrainFloat(_ekf->rngMea,0.5f,100.0f);

        // calculate relative velocity in sensor frame
        Vector3f temp;
        temp.x=_ekf->velNED[0];
        temp.y=_ekf->velNED[1];
        temp.z=_ekf->velNED[2];
        relVelSensor = Tns*temp;

        // divide velocity by range  and include angular rate effects to get predicted angular LOS rates relative to X and Y axes
        _ekf->losData[0] =  relVelSensor.y/range;
        _ekf->losData[1] = -relVelSensor.x/range;

        newOptFlowData = true;
    } else {
        newOptFlowData = false;
    }
}

void readMagData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (MAGtimestamp <= IMUtimestamp && !endOfData)
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
            _ekf->magData.x =  0.001f*(tempMagPrev[2] - tempMagPrev[5]);
            _ekf->magBias.x = -0.001f*tempMagPrev[5];
            _ekf->magData.y =  0.001f*(tempMagPrev[3] - tempMagPrev[6]);
            _ekf->magBias.y = -0.001f*tempMagPrev[6];
            _ekf->magData.z =  0.001f*(tempMagPrev[4] - tempMagPrev[7]);
            _ekf->magBias.z = -0.001f*tempMagPrev[7];
        } else {
            break;
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

void readAirData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    // Currently synthesise a terrain measurement that is 5 m below the baro alt
    while (ADStimestamp <= IMUtimestamp && !endOfData)
    {
        for (uint8_t j=0; j<=9; j++)
        {
            tempAdsPrev[j] = tempAds[j];
            if (fscanf(pAdsFile, "%f", &adsIn) != EOF) {
                tempAds[j] = adsIn;
            } else {
                endOfData = true;
                break;
            }
        }
        if (!endOfData)
        {
            ADStimestamp  = tempAds[0];
            ADSmsec = tempAdsPrev[1];
            _ekf->VtasMeas = _ekf->EAS2TAS*tempAdsPrev[7];
            _ekf->baroHgt = tempAdsPrev[8];
        } else {
            break;
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

void readRngData()
{
    // Currently synthesise a terrain measurement that is 5 m below the baro alt
    if (newAdsData) {
        _ekf->rngMea = (_ekf->baroHgt  - _ekf->hgtRef - _ekf->baroHgtOffset + 5.0f) / _ekf->Tbn.z.z;
        newRngData = true;
    } else {
        newRngData = false;
    }
}

void readOnboardData()
{
    if (pOnboardFile <= 0)
        return;

    float tempOnboard[7];

    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (onboardTimestamp <= IMUtimestamp && !endOfData)
    {
        for (uint8_t j = 0; j < 7; j++)
        {
            float onboardIn;
            if (fscanf(pOnboardFile, "%f", &onboardIn) != EOF) tempOnboard[j] = onboardIn;
            else endOfData = true;
        }
        if (!endOfData)
        {
            onboardTimestamp  = tempOnboard[0];
            onboardLat = deg2rad*tempOnboard[1];
            onboardLon = deg2rad*tempOnboard[2] - pi;
            onboardHgt = tempOnboard[3];
            onboardVelNED[0] = tempOnboard[4];
            onboardVelNED[1] = tempOnboard[5];
            onboardVelNED[2] = tempOnboard[6];
            //printf("velned onboard: %e %e %e %e %e %e\n", onboardLat, onboardLon, onboardHgt, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
        } else {
            break;
        }
    }
    if (onboardMsec > lastOnboardMsec)
    {
        lastOnboardMsec = onboardMsec;
        newOnboardData = true;
    }
    else
    {
        newOnboardData = false;
    }
}

void readAhrsData()
{
    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (AHRStimestamp <= IMUtimestamp && !endOfData)
    {
        for (uint8_t j=0; j<=6; j++)
        {
            tempAhrsPrev[j] = tempAhrs[j];
            if (fscanf(pAhrsFile, "%f", &ahrsIn) != EOF) tempAhrs[j] = ahrsIn;
            else endOfData = true;
        }
        if (!endOfData)
        {
            AHRStimestamp  = tempAhrs[0];
            AHRSmsec = tempAhrsPrev[1];
            for (uint8_t j=0; j<=2; j++)
            {
                ahrsEul[j] = deg2rad*tempAhrsPrev[j+2];
            }
            ahrsErrorRP = tempAhrs[5];
            ahrsErrorYaw = tempAhrs[6];
        } else {
            break;
        }
    }
}

void WriteFilterOutput()
{

    float tempQuat[4];
    for (uint8_t j=0; j<4; j++) tempQuat[j] = _ekf->states[j];
    _ekf->quat2eul(eulerEst, tempQuat);

    // filter states
    fprintf(pStateOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<n_states; i++)
    {
        fprintf(pStateOutFile," %e", _ekf->states[i]);
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
    for (uint8_t i=0; i<n_states; i++)
    {
        fprintf(pCovOutFile," %e", _ekf->P[i][i]);
    }
    fprintf(pCovOutFile,"\n");
    // velocity, position and height observations used by the filter
    fprintf(pRefPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pRefPosVelOutFile," %e %e %e %e %e %e", _ekf->velNED[0], _ekf->velNED[1], _ekf->velNED[2], _ekf->posNE[0], _ekf->posNE[1], _ekf->hgtMea);
    fprintf(pRefPosVelOutFile,"\n");

    fprintf(pOnboardPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pOnboardPosVelOutFile," %e %e %e %e %e %e", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + _ekf->hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    // printf("velned onboard out: %e %e %e %e %e %e\n", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + _ekf->hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    fprintf(pOnboardPosVelOutFile,"\n");

    // raw GPS outputs
    fprintf(pGpsRawOUTFile," %e", float(IMUmsec*0.001f));
    fprintf(pGpsRawOUTFile," %e %e %e %e %e %e", gpsRaw[1], gpsRaw[2], gpsRaw[3], gpsRaw[4], gpsRaw[5], gpsRaw[6]);
    fprintf(pGpsRawOUTFile,"\n");

    // raw GPS put into local frame and integrated gyros
    fprintf(validationOutFile," %e", float(IMUmsec*0.001f));
    fprintf(validationOutFile," %e %e %e %e %e %e", _ekf->delAngTotal.x, _ekf->delAngTotal.y, _ekf->delAngTotal.z, _ekf->posNE[0], _ekf->posNE[1], _ekf->hgtMea);
    fprintf(validationOutFile,"\n");

    // velocity and position innovations and innovation variances
    fprintf(pVelPosFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=5; i++)
    {
        fprintf(pVelPosFuseFile," %e %e", _ekf->innovVelPos[i], _ekf->varInnovVelPos[i]);
    }
    fprintf(pVelPosFuseFile,"\n");

    // magnetometer innovations and innovation variances
    fprintf(pMagFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=2; i++)
    {
        fprintf(pMagFuseFile," %e %e", _ekf->innovMag[i], _ekf->varInnovMag[i]);
    }
    fprintf(pMagFuseFile,"\n");

    // airspeed innovation and innovation variance
    fprintf(pTasFuseFile," %e", float(IMUmsec*0.001f));
    fprintf(pTasFuseFile," %e %e", _ekf->innovVtas, _ekf->varInnovVtas);
    fprintf(pTasFuseFile,"\n");

    // range finder innovation and innovation variance
    fprintf(pRngFuseFile," %e", float(IMUmsec*0.001f));
    fprintf(pRngFuseFile," %e %e", _ekf->innovRng, _ekf->varInnovRng);
    fprintf(pRngFuseFile,"\n");

    // optical flow innovation and innovation variance
    fprintf(pOptFlowFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=1; i++)
    {
        fprintf(pOptFlowFuseFile," %e %e", _ekf->innovOptFlow [i], _ekf->varInnovOptFlow[i]);
    }
    fprintf(pOptFlowFuseFile,"\n");
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
    msecAlignTime = 1000*timeArray[0];
    msecStartTime = 1000*timeArray[1];
    msecEndTime   = 1000*timeArray[2];
    msecVelDelay  = timeArray[3];
    msecPosDelay  = timeArray[4];
    msecHgtDelay  = timeArray[5];
    msecMagDelay  = timeArray[6];
    msecTasDelay  = timeArray[7];
    _ekf->EAS2TAS       = timeArray[8];
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
    fclose (pGpsRawINFile);
    fclose (pGpsRawOUTFile);
    fclose (validationOutFile);
    fclose (pOnboardPosVelOutFile);
    fclose (pOnboardFile);
}

float ConstrainFloat(float val, float min, float max)
{
    float ret;
    if (val > max) {
        ret = max;
        ekf_debug("> max: %8.4f, val: %8.4f", (double)max, (double)val);
    } else if (val < min) {
        ret = min;
        ekf_debug("< min: %8.4f, val: %8.4f", (double)min, (double)val);
    } else {
        ret = val;
    }

    if (!isfinite(val)) {
        ekf_debug("constrain: non-finite!");
    }

    return ret;
}
