#include "estimator_22states.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readIMUData();

void readGpsData();

void readMagData();

void readAirData();

void readAhrsData();

void readTimingData();

void readOnboardData();

void readFlowData();

void readDistData();

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
uint32_t msecOptFlowDelay = 55;

// IMU input data variables
float imuIn;
float tempImu[8];
float IMUtimestamp;
static uint32_t IMUmsec = 0;
static uint64_t IMUusec = 0;

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

uint32_t flowMsec = 0;
uint32_t lastFlowMsec = 0;
bool newFlowData = false;

float flowTimestamp;      // in seconds
float flowRawPixelX;       // in pixels
float flowRawPixelY;       // in pixels
float flowDistance;        // in meters
float flowQuality;   // distance normalized between 0 and 1
float flowSensorId;        // sensor ID
float flowGyroX = 0.0f;
float flowGyroY = 0.0f;
float flowGyroBiasX = 0.0f;
float flowGyroBiasY = 0.0f;

float flowRadX;
float flowRadY;

float flowRawGroundSpeedX;
float flowRawGroundSpeedY;

uint32_t distMsec = 0;
uint32_t lastDistMsec = 0;
bool newDistData = false;

float distTimestamp = 0.0f;
bool distValid = false;
float distGroundDistance;
float distLastValidReading;

// input data timing
uint64_t msecAlignTime;
uint64_t msecStartTime;
uint64_t msecEndTime;

float gpsGndSpd;
float gpsCourse;
float gpsVelD;

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
FILE * pInFlowFile;
FILE * pInDistFile;
FILE * pOutFlowFile;

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

    pInFlowFile = fopen ("FLOW.txt","r");
    pInDistFile = fopen ("DIST.txt","r");
    pOutFlowFile = fopen ("FlowRawOut.txt","w");

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
    bool hgtTested = false;

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

            // Trigger a NaN altitude at 12.5% of the log
            if (!hgtTested && (IMUmsec > (msecEndTime - msecStartTime) / 8)) {
                _ekf->hgtMea = 0.0f / 0.0f;
                hgtTested = true;
                printf("WARNING: TRIGGERING NaN ALT MEASUREMENT ON PURPOSE!\n");
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
            readMagData();
            readAirData();
            readAhrsData();
            readOnboardData();
            readFlowData();
            readDistData();
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
                    _ekf->calcvelNED(_ekf->velNED, gpsCourse, gpsGndSpd, gpsVelD);
                }
                _ekf->InitialiseFilter(_ekf->velNED, _ekf->gpsLat, _ekf->gpsLon, _ekf->gpsHgt, 0.0f);

            } else if ((IMUmsec > msecAlignTime) && !_ekf->statesInitialised) {

                float initVelNED[3];

                initVelNED[0] = 0.0f;
                initVelNED[1] = 0.0f;
                initVelNED[2] = 0.0f;
                _ekf->posNE[0] = posNED[0];
                _ekf->posNE[1] = posNED[1];

                _ekf->InitialiseFilter(initVelNED, 0.0, 0.0, 0.0f, 0.0f);

            } else if (_ekf->statesInitialised) {

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

                // Set global time stamp used by EKF processes
                _ekf->globalTimeStamp_ms = IMUmsec;
            }

            // Fuse optical flow measurements
            if (newFlowData && _ekf->statesInitialised && _ekf->useOpticalFlow && flowQuality > 0.5 && _ekf->Tnb.z.z > 0.71f)
            {
                // recall states and angular rates stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtFlowTime, (IMUmsec - msecOptFlowDelay));
                _ekf->RecallOmega(_ekf->omegaAcrossFlowTime, (IMUmsec - 2*msecOptFlowDelay));

                // Calculate bias errorsfor flow sensor internal gyro
                flowGyroBiasX = 0.999f * flowGyroBiasX + 0.001f * (flowGyroX - _ekf->omegaAcrossFlowTime[0]);
                flowGyroBiasY = 0.999f * flowGyroBiasY + 0.001f * (flowGyroY - _ekf->omegaAcrossFlowTime[1]);

                //use sensor internal rates corrected for bias errors
                _ekf->omegaAcrossFlowTime[0] = flowGyroX - flowGyroBiasX;
                _ekf->omegaAcrossFlowTime[1] = flowGyroY - flowGyroBiasY;

                // calculate rotation matrix
                // Copy required states to local variable names
                float q0 = _ekf->statesAtFlowTime[0];
                float q1 = _ekf->statesAtFlowTime[1];
                float q2 = _ekf->statesAtFlowTime[2];
                float q3 = _ekf->statesAtFlowTime[3];
                float q00 = _ekf->sq(q0);
                float q11 = _ekf->sq(q1);
                float q22 = _ekf->sq(q2);
                float q33 = _ekf->sq(q3);
                float q01 = q0 * q1;
                float q02 = q0 * q2;
                float q03 = q0 * q3;
                float q12 = q1 * q2;
                float q13 = q1 * q3;
                float q23 = q2 * q3;
                _ekf->Tnb_flow.x.x = q00 + q11 - q22 - q33;
                _ekf->Tnb_flow.y.y = q00 - q11 + q22 - q33;
                _ekf->Tnb_flow.z.z = q00 - q11 - q22 + q33;
                _ekf->Tnb_flow.y.x = 2*(q12 - q03);
                _ekf->Tnb_flow.z.x = 2*(q13 + q02);
                _ekf->Tnb_flow.x.y = 2*(q12 + q03);
                _ekf->Tnb_flow.z.y = 2*(q23 - q01);
                _ekf->Tnb_flow.x.z = 2*(q13 - q02);
                _ekf->Tnb_flow.y.z = 2*(q23 + q01);

                // scale from raw pixel flow rate to radians/second
                //float scaleFactor = 0.03f; // best value for quad106.zip data using the 16 mm lens
                //float scaleFactor = 0.06f; // best value for InputFilesPX4_flow.zip data
                //float scaleFactor = 0.882f; // best value for quad123.zip data which outputs flow rates that have already been scaled to rad/sec
                float scaleFactor = 1.000f; // best value for quad-124.zip data which outputs flow rates that have already been scaled to rad/sec
                flowRadX = flowRawPixelX * scaleFactor;
                flowRadY = flowRawPixelY * scaleFactor;

                // calculate motion compensated angular flow rates used for fusion in the main nav filter
                _ekf->flowRadXYcomp[0] = flowRadX/_ekf->flowStates[0] + _ekf->omegaAcrossFlowTime[0];
                _ekf->flowRadXYcomp[1] = flowRadY/_ekf->flowStates[0] + _ekf->omegaAcrossFlowTime[1];

                // these flow rates are not motion compensated and are used for focal length scale factor estimation
                _ekf->flowRadXY[0] = flowRadX;
                _ekf->flowRadXY[1] = flowRadY;

                // perform optical flow fusion
                _ekf->fuseOptFlowData = true;
                _ekf->fuseRngData = false;

                // don't try to estimate focal length scale factor if GPS is not being used or there is no range finder.
                if (_ekf->useGPS && _ekf->useRangeFinder) {
                    _ekf->OpticalFlowEKF();
                }
                _ekf->FuseOptFlow();
                _ekf->fuseOptFlowData = false;

                // estimate speed over ground for cross-check of data (debugging only)
                float tempQuat[4];
                float euler[3];
                for (uint8_t j=0; j<4; j++) tempQuat[j] = _ekf->states[j];
                _ekf->quat2eul(euler, tempQuat);
                float bx = (flowRadX - _ekf->angRate.x) * distLastValidReading;
                float by = (flowRadY - _ekf->angRate.y) * distLastValidReading;
                flowRawGroundSpeedY = cos(euler[2]) * bx + -sin(euler[2]) * by;
                flowRawGroundSpeedX = -(sin(euler[2]) * bx + cos(euler[2]) * by);
            } else {
                _ekf->fuseOptFlowData = false;
            }

            // Fuse Ground distance Measurements
            if (newDistData && _ekf->statesInitialised && _ekf->useRangeFinder)
            {
                if (distValid > 0.0f && _ekf->Tnb.z.z > 0.9f) {
                    distLastValidReading = distGroundDistance;
                    _ekf->rngMea = distGroundDistance;
                    _ekf->fuseRngData = true;
                    _ekf->fuseOptFlowData = false;
                    _ekf->RecallStates(_ekf->statesAtRngTime, (IMUmsec - msecRngDelay));
                    _ekf->OpticalFlowEKF();
                    _ekf->fuseRngData = false;
                }
            }

                // Fuse GPS Measurements
                if (newDataGps)
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
                        _ekf->calcvelNED(_ekf->velNED, gpsCourse, gpsGndSpd, gpsVelD);
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
                _ekf->calcposNED(posNED, _ekf->gpsLat, _ekf->gpsLon, _ekf->gpsHgt, _ekf->latRef, _ekf->lonRef, _ekf->hgtRef);

                if (pOnboardFile > 0) {
                    _ekf->calcposNED(onboardPosNED, onboardLat, onboardLon, onboardHgt, _ekf->latRef, _ekf->lonRef, _ekf->hgtRef);
                }

                _ekf->posNE[0] = posNED[0];
                _ekf->posNE[1] = posNED[1];

                // fuse GPS
                if (_ekf->useGPS && IMUmsec < 1000) {
                    _ekf->fuseVelData = true;
                    _ekf->fusePosData = true;
                    _ekf->fuseHgtData = false;
                    // recall states stored at time of measurement after adjusting for delays
                    _ekf->RecallStates(_ekf->statesAtVelTime, (IMUmsec - msecVelDelay));
                    _ekf->RecallStates(_ekf->statesAtPosTime, (IMUmsec - msecPosDelay));
                    // record the last fix time
                    _ekf->lastFixTime_ms = IMUmsec;
                    // run the fusion step
                    _ekf->FuseVelposNED();
                } else {
                    _ekf->fuseVelData = false;
                    _ekf->fusePosData = false;
                    _ekf->fuseHgtData = false;
                }
            }
            else
            {
                _ekf->fuseVelData = false;
                _ekf->fusePosData = false;
                _ekf->fuseHgtData = false;
            }

            if (newAdsData && _ekf->statesInitialised)
            {
                // Could use a blend of GPS and baro alt data if desired
                _ekf->hgtMea = 1.0f*_ekf->baroHgt + 0.0f*_ekf->gpsHgt - _ekf->hgtRef - _ekf->baroHgtOffset;
                _ekf->fuseVelData = false;
                _ekf->fusePosData = false;
                _ekf->fuseHgtData = true;
                // recall states stored at time of measurement after adjusting for delays
                _ekf->RecallStates(_ekf->statesAtHgtTime, (IMUmsec - msecHgtDelay));
                // run the fusion step
                _ekf->FuseVelposNED();
                printf("time = %e \n", IMUtimestamp);
            }
            else
            {
                _ekf->fuseVelData = false;
                _ekf->fusePosData = false;
                _ekf->fuseHgtData = false;
            }

            // Fuse Magnetometer Measurements
            if (newDataMag && _ekf->statesInitialised && _ekf->useCompass)
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
            if (newAdsData && _ekf->statesInitialised && _ekf->VtasMeas > 8.0f && _ekf->useAirspeed)
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

            }

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

    printf("\n\nSuccess: Finished processing complete dataset. Text files written.\n");
}

uint32_t millis()
{
    return IMUmsec;
}

uint64_t getMicros()
{
    return IMUusec;
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
        IMUusec = tempImu[0] * 1e6;
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

        if (!endOfData && (tempGps[1] > 2) /* 3 or more */)
        {
            float now = tempGps[0];
            float gpsDt = (now - GPStimestamp);
            _ekf->updateDtGpsFilt(gpsDt);

            GPStimestamp  = tempGps[0];
            GPSmsec = tempGpsPrev[2];
            _ekf->GPSstatus = tempGpsPrev[1];
            gpsCourse = deg2rad*tempGpsPrev[11];
            gpsGndSpd = tempGpsPrev[10];
            gpsVelD = tempGpsPrev[12];
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
            float now = tempAds[0];
            float hgtDt = (now - ADStimestamp);
            _ekf->updateDtHgtFilt(hgtDt);

            ADStimestamp  = now;
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
            onboardMsec = tempOnboard[0];
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

void readFlowData()
{
    if (pInFlowFile <= 0)
        return;

    float temp[8];

    // read in current value
    while (flowTimestamp <= IMUtimestamp)
    {
        for (unsigned j = 0; j < (sizeof(temp) / sizeof(temp[0])); j++)
        {
            float in;
            if (fscanf(pInFlowFile, "%f", &in) != EOF) temp[j] = in;
            else endOfData = true;
        }
        if (!endOfData)
        {
            // timestamp, rawx, rawy, distance, quality, sensor id, flowGyroX, flowGyroY
            flowTimestamp = temp[0];       // in milliseconds
            flowRawPixelX = temp[1];        // in pixels
            flowRawPixelY = temp[2];        // in pixels
            flowDistance  = temp[3];         // in meters
            // catch glitches in logged data
            if (flowRawPixelX > 200 || flowRawPixelY > 200 || flowRawPixelX < -200 || flowRawPixelY < -200) {
                flowQuality = 0.0f;    // quality normalized between 0 and 1
            } else {
                flowQuality = temp[4] / 255;    // quality normalized between 0 and 1
            }
            flowSensorId = temp[5];         // sensor ID
            flowGyroX = temp[6];
            flowGyroY = temp[7];

            flowMsec = temp[0];
        }
    }
    if (flowMsec > lastFlowMsec)
    {
        // assume 1/2 interval effective delay associated with averaging inside the sensor
        msecOptFlowDelay = (flowMsec - lastFlowMsec)/2;
        lastFlowMsec = flowMsec;
        newFlowData = true;
    }
    else
    {
        newFlowData = false;
    }
}

void readDistData()
{
    if (pInDistFile <= 0)
        return;

    float temp[3];

    // read in current value
    while (distTimestamp <= IMUtimestamp)
    {
        for (unsigned j = 0; j < (sizeof(temp) / sizeof(temp[0])); j++)
        {
            float in;
            if (fscanf(pInDistFile, "%f", &in) != EOF) temp[j] = in;
            else endOfData = true;
        }
        if (!endOfData)
        {
            // timestamp, distance, flags
            distTimestamp  = temp[0];       // in milliseconds
            distGroundDistance = temp[1];   // in meters
            distValid = (temp[2] > 0.0f);   // reading is valid

            distMsec = temp[0];
        }
    }
    if (distMsec > lastDistMsec)
    {
        lastDistMsec = distMsec;
        newDistData = true;
    }
    else
    {
        newDistData = false;
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
    fprintf(pRefPosVelOutFile," %e %e %e %e %e %e", _ekf->velNED[0], _ekf->velNED[1], _ekf->velNED[2], _ekf->posNE[0], _ekf->posNE[1], -_ekf->hgtMea);
    fprintf(pRefPosVelOutFile,"\n");

    fprintf(pOnboardPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pOnboardPosVelOutFile," %e %e %e %e %e %e", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + _ekf->hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    // printf("velned onboard out: %e %e %e %e %e %e\n", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + _ekf->hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    fprintf(pOnboardPosVelOutFile,"\n");

    fprintf(pOutFlowFile, " %e", float(IMUmsec*0.001f));
    fprintf(pOutFlowFile, " %e %e %e %e %e %e %e %e %e\n", flowRadX, flowRadY, _ekf->angRate.x, _ekf->angRate.y, flowRawGroundSpeedX, flowRawGroundSpeedY, gpsRaw[4], gpsRaw[5], distLastValidReading);

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
    fprintf(pRngFuseFile," %e %e %e %e %e %e %e %e", _ekf->innovRng, _ekf->varInnovRng, _ekf->flowStates[1], _ekf->rngMea, _ekf->hgtMea, -_ekf->statesAtRngTime[9], -_ekf->flowStates[1], -_ekf->states[9] - distGroundDistance);
    fprintf(pRngFuseFile,"\n");

    // optical flow innovation and innovation variance
    fprintf(pOptFlowFuseFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=1; i++)
    {
        fprintf(pOptFlowFuseFile," %e %e", _ekf->innovOptFlow [i], _ekf->varInnovOptFlow[i]);
    }
    // focal length scale factor and height above ground estimate and innovations from optical flow rates
    fprintf(pOptFlowFuseFile," %e %e %e %e %e %e", _ekf->flowStates[0], _ekf->auxFlowObsInnov[0], _ekf->auxFlowObsInnov[1], _ekf->flowStates[1] - _ekf->states[9], distGroundDistance*_ekf->Tbn.z.z, - _ekf->states[9]);
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
