
#include "estimator.h"

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

void WriteFilterOutput();

void CloseFiles();

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
uint32_t msecAlignTime;
uint32_t msecStartTime;
uint32_t msecEndTime;

float gpsGndSpd;

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
FILE * pGpsRawOUTFile;
FILE * pGpsRawINFile;
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
        printf(" %e", states[i]);
    }
    printf("\n");
    for (; i<4+6; i++)
    {
        printf(" %e", states[i]);
    }
    printf("\n");
    for (; i<4+6+6; i++)
    {
        printf(" %e", states[i]);
    }
    printf("\n");
    for (; i<n_states; i++)
    {
        printf(" %e", states[i]);
    }
    printf("\n");
}

int main()
{

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
    pGpsRawINFile = fopen ("GPSraw.txt","r");
    pGpsRawOUTFile = open_with_exit ("GPSrawOut.txt","w");
    pOnboardFile = fopen ("GPOSonboard.txt","r");
    pOnboardPosVelOutFile = open_with_exit ("OnboardVelPosDataOut.txt","w");

    printf("Filter start\n");

    memset(gpsRaw, 0, sizeof(gpsRaw));

    readTimingData();

    printf("First data instances loaded\n");

    float dt = 0.0f; // time lapsed since last covariance prediction

    while (true) {

        // read test data from files for next timestamp
        readIMUData();
        readGpsData();
        readMagData();
        readAirData();
        readAhrsData();
        readOnboardData();
        if (IMUmsec > msecEndTime)
        {
            CloseFiles();
            break;
        }

        if ((IMUmsec >= msecStartTime) && (IMUmsec <= msecEndTime))
        {
            // Initialise states, covariance and other data
            if ((IMUmsec > msecAlignTime) && !statesInitialised && (GPSstatus == 3))
            {
                if (pGpsRawINFile > 0)
                {
                    velNED[0] = gpsRaw[4];
                    velNED[1] = gpsRaw[5];
                    velNED[2] = gpsRaw[6];
                }
                else
                {
                    calcvelNED(velNED, gpsCourse, gpsGndSpd, gpsVelD);
                }
                InitialiseFilter(velNED);
            }

            // If valid IMU data and states initialised, predict states and covariances
            if (statesInitialised)
            {
                // Run the strapdown INS equations every IMU update
                UpdateStrapdownEquationsNED();
                #if 1
                // debug code - could be turned into a filter mnitoring/watchdog function
                float tempQuat[4];
                for (uint8_t j=0; j<4; j++) tempQuat[j] = states[j];
                quat2eul(eulerEst, tempQuat);
                for (uint8_t j=0; j<=2; j++) eulerDif[j] = eulerEst[j] - ahrsEul[j];
                if (eulerDif[2] > pi) eulerDif[2] -= 2*pi;
                if (eulerDif[2] < -pi) eulerDif[2] += 2*pi;
                #endif
                // store the predicted states for subsequent use by measurement fusion
                StoreStates(IMUmsec);
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
                    CovariancePrediction(dt);
                    summedDelAng = summedDelAng.zero();
                    summedDelVel = summedDelVel.zero();
                    dt = 0.0f;
                }
            }

            // Fuse GPS Measurements
            if (newDataGps && statesInitialised)
            {
                // Convert GPS measurements to Pos NE, hgt and Vel NED
                if (pGpsRawINFile > 0)
                {
                    velNED[0] = gpsRaw[4];
                    velNED[1] = gpsRaw[5];
                    velNED[2] = gpsRaw[6];
                }
                else
                {
                    calcvelNED(velNED, gpsCourse, gpsGndSpd, gpsVelD);
                }
                calcposNED(posNED, gpsLat, gpsLon, gpsHgt, latRef, lonRef, hgtRef);

                if (pOnboardFile > 0) {
                    calcposNED(onboardPosNED, onboardLat, onboardLon, onboardHgt, latRef, lonRef, hgtRef);
                    printf("in: %e %e %e posned 1: %e\n", onboardLat, onboardLon, onboardHgt, onboardPosNED[1]);
                }

                posNE[0] = posNED[0];
                posNE[1] = posNED[1];
                 // set fusion flags
                fuseVelData = true;
                fusePosData = true;
                // recall states stored at time of measurement after adjusting for delays
                RecallStates(statesAtVelTime, (IMUmsec - msecVelDelay));
                RecallStates(statesAtPosTime, (IMUmsec - msecPosDelay));
                // run the fusion step
                FuseVelposNED();
            }
            else
            {
                fuseVelData = false;
                fusePosData = false;
            }

            if (newAdsData && statesInitialised)
            {
                // Could use a blend of GPS and baro alt data if desired
                hgtMea = 1.0f*baroHgt + 0.0f*gpsHgt;
                fuseHgtData = true;
                // recall states stored at time of measurement after adjusting for delays
                RecallStates(statesAtHgtTime, (IMUmsec - msecHgtDelay));
                // run the fusion step
                FuseVelposNED();
            }
            else
            {
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

            onGround = false;


            // State vector:
            // 0-3: quaternions (q0, q1, q2, q3)
            // 4-6: Velocity - m/sec (North, East, Down)
            // 7-9: Position - m (North, East, Down)
            // 10-12: Delta Angle bias - rad (X,Y,Z)
            // 13-14: Wind Vector  - m/sec (North,East)
            // 15-17: Earth Magnetic Field Vector - milligauss (North, East, Down)
            // 18-20: Body Magnetic Field Vector - milligauss (X,Y,Z)
            printf("\n");
            printf("dtIMU: %8.6f, dt: %8.6f, imuMsec: %lld\n", dtIMU, dt, IMUmsec);
            printf("posNED: %8.4f, %8.4f, %8.4f, velNED: %8.4f, %8.4f, %8.4f\n", (double)posNED[0], (double)posNED[1], (double)posNED[2],
                (double)velNED[0], (double)velNED[1], (double)velNED[2]);
            printf("vTAS: %8.4f baro alt: %8.4f\n", VtasMeas, hgtMea);
            printf("mag: %8.4f, %8.4f, %8.4f\n", (double)magData.x, (double)magData.y, (double)magData.z);
            printf("states (quat)        [1-4]: %8.4f, %8.4f, %8.4f, %8.4f\n", (double)states[0], (double)states[1], (double)states[2], (double)states[3]);
            printf("states (vel m/s)     [5-7]: %8.4f, %8.4f, %8.4f\n", (double)states[4], (double)states[5], (double)states[6]);
            printf("states (pos m)      [8-10]: %8.4f, %8.4f, %8.4f\n", (double)states[7], (double)states[8], (double)states[9]);
            printf("states (delta ang) [11-13]: %8.4f, %8.4f, %8.4f\n", (double)states[10], (double)states[11], (double)states[12]);
            printf("states (wind)      [14-15]: %8.4f, %8.4f\n", (double)states[13], (double)states[14]);
            printf("states (earth mag) [16-18]: %8.4f, %8.4f, %8.4f\n", (double)states[15], (double)states[16], (double)states[17]);
            printf("states (body mag)  [19-21]: %8.4f, %8.4f, %8.4f\n", (double)states[18], (double)states[19], (double)states[20]);
            printf("states: %s %s %s %s %s %s %s %s %s\n",
                (statesInitialised) ? "INITIALIZED" : "NON_INIT",
                (onGround) ? "ON_GROUND" : "AIRBORNE",
                (fuseVelData) ? "FUSE_VEL" : "INH_VEL",
                (fusePosData) ? "FUSE_POS" : "INH_POS",
                (fuseHgtData) ? "FUSE_HGT" : "INH_HGT",
                (fuseMagData) ? "FUSE_MAG" : "INH_MAG",
                (fuseVtasData) ? "FUSE_VTAS" : "INH_VTAS",
                (useAirspeed) ? "USE_AIRSPD" : "IGN_AIRSPD",
                (useCompass) ? "USE_COMPASS" : "IGN_COMPASS");

        }
    }
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
                }
            }
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

void readAirData()
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
            baroHgt = 0;//tempAdsPrev[8];
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
    float tempOnboard[7];

    // wind data forward to one update past current IMU data time
    // and then take data from previous update
    while (onboardTimestamp <= IMUtimestamp)
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
            printf("velned onboard: %e %e %e %e %e %e\n", onboardLat, onboardLon, onboardHgt, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
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
    while (AHRStimestamp <= IMUtimestamp)
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
        }
    }
}

void WriteFilterOutput()
{

    float tempQuat[4];
    for (uint8_t j=0; j<4; j++) tempQuat[j] = states[j];
    quat2eul(eulerEst, tempQuat);

    // filter states
    fprintf(pStateOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<n_states; i++)
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
    for (uint8_t i=0; i<n_states; i++)
    {
        fprintf(pCovOutFile," %e", P[i][i]);
    }
    fprintf(pCovOutFile,"\n");
    // velocity, position and height observations used by the filter
    fprintf(pRefPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pRefPosVelOutFile," %e %e %e %e %e %e", velNED[0], velNED[1], velNED[2], posNE[0], posNE[1], hgtMea);
    fprintf(pRefPosVelOutFile,"\n");

    fprintf(pOnboardPosVelOutFile," %e", float(IMUmsec*0.001f));
    fprintf(pOnboardPosVelOutFile," %e %e %e %e %e %e", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    printf("velned onboard out: %e %e %e %e %e %e\n", onboardPosNED[0], onboardPosNED[1], -onboardPosNED[2] + hgtRef, onboardVelNED[0], onboardVelNED[1], onboardVelNED[2]);
    fprintf(pOnboardPosVelOutFile,"\n");

    // raw GPS outputs
    fprintf(pGpsRawOUTFile," %e", float(IMUmsec*0.001f));
    fprintf(pGpsRawOUTFile," %e %e %e %e %e %e", gpsRaw[1], gpsRaw[2], gpsRaw[3], gpsRaw[4], gpsRaw[5], gpsRaw[6]);
    fprintf(pGpsRawOUTFile,"\n");

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
    fclose (pGpsRawINFile);
    fclose (pGpsRawOUTFile);
    fclose (pOnboardPosVelOutFile);
    fclose (pOnboardFile);
}
