
#include "estimator.h"

#include <stdint.h>
#include <stdio.h>

void readIMUData();

void readGpsData();

void readMagData();

void readAirSpdData();

void readAhrsData();

void readTimingData();

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
                // // debug code - could be tunred into a filter mnitoring/watchdog function
                // float tempQuat[4];
                // for (uint8_t j=0; j<=3; j++) tempQuat[j] = states[j];
                // quat2eul(eulerEst, tempQuat);
                // for (uint8_t j=0; j<=2; j++) eulerDif[j] = eulerEst[j] - ahrsEul[j];
                // if (eulerDif[2] > pi) eulerDif[2] -= 2*pi;
                // if (eulerDif[2] < -pi) eulerDif[2] += 2*pi;
                // store the predicted states for subsequent use by measurement fusion
                StoreStates();
                // Check if on ground - status is used by covariance prediction
                OnGroundCheck();
                // sum delta angles and time used by covariance prediction
                summedDelAng = summedDelAng + correctedDelAng;
                summedDelVel = summedDelVel + correctedDelVel;
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
            // XXX the above comparison is always true - right fix?
            // !endOfData != EOF
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
    // filter states
    fprintf(pStateOutFile," %e", float(IMUmsec*0.001f));
    for (uint8_t i=0; i<=23; i++)
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
    for (uint8_t i=0; i<=23; i++)
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
