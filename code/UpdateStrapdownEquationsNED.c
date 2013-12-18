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
    accNavMag = sqrt(sq(delVelNav.x) + sq(delVelNav.y) + sq(delVelNav.z)) / dt;

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
