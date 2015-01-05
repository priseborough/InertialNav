dt = 0.1;
Popt = 0.25;
navPosD = -100;
q0 = 1.0;
q1 = 0.0;
q2 = 0.0;
q3 = 0.0;
terrainState = 0.0;
vel_x = 20.0;
vel_y = 0.0;
vel_z = 0.0;
R_LOS = 0.3^2;
flowX = 0.0;
flowY = -0.25;
for index = 1:1000
    % predict covariance
    Popt = Popt + dt*2.0;
    % calculate predicted measurement
    predFlow = calcPredFlow(navPosD,q0,q1,q2,q3,terrainState,vel_x,vel_y,vel_z);
    % calculate innovation
    innovFlow = predFlow - [flowX;flowY];
    % calculate measurement Jacobian
    t2 = sq(q0);
    t3 = sq(q1);
    t4 = sq(q2);
    t5 = sq(q3);
    t6 = navPosD-terrainState;
    t7 = 1.0/sq(t6);
    t8 = q0*q3*2.0;
    t9 = t2-t3-t4+t5;
    H_OPT(1) = -t7*t9*( vel_z*(q0*q1*2.0+q2*q3*2.0)+vel_y*(t2-t3+t4-t5)-vel_x*(t8-q1*q2*2.0));
    H_OPT(2) =  t7*t9*(-vel_z*(q0*q2*2.0-q1*q3*2.0)+vel_x*(t2+t3-t4-t5)+vel_y*(t8+q1*q2*2.0));
    for index = 1:2
    % calculate innovation variances
        auxFlowObsInnovVar(index) = H_OPT(index)*Popt*H_OPT(index) + R_LOS;
    % calculate Kalman gains
        K_OPT(index) = Popt*H_OPT(index)/auxFlowObsInnovVar(index);
    % calculate state update
        terrainState = terrainState - K_OPT(obsIndex)*innovFlow(obsIndex);
    % calculate covariance update
        Popt = Popt - K_OPT(obsIndex)*H_OPT(obsIndex)*Popt;
    end
    tout(index) = dt*index;
    xout(index) = terrainState;
    pout(index) = Popt;
    varout(index,:) = auxFlowObsInnovVar;
    
end