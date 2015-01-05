function H_OPT = calcH_OPT(navPosD,q0,q1,q2,q3,terrainState,vel_x,vel_y,vel_z)
%CALCH_OPT
%    H_OPT = CALCH_OPT(NAVPOSD,Q0,Q1,Q2,Q3,TERRAINSTATE,VEL_X,VEL_Y,VEL_Z)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    05-Jan-2015 13:37:33

t2 = q0.^2;
t3 = q1.^2;
t4 = q2.^2;
t5 = q3.^2;
t6 = navPosD-terrainState;
t7 = 1.0./t6.^2;
t8 = q0.*q3.*2.0;
t9 = t2-t3-t4+t5;
H_OPT = [-t7.*t9.*(vel_z.*(q0.*q1.*2.0+q2.*q3.*2.0)+vel_y.*(t2-t3+t4-t5)-vel_x.*(t8-q1.*q2.*2.0));t7.*t9.*(-vel_z.*(q0.*q2.*2.0-q1.*q3.*2.0)+vel_x.*(t2+t3-t4-t5)+vel_y.*(t8+q1.*q2.*2.0))];