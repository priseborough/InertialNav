function [...
    quat, ... % quaternion state vector after fusion of measurements
    states, ... % state vector after fusion of measurements
    P, ... % state covariance matrix after fusion of corrections
    innovation, ... % Declination innovation - rad
    varInnov] ... %
    = FuseMagnetometer( ...
    quat, ... % predicted quaternion states
    states, ... % predicted states
    P, ... % predicted covariance
    magMea, ... % body frame magnetic flux measurements
    inhibitMag, ... % inhibit magnetic field state and covariance corrections
    Tbn)  % Estimated coordinate transformation matrix from body to NED frame

q0 = quat(1);
q1 = quat(2);
q2 = quat(3);
q3 = quat(4);

magXbias = states(13);
magYbias = states(14);
magZbias = states(15);

magN = states(10);
magE = states(11);
magD = states(12);

R_MAG = 0.05^2;
innovation = zeros(1,3);
varInnov = zeros(1,3);

% Calculate the predicted magnetic declination
magPred = transpose(Tbn)*[magN;magE;magD] + [magXbias;magYbias;magZbias];

for obsIndex = 1:3
    
    % Calculate corrections using X component
    if (obsIndex == 1)
        H = calcH_MAGX(magD,magE,magN,q0,q1,q2,q3);
    elseif (obsIndex == 2)
        H = calcH_MAGY(magD,magE,magN,q0,q1,q2,q3);
    elseif (obsIndex == 3)
        H = calcH_MAGZ(magD,magE,magN,q0,q1,q2,q3);
    end
    if (inhibitMag)
        H(10:15) = 0;
    end
    varInnov(obsIndex) = (H*P*transpose(H) + R_MAG);
    Kfusion = (P*transpose(H))/varInnov(obsIndex);
    innovation(obsIndex) = magPred(obsIndex) - magMea(obsIndex);
    
    % correct the state vector
    states(1:3) = 0; % rotation error vector does not accumulate and is always zeroed before updates are performed
    states = states - Kfusion * innovation(obsIndex);
    
    % the first 3 states represent the angular misalignment vector. This is
    % is used to correct the estimate quaternion
    % Convert the error rotation vector to its equivalent quaternion
    % error = truth - estimate
    rotationMag = sqrt(states(1)^2 + states(2)^2 + states(3)^2);
    if rotationMag<1e-6
        deltaQuat = single([1;0;0;0]);
    else
        deltaQuat = [cos(0.5*rotationMag); [states(1);states(2);states(3)]/rotationMag*sin(0.5*rotationMag)];
    end
    
    % Update the quaternion states by rotating from the previous attitude through
    % the delta angle rotation quaternion
    quat = [quat(1)*deltaQuat(1)-transpose(quat(2:4))*deltaQuat(2:4); quat(1)*deltaQuat(2:4) + deltaQuat(1)*quat(2:4) + cross(quat(2:4),deltaQuat(2:4))];
    
    % normalise the updated quaternion states
    quatMag = sqrt(quat(1)^2 + quat(2)^2 + quat(3)^2 + quat(4)^2);
    if (quatMag > 1e-12)
        quat = quat / quatMag;
    end
    
    % correct the covariance P = P - K*H*P
    P = P - Kfusion*H*P;
    
    % Force symmetry on the covariance matrix to prevent ill-conditioning
    % of the matrix which would cause the filter to blow-up
    P = 0.5*(P + transpose(P));
    
    % ensure diagonals are positive
    for i=1:15
        if P(i,i) < 0
            P(i,i) = 0;
        end
    end
    
end

end