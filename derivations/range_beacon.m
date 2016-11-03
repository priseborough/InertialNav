%% Derive and test an EKF to estinate initial posiiton from a series of range beacon measurements taken over time
clear all;
close all;

% define beacon position in NED frame
pos_bcn(:,1) = [10;10;0];
pos_bcn(:,2) = [-10;10;-10];
pos_bcn(:,3) = [-10;-10;0];
pos_bcn(:,4) = [10;-10;-10];

% define vehicle position in NED frame
veh_pos_truth = [20;20;-10];

% define truth range from vehicle to each beacon
for i=1:4
    range_truth(i) = sqrt(dot((veh_pos_truth - pos_bcn(:,i)),(veh_pos_truth - pos_bcn(:,i))));
end

% assume 0.3m measurement accuracy
rngMeasError = 0.3;

% for initial testing set range measurement to truth
range_meas = range_truth + rngMeasError*randn(size(range_truth));

% derive equations for a simple 3 state estimator
syms bcn_pos_n bcn_pos_e bcn_pos_d real;
syms veh_pos_n veh_pos_e veh_pos_d real;
syms R_RNG real;

stateVector = [veh_pos_n;veh_pos_e;veh_pos_d];
newStateVector = [veh_pos_n;veh_pos_e;veh_pos_d];


% derive the state transition matrix
F = jacobian(newStateVector, stateVector);
f = matlabFunction(F,'file','calcF.m');

% define a symbolic covariance matrix using strings to represent
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')'
% these can be substituted later to create executable code
for rowIndex = 1:3
    for colIndex = 1:3
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

% derive equations for fusion of range measurements
rng_pred = sqrt((veh_pos_n - bcn_pos_n)^2 + (veh_pos_e - bcn_pos_e)^2 + (veh_pos_d - bcn_pos_d)^2);

H_RNG = jacobian(rng_pred,stateVector); % measurement Jacobian
f = matlabFunction(H_RNG,'file','calcH_rng.m');

%% test EKF equations

% initialise state to a position in the middle of the observed beacons
stateEstimate = [0;0;0];
for bcnIndex=1:4
    stateEstimate = stateEstimate + pos_bcn(:,bcnIndex);
end
stateEstimate = stateEstimate/bcnIndex;

% calculate the maximum range innovation and use it to set the initial
% position uncertainty
maxRngErr = 0;
for bcnIndex=1:4
    range_pred = sqrt(dot((stateEstimate - pos_bcn(:,bcnIndex)),(stateEstimate - pos_bcn(:,bcnIndex))));
    innovation = range_pred -  range_meas(bcnIndex);
    if (abs(innovation) > maxRngErr)
        maxRngErr = abs(innovation);
    end
end
covariance = eye(3) * maxRngErr^2;

for timeIndex = 1:10
    % save estimates
    stateError(:,timeIndex) = stateEstimate - veh_pos_truth;
    stateVariance(1,timeIndex) = covariance(1,1);
    stateVariance(2,timeIndex) = covariance(2,2);
    stateVariance(3,timeIndex) = covariance(3,3);
    
    % get new measurement
    range_meas = range_truth + 0.3*randn(size(range_truth));
    
    % state prediction - time invariant process modeul suitable for static
    % alignment application
    stateEstimate = stateEstimate;
    
    % covariance prediction
    covariance = covariance + 0.1*eye(3);
    
    % range fusion
    for measIndex=1:4
        % calculate the observation jacobian
        H_RNG = calcH_rng(pos_bcn(3,measIndex),pos_bcn(2,measIndex),pos_bcn(1,measIndex),stateEstimate(3),stateEstimate(2),stateEstimate(1));
        
        % calculate the innovation
        range_pred = sqrt(dot((stateEstimate - pos_bcn(:,measIndex)),(stateEstimate - pos_bcn(:,measIndex))));
        innovation = range_pred -  range_meas(measIndex);
        
        % calculate the Kalman gain
        K_RNG = (covariance*transpose(H_RNG))/(H_RNG*covariance*transpose(H_RNG) + rngMeasError^2);
        
        % update the state vector
        stateEstimate = stateEstimate - K_RNG * innovation;
        
        % update the covariance matrix
        covariance = covariance - K_RNG*H_RNG*covariance;
        
        % Force symmetry on the covariance matrix to prevent ill-conditioning
        covariance = 0.5*(covariance + transpose(covariance));
        
        % ensure diagonals are positive
        for stateIndex=1:3
            if covariance(stateIndex,stateIndex) < 0
                covariance(stateIndex,stateIndex) = 0;
            end
        end
        
    end
    
end

% save final estimates
stateError(:,timeIndex+1) = stateEstimate - veh_pos_truth;
stateVariance(1,timeIndex+1) = covariance(1,1);
stateVariance(2,timeIndex+1) = covariance(2,2);
stateVariance(3,timeIndex+1) = covariance(3,3);

% plot the state estimation error and uncertainty
close all;
figure;
for index = 1:3
    subplot(3,1,index);
    plot([stateError(index,:);sqrt(stateVariance(index,:));-sqrt(stateVariance(index,:))]');
    grid on;
end