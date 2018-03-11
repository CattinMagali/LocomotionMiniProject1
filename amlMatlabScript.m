
function main()

    clc; close all; clear all; % a bit of cleaning
    %% ======================================================================
    % 1)  ASSIGNEMENT 0: SETTINGS
    %======================================================================
    % replace by your own information
    LAST_NAME = 'Cattin'; % your last name
    FIRST_NAME = 'Magali'; % your first name
    SCIPER = '236251'; % your sciper number
    DATASET_NUMBER = 6; % your dataset number
    
    % settings
    DATASET_NAME = ['amlWalkingDataStruct' num2str(DATASET_NUMBER) '.mat'];
    load(DATASET_NAME); % load data structure
    scriptOutResults = [];
     
    %% ======================================================================
    % 2)  ASSIGNEMENT 1: TEMPORAL AND SPATIAL GAIT ANALYSIS
    %======================================================================
    %% Exercice 2.A.1
    %---------------------
    % align gyroscope TF with the foot AF using function
    % alignGyroscope(data) and plot the results
    % <<< ENTER YOUR CODE HERE >>>
    
    [leftGyro_AlignedwithFoot, rightGyro_AlignedwithFoot] = alignGyroscopeTF2AF(data);
    
    figure(1);
    plot(data.imu.left.time, leftGyro_AlignedwithFoot);
    xlabel('Time [s]'); ylabel('Angular velocity [deg/s]'); legend('x', 'y', 'z'); 
    title('Gyroscope measurement, left foot');
   
    figure(2);
    plot(data.imu.right.time, rightGyro_AlignedwithFoot);
    xlabel('Time [s]'); ylabel('Angular velocity [deg/s]'); legend('x', 'y', 'z'); 
    title('Gyroscope measurement, right foot');
    
    %% Exercice 2.A.2
    %---------------------
    % detect IC, TC, FF for all midswing-to-midswing cycles (Left foot).
    % <<< ENTER YOUR CODE HERE >>>
    
    a = zeros(length(data.imu.left.midswings)-1,1);
    FFL_idx = a; FFL = a; ICL = a; ICL_idx = a; TCL = a; TCL_idx = a; 
       
    for i = 1:(length(data.imu.left.midswings)-1)
        %Determine the index which marks the beginning of the flat-foot
        %region: beginning of the midswing cycle + distance to reach the
        %plateau (40% of the cycle corresponding to swing)
        halfCycleL(i) = round(data.imu.left.midswings(i) + 0.5*length(data.imu.left.midswings(i):data.imu.left.midswings(i+1)));
        %Considering the minimum between the beginning of the midswing
        %cycle and half of the cycle --> get the Initial Contact
        [ICL(i), ICL_idx(i)] = min(leftGyro_AlignedwithFoot(data.imu.left.midswings(i):halfCycleL(i),3));
        %But the index found corresponds to a "local" index starting at the
        %beginning of the midswing cycle, thus we need to add all previous indices
        ICL_idx(i) = ICL_idx(i) + data.imu.left.midswings(i);
        %Considering the minimum between the second half of the cycle and 
        %the end of the midswing cycle --> get the Terminal Contact
        [TCL(i), TCL_idx(i)] = min(leftGyro_AlignedwithFoot(halfCycleL(i):data.imu.left.midswings(i+1),3));
        TCL_idx(i) = TCL_idx(i) + halfCycleL(i);
        %Mid-stance is half way between the IC and the TC
        FFL_idx(i) = round(ICL_idx(i) + (TCL_idx(i) -  ICL_idx(i))/2);
        %Consider the 3rd component, z
        FFL(i) = leftGyro_AlignedwithFoot(FFL_idx(i), 3);
    end
    
    % add results to the output structure
    scriptOutResults.imu.leftIC = ICL_idx; % insert your detection of IC
    scriptOutResults.imu.leftTC = TCL_idx; % insert your detection of TC
    scriptOutResults.imu.leftFF = FFL_idx; % insert your detection of FF
    % detect IC, TC, FF for all midswing-to-midswing cycles (Right foot).
    % <<< ENTER YOUR CODE HERE >>>
    
    b = zeros(length(data.imu.right.midswings)-1,1);
    FFR_idx = b; FFR = b; ICR = b; ICR_idx = b; TCR = b; TCR_idx = b; 
   
    for i = 1:(length(data.imu.right.midswings)-1)
        halfCycleR(i) = round(data.imu.right.midswings(i) + 0.5*length(data.imu.right.midswings(i):data.imu.right.midswings(i+1)));
        [ICR(i), ICR_idx(i)] = min(rightGyro_AlignedwithFoot(data.imu.right.midswings(i):halfCycleR(i),3));
        ICR_idx(i) = ICR_idx(i) + data.imu.right.midswings(i);
        [TCR(i), TCR_idx(i)] = min(rightGyro_AlignedwithFoot(halfCycleR(i):data.imu.right.midswings(i+1),3));
        TCR_idx(i) = TCR_idx(i) + halfCycleR(i);
        FFR_idx(i) = round(ICR_idx(i) + (TCR_idx(i) -  ICR_idx(i))/2);
        FFR(i) = rightGyro_AlignedwithFoot(FFR_idx(i), 3);
    end
       
    % add results to the output structure
    scriptOutResults.imu.rightIC = ICR_idx; % insert your detection of IC
    scriptOutResults.imu.rightTC = TCR_idx; % insert your detection of TC
    scriptOutResults.imu.rightFF = FFR_idx; % insert your detection of FF
    
    %% Exercice 2.A.3
    %---------------------
    % plot detection results for right foot
    % <<< ENTER YOUR CODE HERE >>>
    
    %2 Gait Cycles of right foot
    figure(3);
       
    plot(data.imu.right.time, rightGyro_AlignedwithFoot(:, 3)); 
    hold on;
    plot(data.imu.right.time(ICR_idx),  rightGyro_AlignedwithFoot(ICR_idx, 3), 'rx');
    plot(data.imu.right.time(TCR_idx),  rightGyro_AlignedwithFoot(TCR_idx, 3), 'gx');
    plot(data.imu.right.time(FFR_idx),  rightGyro_AlignedwithFoot(FFR_idx, 3), 'bx'); 
    plot(data.imu.right.time(data.imu.right.midswings),  rightGyro_AlignedwithFoot(data.imu.right.midswings, 3), 'mx');   
    
    title('Gait Cycles, right foot');
    xlabel('Time [s]'); ylabel('Angular velocity [deg/s]');
    legend('Angular velocity [deg/s]', 'Initial Contact', 'Terminal Contact', 'MidStance', 'MidSwing');
    
    %left foot
    figure(4);
       
    plot(data.imu.left.time, leftGyro_AlignedwithFoot(:, 3)); 
    hold on;
    plot(data.imu.left.time(ICL_idx),  leftGyro_AlignedwithFoot(ICL_idx, 3), 'rx');
    plot(data.imu.left.time(TCL_idx),  leftGyro_AlignedwithFoot(TCL_idx, 3), 'gx');
    plot(data.imu.left.time(FFL_idx),  leftGyro_AlignedwithFoot(FFL_idx, 3), 'bx'); 
    plot(data.imu.left.time(data.imu.left.midswings),  leftGyro_AlignedwithFoot(data.imu.left.midswings, 3), 'mx');   
    
    title('Gait Cycles, left foot');
    xlabel('Time [s]'); ylabel('Angular velocity [deg/s]');
    legend('Angular velocity [deg/s]', 'Initial Contact', 'Terminal Contact', 'MidStance', 'MidSwing');
    
    %% Exercice 2.A.4
    %---------------------
    % compute the stance phase duration, gait cycle time and cadence for
    % left and right leg.
    % <<< ENTER YOUR CODE HERE >>>
    
    %Left foot
    %Creation of a tab storing stance & gait cycle duration as well as
    %cadence in different columns for each midswing cycle
    stance_gaitCycle_cadenceL = zeros(length(data.imu.left.midswings)-2, 3);
    
    %the first & last gait cycles are cut in the middle because of midswing segmentation
    for i = 1:(length(data.imu.left.midswings)-2)
        stance_gaitCycle_cadenceL(i,:) = [data.imu.left.time(TCL_idx(i)) - data.imu.left.time(ICL_idx(i)), 
                                         data.imu.left.time(ICL_idx(i+1)) - data.imu.left.time(ICL_idx(i)),
                                         120/(data.imu.left.time(ICL_idx(i+1)) - data.imu.left.time(ICL_idx(i)))];
    end
    
    meanL = mean(stance_gaitCycle_cadenceL, 1);
    stdL = std(stance_gaitCycle_cadenceL,0, 1);
  
    %Same for the right foot
    stance_gaitCycle_CadenceR = zeros(length(data.imu.right.midswings)-2, 3);
    for i = 1:(length(data.imu.right.midswings)-2)
        stance_gaitCycle_CadenceR(i,:) = [data.imu.right.time(TCR_idx(i)) - data.imu.right.time(ICR_idx(i)), 
                                         data.imu.right.time(ICR_idx(i+1)) - data.imu.right.time(ICR_idx(i)),
                                         120/(data.imu.right.time(ICR_idx(i+1)) - data.imu.right.time(ICR_idx(i)))]; 
    end
    
    meanR = mean(stance_gaitCycle_CadenceR, 1);
    stdR = std(stance_gaitCycle_CadenceR,0, 1);
    
    scriptOutResults.imu.leftMeanStance = meanR(1); % insert the mean stance phase duration of left foot
    scriptOutResults.imu.leftSTDStance = stdR(1); % insert the stance phase duration STD of left foot
    scriptOutResults.imu.rightMeanStance = meanL(1); % insert the mean stance phase duration of right foot
    scriptOutResults.imu.rightSTDStance = stdR(1); % insert the stance phase duration STD of right foot
    scriptOutResults.imu.leftMeanGaitCycleTime = meanL(2); % insert the mean gait cycle time of the left foot
    scriptOutResults.imu.leftSTDGaitCycleTime = stdL(2); % insert the gait cycle time STD of the left foot
    scriptOutResults.imu.rightMeanGaitCycleTime = meanR(2); % insert the mean gait cycle time of the right foot
    scriptOutResults.imu.rightSTDGaitCycleTime = stdR(2); % insert the gait cycle time STD of the right foot
    scriptOutResults.imu.leftMeanCadence = meanL(3); % insert the mean cadence of the left foot
    scriptOutResults.imu.leftSTDCadence = stdL(3); % insert the cadence STD of the left foot
    scriptOutResults.imu.rightMeanCadence = meanR(3); % insert the mean cadence of the left foot
    scriptOutResults.imu.rightSTDCadence = stdR(3); % insert the cadence STD of the left foot
    
    %% Exercice 2.A.5
    %---------------------
    % Compare the mean cadence of the right leg to the left leg 
    % <<< ENTER YOUR CODE HERE >>>
    
    diffCadence = meanL(3)-meanR(3);
    diffPropCadence = (meanL(3)-meanR(3))/meanR(3);
    
    %% Exercice 2.A.6
    %---------------------
    % Estimate the coefficient of variation (in %) of the gait cycle time 
    % obtained from of the right foot
    % <<< ENTER YOUR CODE HERE >>>
    
    coeffVar_gaitCycleR = stdR(2)/meanR(2);
    
    scriptOutResults.imu.cvGCT = coeffVar_gaitCycleR; % insert CV GCT right foot
    
    
    %% Exercice 2.B.1 (Bonus)
    %---------------------
    % compute the pitch angle from the gyroscope pitch angular velocity.
    % <<< ENTER YOUR CODE HERE >>>
    
    gaitCycle = data.imu.left.time(FFL_idx(8):FFL_idx(9));
    length_gaitCycle = length(FFL_idx(8):FFL_idx(9));
    w_pitch = zeros(length_gaitCycle-1, 1);
      
    for i = 1:length_gaitCycle-1
        w_pitch(i) = trapz(gaitCycle(1:i+1), leftGyro_AlignedwithFoot(FFL_idx(8):(FFL_idx(8)+i), 3));  
    end
    
    %% Exercice 2.B.2 (Bonus)
    %---------------------
    % correct the drift on the pitch angle signal
    % <<< ENTER YOUR CODE HERE >>>
    
    error_FFL = w_pitch(end) - w_pitch(1);
    linearError = zeros(length_gaitCycle-1, 1);
    w_pitch_corr = zeros(length_gaitCycle-1, 1);
    
    for i = 1:length_gaitCycle-1
        linearError(i) = error_FFL*i/length_gaitCycle;
        w_pitch_corr(i) = w_pitch(i)-linearError(i);
    end
    
    %% Exercice 2.B.3 (Bonus)
    %---------------------
    % plot gyroscope pitch angular velocity, pitch angle with and without
    % drift
    % <<< ENTER YOUR CODE HERE >>>
    
    figure(5)
    plot(data.imu.left.time(FFL_idx(8):FFL_idx(9)-1), w_pitch, '-r');
    hold on; plot(data.imu.left.time(FFL_idx(8):FFL_idx(9)-1), w_pitch_corr, '--r');
    ylabel('Pitch angle [deg]');
    
    title('Foot-flat to foot-flat Gait cycle'); 
    xlabel('Time [s]'); 
    legend('Pitch angle [deg]', 'Corrected pitch angle [deg]');
    axis([9.4640 10.5020 -100 40]);
    
%     yyaxis left;
%     plot(gaitCycle, leftGyro_AlignedwithFoot(FFL_idx(8):(FFL_idx(9)), 3));
%     ylabel('Angular velocity [deg/s]');
%    
%     yyaxis right;
%     plot(data.imu.left.time(FFL_idx(8):FFL_idx(9)-1), w_pitch);
%     hold on; plot(data.imu.left.time(FFL_idx(8):FFL_idx(9)-1), w_pitch_corr);
%     ylabel('Pitch angle [deg/s]');
%     
%     title('Foot-flat to foot-flat Gait cycle'); 
%     xlabel('Time [s]'); 
%     legend('Angular velocity [deg/s]', 'Pitch angle [deg]', 'Corrected pitch angle [deg]');
    
    %% ======================================================================
    %  3) ASSIGNEMENT 2: FRAME & ORIENTATIONS
    %======================================================================
    %% Exercice 3.A.1
    %---------------------
    % gravity vector in the right foot IMU TF
    % <<< ENTER YOUR CODE HERE >>>
    
    TFgR = mean(data.imu.right.accelstatic, 1);
    
    scriptOutResults.imu.rightGravityTF = TFgR; % insert right foot TFg here
    
    
    %% Exercice 3.A.3
    %---------------------
    % find R_TFg_Y_AF between TFg and Y_AF
    % <<< ENTER YOUR CODE HERE >>>
    
    Y_AF = [0, 1, 0];    
    R_TFg_Y_AF = get3DRotationMatrixA2B(TFgR,Y_AF) 
    
    scriptOutResults.imu.rightRotationYAF = R_TFg_Y_AF; % insert R_TFg_Y_AF
    
    
    %% Exercice 3.A.4
    %---------------------
    % plot the static signals before and after the rotation
    % <<< ENTER YOUR CODE HERE >>>
    
    figure(6);
    subplot(2,1,1);
    plot(data.imu.right.accelstatic);
    title('Measure of the accelerometer in TF'); legend('x','y','z');
    xlabel('Samples'); ylabel('Static accceleration [g]');
    
    statAccel_AF = data.imu.right.accelstatic*R_TFg_Y_AF'; 
    subplot(2,1,2);
    plot(statAccel_AF);
    title('Measure of the accelerometer in AF'); legend('x','y','z');
    xlabel('Samples'); ylabel('Static accceleration [g]');
    
    
    
    %% Exercice 3.B.1
    %---------------------
    % construct the technical frame of the left foot
    % <<< ENTER YOUR CODE HERE >>>
    
    %TF is defined by the axis pointing from M1 (leftCenterFoot) towards,
    %M2 (leftLateralFoot) for the x-axis, M3 (leftMedialFoot) for the
    %z-axis and upwards for the y-axis. 
    
    %x-axis: M1 --> M2, then normalization to get a orthonormal dimensional
    %space
    M12 = mean(data.motioncameras.static.leftLateralFoot) - mean(data.motioncameras.static.leftCenterFoot);
    TFx = M12/norm(M12)
    %z-axis: M1 --> M3
    M13 = mean(data.motioncameras.static.leftMedialFoot) - mean(data.motioncameras.static.leftCenterFoot);
    TFz = M13/norm(M13)
    %y-axis corresponds to the cross product of the x- and z-axis since it
    %is orthogonal to them
    M1u = cross(TFz, TFx); 
    TFy = M1u/norm(M1u)
   
    scriptOutResults.motioncameras.tfX = TFx; % insert TF x-axis
    scriptOutResults.motioncameras.tfY = TFy; % insert TF y-axis
    scriptOutResults.motioncameras.tfZ = TFz; % insert TF z-axis
    
    
    %% Exercice 3.B.2
    %---------------------
    % compute R_TF_GF
    % <<< ENTER YOUR CODE HERE >>>
    
    R_TF_GF = [TFx', TFy', TFz']
    
    %Rotation matrix 
    scriptOutResults.motioncameras.R_TF_GF = R_TF_GF; % insert R_TF_GF
    
    
    %% Exercice 3.B.3
    %---------------------
    % construct the anatomical frame of the left foot
    % <<< ENTER YOUR CODE HERE >>>
    
    %AF is defined by the vector pointing from M4 to M5 but beginning at
    %M1, the vector pointing upward from M1 & // to Yo and the vector
    %orthogonal to these 2 latter, pointing forward
    
    %z-axis: M4 --> M5
    M45 = mean(data.motioncameras.static.leftMedialMalleolus) - mean(data.motioncameras.static.leftLateralMalleolus);
    AFz = M45/norm(M45)
    %y-axis: M1 --> up
    AFy = [0, 1, 0];
    %x-axis corresponds to the cross product of the y- and z-axis since it
    %is orthogonal to them
    M1f = cross(AFy, AFz); 
    AFx = M1f/norm(M1f)
    
    scriptOutResults.motioncameras.afX = AFx; % insert AF x-axis
    scriptOutResults.motioncameras.afY = AFy; % insert AF y-axis
    scriptOutResults.motioncameras.afZ = AFz; % insert AF z-axis
    
    
    %% Exercice 3.B.4
    %---------------------
    % compute R_AF_GF
    % <<< ENTER YOUR CODE HERE >>>
    R_AF_GF = [AFx', AFy', AFz']
    
    scriptOutResults.motioncameras.R_AF_GF = R_AF_GF; % insert R_AF_GF
    
    
    %% Exercice 3.B.5
    %---------------------
    % compute R_TF_AF
    % <<< ENTER YOUR CODE HERE >>>
    
    R_TF_AF = inv(R_AF_GF)*R_TF_GF 
    
    scriptOutResults.motioncameras.R_TF_AF = R_TF_AF; % insert R_TF_AF
    
    
    %% Exercice 3.C
    %---------------------
    % (1) compute TF for walking
    % <<< ENTER YOUR CODE HERE >>>
    
   for i = 1:length(data.motioncameras.walking.time)
        M12w = data.motioncameras.walking.leftLateralFoot(i, :) - data.motioncameras.walking.leftCenterFoot(i,:);
        TFwx = M12w./norm(M12w);
        %z-axis: M1 --> M3
        M13w = data.motioncameras.walking.leftMedialFoot(i,:) - data.motioncameras.walking.leftCenterFoot(i,:);
        TFwz = M13w./norm(M13w);
        %y-axis corresponds to the cross product of the x- and z-axis since it
        %is orthogonal to them
        M1wu = cross(TFwz, TFwx); 
        TFwy = M1wu./norm(M1wu);
    end
    
    % (2) compute AF for walking
    % <<< ENTER YOUR CODE HERE >>>
    
   for i = 1:length(data.motioncameras.walking.time)
        %z-axis: M4 --> M5
        M45w = data.motioncameras.walking.leftMedialMalleolus(i,:) - data.motioncameras.dynamic.leftLateralMalleolus(i,:);
        AFwz = M45w./norm(M45w);
        %y-axis: M1 --> up
        AFwy = [0, 1, 0];
        %x-axis corresponds to the cross product of the y- and z-axis since it
        %is orthogonal to them
        M1wf = cross(AFwy, AFwz); 
        AFwx = M1wf/norm(M1wf)
   end
    
    % (3) compute the pitch angle
    % <<< ENTER YOUR CODE HERE >>>
    
    % (4) plot the pitch angle
    % <<< ENTER YOUR CODE HERE >>>
    
    
    
    %% ======================================================================
    %  4) ASSIGNEMENT 3: KINETIC ANALYSIS
    %======================================================================
    %% Exercice 4.A.1
    %---------------------
    % compute the force signal for all cell (transform pressure into force)
    % <<< ENTER YOUR CODE HERE >>>
    
    % detect the sample index of the multiple IC, TS, HO, TO
    % <<< ENTER YOUR CODE HERE >>>
    
    % store IC, TS, HO and TO detection index
    scriptOutResults.insole.rightIC = []; % insert the index of the right foot IC events
    scriptOutResults.insole.rightTS = []; % insert the index of the right foot TS events
    scriptOutResults.insole.rightHO = []; % insert the index of the right foot HO events
    scriptOutResults.insole.rightTO = []; % insert the index of the right foot TO events
    
    
    % Exercice 4.A.2
    %---------------------
    % Add in the frame below a graph showing F_rear and F_Fore at least two
    % strides of the right foot where HS, TS, HO and TO events are 
    % correctly detected and show these event in your plot. Do not forget 
    % to add labels on each axis and a legend for all signals. 
    % <<< ENTER YOUR CODE HERE >>>
    
    
    % Exercice 4.A.3
    %---------------------
    % for the two cycles above, estimate the foot-flat duration
    % <<< ENTER YOUR CODE HERE >>>
    
    
    % Exercice 4.A.4
    %---------------------
    % estimate the total vertical force signal recorded by the insole 
    % during one foot-flat period and plot the result
    % <<< ENTER YOUR CODE HERE >>>

     
    % Exercice 4.B.2
    %---------------------
    % compute the net force at the ankle (F_A) and the net moment at the
    % ankle (M_A) for every timesample during one footflat period
    % <<< ENTER YOUR CODE HERE >>>
    
     scriptOutResults.insole.F_A = [];
     scriptOutResults.insole.M_A = [];
    
     
    % Exercice 4.C.3
    %---------------------
    % compute the mean value of F_A and M_A
    % <<< ENTER YOUR CODE HERE >>>
    
     scriptOutResults.insole.MeanF_A = [];
     scriptOutResults.insole.MeanM_A = [];

    
    %======================================================================
    %  5) ENDING TASKS
    %======================================================================
    % Save the output structure
    %---------------------
    save([LAST_NAME '_' FIRST_NAME '_outStruct.mat'],'scriptOutResults');
    
end %function main


%==========================================================================
%   AML LIBRARY
%==========================================================================
function [alignedLeftGyro,alignedRightGyro] = alignGyroscopeTF2AF(data)
% ALIGNGYROSCOPETF2AF aligns the gyroscope TF with the foot AF
%   [A, B] = alignGyroscopeTF2AF(D) returns the angular velocity measured
%   by the gyroscope, but expressed in the anatomical frame of the foot.
%   Here D is the complete data structure given in the project, A is the
%   Nx3 matrix with the left foot angular velocity, and B is the Nx3 matrix
%   with the right foot angular velocity.
    alignedLeftGyro = data.imu.left.gyro * data.imu.left.calibmatrix.';
    alignedRightGyro = data.imu.right.gyro * data.imu.right.calibmatrix.';
end % function

function [R] = get3DRotationMatrixA2B(A,B)
%GET3DROTATIONMATRIXA2B returns the rotation matrix which rotates vector A onto vector B. 
%   This function must be used as such B = R * A with 
%   R = get3DRotationMatrixA2B(A,B).
%
%   INPUTS:
%       - A: vector in 3D space
%       - B: vector in 3D space
%   OUTPUTS:
%       - R: 3x3 rotation matrix

    % The formula used can be found on: 
    % http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
    A = A(:)/norm(A); B = B(:)/norm(B); % force Nx1 format and normalize
    v = cross(A,B);
    s = norm(v);
    c = dot(A,B); %scalar product
    Vskew = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    R = eye(3) + Vskew + Vskew^2 * ((1-c)/s^2);
end % function
