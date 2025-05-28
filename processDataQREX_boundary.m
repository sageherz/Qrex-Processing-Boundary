clear
close all
clc

addpath ("library\")

%% Selecting boundary test type
testTypeList = {'Ground Effect (Hover)','Ground Effect (Forward Flight)',...
    'Wall/Corner Effect (Hover)','Partial Boundary (z/R constant)', 'Partial Boundary (x/R constant)'};
[testType, ~] = listdlg('PromptString',{'Select test type.'},'SelectionMode','single',...
    'ListString',testTypeList);

%% TURN ON/OFF FOR TEST W/ VS. W/O KIEL PRESSURE RAKE
qPressure = 1;
if qPressure
    probe_angles = input("ENTER PROBE ANGLES (deg): ");
    probe_locs = round(([5, 6.95, 8.8, 10.6, 12.6, 14.4, 16.3, 18.3])./19.05,2);
end

%% ENTER AMBIENT AIR CONDITIONS FROM TIME OF FLIGHT TEST
T = 22; % room temperature (C)
hr = 16; % humidity (%)
P = 99200; % pressure (Pa)

%% Getting ambient room density
rho = air_density(T,hr,P);

%% Rotor geometry and vehicle orientation (where applicable)
R = 0.1905; % rotor radius (m)
A = pi*R^2; % rotor disk area (m^2)

if testType ~= 1
    skew = input("ENTER SKEW ANGLE (CCW = + and CW = -): ");
    L = 13.9375*0.0254; % arm length relative to vehicle center (m)
    delta = 2/100; % offset of lidar sensor for vehicle middle (m)
    X = L*cosd(45-abs(skew));
    x = X - delta; % adjusting for lidar offset
end

%% Save data from ROS (1 or 2) file to .mat file
[rosfile, rospath] = uigetfile({'*.db3';'*.bag'}, 'Select Bag File', 'multiselect', 'on' ); % getting file and path of .bag or .db3 file

rosfile_type = rosfile(end-2:end); % db3 (ROS2) or bag (ROS1) file extension
if strcmp(rosfile_type, "bag") == 1
    ROS_ver = 1; % ROS version (ROS1)
    rawData = saveData_indoor(rosfile,rospath); % returns name of file with extracted data from bag file (ROS1)
    matname = strcat(rosfile(1:19));
elseif strcmp(rosfile_type, "db3") == 1
    ROS_ver = 0; % ROS version (ROS2)
    rawData = saveData_ROS2(rosfile,rospath); % returns name of file with extracted data from db3 file (ROS2)
    matname = strcat(rosfile(1:29));
    rospath = rospath(1:length(rospath)-28);
end

%% Save data from log file to .mat file
[logfile, logpath] = uigetfile( '*.ulg', 'Select Log File', 'multiselect', 'on' ); % getting file and path of .ulg file
logData = Ulog_Parser_Function(logfile,logpath);

%% Process thrust data (unbias, filter, scale)
Thrust_filt = medfilt1(rawData.data.Thrust,100); % filtering raw thrust
Thrust_filt(1,:) = [];

% plotting filtered thrust to select takeoff and landing points
figure
plot(Thrust_filt)
title("Select takeoff and landing points")
[flight_pts,~] = getpts; flight_pts = round(flight_pts,0);
if length(flight_pts) < 2 % constant unbias (subtracting off start/end point, whichever is provided)
    Thrust_unbias = unbiasThrust_constant(Thrust_filt, flight_pts);
else % linear unbias
    Thrust_unbias = unbiasThrust_linear(Thrust_filt, flight_pts);
end

save(fullfile(rospath, matname),"flight_pts","-append") % saving takeoff/landing points to mat file

kThrust = -[0.00093021, -0.00080045, -0.00098768, 0.00082721]; % load cell scale factors (updated 07/08/24)
Thrust_scale = kThrust.*Thrust_unbias;

%%  Aligning and scaling RPM data (from ESC)
RPM_filt = medfilt1(rawData.data.ESC.RPM,15); % filtering raw RPM data
tESC = rawData.times.tESC;

if ROS_ver % checking ROS version...
    tQrex = rawData.times.tData;
else
    tQrex = rawData.times.tQrex;
end

for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tESC-tQrex(i))); %#ok<*SAGROW>
end
RPM_align = RPM_filt(corr_idx,:); RPM_align(1,:) = []; % aligning RPM to thrust data
RPM_scale = RPM_align; % scaling RPM data

%% Computing CT
Omega = RPM_scale.*((2*pi)/60); % converting RPM to rad/s
CT = Thrust_scale./(rho*((Omega.*R).^2)*A);

%% Calculating power (from ESC or Pi)
power_source = "Pi"; % ENTER SELECTION FOR SOURCE OF POWER DATA (either from ESC ("ESC") or from Pi ("Pi"))

% filtering current and voltage
if strcmp(power_source,"Pi") == 1
    Current_filt = medfilt1(rawData.data.Current,100); Current_filt(1,:) = [];
    Voltage_filt = medfilt1(rawData.data.Voltage(:,end),100); Voltage_filt(1,:) = [];

    kCurrent = [0.02082, 0.02002, 0.02064, 0.02036];
    Current_filt = Current_filt - Current_filt(1,:);
    Current_align = kCurrent.*Current_filt;

    kVoltage = 0.006375425; cVoltage = -0.0882915;
    Voltage_align = kVoltage*Voltage_filt + cVoltage;
elseif strcmp(power_source,"ESC") == 1
    Current_filt = medfilt1(rawData.data.ESC.Current,10); Current_filt(1,:) = [];
    Voltage_filt = medfilt1(rawData.data.ESC.Voltage,10); Voltage_filt(1,:) = [];

    Current_align = Current_filt(corr_idx,:); % aligning Current to thrust data
    Voltage_align = Voltage_filt(corr_idx,:); % aligning Voltage to thrust data
end

Power = Current_align.*Voltage_align; % individual arm power (W)
totPower = sum(Power,2); % total power (W)

%% Aligning ground and wall distance data from log file to data from pi
groundDist_pi = rawData.data.px4FlowDist; % ARK flow distance from pi
N_dist_sensor = length(fieldnames(logData.distance_sensor)); % checking number of distance sensor topics
if N_dist_sensor == 1 % if only ARK flow sensor...
    ground_ID = string(fieldnames(logData.distance_sensor));
elseif N_dist_sensor == 2
    dist_fieldnames = string(fieldnames(logData.distance_sensor));
    for n = 1:N_dist_sensor
        dist_sensor_length(n) = length(logData.distance_sensor.(dist_fieldnames(n)).current_distance);
    end
    [~, ground_ID] = max(dist_sensor_length); ground_ID = dist_fieldnames(ground_ID);
    wall_ID = dist_fieldnames(ground_ID ~= dist_fieldnames);
    wall_ID = wall_ID(1);
end
groundDist_log = double(logData.distance_sensor.(string(ground_ID)).current_distance); % ARK flow distance from log file
tLog = seconds(logData.distance_sensor.(string(ground_ID)).timestamp);
tPi = rawData.times.tPX4flow;

% finding take off points:
figure
plot(groundDist_pi)
title("Ground distance data from pi")
figure
plot(groundDist_log)
title("Ground distance data from log")

% RUN UNTIL HERE AND THEN SELECT ALIGNING POINTS FOR PI AND LOG DATA BELOW

takeOffPt_pi = 419; takeOffPt_log = 735; % ENTER ALIGNING POINTS
t_delta = tPi(takeOffPt_pi) - tLog(takeOffPt_log); % time offset between pi and pixhawk

% sanity check plot
figure
plot(tLog+t_delta, groundDist_log)
hold on
plot(tPi, groundDist_pi)

% RUN UNTIL HERE TO MAKE SURE YOUR DATA IS PROPERLY ALIGNED (W/ SANITY
% CHECK PLOT)

tGroundDist = tLog + t_delta; % time vector aligned to pi
clear("corr_idx");
for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tGroundDist-tQrex(i)));
end
groundDist_align = groundDist_log(corr_idx,:);

if N_dist_sensor ~= 1 % if only ARK flow sensor...
    % aligning wall distance to thrust data:
    tWallDist = seconds(logData.distance_sensor.(string(wall_ID)).timestamp) + t_delta; % time vector aligned to pi
    wallDist_log = double(logData.distance_sensor.(string(wall_ID)).current_distance);

    clear("corr_idx"); % aligning wall distance
    for i = 1:length(tQrex)
        [~,corr_idx(i)] = min(abs(tWallDist-tQrex(i)));
    end
    wallDist_align = wallDist_log(corr_idx,:);
else
    wallDist_align = nan(length(groundDist_align),1);
end

%% Aligning local position (x, y) from ROS file to thrust data
pos = rawData.data.pos;
tPos = rawData.times.tPos;
clear("corr_idx");
for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tPos-tQrex(i)));
end
pos_align = pos(corr_idx,:);

%% Aligning angular rate data from ROS file to thrust data (NO ANGULAR VELOCITY DATA FOR ROS2)
angVel = rawData.data.angVel;
tVel = rawData.times.tVel;
clear("corr_idx");
for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tVel-tQrex(i)));
end
angVel_align = angVel(corr_idx,:);

%% Aligning linear velocity data from ROS file to thrust data
linVel = rawData.data.linVel;
clear("corr_idx");
for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tVel-tQrex(i)));
end
linVel_align = linVel(corr_idx,:);

%% Aligning linear velocity data from ROS file to thrust data
tIMU = rawData.times.tIMU;
quat = rawData.data.quat;
clear("corr_idx");
for i = 1:length(tQrex)
    [~,corr_idx(i)] = min(abs(tIMU-tQrex(i)));
end
quat_align = quat(corr_idx,:);

[yaw, pitch, roll] = quat2angle([quat_align(:,4), quat_align(:,1:3)]);
yaw = rad2deg(yaw); pitch = rad2deg(pitch); roll = rad2deg(roll);

%% Aligning, filtering, and de-biasing rotor pressure data from ROS file to thrust data
% (ONLY FOR TESTS WITH ROTOR PRESSURE PROBE)
if qPressure
    n_probes = 4; % ENTER NUMBER OF KIEL PROBES (1, 2, 3, or 4)

    Pressure = rawData.data.Pressure(:,1:n_probes*8);
    Pressure(:,1) = 1.213*Pressure(:,1); Pressure(:,2) = 1.109*Pressure(:,2);
    tPressure = rawData.times.tPressure;
    clear("corr_idx");
    for i = 1:length(tQrex)
        [~,corr_idx(i)] = min(abs(tPressure-tQrex(i)));
    end
    Pressure_align = Pressure(corr_idx,:); % aligning
    % Pressure_align(:,1) = medfilt1(Pressure_align(:,1),25);
    Pressure_align = medfilt1(Pressure_align,100); Pressure_align(1,:) = []; % filtering
    Pressure_align = Pressure_align - mean(Pressure_align(1:flight_pts(1),:)); % debiasing

    % plotting:
    close all
    for n = 1:n_probes
        figure(n)
        plot(Pressure_align(:,(n-1)*8+1), "k")
        hold on
        plot(Pressure_align(:,(n-1)*8+2:n*8))
        title(strcat("Probes ", string((n-1)*8+1), "-", string(n*8)))
        legend(string((n-1)*8+1:n*8),"Location","eastoutside")
        ylim([-200 200])
    end
end

%% Saving data to main .mat file
save(fullfile(rospath, matname))

%% Finding beginning and end of each test
figure % Thrust and ground distance plot (for ALL boundary test types)
plot(Thrust_scale)
hold on
plot(groundDist_align)
title("Ground distance and thrust measurements (find start and end points of test)")

if testType == 3 % if wall effect (hover) test...
    figure % Thrust and wall distance plot
    plot(Thrust_scale)
    hold on
    plot(wallDist_align)
    title("Wall distance and thrust measurements (find start and end points of test)")
end

if testType == 5 % if partial boundary (x/R constant) test...
    figure % x, y position and ground distance plot
    plot(pos_align(:,[1,2]))
    hold on
    plot(groundDist_align)
    title("Ground distance and local position (x,y) measurements (find start and end points of test)")
end

% RUN UNTIL HERE AND SELECT START AND END TEST POINTS (ENTER BELOW):

start_test = 9378; end_test = 50936; % ENTER START AND END TEST POINTS BASED ON PLOTS

%% Eliminating "dynamic" (or non-hover points) based on linear velocities
hover_threshold = 0.15;

hover_pts = find(abs(linVel_align(start_test:end_test,1)) < hover_threshold & ...
    abs(linVel_align(start_test:end_test,2)) < hover_threshold & ...
    abs(linVel_align(start_test:end_test,3)) < hover_threshold)+start_test;

groundDist_hover = groundDist_align(hover_pts); % ground distance measurements at steady-state hover
wallDist_hover = wallDist_align(hover_pts); % wall distance measurements at steady-state hover
Thrust_hover = Thrust_scale(hover_pts,:); % thrust measurements at steady-state hover
RPM_hover = RPM_scale(hover_pts,:); % rpm measurements at steady-state hover
CT_hover = CT(hover_pts,:); % CT measurements at steady-state hover
Power_hover = Power(hover_pts,:); % power measurements at steady-state hover
totPower_hover = totPower(hover_pts); % total power measurements at steady-state hover
angVel_hover = angVel_align(hover_pts,:); % angular rate measurements at steady-state hover

if exist("Pressure", "var") % checking if test included rotor pressure data...
    Pressure_hover = Pressure_align(hover_pts,:); % rotor pressure measurements at steady-state hover
end

%% Eliminating "dynamic" (or non-hover points) based on angular velocities
clear("hover_pts")

hover_threshold = 0.15;

hover_pts = find(abs(angVel_hover(:,1)) < hover_threshold & ...
    abs(angVel_hover(:,2)) < hover_threshold & ...
    abs(angVel_hover(:,3)) < hover_threshold);

groundDist_hover = groundDist_hover(hover_pts); % ground distance measurements at steady-state hover
wallDist_hover = wallDist_hover(hover_pts); % wall distance measurements at steady-state hover
Thrust_hover = Thrust_hover(hover_pts,:); % thrust measurements at steady-state hover
RPM_hover = RPM_hover(hover_pts,:); % rpm measurements at steady-state hover
CT_hover = CT_hover(hover_pts,:); % CT measurements at steady-state hover
Power_hover = Power_hover(hover_pts,:); % power measurements at steady-state hover
totPower_hover = totPower_hover(hover_pts); % total power measurements at steady-state hover
angVel_hover = angVel_hover(hover_pts,:); % angular rate measurements at steady-state hover

if exist("Pressure", "var") % checking if test included rotor pressure data...
    Pressure_hover = Pressure_hover(hover_pts,:); % rotor pressure measurements at steady-state hover
end

%% Taking out points not at correct ground height (only applicable to wall effect and partial boundary constant z/R tests)
if testType == 3 || testType == 4
    zR_target = 6; % ENTER TARGET HEIGHT ABOVE GROUND BEFORE RUNNING SECTION (z/R)
    zR_target = zR_target*R - 0.095; % converting to be wrt ARK flow sensor
    zR_threshold = 0.1; % below/above limit of 20cm from target altitude
    [idx_rm, ~] = find(groundDist_hover < zR_target - zR_threshold | groundDist_hover > zR_target + zR_threshold);

    % eliminating points below GE limit:
    groundDist_hover(idx_rm) = [];
    wallDist_hover(idx_rm) = [];
    Thrust_hover(idx_rm,:) = [];
    RPM_hover(idx_rm,:) = [];
    CT_hover(idx_rm,:) = [];
    Power_hover(idx_rm,:) = [];
    totPower_hover(idx_rm) = [];
    angVel_hover(idx_rm,:) = [];
    if exist("Pressure", "var") % checking if test included rotor pressure data...
        Pressure_hover(idx_rm,:) = [];
    end
end

%% Plotting Pressure in hover
close all
for n = 1:n_probes
    figure(n)
    plot(Pressure_hover(:,(n-1)*8+1), "k.")
    hold on
    plot(Pressure_hover(:,(n-1)*8+2:n*8),".")
    title(strcat("Probes ", string((n-1)*8+1), "-", string(n*8)))
    legend(string((n-1)*8+1:n*8),"Location","eastoutside")
    ylim([-200 200])
end

%% Selecting x/R windows (only applicable to partial boundary constant z/R tests)
if testType == 4
    figure
    plot(wallDist_hover, ".")
    title("Select x/R test windows")
    xlabel("Index"); ylabel("Wall Distance (m)")
    [xR_window,~] = getpts; xR_window = ceil(xR_window);
end

%% Sorting data based on distance (ground or wall depending on test type)
if testType == 1 || testType == 5 % if ground effect (hover) or partial boundary (x/R constant) test...
    [groundDist_sort, idx_sort] = sort(groundDist_hover); % sorted distances and corresponding vector indices (low --> high)
    wallDist_sort = wallDist_hover(idx_sort);
elseif testType == 3 % if wall effect (hover) test...
    [wallDist_sort, idx_sort] = sort(wallDist_hover); % sorted distances and corresponding vector indices (low --> high)
    groundDist_sort = groundDist_hover(idx_sort);
end
if testType ~= 4
    % sorting thrust, rpm, CT, power, and total power data:
    Thrust_sort = Thrust_hover(idx_sort,:);
    RPM_sort = RPM_hover(idx_sort,:);
    CT_sort = CT_hover(idx_sort,:);
    Power_sort = Power_hover(idx_sort,:);
    totPower_sort = totPower_hover(idx_sort);
    angVel_sort = angVel_hover(idx_sort,:);
    if exist("Pressure", "var") % checking if test included rotor pressure data...
        Pressure_sort = Pressure_hover(idx_sort,:);
    end
end

%% Normalizing wall and ground distances (where applicable)
if testType ~= 4
    groundDist_rotor = groundDist_sort + 0.095; % adding distance between rotor plane and ARK flow sensor (+ 9.5 cm)
    groundDist_norm = groundDist_rotor./R; % normalizing by rotor radius
end
if testType == 3
    wallDist_rotor = wallDist_sort - x - R; % adjusting to be distance between rotor tip and wall
    wallDist_norm = wallDist_rotor./R; % normalizing by rotor radius
end

%% Saving all data to .mat file
if testType == 1 % for ground effect (hover) test...
    savefilename = fullfile(rospath, matname);
    save(savefilename)
elseif testType == 3 || testType == 4 % for wall effect (hover) or for partial boundary (z/R constant) test...
    zR = (zR_target+0.095)/R;
    if testType == 4 % for partial boundary (z/R constant) test...
        % xR = [6, 4, 3, 2, 1]; % ENTER x/R points (in order) for test
    end
    savefilename = strcat(fullfile(rospath, matname),"_zR",string(zR));
elseif testType == 5 % for partial boundary (x/R constant) test...
    xR = -2; % ENTER x/R SETPOINT BEFORE RUNNING SECTION FOR testType = 5
    if xR < -2
        groundDist_norm = groundDist_norm + 6;
    end
    savefilename = strcat(fullfile(rospath, matname),"_xR",string(xR));
end

all_files = dir(rospath); 
if testType ~= 1
    if ~isfile(strcat(savefilename,".mat"))
        save(strcat(savefilename,".mat"))
    elseif ~isfile(strcat(savefilename,"_2.mat"))
        savefilename = strcat(fullfile(rospath, matname),"_zR",string(zR),"_2.mat");
        save(savefilename)
    else
        savefilename = strcat(fullfile(rospath, matname),"_zR",string(zR),"_3.mat");
        save(savefilename)
    end
end
