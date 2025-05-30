clear
close all
clc

%% Loading processed data
[matfile, matpath] = uigetfile( '*.mat', 'Select mat File', 'multiselect', 'on' );
load(fullfile(matpath, matfile));

%% Finding averaging windows 
close all;

if testType ~= 4
    figure(1);
    set(gcf,"Color","white")
    if testType == 3
        plot(wallDist_norm,Thrust_sort,".")
        xlabel("x/R")
        title("Thrust as a function of wall distance")
    else
        plot(groundDist_norm,Thrust_sort,".")
        xlabel("z/R")
        title("Thrust as a function of ground distance")
    end
    ylabel("Thrust (N)")
    grid on
    legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location","best")
end

% RUN UNTIL HERE AND SELECT WINDOWS FOR AVERAGING

% ENTER WINDOWS FOR AVERAGING BELOW:
Dist_norm_low = [0.7, 0.9, 1.4, 2, 3, 4.25, 4.85, 6.05, 7.3, 10.1, 11.9]; % low x/R or z/R value for each averaging window
Dist_norm_high = [0.85, 1.3, 1.65, 2.3, 3.15, 4.35, 4.95, 6.3, 7.5, 10.3, 12.1]; % high x/R or z/R value for each averaging window

%% Averaging (and getting standard deviation) data within defined windows
if testType ~= 4
    for n = 1:length(Dist_norm_high) % for each averaging window
        if testType == 3 % for wall effect (hover) test...
            window = find(wallDist_norm > Dist_norm_low(n) & wallDist_norm < Dist_norm_high(n));
            wallDist_mean(n) = mean(wallDist_norm(window));
        else % for ground effect (hover) test and for partial boundary (x/R constant) test...
            window = find(groundDist_norm > Dist_norm_low(n) & groundDist_norm < Dist_norm_high(n));
            groundDist_mean(n) = mean(groundDist_norm(window));
        end

        Thrust_mean(n,:) = mean(Thrust_sort(window,:)); %#ok<*SAGROW>
        Thrust_std(n,:) = std(Thrust_sort(window,:));

        RPM_mean(n,:) = mean(RPM_sort(window,:));
        RPM_std(n,:) = std(RPM_sort(window,:));

        CT_mean(n,:) = mean(CT_sort(window,:));
        CT_std(n,:) = std(CT_sort(window,:));

        Power_mean(n,:) = mean(Power_sort(window,:));
        Power_std(n,:) = std(Power_sort(window,:));

        totPower_mean(n,:) = mean(totPower_sort(window));
        totPower_std(n,:) = std(totPower_sort(window));

        angVel_mean(n,:) = mean(angVel_sort(window,:));
        angVel_std(n,:) = std(angVel_sort(window,:));

        if exist("Pressure", "var") % checking if test included rotor pressure data...
            Pressure_mean(n,:) = mean(Pressure_sort(window,:));
            Pressure_std(n,:) = std(Pressure_sort(window,:));
        end
    end

else % for partial boundary (z/R constant) test...
    for n = 1:2:length(xR_window)-1 % for each averaging window
        groundDist_mean(ceil(n/2)) = mean(groundDist_hover(xR_window(n):xR_window(n+1)));
        wallDist_mean(ceil(n/2)) = mean(wallDist_hover(xR_window(n):xR_window(n+1)));

        Thrust_mean(ceil(n/2),:) = mean(Thrust_hover(xR_window(n):xR_window(n+1),:)); %#ok<*SAGROW>
        Thrust_std(ceil(n/2),:) = std(Thrust_hover(xR_window(n):xR_window(n+1),:));

        RPM_mean(ceil(n/2),:) = mean(RPM_hover(xR_window(n):xR_window(n+1),:));
        RPM_std(ceil(n/2),:) = std(RPM_hover(xR_window(n):xR_window(n+1),:));

        CT_mean(ceil(n/2),:) = mean(CT_hover(xR_window(n):xR_window(n+1),:));
        CT_std(ceil(n/2),:) = std(CT_hover(xR_window(n):xR_window(n+1),:));

        Power_mean(ceil(n/2),:) = mean(Power_hover(xR_window(n):xR_window(n+1),:));
        Power_std(ceil(n/2),:) = std(Power_hover(xR_window(n):xR_window(n+1),:));

        totPower_mean(ceil(n/2),:) = mean(totPower_hover(xR_window(n):xR_window(n+1)));
        totPower_std(ceil(n/2),:) = std(totPower_hover(xR_window(n):xR_window(n+1)));

        angVel_mean(ceil(n/2),:) = mean(angVel_hover(xR_window(n):xR_window(n+1),:));
        angVel_std(ceil(n/2),:) = std(angVel_hover(xR_window(n):xR_window(n+1),:));
    end

end

%% HOBE Points
HOBE_start = [102884, 148095]; HOBE_end = [104921, 149891];

for n = 1:length(HOBE_start)
    T_HOBE(n,:) = mean(Thrust_scale(HOBE_start(n):HOBE_end(n),:));
    RPM_HOBE(n,:) = mean(RPM_scale(HOBE_start(n):HOBE_end(n),:));
    CT_HOBE(n,:) = mean(CT(HOBE_start(n):HOBE_end(n),:));
    Power_HOBE(n,:) = mean(Power(HOBE_start(n):HOBE_end(n),:));
    totPower_HOBE(n,1) = mean(totPower(HOBE_start(n):HOBE_end(n),:));
    if qPressure == 1
        Pressure_HOBE(n,:) = mean(Pressure_align(HOBE_start(n):HOBE_end(n),:));
    end
end

%% Finding normalization point (farthest away from the boundary or hover out of ALL boundary effect - HOBE)
norm_point = "HOBE"; % options: "far boundary" or "HOBE"

if strcmp(norm_point, "far boundary") == 1 % normalizing by point farthest from boundary...
    Thrust_norm_pt = Thrust_mean(end,:);
    RPM_norm_pt = RPM_mean(end,:);
    CT_norm_pt = CT_mean(end,:);
    Power_norm_pt = Power_mean(end,:);
    totPower_norm_pt = totPower_mean(end);
    if exist("Pressure", "var") % checking if test included rotor pressure data...
        Pressure_norm_pt = Pressure_mean(end,:);
    end

elseif strcmp(norm_point, "HOBE") == 1
    Thrust_norm_pt = mean(T_HOBE, 1);
    RPM_norm_pt = mean(RPM_HOBE, 1);
    CT_norm_pt = mean(CT_HOBE, 1);
    Power_norm_pt = mean(Power_HOBE, 1);
    totPower_norm_pt = mean(totPower_HOBE, 1);
    if exist("Pressure", "var") % checking if test included rotor pressure data...
        Pressure_norm_pt = mean(Pressure_HOBE, 1);
    end
end

%% Normalizing data and finding standard deviation/error estimation
sigma = 2;

Thrust_norm = Thrust_mean./Thrust_norm_pt;
Thrust_err = sigma*(Thrust_std./Thrust_mean);

RPM_norm = RPM_mean./RPM_norm_pt;
RPM_err = sigma*(RPM_std./RPM_mean);

CT_norm = CT_mean./CT_norm_pt;
CT_err = sigma*(CT_std./CT_mean);

Power_norm = Power_mean./Power_norm_pt;
Power_err = sigma*(Power_std./Power_mean);

totPower_norm = totPower_mean./totPower_norm_pt;
totPower_err = sigma*(totPower_std./totPower_mean);

if exist("Pressure", "var") % checking if test included rotor pressure data...
    Pressure_norm = Pressure_mean./Pressure_norm_pt;
    Pressure_err = sigma*(Pressure_std./Pressure_mean);
end

%% Plotting normalized data wrt either ground or wall distance
close all

if testType == 1 || testType == 5
    xplot = groundDist_mean;
    x_label = "z/R";
elseif testType == 3
    xplot = wallDist_mean;
    x_label = "x/R";
elseif testType == 4
    xplot = xR;
    x_label = "x/R";
end

figure % normalized thrust
set(gcf,"Color","white")
errorbar(xplot, Thrust_norm, 2*(Thrust_std./Thrust_mean),".","MarkerSize",14)
xlabel(x_label)
ylabel("$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
title("Normalized thrust as a function of boundary distance")
grid on
legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location","best")
ylim([0.7 1.1])

figure % normalized RPM
set(gcf,"Color","white")
errorbar(xplot, RPM_norm, 2*(RPM_std./RPM_mean),".","MarkerSize",14)
xlabel(x_label)
ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
title("Normalized RPM as a function of boundary distance")
grid on
legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location","best")
ylim([0.7 1.1])

figure % normalized CT
set(gcf,"Color","white")
errorbar(xplot, CT_norm, 2*(CT_std./CT_mean),".","MarkerSize",14)
xlabel(x_label)
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
title("Normalized C_T as a function of boundary distance")
grid on
legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location","best")
ylim([0.7 1.1])

figure % normalized power
set(gcf,"Color","white")
errorbar(xplot, Power_norm, 2*(Power_std./Power_mean),".","MarkerSize",14)
xlabel(x_label)
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
title("Normalized power as a function of boundary distance")
grid on
legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location","best")
ylim([0.7 1.1])

figure % normalized total power
set(gcf,"Color","white")
errorbar(xplot, totPower_norm, 2*(totPower_std./totPower_mean),".","MarkerSize",14)
xlabel(x_label)
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
title("Normalized total power as a function of boundary distance")
grid on
ylim([0.7 1.1])

if exist("Pressure", "var") % checking if test included rotor pressure data...
    if ~exist("skew","var")
        if probe_angles(1) == 0
            Pressure_norm_deg = [Pressure_norm(:,25:32), Pressure_norm(:,1:24)];
            Pressure_mean_deg = [Pressure_mean(:,25:32), Pressure_mean(:,1:24)];
            Pressure_std_deg = [Pressure_std(:,25:32), Pressure_std(:,1:24)];
        elseif probe_angles(1) == 45 && skew == 0 || probe_angles(1) == 0 && skew == 90 % for kiel rakes at 45, 135, 225, and 315 degrees
            Pressure_norm_deg = Pressure_norm;
            Pressure_mean_deg = Pressure_mean;
            Pressure_std_deg = Pressure_std;
        end
    else

        if probe_angles(1) == 0 && skew == 0 % for kiel rakes at 0, 90, 180, and 270 degrees
            Pressure_norm_deg = [Pressure_norm(:,25:32), Pressure_norm(:,1:24)];
            Pressure_mean_deg = [Pressure_mean(:,25:32), Pressure_mean(:,1:24)];
            Pressure_std_deg = [Pressure_std(:,25:32), Pressure_std(:,1:24)];
        elseif probe_angles(1) == 45 && skew == 0 || probe_angles(1) == 0 && skew == 90 % for kiel rakes at 45, 135, 225, and 315 degrees
            Pressure_norm_deg = Pressure_norm;
            Pressure_mean_deg = Pressure_mean;
            Pressure_std_deg = Pressure_std;
        elseif probe_angles(1) == 45 && skew == 90 % for kiel rakes at 45, 135, 225, and 315 degrees
            Pressure_norm_deg = [Pressure_norm(:,9:32), Pressure_norm(:,1:8)];
            Pressure_mean_deg =  [Pressure_mean(:,9:32), Pressure_mean(:,1:8)];
            Pressure_std_deg =  [Pressure_std(:,9:32), Pressure_std(:,1:8)];
        end

    end
    
    for n = 1:n_probes
        figure(5+n) % normalized rotor pressure
        set(gcf,"Color","white")
        plot(xplot, Pressure_norm_deg(:,(n-1)*8+1), "k", "Marker", ".",...
                "MarkerSize",14,"LineStyle","none","LineWidth",1.5)
            hold on
        if testType == 1
            plot(xplot, Pressure_norm_deg(:,(n-1)*8+2:(n*8)-1), "Marker", ".",...
                "MarkerSize",14,"LineStyle","none","LineWidth",1.5)
            ylim([0.1 1.5])
        elseif testType == 3
            plot(xplot, Pressure_norm_deg(:,(n-1)*8+2:(n*8)-1), "Marker", ".",...
                "MarkerSize",14,"LineStyle","none","LineWidth",1.5)
            ylim([0.8 1.1])
        end
        xlabel(x_label)
        ylabel("p (Pa)","Interpreter","latex","Rotation",0,...
            "FontSize",8,"LineStyle","none","LineWidth",1.5)
        title(strcat("Rotor pressure ratio as a function of boundary distance (",...
            string(probe_angles(n))," deg)"))
        grid on
        % legend_strings = strcat(repmat("r/R = ", 1, 8), string(probe_locs));
        % legend(legend_strings, "Location","best")
    end
end

% figure 
% set(gcf,"Color","white")
% yyaxis left
% plot(xplot(1:8), Thrust_std(1:8,1), ".-","MarkerSize",14,"LineWidth",1.5)
% ylabel("T (N)","Interpreter","latex","Rotation",0,...
%     "FontSize",14)
% hold on
% yyaxis right
% plot(xplot(1:8), Power_std(1:8,1), ".-","MarkerSize",14,"LineWidth",1.5)
% ylabel("P (W)","Interpreter","latex","Rotation",0,...
%     "FontSize",14)
% xlabel(x_label)
% xlim([0 7])
% ax = gca; ax.FontSize = 14;
% grid on; grid minor;

%% 3D pressure plots
for n = 1:n_probes % mean pressure ratio
    figure
    set(gcf,"Color","white")
    plot3(xplot, probe_locs(1)*ones(1,length(xplot)), Pressure_mean_deg(:,(n-1)*8+1), "k", "Marker",...
        ".","MarkerSize",14,"LineWidth",1.5)
    hold on
   
    last_probe = 7;
    zlim([0 110])
    
    for i = 2:last_probe
        plot3(xplot, probe_locs(i)*ones(1,length(xplot)), Pressure_mean_deg(:,(n-1)*8+i), "Marker",...
            ".","MarkerSize",14,"LineWidth",1.5)
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel(x_label,"FontSize",14,"Interpreter","latex")
    zlabel("$\frac{p_{IGE}}{p_{OGE}}$","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor pressure as a function of boundary distance (",...
            string(probe_angles(n))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
end

for n = 1:n_probes % standard deviation (pressure fluctuations)
    figure
    set(gcf,"Color","white")
    plot3(xplot, probe_locs(1)*ones(1,length(xplot)), Pressure_std_deg(:,(n-1)*8+1), "k", "Marker",...
        ".","MarkerSize",14,"LineWidth",1.5)
    hold on
    for i = 2:7
        plot3(xplot, probe_locs(i)*ones(1,length(xplot)), Pressure_std_deg(:,(n-1)*8+i), "Marker",...
            ".","MarkerSize",14,"LineWidth",1.5)
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel(x_label,"FontSize",14,"Interpreter","latex")
    zlabel("$\sigma$ (Pa)","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    zlim([0 50])
    title(strcat("Rotor pressure fluctuations as a function of boundary distance (",...
            string(probe_angles(n))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
end

%% Converting ground distances to be wrt distance to rotor plane and normalizing (for partial boundary constant z/R test only)
if testType == 4
    groundDist_rotor = groundDist_mean + 0.095; % adding distance between rotor plane and ARK flow sensor (+ 9.5 cm)
    groundDist_norm = groundDist_rotor./R; % normalizing by rotor radius
    groundDist_mean = groundDist_norm; % resetting variable to remain consistent with past implementation
end

%% Saving all data to .mat file
save(fullfile(matpath, matfile));





