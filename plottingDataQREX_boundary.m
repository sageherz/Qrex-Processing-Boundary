clear
close all
clc

%% Selecting boundary test type to process all "good" files
processTestTypeList = {'Ground', 'Ground w Pressure Probes','Wall',...
    'Wall w Pressure Probes', 'Partial'};
[processTestType, ~] = listdlg('PromptString',{'Select test type.'},'SelectionMode','single',...
    'ListString',processTestTypeList);

if processTestType ~= 1 && processTestType ~= 2 
    skew = input("ENTER SKEW ANGLE (CCW = + and CW = -): ");
    if processTestType == 3 || processTestType == 4
        zR = input("ENTER z/R TARGET HEIGHT: ");
    end
end

if processTestType== 2 || processTestType == 4
    kiel_az = input("ENTER KIEL PROBE AZIMUTH (90 or 45):");
end

%% Getting all filepaths associated with desired test
path_allfiles = "\GOOD_BOUNDARY_TESTS.xlsx"; 
files_all = readtable(path_allfiles, "Sheet", string(processTestTypeList(processTestType)));

if processTestType ~= 1 && processTestType ~= 2 
    files_all = files_all(files_all.SkewAngle==skew, :);
    if processTestType == 3 || processTestType == 4
        zR_ext = strcat("zR",string(zR));
        files_all = files_all(contains(files_all.Extension, zR_ext), :);
    end
end

if processTestType == 2 || processTestType == 4
    files_all = files_all(contains(string(files_all.KielAzimuth), string(kiel_az)), :);
end

for n_test = 1:height(files_all)
    if ~isempty(files_all.Extension(n_test))
        filenames(n_test) = strcat(string(files_all.TestType(n_test)),"/",...
            string(files_all.Date(n_test)),"/", string(files_all.FlightNumber(n_test))...
            ,"/",string(files_all.Filename(n_test)),string(files_all.Extension(n_test)),".mat"); %#ok<*SAGROW>
    else
        filenames(n_test) = strcat(string(files_all.TestType(n_test)),"/",...
            string(files_all.Date(n_test)),"/", string(files_all.FlightNumber(n_test))...
            ,"/",string(files_all.Filename(n_test)),".mat");
    end
end

%% Combining data from all trials 
if processTestType ~= 5
    n_trials = length(filenames);
else
    n_trials = files_all.TrialNumber(end);
end

for n_test = 1:length(filenames)
    load(filenames(n_test))

    if processTestType == 1 || processTestType == 2 % if GROUND test...
        xR_target = 0; % setting null x/R target value
        zR_target = fliplr([14,12,10,7.4,6,5,4,3,2,1.5,1,0.75]); % target z/R test points for GROUND test
        [~,corr_idx] = min(abs(groundDist_mean-zR_target.')); % finding correlation indices between target and actual distance values
    elseif processTestType == 3 || processTestType == 4
        xR_target = fliplr([14,12,10,7.4,6,5,4,3,2,1.5,1,0.75,0.5]); % target x/R test points for WALL test
        zR_target = zR; % target z/R test point for WALL test
        [~,corr_idx] = min(abs(wallDist_mean-xR_target.'));
    elseif processTestType == 5
        xR_target = [6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6]; % target x/R test points for PARTIAL boundary test
        zR_target = fliplr([12,11,10,9,8,7,6,5,4,3,2,1]); % target z/R test points for PARTIAL boundary test
        if contains(files_all.Extension(n_test),"xR")
            [~,corr_idx_zR] = min(abs(groundDist_mean-zR_target.')); % finding correlation indices between target and actual distance values
            corr_idx_xR = find(xR == xR_target); 
        else
            corr_idx_zR = find(round(zR(1)) == zR_target); 
            [~,corr_idx_xR] = min(abs(xR-xR_target.'));
        end
    end

    % initializing data matrices (z/R, x/R, arm#/dimension, trial#):
    if n_test == 1
        if processTestType == 1 || processTestType == 2 || processTestType == 5
            groundDist_all = nan(length(zR_target),length(xR_target),1,n_trials);
        end
        if processTestType == 3 || processTestType == 4 || processTestType == 5
            wallDist_all = nan(length(zR_target),length(xR_target),1,n_trials);
        end
        if processTestType == 2 || processTestType == 4
            Pressure_all = nan(length(zR_target), length(xR_target), 32, n_trials);
            Pressure_norm_all = nan(length(zR_target), length(xR_target), 32, n_trials);
            Pressure_std_all = nan(length(zR_target), length(xR_target), 32, n_trials);
        
            % v_all = nan(length(zR_target), length(xR_target), 32, n_trials);
            % v_norm_all = nan(length(zR_target), length(xR_target), 32, n_trials);
            % v_std_all = nan(length(zR_target), length(xR_target), 32, n_trials);
        end
        Thrust_all = nan(length(zR_target),length(xR_target),4,n_trials);
        Thrust_std_all = nan(length(zR_target),length(xR_target),4,n_trials);
        RPM_all = nan(length(zR_target),length(xR_target),4,n_trials);
        CT_all = nan(length(zR_target),length(xR_target),4,n_trials);
        Power_all = nan(length(zR_target),length(xR_target),4,n_trials);
        Power_std_all = nan(length(zR_target),length(xR_target),4,n_trials);
        totPower_all = nan(length(zR_target),length(xR_target),1,n_trials);
        angVel_all = nan(length(zR_target),length(xR_target),3,n_trials);
    end

    % populating normalized data matrices:
    if processTestType == 1 || processTestType == 2
        groundDist_all(corr_idx,:,:,n_test) = groundDist_mean;
        Thrust_all(corr_idx,:,:,n_test) = Thrust_mean;
        Thrust_std_all(corr_idx,:,:,n_test) = Thrust_std;
        RPM_all(corr_idx,:,:,n_test) = RPM_mean;
        CT_all(corr_idx,:,:,n_test) = CT_mean;
        Power_all(corr_idx,:,:,n_test) = Power_mean;
        Power_std_all(corr_idx,:,:,n_test) = Power_std;
        totPower_all(corr_idx,:,:,n_test) = totPower_mean;
        angVel_all(corr_idx,:,:,n_test) = angVel_mean;
        if processTestType == 2
            Pressure_all(corr_idx,:,:,n_test) = Pressure_mean_deg;
            Pressure_norm_all(corr_idx,:,:,n_test) = Pressure_norm_deg;
            Pressure_std_all(corr_idx,:,:,n_test) = Pressure_std_deg;
        
            % v_all(corr_idx,:,:,n_test) = v_mean;
            % v_norm_all(corr_idx,:,:,n_test) = v_norm;
            % v_std_all(corr_idx,:,:,n_test) = v_std;
        end
    elseif processTestType == 3 || processTestType == 4
        wallDist_all(:,corr_idx,:,n_test) = wallDist_mean;
        Thrust_all(:,corr_idx,:,n_test) = Thrust_norm;
        RPM_all(:,corr_idx,:,n_test) = RPM_norm;
        CT_all(:,corr_idx,:,n_test) = CT_norm;
        Power_all(:,corr_idx,:,n_test) = Power_norm;
        totPower_all(:,corr_idx,:,n_test) = totPower_norm;
        % angVel_all(:,corr_idx,:,n_test) = angVel_norm;
        if processTestType == 4
            Pressure_all(:,corr_idx,:,n_test) = Pressure_mean_deg;
            Pressure_norm_all(:,corr_idx,:,n_test) = Pressure_norm_deg;
            Pressure_std_all(:,corr_idx,:,n_test) = Pressure_std_deg;
        end
    elseif processTestType == 5
        groundDist_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = groundDist_mean;
        Thrust_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = Thrust_mean;
        Thrust_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = Thrust_std;
        Thrust_std_all(Thrust_std_all(:) == 0) = nan;
        RPM_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = RPM_mean;
        RPM_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = RPM_std;
        RPM_std_all(RPM_std_all(:) == 0) = nan;
        CT_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = CT_mean;
        CT_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = CT_std;
        CT_std_all(CT_std_all(:) == 0) = nan;
        Power_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = Power_mean;
        Power_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = Power_std;
        Power_std_all(Power_std_all(:) == 0) = nan;
        totPower_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = totPower_mean;
        totPower_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = totPower_std;
        totPower_std_all(totPower_std_all(:) == 0) = nan;
        angVel_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = angVel_mean;
        angVel_std_all(corr_idx_zR,corr_idx_xR,:,files_all.TrialNumber(n_test)) = angVel_std;
        angVel_std_all(angVel_std_all(:) == 0) = nan;
    end
end

%% Saving all combined data
% % FOR GROUND TEST:
% save("Ground/ALL_PROCESSED_TESTS.mat","groundDist_all","Thrust_all","RPM_all",...
%     "CT_all","Power_all","totPower_all","angVel_all","processTestType","filenames")

% % FOR GROUND TEST w/ KIEL PROBES:
% save(strcat("Ground/ALL_PROCESSED_TESTS_P",string(kiel_az),"deg.mat"),...
%     "groundDist_all", "Thrust_all","RPM_all","CT_all","Power_all",...
%     "totPower_all","angVel_all","Pressure_all","Pressure_norm_all",...
%     "Pressure_std_all","processTestType",...
%     "filenames")

% % FOR WALL TEST:
% save(strcat("Wall\PROCESSED TESTS\Cross 0\ALL_PROCESSED_TESTS_zR",string(zR),".mat"),...
%     "wallDist_all","Thrust_all","RPM_all","CT_all","Power_all","totPower_all",...
%     "angVel_all","processTestType","zR","skew","filenames")

% % FOR WALL TEST w/ KIEL PROBES:
% save(strcat("Wall\PROCESSED TESTS\Cross ",string(skew),"\ALL_PROCESSED_TESTS_zR",string(zR),...
%     "_P",string(kiel_az),"deg.mat"),"wallDist_all","Thrust_all",...
%     "RPM_all","CT_all","Power_all","totPower_all","angVel_all",...
%     "Pressure_all","Pressure_norm_all","Pressure_std_all",...
%     "processTestType","zR","skew","filenames")

% FOR PARTIAL TEST:
save("Partial\PROCESSED TESTS\ALL_PROCESSED_TESTS_Cross0.mat","groundDist_all","Thrust_all","RPM_all",...
    "CT_all","Power_all","totPower_all","angVel_all","processTestType","filenames",...
    "xR_target","zR_target","skew","Thrust_std_all","RPM_std_all",...
    "CT_std_all","Power_std_all","totPower_std_all","angVel_std_all")


%% Defining fore/aft/side rotors based on skew angle (WALL or PARTIAL boundary test)
if processTestType ~= 1 && processTestType ~= 2
    if skew == 0
        fore_rotors = [1,3]; aft_rotors = [2,4];
    elseif skew == 45
        fore_rotors = 3; aft_rotors = 4; side_rotors = [1,2];
    elseif skew == 90
        fore_rotors = [1,4]; aft_rotors = [2,3];
    end
end

%% Comparison plotting for all trials (Ground + Wall)
close all

defaultColors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; ...
    0.6350, 0.0780, 0.1840];
all_marks = {'o','s','d','x','*','^','v','>','<','p','h'};

for n_test = 1:length(filenames)
    figure(1) % thrust comp plot
    set(gcf,"Color","white")
    hold on
    if processTestType == 1
        for ARM = 1:4
        plot(groundDist_all(:,:,1,n_test),Thrust_all(:,:,ARM,n_test),...
            "Marker",all_marks(n_test),"Color",defaultColors(ARM,:),"LineStyle","none")
        end
    elseif processTestType == 2
        plot(wallDist_all(:,:,1,n_test), mean(Thrust_all(:,:,fore_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 0], "LineStyle","none")
        plot(wallDist_all(:,:,1,n_test), mean(Thrust_all(:,:,aft_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 1], "LineStyle","none")
        if exist("side_rotors", "var")
            plot(wallDist_all(:,:,1,n_test), mean(Thrust_all(:,:,side_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[1, 0, 0], "LineStyle","none")
        end
    end

    figure(2) % rpm comp plot
    set(gcf,"Color","white")
    hold on
    if processTestType == 1
        for ARM = 1:4
        plot(groundDist_all(:,:,1,n_test),RPM_all(:,:,ARM,n_test),...
            "Marker",all_marks(n_test),"Color",defaultColors(ARM,:),"LineStyle","none")
        end
    elseif processTestType == 2
        plot(wallDist_all(:,:,1,n_test), mean(RPM_all(:,:,fore_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 0], "LineStyle","none")
        plot(wallDist_all(:,:,1,n_test), mean(RPM_all(:,:,aft_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 1], "LineStyle","none")
        if exist("side_rotors", "var")
            plot(wallDist_all(:,:,1,n_test), mean(RPM_all(:,:,side_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[1, 0, 0], "LineStyle","none")
        end
    end

    figure(3) % CT comp plot
    set(gcf,"Color","white")
    hold on
    if processTestType == 1
        for ARM = 1:4
        plot(groundDist_all(:,:,1,n_test),CT_all(:,:,ARM,n_test),...
            "Marker",all_marks(n_test),"Color",defaultColors(ARM,:),"LineStyle","none")
        end
    elseif processTestType == 2
        plot(wallDist_all(:,:,1,n_test), mean(CT_all(:,:,fore_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 0], "LineStyle","none")
        plot(wallDist_all(:,:,1,n_test), mean(CT_all(:,:,aft_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 1], "LineStyle","none")
        if exist("side_rotors", "var")
            plot(wallDist_all(:,:,1,n_test), mean(CT_all(:,:,side_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[1, 0, 0], "LineStyle","none")
        end
    end

    figure(4) % power comp plot
    set(gcf,"Color","white")
    hold on
    if processTestType == 1
        for ARM = 1:4
        plot(groundDist_all(:,:,1,n_test),Power_all(:,:,ARM,n_test),...
            "Marker",all_marks(n_test),"Color",defaultColors(ARM,:),"LineStyle","none")
        end
    elseif processTestType == 2
        plot(wallDist_all(:,:,1,n_test), mean(Power_all(:,:,fore_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 0], "LineStyle","none")
        plot(wallDist_all(:,:,1,n_test), mean(Power_all(:,:,aft_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[0, 0, 1], "LineStyle","none")
        if exist("side_rotors", "var")
            plot(wallDist_all(:,:,1,n_test), mean(Power_all(:,:,side_rotors,n_test),3,"omitnan"),...
            "Marker",all_marks(n_test),"Color",[1, 0, 0], "LineStyle","none")
        end
    end

    figure(5) % total power comp plot
    set(gcf,"Color","white")
    hold on
    if processTestType == 1
        plot(groundDist_all(:,:,1,n_test),totPower_all(:,:,1,n_test),"Marker",...
            all_marks(n_test),"Color",[0, 0, 0],"LineStyle","none")
    elseif processTestType == 2
        plot(wallDist_all(:,:,1,n_test), totPower_all(:,:,1,n_test),...
            "Marker",all_marks(n_test),"Color",[0, 0, 0], "LineStyle","none")
    end
 
end

% ADDING AXES TITLES, LEGENDS,...
if processTestType == 1
    x_label = "z/R";
    legend_strings = [];
    for n_test = 1:length(filenames)
        legend_strings = cat(2, legend_strings, [strcat("test ",string(n_test)," - Arm 1"),...
            strcat("test ",string(n_test)," - Arm 2"), strcat("test ",string(n_test)," - Arm 3"),...
            strcat("test ",string(n_test)," - Arm 4")]);
    end
    title_ext = "GE";
elseif processTestType == 2
    x_label = "x/R";
    legend_strings = [];
    for n_test = 1:length(filenames)
        legend_strings = cat(2, legend_strings, [strcat("test ",string(n_test)," - fore"), strcat("test ",string(n_test)," - aft")]);
        if exist("side_rotors", "var")
            legend_strings = cat(2, legend_strings, [strcat("test ",string(n_test)," - side")]);
        end
    end
    title_ext = strcat("WE ($\psi$ = ", string(skew), "$^{\circ}$, z/R = ", string(zR));
end

figure(1)
xlabel(x_label)
ylabel("$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
legend(legend_strings,"Location","south","NumColumns",length(filenames))
title(strcat("Thrust comparison - ", title_ext,")"),"Interpreter","latex")
grid on
ylim([0.7 1.1])

figure(2)
xlabel(x_label)
ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
legend(legend_strings,"Location","south","NumColumns",length(filenames))
title(strcat("RPM comparison - ", title_ext,")"),"Interpreter","latex")
grid on
ylim([0.7 1.1])

figure(3)
xlabel(x_label)
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
legend(legend_strings,"Location","south","NumColumns",length(filenames))
title(strcat("$C_T$ comparison - ", title_ext,")"),"Interpreter","latex")
grid on
ylim([0.7 1.1])


figure(4)
xlabel(x_label)
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
legend(legend_strings,"Location","south","NumColumns",length(filenames))
title(strcat("Power comparison - ", title_ext,")"),"Interpreter","latex")
grid on
ylim([0.7 1.2])

legend_strings = [];
for n_test = 1:length(filenames)
    legend_strings = cat(2, legend_strings, [strcat("test ",string(n_test))]);
end

figure(5)
xlabel(x_label)
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",14)
legend(legend_strings,"Location","southeast")
title(strcat("Total power comparison - ", title_ext,")"),"Interpreter","latex")
grid on
ylim([0.7 1.1])

%% Individual (per trial) THRUST plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % thrust surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Thrust_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Thrust (fore)")
view(2)

figure % thrust surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Thrust_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Thrust (aft)")
view(2)

figure % thrust regular plot (fore)
set(gcf,"Color","white")
plot(zR_target, mean(Thrust_all(:,1:6,fore_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(Thrust_all(:,7,fore_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(Thrust_all(:,8:13,fore_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("Thrust (N)","FontSize",font_size,"Interpreter","latex")
title("Thrust (fore)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % thrust regular plot (aft)
set(gcf,"Color","white")
plot(zR_target, mean(Thrust_all(:,1:6,aft_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(Thrust_all(:,7,aft_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(Thrust_all(:,8:13,aft_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("Thrust (N)","FontSize",font_size,"Interpreter","latex")
title("Thrust (aft)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % thrust std surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Thrust_std_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Thrust standard deviation (fore)")
view(2)

figure % thrust std surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Thrust_std_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Thrust standard deviation (aft)")
view(2)

%% Individual (per trial) RPM plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % RPM surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(RPM_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("RPM (fore)")
view(2)

figure % RPM surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(RPM_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("RPM (aft)")
view(2)

figure % RPM regular plot (fore)
set(gcf,"Color","white")
plot(zR_target, mean(RPM_all(:,1:6,fore_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(RPM_all(:,7,fore_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(RPM_all(:,8:13,fore_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("RPM","FontSize",font_size,"Interpreter","latex")
title("RPM (fore)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % RPM regular plot (aft)
set(gcf,"Color","white")
plot(zR_target, mean(RPM_all(:,1:6,aft_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(RPM_all(:,7,aft_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(RPM_all(:,8:13,aft_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("RPM","FontSize",font_size,"Interpreter","latex")
title("RPM (aft)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % RPM std surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(RPM_std_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("RPM standard deviation (fore)")
view(2)

figure % RPM std surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(RPM_std_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("RPM standard deviation (aft)")
view(2)

%% Individual (per trial) CT plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % CT surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(CT_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("C_T (fore)")
view(2)

figure % CT surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(CT_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("C_T (aft)")
view(2)

figure % CT regular plot (fore)
set(gcf,"Color","white")
plot(zR_target, mean(CT_all(:,1:6,fore_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(CT_all(:,7,fore_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(CT_all(:,8:13,fore_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("$C_T$","FontSize",font_size,"Interpreter","latex")
title("C_T (fore)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % CT regular plot (aft)
set(gcf,"Color","white")
plot(zR_target, mean(CT_all(:,1:6,aft_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(CT_all(:,7,aft_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(CT_all(:,8:13,aft_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("$C_T$","FontSize",font_size,"Interpreter","latex")
title("C_T (aft)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % CT std surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(CT_std_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("C_T standard deviation (fore)")
view(2)

figure % CT std surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(CT_std_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("C_T standard deviation (aft)")
view(2)

%% Individual (per trial) POWER plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % power surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Power_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Power (fore)")
view(2)

figure % power surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Power_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Power (aft)")
view(2)

figure % power regular plot (fore)
set(gcf,"Color","white")
plot(zR_target, mean(Power_all(:,1:6,fore_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(Power_all(:,7,fore_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(Power_all(:,8:13,fore_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("Power (W)","FontSize",font_size,"Interpreter","latex")
title("Power (fore)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % power regular plot (aft)
set(gcf,"Color","white")
plot(zR_target, mean(Power_all(:,1:6,aft_rotors,trial_num),3), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, mean(Power_all(:,7,aft_rotors,trial_num),3), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, mean(Power_all(:,8:13,aft_rotors,trial_num),3), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("Power (W)","FontSize",font_size,"Interpreter","latex")
title("Power (aft)")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % power std surface plot (fore)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Power_std_all(:,:,fore_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Power standard deviation (fore)")
view(2)

figure % power std surface plot (aft)
set(gcf,"Color","white")
contourf(xR_target, zR_target, mean(Power_std_all(:,:,aft_rotors,trial_num),3),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Power standard deviation (aft)")
view(2)

%% Individual (per trial) TOTAL POWER plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % total power surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, totPower_all(:,:,1,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Total Power")
view(2)

figure % total power regular plot 
set(gcf,"Color","white")
plot(zR_target, totPower_all(:,1:6,1,trial_num), "o", "MarkerSize", 7, "LineWidth",1.5)
hold on
plot(zR_target, totPower_all(:,7,1,trial_num), "x", "MarkerSize", 7, "LineWidth",1.5)
plot(zR_target, totPower_all(:,8:13,1,trial_num), "^", "MarkerSize", 7, "LineWidth",1.5)
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
ylabel("Power (W)","FontSize",font_size,"Interpreter","latex")
title("Total Power")
legend("x/R = 6", "x/R = 5", "x/R = 4", "x/R = 3", "x/R = 2", "x/R = 1", ...
    "x/R = 0", "x/R = -1", "x/R = -2", "x/R = -3", "x/R = -4", "x/R = -5",...
    "x/R = -6")

figure % total power std surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, totPower_std_all(:,:,1,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Total Power standard deviation")
view(2)

%% Individual (per trial) ROLL RATE plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % roll rate surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, angVel_all(:,:,1,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{p_{IBE}}{p_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Roll rate")
view(2)

figure % roll rate std surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, angVel_std_all(:,:,1,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{p_{IBE}}{p_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Roll rate standard deviation")
view(2)

%% Individual (per trial) PITCH RATE plots (Partial boundary tests ONLY)
close all
addpath("library\")
trial_num = 1; % ENTER TRIAL NUMBER
font_size = 16;
txt = "BOX";

figure % pitch rate surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, angVel_all(:,:,2,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{q_{IBE}}{q_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Pitch rate")
view(2)

figure % pitch rate std surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, angVel_std_all(:,:,2,trial_num),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(50))
cbar = colorbar;
xlabel("x/R","FontSize",font_size,"Interpreter","latex")
xlim([-6 6])
ylabel("z/R","FontSize",font_size,"Interpreter","latex")
ylabel(cbar,"$\frac{q_{IBE}}{q_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
title("Pitch rate standard deviation")
view(2)













