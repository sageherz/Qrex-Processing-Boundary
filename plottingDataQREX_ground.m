clear
close all
clc

%% Loading file containing processed data
qKiel = 1; % 1 = test w/ pressure measurements, 0 = test w/o pressure measurements

processed_filepath = "Ground\ALL_PROCESSED_TESTS";
if qKiel
    kiel_az = input("ENTER KIEL PROBE AZIMUTH (90 or 45):");
    processed_filepath = strcat(processed_filepath, "_P",string(kiel_az),"deg");
end
load(processed_filepath)

%% Finding total mean values and computing error
sigma = 2; % # of standard deviations for uncertainty estimation

N = sum(~isnan(groundDist_all),4);

groundDist.mean = mean(groundDist_all,4,"omitnan");
groundDist.std = std(groundDist_all,[],4,"omitnan");
groundDist.err = (sigma.*groundDist.std)./(N-1);

for ARM = 1:4
    fieldname = strcat("arm",string(ARM));

    Thrust.(fieldname).mean = mean(Thrust_all(:,:,ARM,:),4,"omitnan");
    Thrust.(fieldname).std = std(Thrust_all(:,:,ARM,:),[],4,"omitnan");
    Thrust.(fieldname).err = (sigma.*Thrust.(fieldname).std)./(N-1);

    RPM.(fieldname).mean = mean(RPM_all(:,:,ARM,:),4,"omitnan");
    RPM.(fieldname).std = std(RPM_all(:,:,ARM,:),[],4,"omitnan");
    RPM.(fieldname).err = (sigma.*RPM.(fieldname).std)./(N-1);

    CT.(fieldname).mean = mean(CT_all(:,:,ARM,:),4,"omitnan");
    CT.(fieldname).std = std(CT_all(:,:,ARM,:),[],4,"omitnan");
    CT.(fieldname).err = (sigma.*CT.(fieldname).std)./(N-1);

    Power.(fieldname).mean = mean(Power_all(:,:,ARM,:),4,"omitnan");
    Power.(fieldname).std = std(Power_all(:,:,ARM,:),[],4,"omitnan");
    Power.(fieldname).err = (sigma.*Power.(fieldname).std)./(N-1);
end

totPower.mean = mean(totPower_all,4,"omitnan");
totPower.std = std(totPower_all,[],4,"omitnan");
totPower.err = (sigma.*totPower.std)./(N-1);

if qKiel
    if kiel_az == 90
        Pressure_90.mean = mean(Pressure_all,4,"omitnan");
        Pressure_90.std = std(Pressure_all,[],4,"omitnan");
        Pressure_90.err = (sigma.*Pressure_90.std)./(N-1);

        Pressure_norm_90.mean = mean(Pressure_norm_all,4,"omitnan");
        Pressure_norm_90.std = std(Pressure_norm_all,[],4,"omitnan");
        Pressure_norm_90.err = (sigma.*Pressure_norm_90.std)./(N-1);

        Pressure_90.fluct = mean(Pressure_std_all,4,"omitnan");

    elseif kiel_az == 45
        Pressure_45.mean = mean(Pressure_all,4,"omitnan");
        Pressure_45.std = std(Pressure_all,[],4,"omitnan");
        Pressure_45.err = (sigma.*Pressure_45.std)./(N-1);

        Pressure_norm_45.mean = mean(Pressure_norm_all,4,"omitnan");
        Pressure_norm_45.std = std(Pressure_norm_all,[],4,"omitnan");
        Pressure_norm_45.err = (sigma.*Pressure_norm_45.std)./(N-1);

        Pressure_45.fluct = mean(Pressure_std_all,4,"omitnan");
    end
end

%% Plotting final, averaged data with error estimations (GROUND + WALL tests)
qNorm = 1; % 1 = normalize data w/ far boundary point, 0 = non-normalized

if qNorm
    for ARM = 1:4
    fieldname = strcat("arm",string(ARM));
    
    Thrust.(fieldname).err = Thrust.(fieldname).err./Thrust.(fieldname).mean(end);
    Thrust.(fieldname).mean = Thrust.(fieldname).mean./Thrust.(fieldname).mean(end);
    
    RPM.(fieldname).err = RPM.(fieldname).err./RPM.(fieldname).mean(end);
    RPM.(fieldname).mean = RPM.(fieldname).mean./RPM.(fieldname).mean(end);
    
    CT.(fieldname).err = CT.(fieldname).err./CT.(fieldname).mean(end);
    CT.(fieldname).mean = CT.(fieldname).mean./CT.(fieldname).mean(end);
    
    Power.(fieldname).err = Power.(fieldname).err./Power.(fieldname).mean(end);
    Power.(fieldname).mean = Power.(fieldname).mean./Power.(fieldname).mean(end);
    
    end
    
    totPower.err = totPower.err./totPower.mean(end);
    totPower.mean = totPower.mean./totPower.mean(end);
   

end

close all
font_size = 16;

figure % thrust
set(gcf,"Color","white")
hold on
for ARM = 1:4
    fieldname = strcat("arm",string(ARM));
    errorbar(groundDist.mean, Thrust.(fieldname).mean, Thrust.(fieldname).err,...
        Thrust.(fieldname).err, groundDist.err, groundDist.err,'o',...
        "LineWidth", 1.5)
end
xlabel("z/R","Interpreter","latex","FontSize",font_size)
title("Mean thrust as a function of ground distance")
legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location", "southeast")
ylabel("$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % rpm
set(gcf,"Color","white")
hold on
    for ARM = 1:4
        fieldname = strcat("arm",string(ARM));
        errorbar(groundDist.mean, RPM.(fieldname).mean, RPM.(fieldname).err,...
            RPM.(fieldname).err, groundDist.err, groundDist.err,'o',...
            "LineWidth", 1.5)
    end
    xlabel("z/R","Interpreter","latex","FontSize",font_size)
    title("Mean RPM as a function of ground distance")
    legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location", "southeast")
ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % CT
set(gcf,"Color","white")
hold on
    for ARM = 1:4
        fieldname = strcat("arm",string(ARM));
        errorbar(groundDist.mean, CT.(fieldname).mean, CT.(fieldname).err,...
            CT.(fieldname).err, groundDist.err, groundDist.err,'o',...
            "LineWidth", 1.5)
    end
    xlabel("z/R","Interpreter","latex","FontSize",font_size)
    title("Mean C_T as a function of ground distance")
    legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location", "southeast")
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % individual power
set(gcf,"Color","white")
hold on
    for ARM = 1:4
        fieldname = strcat("arm",string(ARM));
        errorbar(groundDist.mean, Power.(fieldname).mean, Power.(fieldname).err,...
            Power.(fieldname).err, groundDist.err, groundDist.err,'o',...
            "LineWidth", 1.5)
    end
    xlabel("z/R","Interpreter","latex","FontSize",font_size)
    title("Mean power as a function of ground distance")
    legend("Arm 1", "Arm 2", "Arm 3", "Arm 4", "Location", "southeast")
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % total power
set(gcf,"Color","white")
    errorbar(groundDist.mean, totPower.mean, totPower.err, totPower.err,...
        groundDist.err, groundDist.err,'o',"LineWidth", 1.5)
    xlabel("z/R")
    title("Mean total power as a function of ground distance")
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % total (combined) thrust
totThrust.mean = zeros(size(Thrust.arm1.mean));
totThrust.err = zeros(size(Thrust.arm1.err));
set(gcf,"Color","white")
    for ARM = 1:4
        fieldname = strcat("arm",string(ARM));
        totThrust.mean = totThrust.mean + Thrust.(fieldname).mean;
        totThrust.err = totThrust.err + Thrust.(fieldname).err;
    end
    totThrust.mean = totThrust.mean./4; totThrust.err = totThrust.err./4;

    errorbar(groundDist.mean, totThrust.mean, totThrust.err, totThrust.err,...
            groundDist.err, groundDist.err,'o',"LineWidth", 1.5)
    xlabel("z/R","FontSize",font_size,"Interpreter","latex")
    title("Mean total thrust as a function of ground distance")
ylabel("$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

figure % total (combined) ct
totCT.mean = zeros(size(CT.arm1.mean));
totCT.err = zeros(size(CT.arm1.err));
set(gcf,"Color","white")
    for ARM = 1:4
        fieldname = strcat("arm",string(ARM));
        totCT.mean = totCT.mean + CT.(fieldname).mean;
        totCT.err = totCT.err + CT.(fieldname).err;
    end
    totCT.mean = totCT.mean./4; totCT.err = totCT.err./4;

    errorbar(groundDist.mean, totCT.mean, totCT.err, totCT.err,...
            groundDist.err, groundDist.err,'o',"LineWidth", 1.5)
    xlabel("z/R","FontSize",font_size,"Interpreter","latex")
    title("Mean total C_T as a function of ground distance")
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
if qNorm
    ylim([0.75 1.05])
end

%% Plotting PRESSURE (non-normalized)
close all;

pressure_plot_colors = [0 0 0; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840; 0 0.4470 0.7410];
probe_locs = round(([5, 6.95, 8.8, 10.6, 12.6, 14.4, 16.3, 18.3])./19.05,2);
if kiel_az == 90
    probe_angles = [0, 90, 180, 270];
elseif kiel_az == 45
    probe_angles = [45, 135, 225, 315];
end

for n_probe = 1:4 % mean pressure ratio
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = Pressure_90.mean(:,1,(n_probe-1)*8+1);
        err_z = Pressure_90.err(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_45.mean(:,1,(n_probe-1)*8+1);
        err_z = Pressure_45.err(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = Pressure_90.mean(:,1,(n_probe-1)*8+i);
            err_z = Pressure_90.err(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_45.mean(:,1,(n_probe-1)*8+i);
            err_z = Pressure_45.err(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
        plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
            pressure_plot_colors(i,:),"LineWidth",0.5) % error bar
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("z/R","FontSize",14,"Interpreter","latex")
    zlabel("P (Pa)","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor pressure as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 110])
end

for n_probe = 1:4 % standard deviation (pressure fluctuations)
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = Pressure_90.fluct(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_45.fluct(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) 
    hold on
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = Pressure_90.fluct(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_45.fluct(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) 
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("z/R","FontSize",14,"Interpreter","latex")
    zlabel("P (Pa)","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor pressure fluctuations as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 50])
end

% Plotting PRESSURE (normalized)
for n_probe = 1:4 % mean normalized pressure ratio
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = Pressure_norm_90.mean(:,1,(n_probe-1)*8+1);
        err_z = Pressure_norm_90.err(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_norm_45.mean(:,1,(n_probe-1)*8+1);
        err_z = Pressure_norm_45.err(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = Pressure_norm_90.mean(:,1,(n_probe-1)*8+i);
            err_z = Pressure_norm_90.err(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_norm_45.mean(:,1,(n_probe-1)*8+i);
            err_z = Pressure_norm_45.err(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
        plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
            pressure_plot_colors(i,:),"LineWidth",0.5) % error bar
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("x/R","FontSize",14,"Interpreter","latex")
    zlabel("$\frac{p_{IGE}}{p_{OGE}}$","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor pressure ratio as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 1.5])
end

%% Plotting VELOCITY 
close all;

pressure_plot_colors = [0 0 0; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840; 0 0.4470 0.7410];
probe_locs = round(([5, 6.95, 8.8, 10.6, 12.6, 14.4, 16.3, 18.3])./19.05,2);
if kiel_az == 90
    probe_angles = [0, 90, 180, 270];
elseif kiel_az == 45
    probe_angles = [45, 135, 225, 315];
end

for n_probe = 1:4 % mean pressure ratio
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = v_90.mean(:,1,(n_probe-1)*8+1);
        err_z = v_90.err(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = v_45.mean(:,1,(n_probe-1)*8+1);
        err_z = v_45.err(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = v_90.mean(:,1,(n_probe-1)*8+i);
            err_z = v_90.err(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = v_45.mean(:,1,(n_probe-1)*8+i);
            err_z = v_45.err(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
        plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
            pressure_plot_colors(i,:),"LineWidth",0.5) % error bar
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("z/R","FontSize",14,"Interpreter","latex")
    zlabel("v (m/s)","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor velocity as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 15])
end

for n_probe = 1:4 % standard deviation (pressure fluctuations)
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = v_90.fluct(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = v_45.fluct(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) 
    hold on
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = v_90.fluct(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = v_45.fluct(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) 
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("z/R","FontSize",14,"Interpreter","latex")
    zlabel("v (m/s)","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor velocity fluctuations as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 10])
end

% Plotting PRESSURE (normalized)
for n_probe = 1:4 % mean normalized pressure ratio
    figure
    set(gcf,"Color","white")
    x = groundDist.mean;
    y = probe_locs(1)*ones(size(groundDist.mean));
    if kiel_az == 90
        z = v_norm_90.mean(:,1,(n_probe-1)*8+1);
        err_z = v_norm_90.err(:,1,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = v_norm_45.mean(:,1,(n_probe-1)*8+1);
        err_z = v_norm_45.err(:,1,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = groundDist.mean;
        y = probe_locs(i)*ones(size(groundDist.mean));
        if kiel_az == 90
            z = v_norm_90.mean(:,1,(n_probe-1)*8+i);
            err_z = v_norm_90.err(:,1,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = v_norm_45.mean(:,1,(n_probe-1)*8+i);
            err_z = v_norm_45.err(:,1,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
        plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
            pressure_plot_colors(i,:),"LineWidth",0.5) % error bar
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("x/R","FontSize",14,"Interpreter","latex")
    zlabel("$\frac{v_{IGE}}{v_{OGE}}$","FontSize",14,"Interpreter","latex",...
        "Rotation",0)
    title(strcat("Rotor velocity ratio as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 1.5])
end

%% Saving data (for pressure tests only rn)
save_filepath = "GROUND\Pressure_GE.mat";
if ~isfile(save_filepath)
    if kiel_az == 90
        groundDist_90 = groundDist;
        save(save_filepath, "probe_locs","groundDist_90", "Pressure_90",...
            "Pressure_norm_90")
    elseif kiel_az == 45
        groundDist_45 = groundDist;
        save(save_filepath, "probe_locs","groundDist_45", "Pressure_45",...
            "Pressure_norm_45")
    end
else
    if kiel_az == 90
        groundDist_90 = groundDist;
        save(save_filepath,"probe_locs","groundDist_90", "Pressure_90",...
            "Pressure_norm_90","-append")
    elseif kiel_az == 45
        groundDist_45 = groundDist;
        save(save_filepath,"probe_locs","groundDist_45", "Pressure_45",...
            "Pressure_norm_45", "-append")
    end
end
