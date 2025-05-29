clear
close all
clc

%% Loading file containing processed data
skew = input("ENTER SKEW ANGLE (CCW = + and CW = -): ");

processed_filepath = "Partial\PROCESSED TESTS\";
if skew == 45
    fore_rotors = 3; aft_rotors = 4; side_rotors = [1,2];
elseif skew == 0 || skew == 90
    processed_filepath = strcat(processed_filepath,"\ALL_PROCESSED_TESTS_Cross",...
        string(skew));
    fore_rotors = [1,3]; aft_rotors = [2,4];
end

%% Averaging wrt to fore/aft/side rotors
load(processed_filepath)

groundDist.mean = mean(groundDist_all,4,"omitnan");
% groundDist.std = mean(groundDist_std_all,4,"omitnan");

Thrust.fore.mean = mean(mean(Thrust_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
Thrust.fore.std = mean(mean(Thrust_std_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");

Thrust.aft.mean = mean(mean(Thrust_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
Thrust.aft.std = mean(mean(Thrust_std_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");

RPM.fore.mean = mean(mean(RPM_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
RPM.fore.std = mean(mean(RPM_std_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");

RPM.aft.mean = mean(mean(RPM_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
RPM.aft.std = mean(mean(RPM_std_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");

CT.fore.mean = mean(mean(CT_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
CT.fore.std = mean(mean(CT_std_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");

CT.aft.mean = mean(mean(CT_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
CT.aft.std = mean(mean(CT_std_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");

Power.fore.mean = mean(mean(Power_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
Power.fore.std = mean(mean(Power_std_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");

Power.aft.mean = mean(mean(Power_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
Power.aft.std = mean(mean(Power_std_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");

if exist("side_rotors","var")
    Thrust.side.mean = mean(mean(Thrust_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
    Thrust.side.std = mean(mean(Thrust_std_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");

    RPM.side.mean = mean(mean(RPM_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
    RPM.side.std = mean(mean(RPM_std_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");

    CT.side.mean = mean(mean(CT_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
    CT.side.std = mean(mean(CT_std_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");

    Power.side.mean = mean(mean(Power_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
    Power.side.std = mean(mean(Power_std_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
end

totPower.mean = mean(totPower_all,4,"omitnan");
totPower.std = mean(totPower_std_all,4,"omitnan");

rollRate.std = mean(angVel_std_all(:,:,1,:),4,"omitnan");
pitchRate.std = mean(angVel_std_all(:,:,2,:),4,"omitnan");

%% POWER PLOTS
font_size = 22;
txt = "BOX";
addpath("library\")

figure % fore rotor power surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, Power.fore.mean./Power.fore.mean(end,1),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [0.8 1.15]; clim(clims)
xlabel("x/R")
xlim([-6 6])
ylabel("z/R")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

figure % aft rotor power surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, Power.aft.mean./Power.aft.mean(end,1),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [0.8 1.15]; clim(clims)
xlabel("x/R")
xlim([-6 6])
ylabel("z/R")
ylabel(cbar,"$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

figure % fore rotor power std surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, 2.*Power.fore.std, "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = 2.*[1.5 5]; clim(clims)
xlabel("x/R")
xlim([-6 6])
ylabel("z/R")
ylabel(cbar,"$P_{IBE}$ (W)","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

figure % aft rotor power std surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, 2.*Power.aft.std, "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = 2.*[1.5 5]; clim(clims)
xlabel("x/R")
xlim([-6 6])
ylabel("z/R")
ylabel(cbar,"$P_{IBE}$ (W)","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

%% CT PLOTS
font_size = 16;
txt = "BOX";

figure % fore rotor CT surface plot
set(gcf,"Color","white")
contourf(xR_target, zR_target, CT.fore.mean./CT.fore.mean(end,1),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [0.85 1.05]; clim(clims)
xlabel("x/R")
xlim([-6 10])
ylabel("z/R")
ylabel(cbar,"$\frac{C_{T{IBE}}}{C_{T{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

figure % aft rotor CT surface plot
aft_offset = 19/7.5; 
set(gcf,"Color","white")
contourf(xR_target + aft_offset + 1, zR_target, CT.aft.mean./CT.aft.mean(end,1),...
    "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [0.85 1.05]; clim(clims)
xlabel("x/R")
xlim([-6 10])
ylabel("z/R")
ylabel(cbar,"$\frac{C_{T{IBE}}}{C_{T{OBE}}}$","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

%% ANGULAR RATE PLOTS
font_size = 22;
txt = "BOX";

figure % roll rate std surface plot
set(gcf,"Color","white")
contourf(xR_target + 1, zR_target, 2*rad2deg(rollRate.std), "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [3 14]; clim(clims)
xlabel("x/R")
xlim([-6 10])
ylabel("z/R")
ylabel(cbar,"p ($^{\circ}$/s)","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

figure % pitch rate std surface plot
set(gcf,"Color","white")
contourf(xR_target + 1, zR_target, 2*rad2deg(pitchRate.std), "ShowText","on")
hold on
plot([-6, 0], [6, 6], "k-", "LineWidth", 3.5)
plot([0, 0], [1, 6], "k-", "LineWidth", 3.5)
colormap(redblue(100))
cbar = colorbar;
clims = [3 14]; clim(clims)
xlabel("x/R")
xlim([-6 10])
ylabel("z/R")
ylabel(cbar,"q ($^{\circ}$/s)","Interpreter","latex","Rotation",0);
ax = gca;
ax.FontSize = font_size;
text(-3.75, 3.5, txt,"FontSize",font_size,"FontWeight","bold")
view(2)
pbaspect([1.25 1 1])

