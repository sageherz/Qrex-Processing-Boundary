clear
close all
clc

font_size = 16;
default_colors = [0 0 0; 1 0 1; 0 0 1; 1 0 0];

%% Loading file containing processed data
skew = input("ENTER SKEW ANGLE (CCW = + and CW = -): ");
zR_all = input("ENTER z/R TARGET HEIGHT(S): ");
qKiel = 0; % 0 = no Kiel pressure data, 1 = w/ Kiel pressure data

processed_filepath = "Wall\PROCESSED TESTS\";
if skew == 45
    processed_filepath = strcat(processed_filepath,"Plus ",string(skew),...
        "\ALL_PROCESSED_TESTS_zR");
    fore_rotors = 3; aft_rotors = 4; side_rotors = [1,2];
elseif skew == 0 || skew == 90
    processed_filepath = strcat(processed_filepath,"Cross ",string(skew),...
        "\ALL_PROCESSED_TESTS_zR");
    fore_rotors = [1,3]; aft_rotors = [2,4];
end
processed_filepath = repmat(processed_filepath,1,length(zR_all));
processed_filepath = strcat(processed_filepath, string(zR_all));

if qKiel
    kiel_az = input("ENTER KIEL PROBE AZIMUTH (90 or 45):");
    processed_filepath = strcat(processed_filepath, "_P",string(kiel_az),"deg");
end

% processed_filepath = strcat(processed_filepath, repmat("_VFS81",1,4));

%% Averaging wrt to fore/aft/side rotors
for n = 1:length(processed_filepath) % for all z/R heights...
    load(processed_filepath(n))

    sigma = 2; % # of standard deviations for uncertainty estimation
    N = sum(~isnan(wallDist_all),4);

    wallDist.mean(n,:) = mean(wallDist_all,4,"omitnan");
    wallDist.std(n,:) = std(wallDist_all,[],4,"omitnan");
    wallDist.err(n,:) = (sigma.*wallDist.std(n,:))./(N-1);

    Thrust.fore.mean(n,:) = mean(mean(Thrust_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
    Thrust.fore.std(n,:) = mean(std(Thrust_all(:,:,fore_rotors,:),[],4,"omitnan"),3,"omitnan");
    Thrust.fore.err(n,:) = (sigma.*Thrust.fore.std(n,:))./(N-1);

    Thrust.aft.mean(n,:) = mean(mean(Thrust_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
    Thrust.aft.std(n,:) = mean(std(Thrust_all(:,:,aft_rotors,:),[],4,"omitnan"),3,"omitnan");
    Thrust.aft.err(n,:) = (sigma.*Thrust.aft.std(n,:))./(N-1);

    RPM.fore.mean(n,:) = mean(mean(RPM_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
    RPM.fore.std(n,:) = mean(std(RPM_all(:,:,fore_rotors,:),[],4,"omitnan"),3,"omitnan");
    RPM.fore.err(n,:) = (sigma.*RPM.fore.std(n,:))./(N-1);

    RPM.aft.mean(n,:) = mean(mean(RPM_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
    RPM.aft.std(n,:) = mean(std(RPM_all(:,:,aft_rotors,:),[],4,"omitnan"),3,"omitnan");
    RPM.aft.err(n,:) = (sigma.*RPM.aft.std(n,:))./(N-1);

    CT.fore.mean(n,:) = mean(mean(CT_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
    CT.fore.std(n,:) = mean(std(CT_all(:,:,fore_rotors,:),[],4,"omitnan"),3,"omitnan");
    CT.fore.err(n,:) = (sigma.*CT.fore.std(n,:))./(N-1);

    CT.aft.mean(n,:) = mean(mean(CT_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
    CT.aft.std(n,:) = mean(std(CT_all(:,:,aft_rotors,:),[],4,"omitnan"),3,"omitnan");
    CT.aft.err(n,:) = (sigma.*CT.aft.std(n,:))./(N-1);

    Power.fore.mean(n,:) = mean(mean(Power_all(:,:,fore_rotors,:),4,"omitnan"),3,"omitnan");
    Power.fore.std(n,:) = mean(std(Power_all(:,:,fore_rotors,:),[],4,"omitnan"),3,"omitnan");
    Power.fore.err(n,:) = (sigma.*Power.fore.std(n,:))./(N-1);

    Power.aft.mean(n,:) = mean(mean(Power_all(:,:,aft_rotors,:),4,"omitnan"),3,"omitnan");
    Power.aft.std(n,:) = mean(std(Power_all(:,:,aft_rotors,:),[],4,"omitnan"),3,"omitnan");
    Power.aft.err(n,:) = (sigma.*Power.aft.std(n,:))./(N-1);

    if exist("side_rotors","var")
        Thrust.side.mean(n,:) = mean(mean(Thrust_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
        Thrust.side.std(n,:) = mean(std(Thrust_all(:,:,side_rotors,:),[],4,"omitnan"),3,"omitnan");
        Thrust.side.err(n,:) = (sigma.*Thrust.side.std(n,:))./(N-1);

        RPM.side.mean(n,:) = mean(mean(RPM_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
        RPM.side.std(n,:) = mean(std(RPM_all(:,:,side_rotors,:),[],4,"omitnan"),3,"omitnan");
        RPM.side.err(n,:) = (sigma.*RPM.side.std(n,:))./(N-1);

        CT.side.mean(n,:) = mean(mean(CT_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
        CT.side.std(n,:) = mean(std(CT_all(:,:,side_rotors,:),[],4,"omitnan"),3,"omitnan");
        CT.side.err(n,:) = (sigma.*CT.side.std(n,:))./(N-1);

        Power.side.mean(n,:) = mean(mean(Power_all(:,:,side_rotors,:),4,"omitnan"),3,"omitnan");
        Power.side.std(n,:) = mean(std(Power_all(:,:,side_rotors,:),[],4,"omitnan"),3,"omitnan");
        Power.side.err(n,:) = (sigma.*Power.side.std(n,:))./(N-1);
    end

    totPower.mean(n,:) = mean(totPower_all,4,"omitnan");
    totPower.std(n,:) = std(totPower_all,[],4,"omitnan");
    totPower.err(n,:) = (sigma.*totPower.std(n,:))./(N-1);

    if qKiel
        if kiel_az == 90
            Pressure_90.mean(n,:,:) = mean(Pressure_all,4,"omitnan");
            Pressure_90.std(n,:,:) = std(Pressure_all,[],4,"omitnan");
            Pressure_90.err(n,:,:) = (sigma.*Pressure_90.std(n,:,:))./(N-1);

            Pressure_norm_90.mean(n,:,:) = mean(Pressure_norm_all,4,"omitnan");
            Pressure_norm_90.std(n,:,:) = std(Pressure_norm_all,[],4,"omitnan");
            Pressure_norm_90.err(n,:,:) = (sigma.*Pressure_norm_90.std(n,:,:))./(N-1);
        
            Pressure_90.fluct(n,:,:) = mean(Pressure_std_all,4,"omitnan");
        
        elseif kiel_az == 45
            Pressure_45.mean(n,:,:) = mean(Pressure_all,4,"omitnan");
            Pressure_45.std(n,:,:) = std(Pressure_all,[],4,"omitnan");
            Pressure_45.err(n,:,:) = (sigma.*Pressure_45.std(n,:,:))./(N-1);

            Pressure_norm_45.mean(n,:,:) = mean(Pressure_norm_all,4,"omitnan");
            Pressure_norm_45.std(n,:,:) = std(Pressure_norm_all,[],4,"omitnan");
            Pressure_norm_45.err(n,:,:) = (sigma.*Pressure_norm_45.std(n,:,:))./(N-1);

            Pressure_45.fluct(n,:,:) = mean(Pressure_std_all,4,"omitnan");

        end
    end
end

%% Normalizing (at z/R = 6, x/R = 14)
q_norm = 0;

if q_norm
    Thrust_norm.fore.mean = Thrust.fore.mean./Thrust.fore.mean(1,end);
    Thrust_norm.fore.err = Thrust.fore.err./Thrust.fore.mean(1,end); 
    Thrust_norm.aft.mean = Thrust.aft.mean./Thrust.aft.mean(1,end);
    Thrust_norm.aft.err = Thrust.aft.err./Thrust.aft.mean(1,end);

    RPM_norm.fore.mean = RPM.fore.mean./RPM.fore.mean(1,end);
    RPM_norm.fore.err = RPM.fore.err./RPM.fore.mean(1,end); 
    RPM_norm.aft.mean = RPM.aft.mean./RPM.aft.mean(1,end);
    RPM_norm.aft.err = RPM.aft.err./RPM.aft.mean(1,end);

    CT_norm.fore.mean = CT.fore.mean./CT.fore.mean(1,end);
    CT_norm.fore.err = CT.fore.err./CT.fore.mean(1,end); 
    CT_norm.aft.mean = CT.aft.mean./CT.aft.mean(1,end);
    CT_norm.aft.err = CT.aft.err./CT.aft.mean(1,end);

    Power_norm.fore.mean = Power.fore.mean./Power.fore.mean(1,end);
    Power_norm.fore.err = Power.fore.err./Power.fore.mean(1,end); 
    Power_norm.aft.mean = Power.aft.mean./Power.aft.mean(1,end);
    Power_norm.aft.err = Power.aft.err./Power.aft.mean(1,end);

    totPower_norm.mean = totPower.mean./totPower.mean(1,end);
    totPower_norm.err = totPower.err./totPower.mean(1,end);

    if exist("side_rotors", "var")
        Thrust_norm.side.mean = Thrust.side.mean./Thrust.side.mean(1,end);
        Thrust_norm.side.err = Thrust.side.err./Thrust.side.mean(1,end);

        RPM_norm.side.mean = RPM.side.mean./RPM.side.mean(1,end);
        RPM_norm.side.err = RPM.side.err./RPM.side.mean(1,end);

        CT_norm.side.mean = CT.side.mean./CT.side.mean(1,end);
        CT_norm.side.err = CT.side.err./CT.side.mean(1,end);

        Power_norm.side.mean = Power.side.mean./Power.side.mean(1,end);
        Power_norm.side.err = Power.side.err./Power.side.mean(1,end);
    end 
end

%% Plotting THRUST and TOTAL POWER (all) with or without error estimations 
q_error = 0; % ENTER WHETHER TO INCLUDE ERROR ESTIMATION (1 = yes, 0 = no)

% defining legend:
legend_strings = strcat(repmat("z/R = ", 1, length(zR_all)),string(zR_all));
if zR_all(1) > zR_all(end) % checking if z/R is in ascending or descending order...
    legend_strings = sort(repmat(legend_strings, 1, 4-length(fore_rotors)),"descend");
else
    legend_strings = sort(repmat(legend_strings, 1, 4-length(fore_rotors)),"ascend");
end
if exist("side_rotors","var") % if plus config...
    legend_strings = strcat(legend_strings, repmat([" (fore)", " (aft)", " (side)"],...
        1,length(zR_all)));
else % else cross config...
    legend_strings = strcat(legend_strings, repmat([" (fore)", " (aft)"],...
        1,length(zR_all)));
end

figure % thrust
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        if q_norm == 1
            errorbar(wallDist.mean(n,:), Thrust_norm.fore.mean(n,:), Thrust_norm.fore.err(n,:),...
                Thrust_norm.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
            errorbar(wallDist.mean(n,:), Thrust_norm.aft.mean(n,:), Thrust_norm.aft.err(n,:),...
                Thrust_norm.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
            if exist("side_rotors","var")
                errorbar(wallDist.mean(n,:), Thrust_norm.side.mean(n,:), Thrust_norm.side.err(n,:),...
                    Thrust_norm.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                    "LineWidth", 1.5,"Color",default_colors(n,:))
            end
        else
            errorbar(wallDist.mean(n,:), Thrust.fore.mean(n,:), Thrust.fore.err(n,:),...
                Thrust.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
            errorbar(wallDist.mean(n,:), Thrust.aft.mean(n,:), Thrust.aft.err(n,:),...
                Thrust.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
            if exist("side_rotors","var")
                errorbar(wallDist.mean(n,:), Thrust.side.mean(n,:), Thrust.side.err(n,:),...
                    Thrust.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                    "LineWidth", 1.5,"Color",default_colors(n,:))
            end
        end
    else % excluding error bars...
        if q_norm == 1
            plot(wallDist.mean(n,:), Thrust_norm.fore.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
            plot(wallDist.mean(n,:), Thrust_norm.aft.mean(n,:),'^',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
            if exist("side_rotors","var")
                plot(wallDist.mean(n,:), Thrust_norm.side.mean(n,:),'s',"LineWidth",...
                    1.5,"Color",default_colors(n,:))
            end
        else
            plot(wallDist.mean(n,:), Thrust.fore.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
            plot(wallDist.mean(n,:), Thrust.aft.mean(n,:),'^',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
            if exist("side_rotors","var")
                plot(wallDist.mean(n,:), Thrust.side.mean(n,:),'s',"LineWidth",...
                    1.5,"Color",default_colors(n,:))
            end
        end
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean thrust as a function of wall distance")
legend(legend_strings, "Location", "south","NumColumns",length(zR_all))
ylabel("$\frac{T_{IBE}}{T_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

figure % total power
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        if q_norm == 1
            errorbar(wallDist.mean(n,:), totPower_norm.mean(n,:), totPower_norm.err(n,:),...
                totPower_norm.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        else
            errorbar(wallDist.mean(n,:), totPower.mean(n,:), totPower.err(n,:),...
                totPower.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        end
    else % excluding error bars...
        if q_norm == 1
            plot(wallDist.mean(n,:), totPower_norm.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        else
            plot(wallDist.mean(n,:), totPower.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        end
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean total power as a function of wall distance")
legend(strcat(repmat("z/R = ", 1, length(zR_all)),string(zR_all)),...
    "Location", "south","NumColumns",length(zR_all))
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

%% Plotting RPM (fore/aft/side separated) with or without error estimations 
figure % RPM (fore)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        if q_norm == 1
            errorbar(wallDist.mean(n,:), RPM_norm.fore.mean(n,:), RPM_norm.fore.err(n,:),...
                RPM_norm.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        else
            errorbar(wallDist.mean(n,:), RPM.fore.mean(n,:), RPM.fore.err(n,:),...
                RPM.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        end
    else % excluding error bars...
        if q_norm == 1
            plot(wallDist.mean(n,:), RPM_norm.fore.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        else
            plot(wallDist.mean(n,:), RPM.fore.mean(n,:),'o',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        end
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean RPM (fore) as a function of wall distance")
legend(legend_strings(1:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

figure % RPM (aft)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        if q_norm == 1
            errorbar(wallDist.mean(n,:), RPM_norm.aft.mean(n,:), RPM_norm.aft.err(n,:),...
                RPM_norm.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        else
            errorbar(wallDist.mean(n,:), RPM.aft.mean(n,:), RPM.aft.err(n,:),...
                RPM.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        end
        
    else % excluding error bars...
        if q_norm == 1
            plot(wallDist.mean(n,:), RPM_norm.aft.mean(n,:),'^',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        else
            plot(wallDist.mean(n,:), RPM.aft.mean(n,:),'^',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        end
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean RPM (aft) as a function of wall distance")
legend(legend_strings(2:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

if exist("side_rotors","var")
    figure % RPM (side)
    set(gcf,"Color","white")
    hold on
    for n = 1:length(processed_filepath)
        if q_error == 1 % including error bars...
            if q_norm == 1
                errorbar(wallDist.mean(n,:), RPM_norm.side.mean(n,:), RPM_norm.side.err(n,:),...
                    RPM_norm.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                    "LineWidth", 1.5,"Color",default_colors(n,:))
            else
                errorbar(wallDist.mean(n,:), RPM.side.mean(n,:), RPM.side.err(n,:),...
                    RPM.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                    "LineWidth", 1.5,"Color",default_colors(n,:))
            end
        else % excluding error bars...
            if q_norm == 1
                plot(wallDist.mean(n,:), RPM_norm.side.mean(n,:),'s',"LineWidth", 1.5,...
                    "Color",default_colors(n,:))
            else
                plot(wallDist.mean(n,:), RPM.side.mean(n,:),'s',"LineWidth", 1.5,...
                    "Color",default_colors(n,:))
            end
        end
    end
    xlabel("x/R","Interpreter","latex","FontSize",font_size)
    title("Mean RPM (side) as a function of wall distance")
    legend(legend_strings(3:2:length(legend_strings)), "Location", "south",...
        "NumColumns",length(zR_all))
    ylabel("$\frac{RPM_{IBE}}{RPM_{OBE}}$","Interpreter","latex","Rotation",0,...
        "FontSize",font_size)
    ax = gca;
    ax.FontSize = font_size;
    grid on; grid minor;
    ylim([0.75 1.05])
end

%% Plotting CT (fore/aft/side separated) with or without error estimations 
figure % CT (fore)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        errorbar(wallDist.mean(n,:), CT.fore.mean(n,:), CT.fore.err(n,:),...
            CT.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
            "LineWidth", 1.5,"Color",default_colors(n,:))
    else % excluding error bars...
        plot(wallDist.mean(n,:), CT.fore.mean(n,:),'o',"LineWidth", 1.5,...
            "Color",default_colors(n,:))
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean C_T (fore) as a function of wall distance")
legend(legend_strings(1:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

figure % CT (aft)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        errorbar(wallDist.mean(n,:), CT.aft.mean(n,:), CT.aft.err(n,:),...
            CT.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
            "LineWidth", 1.5,"Color",default_colors(n,:))
    else % excluding error bars...
        plot(wallDist.mean(n,:), CT.aft.mean(n,:),'^',"LineWidth", 1.5,...
            "Color",default_colors(n,:))
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean C_T (aft) as a function of wall distance")
legend(legend_strings(2:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.05])

if exist("side_rotors","var")
    figure % CT (side)
    set(gcf,"Color","white")
    hold on
    for n = 1:length(processed_filepath)
        if q_error == 1 % including error bars...
            errorbar(wallDist.mean(n,:), CT.side.mean(n,:), CT.side.err(n,:),...
                CT.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        else % excluding error bars...
            plot(wallDist.mean(n,:), CT.side.mean(n,:),'s',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        end
    end
    xlabel("x/R","Interpreter","latex","FontSize",font_size)
    title("Mean C_T (side) as a function of wall distance")
    legend(legend_strings(3:2:length(legend_strings)), "Location", "south",...
        "NumColumns",length(zR_all))
    ylabel("$\frac{C_{T_{IBE}}}{C_{T_{OBE}}}$","Interpreter","latex","Rotation",0,...
        "FontSize",font_size)
    ax = gca;
    ax.FontSize = font_size;
    grid on; grid minor;
    ylim([0.75 1.05])
end

%% Plotting POWER (fore/aft/side separated) with or without error estimations 
figure % Power (fore)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        errorbar(wallDist.mean(n,:), Power.fore.mean(n,:), Power.fore.err(n,:),...
            Power.fore.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'o',...
            "LineWidth", 1.5,"Color",default_colors(n,:))
    else % excluding error bars...
        plot(wallDist.mean(n,:), Power.fore.mean(n,:),'o',"LineWidth", 1.5,...
            "Color",default_colors(n,:))
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean power (fore) as a function of wall distance")
legend(legend_strings(1:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.1])

figure % Power (aft)
set(gcf,"Color","white")
hold on
for n = 1:length(processed_filepath)
    if q_error == 1 % including error bars...
        errorbar(wallDist.mean(n,:), Power.aft.mean(n,:), Power.aft.err(n,:),...
            Power.aft.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'^',...
            "LineWidth", 1.5,"Color",default_colors(n,:))
    else % excluding error bars...
        plot(wallDist.mean(n,:), Power.aft.mean(n,:),'^',"LineWidth", 1.5,...
            "Color",default_colors(n,:))
    end
end
xlabel("x/R","Interpreter","latex","FontSize",font_size)
title("Mean power (aft) as a function of wall distance")
legend(legend_strings(2:2:length(legend_strings)), "Location", "south",...
    "NumColumns",length(zR_all))
ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
    "FontSize",font_size)
ax = gca;
ax.FontSize = font_size;
grid on; grid minor;
ylim([0.75 1.1])

if exist("side_rotors","var")
    figure % Power (side)
    set(gcf,"Color","white")
    hold on
    for n = 1:length(processed_filepath)
        if q_error == 1 % including error bars...
            errorbar(wallDist.mean(n,:), Power.side.mean(n,:), Power.side.err(n,:),...
                Power.side.err(n,:), wallDist.err(n,:), wallDist.err(n,:),'s',...
                "LineWidth", 1.5,"Color",default_colors(n,:))
        else % excluding error bars...
            plot(wallDist.mean(n,:), Power.side.mean(n,:),'s',"LineWidth", 1.5,...
                "Color",default_colors(n,:))
        end
    end
    xlabel("x/R","Interpreter","latex","FontSize",font_size)
    title("Mean power (side) as a function of wall distance")
    legend(legend_strings(3:2:length(legend_strings)), "Location", "south",...
        "NumColumns",length(zR_all))
    ylabel("$\frac{P_{IBE}}{P_{OBE}}$","Interpreter","latex","Rotation",0,...
        "FontSize",font_size)
    ax = gca;
    ax.FontSize = font_size;
    grid on; grid minor;
    ylim([0.75 1.1])
end

%% Plotting PRESSURE (non-normalized)
close all;

zR_plot = 6; % ENTER z/R height to plot data for
zR_idx = find(zR_plot == zR_all); 

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
    x = wallDist.mean(zR_idx,:);
    y = probe_locs(1)*ones(1,length(wallDist.mean(zR_idx,:)));
    if kiel_az == 90
        z = Pressure_90.mean(zR_idx,:,(n_probe-1)*8+1);
        err_z = Pressure_90.err(zR_idx,:,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_45.mean(zR_idx,:,(n_probe-1)*8+1);
        err_z = Pressure_45.err(zR_idx,:,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = wallDist.mean(zR_idx,:);
        y = probe_locs(i)*ones(1,length(wallDist.mean(zR_idx,:)));
        if kiel_az == 90
            z = Pressure_90.mean(zR_idx,:,(n_probe-1)*8+i);
            err_z = Pressure_90.err(zR_idx,:,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_45.mean(zR_idx,:,(n_probe-1)*8+i);
            err_z = Pressure_45.err(zR_idx,:,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
        plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
            pressure_plot_colors(i,:),"LineWidth",0.5) % error bar
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("x/R","FontSize",14,"Interpreter","latex")
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
    x = wallDist.mean(zR_idx,:);
    y = probe_locs(1)*ones(1,length(wallDist.mean(zR_idx,:)));
    if kiel_az == 90
        z = Pressure_90.fluct(zR_idx,:,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_45.fluct(zR_idx,:,(n_probe-1)*8+1);
    end

    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    for i = 2:7
        x = wallDist.mean(zR_idx,:);
        y = probe_locs(i)*ones(1,length(wallDist.mean(zR_idx,:)));
        if kiel_az == 90
            z = Pressure_90.fluct(zR_idx,:,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_45.fluct(zR_idx,:,(n_probe-1)*8+i);
        end

        plot3(x, y, z, "Color",pressure_plot_colors(i,:), "Marker",".",...
            "MarkerSize",14,"LineWidth",1) % mean
    end
    ylabel("r/R","FontSize",14,"Interpreter","latex")
    xlabel("x/R","FontSize",14,"Interpreter","latex")
    zlabel("$\sigma$ (Pa)","FontSize",14,"Interpreter","latex",...
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
    x = wallDist.mean(zR_idx,:);
    y = probe_locs(1)*ones(1,length(wallDist.mean(zR_idx,:)));
    if kiel_az == 90
        z = Pressure_norm_90.mean(zR_idx,:,(n_probe-1)*8+1);
        err_z = Pressure_norm_90.err(zR_idx,:,(n_probe-1)*8+1);
    elseif kiel_az == 45
        z = Pressure_norm_45.mean(zR_idx,:,(n_probe-1)*8+1);
        err_z = Pressure_norm_45.err(zR_idx,:,(n_probe-1)*8+1);
    end
    plot3(x, y, z, "Color",pressure_plot_colors(1,:), "Marker",".",...
        "MarkerSize",14,"LineWidth",1) % mean
    hold on
    plot3([x; x], [y; y], [-err_z; err_z] + [z; z],"Color",...
        pressure_plot_colors(1,:),"LineWidth",0.5) % error bar
    for i = 2:7
        x = wallDist.mean(zR_idx,:);
        y = probe_locs(i)*ones(1,length(wallDist.mean(zR_idx,:)));
        if kiel_az == 90
            z = Pressure_norm_90.mean(zR_idx,:,(n_probe-1)*8+i);
            err_z = Pressure_norm_90.err(zR_idx,:,(n_probe-1)*8+i);
        elseif kiel_az == 45
            z = Pressure_norm_45.mean(zR_idx,:,(n_probe-1)*8+i);
            err_z = Pressure_norm_45.err(zR_idx,:,(n_probe-1)*8+i);
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
    title(strcat("Rotor pressure as a function of boundary distance (",...
        string(probe_angles(n_probe))," deg)"))
    grid on; grid minor;
    ax = gca; ax.FontSize = 14;
    zlim([0 1.5])
end

%% Saving data (for pressure tests only rn)
save_filepath = strcat("Wall\MEAN PRESSURE DATA\Pressure_Cross",string(skew),"deg.mat");
if ~isfile(save_filepath)
    if kiel_az == 90
        wallDist_90 = wallDist;
        save(save_filepath, "probe_locs","wallDist_90", "Pressure_90",...
            "Pressure_norm_90")
    elseif kiel_az == 45
        wallDist_45 = wallDist;
        save(save_filepath, "probe_locs","wallDist_45", "Pressure_45",...
            "Pressure_norm_45")
    end
else
    if kiel_az == 90
        wallDist_90 = wallDist;
        save(save_filepath,"probe_locs","wallDist_90", "Pressure_90",...
            "Pressure_norm_90", "-append")
    elseif kiel_az == 45
        wallDist_45 = wallDist;
        save(save_filepath,"probe_locs","wallDist_45", "Pressure_45",...
            "Pressure_norm_45", "-append")
    end
end

