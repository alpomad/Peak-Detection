clear all; close all;

%%%THINGS TO CHECK BEFORE and AFTER THE RUN%%%

% % BEFORE
% % << Voltage value and sign (For negative peaks -> '-' sign)
% % 
% % << Height Filter Value for Freq 1 and Freq 2 -> Decide with the help of videos. 
% % After plots are created check if there are any obvious peak that is not counted.
% % If there are any peak that is not counted -> Adjust the height filter value accordingly and run again.
% % Note: Set Height Filter value to '0' if the data is a PBS measurement without beads. 
% % 
% % <<Set Time Difference (td) value. For 450 sa/s measurement -> td can be '0'
% % 
% % AFTER

% % For Noise:
% % << Save Mean 1-2 (Averages), RMS_Noise1-2 (RMS Noise), Percent_Noise1-2 (RMS %)
% % 
% % For Bead:
% % << Check E1_match_imp and E1_match_imp_repeat matrix size -> If they are different check repeating values and estimate the correct peak
% % Do the same for E2 data also.
% % 
% % << Copy Freq1_Freq2_matched_imp data to an excell file
% % 

%%%
%%% FREQ-1 %%%
%%%

%%% Data import %%%
currentPath = pwd;

inputPath = uigetdir('curentPath')

file1Path = [inputPath, '\Freq1.csv'];

importfilezurich(file1Path);

%%% Impedance calculation %%%
E1_complexCurrent = Freq1(:,2)+i*Freq1(:,3);

E1_current = abs(E1_complexCurrent);

voltage = 0.5; % real reading from 4 point measurements (0.5V was set) For Negative Peak Detection --> Set it to -0.5V

E1_impedance = voltage./(E1_current*sqrt(2)); % correct for differential measurements

E1_complexImp = voltage./E1_complexCurrent;

E1_realImp = real(E1_complexImp);
E1_imagImp = imag(E1_complexImp);

E1_time = Freq1(:,1);

E1_frequency = Freq1(:,4);

time1 = E1_time;
time10 = time1(1);
time1 = time1 - time10; %seconds
time1 = 1000*time1; %miliseconds

disp('Impedance1 calculated!');

%% Baseline removal %%%
E1_impedance_blineRemoved = msbackadj(time1, E1_impedance,...
    'WindowSize', 5000,...
    'StepSize', 300,...
    'PreserveHeights', 'false',...
    'SmoothMethod', 'none',...
    'QuantileValue', 0.5,...
    'ShowPlot', 1);

disp('Baseline1 removed!');



%% Peak determination %%%
[E1_imp_peaks, E1_imp_PFWHH, E1_imp_PExt] = mspeaks(time1, E1_impedance_blineRemoved,...
    'Denoising', 'false',...
    'FWHHFilter', 1,...
    'OverSegmentationFilter', 5,...
    'HeightFilter',1000,...  %Min. Peak Value for Freq1%
    'ShowPlot', 1);

disp('Peaks1 detected!');

%%% Peak width calculation %%%
E1_imp_numPeaks = length(E1_imp_peaks);

for g = 1:E1_imp_numPeaks
    E1_imp_PExt_diff_temp(g) = E1_imp_PExt(g,2) - E1_imp_PExt(g,1);
    E1_imp_PFWHH_diff_temp(g) = E1_imp_PFWHH(g,2) - E1_imp_PFWHH(g,1);
end

E1_imp_PExt_diff  = transpose(E1_imp_PExt_diff_temp);
E1_imp_PFWHH_diff  = transpose(E1_imp_PFWHH_diff_temp);

clear E1_imp_PExt_diff_temp E1_imp_PFWHH_diff_temp

%% Time conversion %%%
time1_export = time1 ; % miliseconds to seconds

%% Electrode impedance data export %%%
E1_impedance_dataexport = [time1_export, E1_impedance, E1_impedance_blineRemoved];

%% Peak data export %%%
E1_imp_peak_time_export = E1_imp_peaks(:,1) / 1000; % miliseconds to seconds

E1_imp_peak_dataexport = [E1_imp_peak_time_export, E1_imp_peaks(:,2), E1_imp_PExt_diff, E1_imp_PFWHH_diff];

%% Data plotting %%%
figure
plot(time1_export, E1_impedance);
title('Freq1 Absolute Impedance')

figure
plot(time1_export, E1_realImp);
title('Freq1 Real Impedance')

figure
plot(time1_export, E1_imagImp);
title('Freq1 Imaginary Impedance')

figure
plot(time1_export,E1_impedance_blineRemoved);
title('Freq1 Baseline Removed Impedance')

%%%
%%% FREQ-2 %%%
%%%

%%% Data import %%%

file2Path = [inputPath, '\Freq2.csv'];

importfilezurich(file2Path);

%%% Impedance calculation %%%
E2_complexCurrent = Freq2(:,2)+i*Freq2(:,3);

E2_current = abs(E2_complexCurrent);


E2_impedance = voltage./(E2_current*sqrt(2)); % correct for differential measurements

E2_complexImp = voltage./E2_complexCurrent;

E2_realImp = real(E2_complexImp);
E2_imagImp = imag(E2_complexImp);

E2_time = Freq2(:,1);

E2_frequency = Freq2(:,4);

time2 = E2_time;
time20 = time2(1);
time2 = time2 - time20; %seconds
time2 = 1000*time2; %miliseconds

disp('Impedance2 calculated!');

%% Baseline removal %%%
E2_impedance_blineRemoved = msbackadj(time2, E2_impedance,...
    'WindowSize', 5000,...
    'StepSize', 300,...
    'PreserveHeights', 'false',...
    'SmoothMethod', 'none',...
    'QuantileValue', 0.5,...
    'ShowPlot', 1);

disp('Baseline2 removed!');

%% Peak determination %%%
[E2_imp_peaks, E2_imp_PFWHH, E2_imp_PExt] = mspeaks(time2, E2_impedance_blineRemoved,...
    'Denoising', 'false',...
    'FWHHFilter', 1,...
    'OverSegmentationFilter', 5,...
    'HeightFilter',1000,...   %Min. Peak Value for Freq2%
    'ShowPlot', 1);

disp('Peaks2 detected!');

%%% Peak width calculation %%%
E2_imp_numPeaks = length(E2_imp_peaks);

for k = 1:E2_imp_numPeaks
    E2_imp_PExt_diff_temp(k) = E2_imp_PExt(k,2) - E2_imp_PExt(k,1);
    E2_imp_PFWHH_diff_temp(k) = E2_imp_PFWHH(k,2) - E2_imp_PFWHH(k,1);
end

E2_imp_PExt_diff  = transpose(E2_imp_PExt_diff_temp);
E2_imp_PFWHH_diff  = transpose(E2_imp_PFWHH_diff_temp);

clear E2_imp_PExt_diff_temp E2_imp_PFWHH_diff_temp

%% Time conversion %%%
time2_export = time2 ; % miliseconds to seconds

%% Electrode impedance data export %%%
E2_impedance_dataexport = [time2_export, E2_impedance, E2_impedance_blineRemoved];

%% Peak data export %%%
E2_imp_peak_time_export = E2_imp_peaks(:,1) / 1000; % miliseconds to seconds

E2_imp_peak_dataexport = [E2_imp_peak_time_export, E2_imp_peaks(:,2), E2_imp_PExt_diff, E2_imp_PFWHH_diff];

%% Data plotting %%%
figure
plot(time2_export, E2_impedance);
title('Freq2 Absolute Impedance')

figure
plot(time2_export, E2_realImp);
title('Freq2 Real Impedance')

figure
plot(time2_export, E2_imagImp);
title('Freq2 Imaginary Impedance')

figure
plot(time2_export,E2_impedance_blineRemoved);
title('Freq2 Baseline Removed Impedance')

% NOISE CALCULATION %%

Mean1 = mean(E1_impedance);
Mean2 = mean(E2_impedance);
OffsetMean1 = E1_impedance - Mean1;
OffsetMean2 = E2_impedance - Mean2;
RMS_Noise1 = sqrt((sum(OffsetMean1.^2))./(length(OffsetMean1)));
RMS_Noise2 = sqrt((sum(OffsetMean2.^2))./(length(OffsetMean2)));
Percent_Noise1 = (RMS_Noise1/Mean1)*100;
Percent_Noise2 = (RMS_Noise2/Mean2)*100;

%% PEAK MATCHING %%

ImpLen1 = length(E1_imp_peaks);
ImpLen2 = length(E2_imp_peaks);
td=0; %Time Difference between peaks
na=1;
nb=1;
k1=1;
k2=1;
for aind = na:ImpLen1
    for bind = nb:ImpLen2
        if abs(E1_imp_peaks(aind,1) - E2_imp_peaks(bind,1)) <= td
            E1_match_imp(k1,1) = E1_imp_peaks(aind,1);
            E1_match_imp(k1,2) = E1_imp_peaks(aind,2);
            E2_match_imp(k2,1) = E2_imp_peaks(bind,1);
            E2_match_imp(k2,2) = E2_imp_peaks(bind,2);
            na = aind+1;
            nb = bind+1;
            k1=k1+1;
            k2=k2+1;
            break
        else
        
        end
    end
end

m1=1;
m2=1;
for cind = 1:ImpLen1
    for dind = 1:ImpLen2
        if abs(E1_imp_peaks(cind,1)-E2_imp_peaks(dind,1)) <= td
            E1_match_imp_repeat(m1,1) = E1_imp_peaks(cind,1);
            E1_match_imp_repeat(m1,2) = E1_imp_peaks(cind,2);
            E2_match_imp_repeat(m2,1) = E2_imp_peaks(dind,1);
            E2_match_imp_repeat(m2,2) = E2_imp_peaks(dind,2);
            m1=m1+1;
            m2=m2+1;
         
        else
        
        end
    end
end

Freq1_Freq2_matched_imp = [E1_match_imp(:,2), E2_match_imp(:,2)];