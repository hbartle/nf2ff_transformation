% 
% NF2FF Conversion using Planar Scanner
%
clear 
close all
clc

disp('************************************************')
disp('          Near-To-Far-Field Conversion')
disp('          Cylindrical Scanner')
disp('************************************************')

addpath('misc_functions')
addpath('plot_functions')
addpath('transformation_functions')
%% Load Data
disp('Load Data...')
f = 5e9;

setup = '../measurement/CylindricalScan_HornAntenna/';
% setup = '../measurement/CylindricalScan_HornAntenna_Rotated/';
% setup = '../measurement/CylindricalScan_PatchAntenna/';

data_ff = readtable([setup 'Farfield/farfield (f=' num2str(f*1e-9) ') [1].txt']);

scans = dir(setup);
scan_names = {scans.name};
scan_names(ismember(scan_names,{'.','..','Farfield','Farfield_old','Scan-PhiSamples-ZSpacing-ZSamples_R300.txt'})) = [];

data_nf= cellfun(@(scan_names) readtable([setup,scan_names,'/NearFieldProbeResults' num2str(f*1e-9) 'GHz.txt']),scan_names,'UniformOutput',false);

% Rearrange farfield table
data_ff = table(data_ff.Var1*pi/180,data_ff.Var2*pi/180,data_ff.Var3,data_ff.Var4,data_ff.Var6);
data_ff.Properties.VariableNames = {'theta' 'phi' 'Eabs' 'Ethetaabs' 'Ephiabs'};

% Rearrange nearfield table
data_nf = cellfun(@rearrangeTables,data_nf,'UniformOutput',false);
data_nf = cellfun(@(data_nf) rotateCylindricalNFData(data_nf,'z'),data_nf,'UniformOutput',false);

% Select measurements to process
% data_nf = data_nf([10,6,7]);
% scan_names = scan_names([10,6,7]);
% data_nf = data_nf([8,9,10,11,12]);
% scan_names = scan_names([8,9,10,11,12]);

disp('Done!')
%% NF2FF transformation
disp('NF2FF Transformation...')
delta_theta=1;
theta_range = (0:delta_theta:180-delta_theta)*pi/180;
delta_phi = 1;
phi_range = (0:delta_phi:360-delta_phi)*pi/180;

% data_nf2ff_rect = cellfun(@(data_nf) nf2ff_cylindrical_fft(data_nf,f,theta_range,phi_range,'none'),data_nf,'Uniformoutput',false);
% data_nf2ff_hamming = cellfun(@(data_nf) nf2ff_cylindrical_fft(data_nf,f,theta_range,phi_range,'hamming'),data_nf,'Uniformoutput',false);
% data_nf2ff_tukey = cellfun(@(data_nf) nf2ff_cylindrical_fft(data_nf,f,theta_range,phi_range,'tukey'),data_nf,'Uniformoutput',false);

data_nf2ff_rect = cellfun(@(data_nf) nf2ff_cylindrical_manual(data_nf,f,theta_range,phi_range,'none'),data_nf,'Uniformoutput',false);
data_nf2ff_hamming = cellfun(@(data_nf) nf2ff_cylindrical_manual(data_nf,f,theta_range,phi_range,'hamming'),data_nf,'Uniformoutput',false);
data_nf2ff_tukey = cellfun(@(data_nf) nf2ff_cylindrical_manual(data_nf,f,theta_range,phi_range,'tukey'),data_nf,'Uniformoutput',false);

data_nf2ff = data_nf2ff_rect;
disp('Done!')
%% Plots
disp('Plotting...')
close all
fontsize = 14;

normalized = true;
logarithmic = false;

% Phi=0 cut
phi_cut = 0;
figure('name','Far-Field Cuts,Phi=0°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plotFFPhiCut(data_ff,phi_cut,normalized,logarithmic)
cellfun(@(data_nf2ff) plotNFPhiCutCylindrical(data_nf2ff,phi_cut,normalized,logarithmic),data_nf2ff)
grid on
xlabel('Theta [°]')
if logarithmic == true
    ylim([-50 0])
    ylabel('E-Field Pattern [dB]')
elseif normalized == true
    ylim([0 1])
    ylabel('E-Field Pattern [-]')
end
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Phi=0°')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=0°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
cellfun(@(data_nf2ff) plotDiffPhiCutCylindrical(data_nf2ff,data_ff,phi_cut,theta_range),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Phi=0°')
legend(scan_names)


% Theta=90 cut
theta_cut = pi/2;
figure('name','Far-Field Cuts,Theta=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plotFFThetaCut(data_ff,theta_cut,normalized,logarithmic)
cellfun(@(data_nf2ff) plotNFThetaCutCylindrical(data_nf2ff,theta_cut,normalized,logarithmic),data_nf2ff)
grid on
xlabel('Phi [°]')
if logarithmic == true
    ylim([-50 0])
    ylabel('E-Field Pattern [dB]')
elseif normalized == true
    ylim([0 1])
    ylabel('E-Field Pattern [-]')
end
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Theta=90°')
legend(['Far-Field',scan_names])
 
 
figure('name','Far-Field Error,Theta=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
cellfun(@(data_nf2ff) plotDiffThetaCutCylindrical(data_nf2ff,data_ff,theta_cut),data_nf2ff)
grid on
xlabel('Phi [°]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Theta=90°')
legend(scan_names)


disp('Done!')

%% Accumulated Pattern Error Analysis
close all
fontsize = 14;

% Calculate APE
error_rect = cellfun(@(data_nf2ff) ErrorAnalysis(data_ff,data_nf2ff,'cylindrical'),data_nf2ff);
error_hamming = cellfun(@(data_nf2ff) ErrorAnalysis(data_ff,data_nf2ff,'planar'),data_nf2ff_hamming);
error_tukey = cellfun(@(data_nf2ff) ErrorAnalysis(data_ff,data_nf2ff,'planar'),data_nf2ff_tukey);

% Increasing Area
samples=[8,9,10,11,12];
xValues = [30,35,40,45,50];
figure('name','Accumulated Pattern Error, Varying Area','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(xValues,100*error_rect(samples),'--*')
hold on
plot(xValues,100*error_hamming(samples),'--*')
plot(xValues,100*error_tukey(samples),'--*')
grid on
xlabel('Number of probes in linear dimension');
ylabel('Accumulated Pattern Error [%]')
title({'Accumulated Pattern Error','Varying Measurement Area, \lambda/2 spacing,'})
legend('Rectangular','Hamming','Tukey')

% Varying Spacing
samples = [10,7,6];
figure('name','Accumulated Pattern Error, Varying Spacing','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(100*error_rect(samples),'--*')
hold on
plot(100*error_hamming(samples),'--*')
plot(100*error_tukey(samples),'--*')
grid on
ylabel('Accumulated Pattern Error [%]')
title({'Accumulated Pattern Error','Varying Measurement Spacing'})
xtickangle(45)
xticks(1:3)
xticklabels({'42x40','42x45','42x50'})
legend('Rectangular','Hamming','Tukey')




%% HPBW Error Analysis
close all
fontsize = 14;

% Calculate HPBW error
[~,~, hpbw_err_rect] = cellfun(@(data_nf2ff) HPBWError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_rect,'UniformOutput',false);
[~,~, hpbw_err_hamming] = cellfun(@(data_nf2ff) HPBWError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_hamming,'UniformOutput',false);
[~,~, hpbw_err_tukey] = cellfun(@(data_nf2ff) HPBWError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_tukey,'UniformOutput',false);

hpbw_err_rect = cell2mat(hpbw_err_rect);
hpbw_err_hamming = cell2mat(hpbw_err_hamming);
hpbw_err_tukey = cell2mat(hpbw_err_tukey);

xValues = [30,35,40,45,50];

% Increasing Area / Phi = 0 Cut
samples = [8,9,10,11,12];
figure('name','HPBW Error, Varying Area','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(xValues,hpbw_err_rect(1,samples),'-*')
hold on
plot(xValues,hpbw_err_hamming(1,samples),'-*')
plot(xValues,hpbw_err_tukey(1,samples),'-*')
grid on
xlabel('Number of probes in linear dimension');
xticks(xValues)
ylabel('HPBW Error [Degrees]')
title({'HPBW Error','Varying Measurement Area, \lambda/2 spacing, Phi=0 Cut'})
legend('Rectangular','Hamming','Tukey')

% Increasing Area / Theta = 90 Cut
samples = [8,9,10,11,12];
figure('name','HPBW Error, Varying Area','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(xValues,hpbw_err_rect(2,samples),'-*')
hold on
plot(xValues,hpbw_err_hamming(2,samples),'-*')
plot(xValues,hpbw_err_tukey(2,samples),'-*')
grid on
xlabel('Number of probes in linear dimension');
xticks(xValues)
ylabel('HPBW Error [Degrees]')
title({'HPBW Error','Varying Measurement Area, \lambda/2 spacing, Theta=90 Cut'})
legend('Rectangular','Hamming','Tukey')

% Varying Spacing Phi=0
samples = [10,7,6];
figure('name','Accumulated Pattern Error, Varying Spacing','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(hpbw_err_rect(1,samples),'-*')
hold on
plot(hpbw_err_hamming(1,samples),'-*')
plot(hpbw_err_tukey(1,samples),'-*')
grid on
ylabel('HPBW Error [Degrees]')
title({'HPBW Error','Varying Measurement Spacing, 30\lambda/2 x 30\lambda/2, Phi=0 Cut'})
xtickangle(45)
xticks(1:3)
xticklabels({'42x40','42x45','42x50'})
legend('Rectangular','Hamming','Tukey')

% Varying Spacing Theta=90
samples = [10,7,6];
figure('name','Accumulated Pattern Error, Varying Spacing','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(hpbw_err_rect(2,samples),'-*')
hold on
plot(hpbw_err_hamming(2,samples),'-*')
plot(hpbw_err_tukey(2,samples),'-*')
grid on
ylabel('HPBW Error [Degrees]')
title({'HPBW Error','Varying Measurement Spacing, 30\lambda/2 x 30\lambda/2, Theta=90 Cut'})
xtickangle(45)
xticks(1:3)
xticklabels({'42x40','42x45','42x50'})
legend('Rectangular','Hamming','Tukey')

%% Side-Lobe Level Error Analysis
close all
fontsize = 14;

% Calculate Side-Lobe Level error
[~,~, sll_err_rect] = cellfun(@(data_nf2ff) SLLError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_rect,'UniformOutput',false);
[~,~, sll_err_hamming] = cellfun(@(data_nf2ff) SLLError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_hamming,'UniformOutput',false);
[~,~, sll_err_tukey] = cellfun(@(data_nf2ff) SLLError(data_ff,data_nf2ff,'cylindrical'),data_nf2ff_tukey,'UniformOutput',false);

sll_err_rect = cell2mat(sll_err_rect);
sll_err_hamming = cell2mat(sll_err_hamming);
sll_err_tukey = cell2mat(sll_err_tukey);

xValues = [30,35,40,45,50];

% Increasing Area / Phi = 0 Cut
samples = [8,9,10,11,12];
figure('name','Side-Lobe-Level Error, Varying Area','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(xValues,sll_err_rect(samples),'-*')
hold on
plot(xValues,sll_err_hamming(samples),'-*')
plot(xValues,sll_err_tukey(samples),'-*')
grid on
xlabel('Number of probes in linear dimension');
xticks(xValues)
ylabel('SLL Error [dB]')
title({'Side-Lobe-Level Error','Varying Measurement Area, \lambda/2 spacing'})
legend('Rectangular','Hamming','Tukey')

% Varying Spacing Phi=0
samples = [10,7,6];
figure('name','Side-Lobe-Level Error, Varying Spacing','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
plot(sll_err_rect(samples),'-*')
hold on
plot(sll_err_hamming(samples),'-*')
plot(sll_err_tukey(samples),'-*')
grid on
ylabel('SLL Error [Degrees]')
title({'Side-Lobe-Level Error','Varying Measurement Spacing, 30\lambda/2 x 30\lambda/2'})
xtickangle(45)
xticks(1:3)
xticklabels({'42x40','42x45','42x50'})
legend('Rectangular','Hamming','Tukey')








