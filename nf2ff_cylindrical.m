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
scan_names(ismember(scan_names,{'.','..','Farfield','Scan-PhiSamples-ZSpacing-ZSamples.txt'})) = [];

data_nf= cellfun(@(scan_names) readtable([setup,scan_names,'/NearFieldProbeResults' num2str(f*1e-9) 'GHz.txt']),scan_names,'UniformOutput',false);

% Rearrange farfield table
data_ff = table(data_ff.Var1*pi/180,data_ff.Var2*pi/180,data_ff.Var3,data_ff.Var4,data_ff.Var6);
data_ff.Properties.VariableNames = {'theta' 'phi' 'Eabs' 'Ethetaabs' 'Ephiabs'};

% Rearrange nearfield table
data_nf = cellfun(@rearrangeTables,data_nf,'UniformOutput',false);
data_nf = cellfun(@(data_nf) rotateCylindricalNFData(data_nf,'z'),data_nf,'UniformOutput',false);

% Select measurements to process
% data_nf = data_nf([1,5]);
% scan_names = scan_names([1,5]);

disp('Done!')
%% NF2FF transformation
disp('NF2FF Transformation...')
delta_theta=1;
theta_range = (0:delta_theta:180-delta_theta)*pi/180;
delta_phi = 1;
phi_range = (0:delta_phi:360-delta_phi)*pi/180;
% NF2FF Algorithm Parameters
window = 'tukey';

% data_nf2ff = cellfun(@(data_nf) nf2ff_cylindrical_fft(data_nf,f,theta_range,window),data_nf,'Uniformoutput',false);
data_nf2ff = cellfun(@(data_nf) nf2ff_cylindrical_manual(data_nf,f,theta_range,phi_range,window),data_nf,'Uniformoutput',false);

disp('Done!')
%% Plots
disp('Plotting...')
close all

normalized = true;
logarithmic = true;

% Phi=0 cut
figure('name','Far-Field Cuts,Phi=0�','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
plotFFPhiCut(data_ff,0,normalized,logarithmic)
cellfun(@(data_nf2ff) plotNFPhiCutCylindrical(data_nf2ff,0,normalized,logarithmic),data_nf2ff)
grid on
xlabel('Theta [�]')
if logarithmic == true
    ylim([-50 0])
    ylabel('E-Field Pattern [dB]')
elseif normalized == true
    ylim([0 1])
    ylabel('E-Field Pattern [-]')
end
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Phi=0�')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=0�','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffPhiCutCylindrical(data_nf2ff,data_ff,0,theta_range),data_nf2ff)
grid on
xlabel('Theta [�]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Phi=0�')
legend(scan_names)


% Theta=90 cut
figure('name','Far-Field Cuts,Theta=90�','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
plotFFThetaCut(data_ff,pi/2,normalized,logarithmic)
cellfun(@(data_nf2ff) plotNFThetaCutCylindrical(data_nf2ff,pi/2,normalized,logarithmic),data_nf2ff)
grid on
xlabel('Phi [�]')
if logarithmic == true
    ylim([-50 0])
    ylabel('E-Field Pattern [dB]')
elseif normalized == true
    ylim([0 1])
    ylabel('E-Field Pattern [-]')
end
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Theta=90�')
legend(['Far-Field',scan_names])
 
 
figure('name','Far-Field Error,Theta=90�','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffThetaCutCylindrical(data_nf2ff,data_ff,pi/2),data_nf2ff)
grid on
xlabel('Phi [�]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Theta=90�')
legend(scan_names)


disp('Done!')





