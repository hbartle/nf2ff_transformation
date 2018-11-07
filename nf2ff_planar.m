% 
% NF2FF Conversion using Planar Scanner
%
clear 
close all
clc

disp('************************************************')
disp('          Near-To-Far-Field Conversion')
disp('          Planar Scanner')
disp('************************************************')
%% Load Data
disp('Load Data...')
f = 6e9;
% setup = '../measurement/PlanarScan_Waveguide/';
setup = '../measurement/PlanarScan_HornAntenna/';
% setup = '../measurement/PlanarScan_PatchAntenna/';

data_ff = readtable([setup 'Farfield/farfield (f=' num2str(f*1e-9) ') [1].txt']);

scans = dir(setup);
scan_names = {scans.name};
scan_names(ismember(scan_names,{'.','..','Farfield'})) = [];

data_nf= cellfun(@(scan_names) readtable([setup,scan_names,'/NearFieldProbeResults' num2str(f*1e-9) 'GHz.txt']),scan_names,'UniformOutput',false);

% Rearrange Farfield Table
data_ff = table(data_ff.Var1*pi/180,data_ff.Var2*pi/180,data_ff.Var3,data_ff.Var4,data_ff.Var6);
data_ff.Properties.VariableNames = {'theta' 'phi' 'Eabs' 'Ethetaabs' 'Ephiabs'};

% Rearrange Nearfield table
data_nf = cellfun(@rearrangeTables,data_nf,'UniformOutput',false);
 
data_nf = data_nf(8:10);
scan_names = scan_names(8:10);
disp('Done!')
%% NF2FF transformation
disp('NF2FF Transformation...')
delta_theta=1;
delta_phi=1;
theta_range = (-90:delta_theta:90-delta_theta)*pi/180;
phi_range= (0:delta_phi:180-delta_phi)*pi/180;
[theta_grid,phi_grid]=meshgrid(theta_range,phi_range);

padding = 4;

data_nf2ff = cellfun(@(data_nf) nf2ffTransformation_fft(data_nf,f,padding,phi_range,theta_range),data_nf,'Uniformoutput',false);
% data_nf2ff = cellfun(@(data_nf) nf2ffTransformation_manual(data_nf,f,phi_range,theta_range),data_nf,'Uniformoutput',false);

disp('Done!')
%% Plots
disp('Plotting...')
close all

normalized = true;




% Phi=0 cut
figure('name','Far-Field Cuts,Phi=0°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
plotFFPhiCut(data_ff,0,normalized)
cellfun(@(data_nf2ff) plotNFPhiCut(data_nf2ff,0,normalized),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Phi=0°')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=0°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffPhiCut(data_nf2ff,data_ff,0,theta_range),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Phi=0°')
legend(scan_names)

% Phi=90 cut
figure('name','Far-Field Cuts,Phi=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
plotFFPhiCut(data_ff,pi/2,normalized)
cellfun(@(data_nf2ff) plotNFPhiCut(data_nf2ff,pi/2,normalized),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('E-Field Pattern [-]')
title('Far-Field Cut Phi=90°')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffPhiCut(data_nf2ff,data_ff,pi/2,theta_range),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('Difference to Reference Far-Field [dB]')
title('Difference to Reference Far-Field, Phi=90°')
legend(scan_names)


disp('Done!')





