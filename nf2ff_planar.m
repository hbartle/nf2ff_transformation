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
f = 5e9;
% setup = 'PlanarScan_Waveguide/';
setup = 'PlanarScan_HornAntenna/';

data_ff = readtable([setup 'Farfield/farfield (f=' num2str(f*1e-9) ') [1].txt']);

scans = dir(setup);
scan_names = {scans.name};
scan_names(ismember(scan_names,{'.','..','Farfield'})) = [];

data_nf= cellfun(@(scan_names) readtable([setup,scan_names,'/NearFieldProbeResults' num2str(f*1e-9) 'GHz.txt']),scan_names,'UniformOutput',false);

% Rearrange Farfield Table
data_ff = table(data_ff.Var1,data_ff.Var2,data_ff.Var3,data_ff.Var4,data_ff.Var6);
data_ff.Properties.VariableNames = {'theta' 'phi' 'Eabs' 'Ethetaabs' 'Ephiabs'};

% Rearrange Nearfield table
data_nf = cellfun(@rearrangeTables,data_nf,'UniformOutput',false);

disp('Done!')
%% NF2FF transformation
disp('NF2FF Transformation...')
phi_range = 0:5:355;
theta_range=0:1:80;

data_nf2ff = cellfun(@(data_nf) nf2ffTransformation_w_expansion(data_nf,f,phi_range,theta_range),data_nf,'Uniformoutput',false);
%data_nf2ff = cellfun(@(data_nf) nf2ffTransformation(data_nf,f,phi_range,theta_range),data_nf,'Uniformoutput',false);

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
title('Phi=0°')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=0°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffPhiCut(data_nf2ff,data_ff,0,theta_range),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('E-Field Pattern [-]')
title('Phi=0°')
legend(scan_names)

% Phi=90 cut
figure('name','Far-Field Cuts,Phi=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
plotFFPhiCut(data_ff,0,normalized)
cellfun(@(data_nf2ff) plotNFPhiCut(data_nf2ff,90,normalized),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('E-Field Pattern [-]')
title('Phi=90°')
legend(['Far-Field',scan_names])


figure('name','Far-Field Error,Phi=90°','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1]);
cellfun(@(data_nf2ff) plotDiffPhiCut(data_nf2ff,data_ff,90,theta_range),data_nf2ff)
grid on
xlabel('Theta [°]')
ylabel('E-Field Pattern [-]')
title('Phi=90°')
legend(scan_names)



disp('Done!')





