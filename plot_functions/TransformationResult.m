function [] = TransformationResult(data_ff,data_nf2ff,scanner_case)

switch scanner_case
    case 'planar'

        % Transformed Near-field
        n_theta_nf2ff = length(unique(data_nf2ff.theta));
        n_phi_nf2ff = length(unique(data_nf2ff.phi));

        X_nf2ff = reshape(data_nf2ff.theta,n_theta_nf2ff,n_phi_nf2ff);
        Y_nf2ff = reshape(data_nf2ff.phi,n_theta_nf2ff,n_phi_nf2ff);
        Z_nf2ff = reshape(data_nf2ff.Eabs,n_theta_nf2ff,n_phi_nf2ff);


        X_nf2ff = [X_nf2ff(:,end/2+1:end);fliplr(-X_nf2ff(:,2:end/2+1)) ]';
        Y_nf2ff = [Y_nf2ff(:,end/2+1:end); pi+fliplr(Y_nf2ff(:,2:end/2+1))]';
        Z_nf2ff = [Z_nf2ff(:,end/2+1:end); fliplr(Z_nf2ff(:,2:end/2+1))]';

        maxTheta = size(Z_nf2ff,1);
        maxPhi = size(Z_nf2ff,2);

        % Far Field
        i =find(data_ff.theta~=pi); % Clean up the duplicate samples
        n_theta_ff = length(unique(data_ff.theta(i)));
        n_phi_ff = length(unique(data_ff.phi(i)));

        X_ff = reshape(data_ff.theta(i),n_theta_ff,n_phi_ff);
        Y_ff = reshape(data_ff.phi(i),n_theta_ff,n_phi_ff);
        Z_ff = reshape(data_ff.Eabs(i),n_theta_ff,n_phi_ff);

        X_ff = [X_ff(1:end/2,:); X_ff(end/2+1:end,:)];
        Y_ff = [Y_ff(1:end/2,:); Y_ff(end/2+1:end,:)];
        Z_ff = [Z_ff(1:end/2,:); Z_ff(end/2+1:end,:)];


        X_ff = X_ff(1:maxTheta,1:maxPhi);
        Y_ff = Y_ff(1:maxTheta,1:maxPhi);
        Z_ff = Z_ff(1:maxTheta,1:maxPhi);

    case 'cylindrical'
        
        % Transformed Near-field
        n_theta_nf2ff = length(unique(data_nf2ff.theta));
        n_phi_nf2ff = length(unique(data_nf2ff.phi));

        X_nf2ff = reshape(data_nf2ff.theta,n_phi_nf2ff,n_theta_nf2ff)';
        Y_nf2ff = reshape(data_nf2ff.phi,n_phi_nf2ff,n_theta_nf2ff)';
        Z_nf2ff = reshape(data_nf2ff.Eabs,n_phi_nf2ff,n_theta_nf2ff)';


        
        % Far Field
        i =find(data_ff.theta~=pi); % Clean up the duplicate samples
        n_theta_ff = length(unique(data_ff.theta(i)));
        n_phi_ff = length(unique(data_ff.phi(i)));

        X_ff = reshape(data_ff.theta(i),n_theta_ff,n_phi_ff);
        Y_ff = reshape(data_ff.phi(i),n_theta_ff,n_phi_ff);
        Z_ff = reshape(data_ff.Eabs(i),n_theta_ff,n_phi_ff);

    case 'spherical'
        
        
        
    
        

%% Normalized Field Difference
Z_ff_norm = Z_ff/max(max(Z_ff));
Z_nf2ff_norm = Z_nf2ff/max(max(Z_nf2ff));

Z_diff = 10*log10(abs(Z_ff_norm - Z_nf2ff_norm));


%% Plot
fontsize = 14;
close all

figure('name','Tranformation Result - Far-Field','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
s = surf(X_ff*180/pi,Y_ff*180/pi,Z_ff);
s.EdgeColor = 'flat';
rotate3d on
xlabel('Theta [°]')
ylabel('Phi [°]')
zlabel('Eabs [V/m]')
colorbar

figure('name','Tranformation Result - Transformed Near-Field','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
s=surf(X_nf2ff,Y_nf2ff,Z_nf2ff);
s.EdgeColor = 'flat';
rotate3d on
xlabel('Theta [°]')
ylabel('Phi [°]')
zlabel('Eabs [V/m]')
colorbar

figure('name','Tranformation Result - Error','numbertitle','off',...
        'units','normalized','outerposition',[0 0 1 1],...
        'DefaultAxesFontSize',fontsize);
s=surf(X_nf2ff,Y_nf2ff,Z_diff);
s.EdgeColor = 'flat';
rotate3d on
xlabel('Theta [°]')
ylabel('Phi [°]')
zlabel('Difference [dB]')
colorbar


end

