function [data_nf2ff] = nf2ffTransformation_fft(data_nf,f,padding,phi_range,theta_range)

% Create Results Table
p = length(phi_range);
t = length(theta_range);

data_nf2ff = table(zeros(p*t,1),zeros(p*t,1),zeros(p*t,3),zeros(p*t,1),zeros(p*t,1),zeros(p*t,1));
data_nf2ff.Properties.VariableNames = {'theta','phi','E','Etheta','Ephi','Eabs'};


% Step size in x,y direction [m]
delta_x = unique(data_nf.x);
number_of_samples_x = numel(delta_x);
delta_x = abs(delta_x(1)-delta_x(2));
delta_y = unique(data_nf.y);
number_of_samples_y = numel(delta_y);
delta_y = abs(delta_y(1)-delta_y(2));

% Scanner Plane Dimensions
length_x = delta_x * (number_of_samples_x-1);
length_y = delta_y * (number_of_samples_y-1);
z0 = mean(data_nf.z);


% Rectangular Wavenumber Vector as according to Balanis
lambda=physconst('lightspeed')/f;
k0=2*pi/lambda;
      

number_of_samples_x_padded = padding*number_of_samples_x;
number_of_samples_y_padded = padding*number_of_samples_y;
m= -number_of_samples_x_padded/2:1:number_of_samples_x_padded/2-1;
n= -number_of_samples_y_padded/2:1:number_of_samples_y_padded/2-1;

kx=2*pi*m/(number_of_samples_x_padded*delta_x);
ky=2*pi*n/(number_of_samples_y_padded*delta_y);
[ky_grid,kx_grid] = meshgrid(ky,kx);
kz_grid = sqrt(k0^2-kx_grid.^2-ky_grid.^2);

% Spherical Far-Field Wavenumber Vector
[theta_grid,phi_grid]=meshgrid(theta_range,phi_range);
kx_grid_spherical = k0*sin(theta_grid).*cos(phi_grid);
ky_grid_spherical = k0*sin(theta_grid).*sin(phi_grid);
kz_grid_spherical = k0*cos(theta_grid);


% Reshape E field to fit grid
Ex_nf = reshape(data_nf.E(:,1),number_of_samples_x,number_of_samples_y);
Ey_nf = reshape(data_nf.E(:,2),number_of_samples_x,number_of_samples_y);
Ez_nf = reshape(data_nf.E(:,3),number_of_samples_x,number_of_samples_y);


% Sampling Window Function
% L = length(fx);
% h = hamming(L);
% H = h*h'/L;


% Retrieve Plane Wave Modes through FFT
fx=ifftshift(ifft2(Ex_nf,number_of_samples_x_padded,number_of_samples_y_padded));
fy=ifftshift(ifft2(Ey_nf,number_of_samples_x_padded,number_of_samples_y_padded));
fz=-(fx.*kx_grid+fy.*ky_grid)./kz_grid;


% Interpolate Modes in spherical coordinates
fx_ff=interp2(kx,ky,abs(fx),kx_grid_spherical,ky_grid_spherical,'spline');
fy_ff=interp2(kx,ky,abs(fy),kx_grid_spherical,ky_grid_spherical,'spline');
fz_ff=interp2(kx,ky,abs(fz),kx_grid_spherical,ky_grid_spherical,'spline');

% Far Field 
r=100;
C=1j*(k0*exp(-1j*k0*r))/(2*pi*r);

Etheta=C*(fx_ff.*cos(phi_grid)+fy_ff.*sin(phi_grid));
Ephi=C*cos(theta_grid).*(-fx_ff.*sin(phi_grid)+fy_ff.*cos(phi_grid));
Ex = C*cos(theta_grid).*fx_ff;
Ey = C*cos(theta_grid).*fy_ff;
Ez = C*cos(theta_grid).*fz_ff;

% Fill Results Table
s = numel(theta_grid);
data_nf2ff.Etheta = reshape(Etheta,s,1);
data_nf2ff.Ephi = reshape(Ephi,s,1);
data_nf2ff.E = [reshape(Ex,s,1),reshape(Ey,s,1),reshape(Ez,s,1)];
data_nf2ff.Eabs = vecnorm(data_nf2ff.E,2,2);
data_nf2ff.phi = reshape(phi_grid,s,1);
data_nf2ff.theta= reshape(theta_grid,s,1);

end


