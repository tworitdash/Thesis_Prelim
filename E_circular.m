close all;
clear;
clc;

%% Waveguide parameters
d = 0.0405319403216; % diameter in meters
r = d/2;             % radius


c0 = 3e8;
er0 = 8.85418782e-12;
mu0 = 1.25663706e-6;

er = 1;                 %relative permittivity
mur = 1;                %relative permeability
epsilon = er * er0;
mu = mu0 * mur;
%% Zeros of Bessel's functions and their derivatives:
m = input('Input the first digit of the mode number: ');
n = input('Input the second digit of the mode number: ');
mode = input('Which mode (TE/TM): ', 's');

xm = linspace(0.1, 100, 100000);

Jm = @(z) besselj(m, z);                                         % Bessel's function

dz = 1e-5;

Jm_der = @(z) (besselj(m, z + dz) - besselj(m, z - dz))./(2 * dz); % Derivative of the Bessel's function

if mode == "TE"
    fun = Jm_der;
else
    fun = Jm;
end

ym = fun(xm);                   % Values of the Bessel's function at the positions defined by xm

chsign = find(diff(sign(ym)));  % Detection of the sign changes 

xmn_all = zeros(size(chsign));

for i = 1:size(chsign, 2)
    xmn_all(i) = fzero(fun, xm(chsign(i)));  % finding the roots near the points where Jm changes sign
end



xmn = xmn_all(n); 


beta_rho = xmn/r;  % wave number along the rho direction (\beta_{\rho})

z = 0.2;

%% Cut of frequency

fc = xmn ./ (2 * pi * r * sqrt(mu .* epsilon));

%% Source Parameters
f = fc + eps;

omega = 2 * pi * f;

lamb = c0./f; % wavelength
 
beta = 2 * pi ./ lamb; % wave number

rho = eps:r/1000:r;

phi = -pi+eps:pi/2000:pi-eps;

[rho_, phi_] = meshgrid(rho, phi);


%% Wave equations

if mode == "TE"
    
    [Erho, Ephi, Ez] = E_TE(epsilon, m, rho_, phi_, beta_rho, z, beta);
    [Hrho, Hphi, Hz] = H_TE(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
   
elseif mode == "TM"
    [Erho, Ephi, Ez] = E_TM(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
    [Hrho, Hphi, Hz] = H_TM(mu, m, rho_, phi_, beta_rho, z, beta);
end



E = sqrt(abs(Erho).^2 + abs(Ephi).^2 + Ez.^2);
H = sqrt(abs(Hrho).^2 + abs(Hphi).^2 + Hz.^2);
   
figure;
   
polarPcolor(rho, phi.*180/pi, abs(Erho)); shading flat;
colormap('jet');

title('E_{\rho}', 'FontSize', 12, 'FontWeight', 'bold');

figure;
   
polarPcolor(rho, phi.*180/pi, abs(Ephi)); shading flat;
colormap('jet');
title('E_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');

figure;
polarPcolor(rho, phi.*180/pi, abs(Hrho)); shading flat;
colormap('jet');
title('H_{\rho}', 'FontSize', 12, 'FontWeight', 'bold');

figure;
   
polarPcolor(rho, phi.*180/pi, abs(Hphi)); shading flat;
colormap('jet');
title('H_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');

% 
% figure;
% polarPcolor(rho, phi.*180/pi, abs(Hz)); shading flat; colormap('jet');
% title('Hz', 'FontSize', 12, 'FontWeight', 'bold');
% colormap('jet');
% figure;
% surface(rho, phi.*180/pi, abs(Erho)); shading flat; colormap('jet');
% 
% figure;
% surface(rho, phi.*180/pi, abs(Hphi)); shading flat; colormap('jet');



% figure;
% surface(rho, phi.*180/pi, abs(E)); shading flat; colormap('jet');
% 
% figure;
% surface(rho, phi.*180/pi, abs(H)); shading flat; colormap('jet');


figure;
polarPcolor(rho, phi.*180/pi, abs(E)); shading flat; colormap('jet');
title('E', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');
figure;
polarPcolor(rho, phi.*180/pi, abs(H)); shading flat; colormap('jet');
title('H', 'FontSize', 12, 'FontWeight', 'bold');
colormap('jet');

% figure;
% contour(rho, phi.*180/pi, abs(E));
% figure;
% contour(rho, phi.*180/pi, abs(H));
% figure;
% polarPcolor(rho, phi.*180/pi, abs(sqrt(abs(Hrho).^2 + abs(Hphi).^2))); shading flat; colormap('jet');
% title('H_{\rho} + H_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');



% figure;
% contour(rho, phi.* 180/pi, abs(E));
%    
% figure;
%    
% polarPcolor(rho, phi.*180/pi, abs(H)); shading flat;
% 
% figure;
% contour(rho, phi.* 180/pi, abs(H));
