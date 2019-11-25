%% Normalization matrix with different modes on the waveguide:
clear;
close all;

m = 1:1:50; % first digit of the mode number
N = 2; % second digit of the mode number
%N = 3;
mode = "TE"; % Waveguide mode polarization
% mode = "TM"

F = 10e9;

r = 0.0405319403216/2; % radius of the waveguide
er = 1; % relative  permittivity
mur = 1; % relative Permeability
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;

drho = r/1000;
dphi = pi/2000;

[rho_, phi_] = meshgrid(eps:drho:r, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide

z = 0; 
% Q = zeros(size(F, 2), N(end), N(end));
Q = zeros(m(end), m(end));
Z = zeros(1, size(F, 2));
Y = zeros(1, size(F, 2));

for k = 1:length(F)

for i = 1:length(m)
    [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z] = E_and_H(rho_, phi_, er, mur, z, r, m(i), N, mode, F(k));
    
    Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho_ * drho .* dphi;
    Qij = sum(sum(Poyn));
    
    Q(i, i) = Qij;
%     if mode == "TE"
%         z = 2 * pi * F(k) * mu./ beta_z;
%     elseif mode == "TM"
%         z = beta_z ./ (2 * pi * F(k) .* epsilon);
%     end
%     Z(k) = z;
%     Y(k) = 1./Z(k);
end
end

figure;
plot(m, db(abs(diag(Q)))/10, 'LineWidth', 2); grid on;

xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalization Constant Q_{1, n}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');

figure;
plot(m, real(diag(Q)), 'LineWidth', 2); grid on;
hold on;
plot(m, imag(diag(Q)), 'LineWidth', 2); grid on;

xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalization Constant Q_{1, n}', 'FontSize', 12, 'FontWeight', 'bold');
title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');
legend({'Re(Q)', 'Im(Q)'}, 'FontSize', 12, 'FontWeight', 'bold');

% figure;
% 
% plot(F * 1e-9, real(Z), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Z), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Impedance Z (\Omega)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Impedance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Z)', 'Im(Z)'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% plot(F * 1e-9, real(Y), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Y), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Admittance Y (Mho)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Admittance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Y)', 'Im(Y)'}, 'FontSize', 12, 'FontWeight', 'bold');
