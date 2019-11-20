%% Normalization matrix with different modes on the waveguide:
%clear;
%close all;

m = 5; % first digit of the mode number
N = 1:1:100; % second digit of the mode number

mode = "TE"; % Waveguide mode polarization
% mode = "TM"

r = 0.0405319403216/2; % radius of the waveguide

drho = r/1000;
dphi = pi/2000;

[rho_, phi_] = meshgrid(eps:drho:r, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide

z = 0; 

Q = zeros(N(end), N(end));

for i = 1:length(N)
    [Erho, Ephi, Ez, Hrho, Hphi, Hz] = E_and_H(rho_, phi_, er, mur, z, r, m, N(i), mode);
    
    Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho * drho .* dphi;
    Qij = sum(sum(Poyn));
    
    Q(i, i) = Qij;
end
hold on;
plot(N, db(abs(diag(Q)))/10, 'LineWidth', 2); grid on;

xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalization Constant Q_{1, n}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Normalization Constant for TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');