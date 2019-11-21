function [Erho, Ephi, Ez] = E_TE(epsilon, m, rho, phi, beta_rho, z, beta)
    A = 1;
   
    C1 = 1; % C1 and D1 are for the rho component
    D1 = 0; 
%     
%     C2 = 0; % C2 and D2 are for the phi component
%     D2 = 1;
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
    Erho = -A .* m./(epsilon .* rho) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
        + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ephi = A .* beta_rho./epsilon .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
        + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ez = 0;
end