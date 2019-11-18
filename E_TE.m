function [Erho, Ephi, Ez] = E_TE(epsilon, m, rho, phi, beta_rho, z, beta)
    A = 1;
   
    C = 0;
    D = 1;
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
    Erho = -A .* m./(epsilon .* rho) .* besselj(m, beta_rho .* rho) .* (-C .* sin(m .* phi)...
        + D .* cos(m .* phi)) * exp(-1j .* beta_z .* z);
    Ephi = A .* beta_rho./epsilon .* besselj_der(m, beta_rho .* rho) .* (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ez = 0;
end