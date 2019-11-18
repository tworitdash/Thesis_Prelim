function [Erho, Ephi, Ez] = E_TM(epsilon, m, rho, phi, beta_rho, z, beta, omega, mu)
    B = 1;
   
    C = 1;
    D = 0;
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
    Erho = -B .* (beta_rho .* beta_z ./ (omega .* mu .* epsilon)) .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
    Ephi = -B .* (m .* beta_z ./ (omega .* mu .* epsilon .* rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
        + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
    Ez = -1j .* B .* (beta_rho.^2 ./ (omega .* mu .* epsilon)) .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
end