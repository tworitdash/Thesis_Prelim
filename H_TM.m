function [Hrho, Hphi, Hz] = H_TM(mu, m, rho, phi, beta_rho, z, beta)
    B = 1;
   
    C = 0;
    D = 1;
    
    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    
    Hrho = B .* (m ./ (mu .* rho)) .* besselj(m, beta_rho .* rho).* (-C .* sin(m .* phi)...
        + D .* cos(m .* phi)) * exp(-1j .* beta_z .* z);
    Hphi = -B .* (beta_rho ./ mu) .* besselj_der(m, beta_rho .* rho) .* (C .* cos(m .* phi)...
        + D .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
    Hz = 0;
end