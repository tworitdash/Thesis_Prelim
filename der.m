z = eps:0.1:10;

m = 0:1:5;

for i = 1:length(m)
    der_ = (besselj(m(i) - 1, z) - besselj(m(i) + 1, z))./2;
    
    dz = 1e-5;
   
    bes = besselj(m(i), z);
    hold on;
    plot(z, der_, 'LineWidth', 2);
%     plot(z, bes, 'LineWidth', 2);
    hold on;
    grid on;
end
