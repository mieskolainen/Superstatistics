% Non-linear differential equation test system
% 
% d^2x/dt^2 + 0.6dx/dt + 3x + x^2 = 0
% 
function output = F(t,x)

output = [ x(2);
          -0.6*x(2) - 3.0*x(1) - x(1)^2];

end