function f = func_lecture3_example3(t,state)

global mu

X  = state(1:3);    % Position
V  = state(4:6);    % Velocity

R3 = (X(1)^2 + X(2)^2 + X(3)^2)^(3/2);

G  = -mu./R3.*X;    % Two-body forcing function

f  = [V; G];        % Output

return