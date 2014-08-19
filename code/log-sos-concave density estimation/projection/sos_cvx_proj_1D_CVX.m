% 
% 1-dimension sos-convex projection (CVX)
%
% Professor John Lafferty's Group
% Wei Hu
% The University of Chicago
%
% ------------------------------------------------------------------------------
% Aug 18, 2014
%
%

echo on;
d = 50; %   d: degree (even)
gamma = 10*randn(d+1, 1); % gamma: a random coefficient vector
%  we want to find the projection from gamma to sos-convex cone

% optimization problem:
%       minimize norm(theta-gamma)
%           s.t. f = dot(theta, vx) is sos-convex
% where vx is the vector of d+1 monomials

% write the optimization
cvx_begin
variable theta(d+1)
variable Q(d, d) semidefinite
expression b(d-1)
for k = 1:d-1
    b(k) = 0;
    for i = max(1, k+1-d/2):min(d/2, k)
        b(k) = b(k) + Q(i, k+1-i);
    end
    b(k) == (k+1)*k*theta(k+2)
end
minimize(norm(theta - gamma))
cvx_end

% obtain the polynomial and plot it
syms x; x = sym(x, 'real');
vx = sym('vx', [1 d+1]);
for i = 0:d
    vx(i+1) = x^i;
end
g = vpa(dot(gamma, vx));
f = vpa(dot(theta, vx));
ezplot(f, [-5, 5])
hold on
ezplot(g, [-5, 5])
hold off
echo off;