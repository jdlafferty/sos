% 
% 1-dimension sos-convex regression (CVX)
% 
% Professor John Lafferty's Group
% Wei Hu
% The University of Chicago
%
% ------------------------------------------------------------------------------
% Aug 15, 2014
%
%

echo on;
n = 100;   % n: number of data points
xx = 10*rand(1, n) - 5;    % randam data  i.i.d.~U[-5, 5]
y = exp(xx) + 0.5*randn(1, n);  % y = exp(x) + Gaussian error
plot(xx, y, 'xr')  % plot data
hold on
d = 10;  % d: degree of polynomial

% optimization problem:
%       minimize norm(y-z)
%           s.t. V*theta = z
%                f = dot(theta, vx) is sos-convex
%
% where vx is the vector of d+1 monomials, V is the matrix of the vx's at
% the n data points

% calculate V
for i = 1:n
    for j = 0:d
        V(i, j+1) = xx(i)^j;
    end
end
V = double(V);

% write the optimization
cvx_begin
variables theta(d+1) z(n)
variable Q(d, d) semidefinite
expression b(d-1)
for k = 1:d-1
    b(k) = 0;
    for i = max(1, k+1-d/2):min(d/2, k)
        b(k) = b(k) + Q(i, k+1-i);
    end
    b(k) == (k+1)*k*theta(k+2)
end
V*theta == z
minimize(norm(y' - z))
cvx_end

% obtain the polynomial and plot it
syms x; x = sym(x, 'real');
vx = sym('vx', [1 d+1]);
for i = 0:d
    vx(i+1) = x^i;
end
f = dot(theta, vx);
ezplot(f, [-5, 5])
hold off;
echo off;