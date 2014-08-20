% 
% 2-dimension sos-convex projection (CVX)
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
d = 10; %   d: degree (even)
s = (d+2)*(d+1)/2;   % s: number of monomials in 2 variables of degree up to d

gamma = 10*rand(s, 1); % gamma: a random coefficient vector
%  we want to find the projection from gamma to sos-convex cone

x = sym('x', [2 1]); x = sym(x, 'real');
vx = monomials(x, [0:d]);  % vx: the vector of all s monomials

% optimization problem:
%       minimize norm(theta-gamma)
%           s.t. f = dot(theta, vx) is sos-convex
% where vx is the vector of all s monomials

% write the optimization
cvx_begin
% two different ways to store the coefficient: 1d-array and 2d-array
% f = dot(theta_vector, x) = sum(theta_matrix(i, j)*x_1^i*x_2^j)
variable theta_matrix(d+1, d+1) % theta_matrix: 2d index
variable theta_vector(s)   % theta_vector: 1d index
% make sure the two ways represent the same polynomial
k = 1;
for i = 0:d
    for j = i:-1:0
        theta_vector(k) == theta_matrix(j+1, i-j+1)
        k = k + 1;
    end
end

% details of writing the sos constraint as semidefinite constraint
r = (d/2+1)*d/4;
variable Q(2*r, 2*r) semidefinite
expression b11(d-1, d-1)
expression b12(d-1, d-1)
expression b22(d-1, d-1)
b11(:, :) = zeros(d-1, d-1)
b12(:, :) = zeros(d-1, d-1)
b22(:, :) = zeros(d-1, d-1)
for i1 = 0:d/2-1
    for j1 = 0:d/2-1-i1
        k1 = (i1+j1)*(i1+j1+1)/2 + (j1+1);
        for i2 = 0:d/2-1
            i = i1 + i2;
            for j2 = 0:d/2-1-i2
                j = j1 + j2;
                k2 = (i2+j2)*(i2+j2+1)/2 + (j2+1);
                b11(i+1, j+1) = b11(i+1, j+1) + Q(k1, k2);
                b12(i+1, j+1) = b12(i+1, j+1) + Q(r+k1, k2) + Q(k1, r+k2);
                b22(i+1, j+1) = b22(i+1, j+1) + Q(r+k1, r+k2);
            end
        end
    end
end
for i = 0:d-2
    for j = 0:d-2-i
        b11(i+1, j+1) == (i+2)*(i+1)*theta_matrix(i+3, j+1)
        b12(i+1, j+1) == 2*(i+1)*(j+1)*theta_matrix(i+2, j+2)
        b22(i+1, j+1) == (j+2)*(j+1)*theta_matrix(i+1, j+3)
    end
end
minimize(norm(theta_vector - gamma))
cvx_end

% obtain the polynomial and plot it
f = vpa(dot(theta_vector, vx));
g = vpa(dot(gamma, vx));
ezmesh(f, [-5, 5])
hold on;
ezmesh(g, [-5, 5])
hold off;
echo off;
