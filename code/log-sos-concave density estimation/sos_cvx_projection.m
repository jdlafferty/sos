% 
% sos-convex projection
% 
% Professor John Lafferty's Group
% Wei Hu
% The University of Chicago
%
% ------------------------------------------------------------------------------
% Aug 15, 2014
% - Using SOSTOOLS to calculate the projection to the sos-convex cone
%

echo on;
p = 2; d = 4;  % p: dimension    d: degree (even)
s = nchoosek(p+d, d);  % s: number of monomials in p variables of degree at most d
gamma = 10*rand(s, 1); % gamma: a random coefficient vector
%  we want to find the projection from gamma to sos-convex cone

x = sym('x', [p 1]); x = sym(x, 'real');  % x: symbolic vector of length p
theta = sym('theta', [s 1]); theta = sym(theta, 'real'); % theta: the projection of gamma that we want to find
vx = monomials(x, [0:d]);  % vx is the symbolic vector of all s monomials
f = dot(theta, vx);   % f is the polynomial with coefficient theta
H = hessian(f, x);    % H is the hessian of x
                      % {H is sos} iff {f is sos-cvx}

% the optimization
%          minimize  norm(gamma-theta)
%              s.t.  H is sos
% is transformed to
%          minimize  t
%              s.t.  M = [I, theta; theta', 2*dot(gamma, theta) - dot(gamma,
%              gamma) + t] is psd (or equivalently, sos)
%                    H is sos

% the t and M below are desision variables that represent the ones in the
% above formulation
syms t; t = sym(t, 'real'); 
M = sym(eye(s+1));
M(1:s, s+1) = theta;
M(s+1, 1:s) = theta';
M(s+1, s+1) = 2*dot(gamma, theta) - dot(gamma, gamma) + t;

% set up the program
prog = sosprogram(x, [theta; t]);  % initialize the sos program
prog = sosmatrixineq(prog, H, 'Mineq');  % add constraint: H is sos
prog = sosmatrixineq(prog, M, 'Mineq');  % add constraint: M is sos
prog = sossetobj(prog, t);  % add objective t
prog = sossolve(prog);   % solve the program
f = sosgetsol(prog, f)   % get f
theta = sosgetsol(prog, theta); % get theta

% if p = 1 or 2, plot the resulting polynomial
if p == 1
    ezplot(f, [-10, 10])
elseif p ==2
    ezmesh(f, [-10, 10])
end

% for verification, check if f is sos-convex
%findsos(hessian(f, x))
echo off;