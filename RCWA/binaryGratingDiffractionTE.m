function [R, T, kIIzi] = binaryGratingDiffractionTE(nI, nGR, nRD, nII, Lambda, f, d, lambda, theta, N)
format short
% Returns the reflection and transmission efficiency of a binary grating
% using the rigurous coupled wave analysis technique given in the paper by
% Moharam et. al. in the paper, "Formulation for stable and efficient
% implementation of rigorous coupled-wave analysis of binary gratings,"
% JOSAA 12, 1068 (1995).
% 
% Inputs: nI - refractive index in incident plane
% nGR - refractive index in ridge
% nRD - refractive index in groove
% nII - refractive index in output plane
% Lambda - period of grating
% f - fraction of period occupied by ridge
% d - height of grating
% lambda - wavelength of incident light (in a vacuum)
% theta - angle down from normal (degrees)
% N - trunction order (odd)
% Outputs: R - reflection efficiencies
% T - transmission efficiencies
% [order, efficiency]

% Convert angle to radians
theta = theta * pi / 180;

% Set order array
if (mod(N,2) == 0)
N = N + 1;
display('Order not odd: Increased by 1 for calculation');
end
n = ((1-N):(N-1))';
m = (-(N-1)/2:(N-1)/2)';
m0 = (N-1)/2 + 1;

% Set up wave-vectors: Eq. 6, 7
k0 = 2 * pi / lambda;
kxi = k0 .* (nI * sin(theta) - m .* (lambda/Lambda));
kIzi = conj(k0 .* sqrt(nI^2 - (kxi/k0).^2));
kIIzi = conj(k0 .* sqrt(nII^2 - (kxi/k0).^2));

% Calculatet Fourier harmonics of dielectric in grating region
% Note: Matlab calculates sin(n*pi) = (+/-)3.6739E-16 for integer n instead
% of 0... could be a source of some error.
epsilon = (nRD^2-nGR^2) .* sin(pi * f .* n) ./ (pi .* n);
epsilon(N) = nRD^2 * f + nGR^2*(1-f);

% Set up the A matrix: Eq. 16
Kx = diag(kxi/k0);
E = zeros(N);
for j = 1:N
E(j, :) = epsilon(j:j+N-1);
end
E = fliplr(E);
A = Kx * Kx - E;

% Find eigenvalues and eigenvectors of A matrix
[W, Q] = eig(A, 'nobalance');
Q = sqrt(Q);
q = diag(Q);

% Set up matrix V
V = W * Q;

% Set up matrix equation for cp and cm: Eq. 19, 20, 22, and 23
YI = diag(kIzi/k0);
YII = diag(kIIzi/k0);
X = diag(exp(-k0*q*d));
L = [i*YI*W + V, (i*YI*W - V) * X; (V - i*YII*W)*X, -V - i*YII*W];
R = [zeros(m0-1, 1); i*(kIzi(m0)/k0 + nI * cos(theta)); zeros(N+m0-1, 1)];

% Solve equations
C = L\R;
% C = pinv(L) * R;

% Find Ri: Eq. 21
LR1 = [eye(N); -i*YI];
RR1 = [W, W * X; V, -V*X] * C - [zeros(m0-1, 1); 1; zeros(m0-1,1); zeros(m0-1,1); i*cos(theta)* nI; zeros(m0-1,1)];
RI = LR1\RR1;
DEri = RI .* conj(RI) .* real(kIzi / (k0 * nI * cos(theta)));
R = [m, DEri];

% Find Ti: Eq. 24
LT1 = [eye(N); i*YII];
RT1 = [W*X, W; V*X, -V] * C;
TI = LT1\RT1;
DEti = TI .* conj(TI) .* real(kIIzi / (k0 * nI * cos(theta)));
T = [m, DEti];