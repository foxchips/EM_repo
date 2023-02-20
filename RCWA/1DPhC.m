%-----Maharam et al Vol12, no.5, May 1995, J. Opt. Soc.b %Am.-------------
%-----Planar diffraction: TE polarization-------------------
clear all
n_1=2.5;  %region 1 where EM is incident and backward diffracted waves (reflected)
n_3=1;    %region 3 containing forward diffracted waves (transmitted)

NR_pitch=1e-6;
Lambda=NR_pitch;    %pitch of grating
K=2*pi/Lambda;  %grating vector
lambda=450e-9;  %free space wavelength of incident light
theta=0;    %angle of incidence relative to normal to grating interface/elevation
theta=theta*pi/180; %convert angle to radians
phi=0;  %"                                                            "/azimuth - if 0 this is planar diffraction - in which can can decompose into TE and TM-polarization, conical otherwise

ngr=n_1;
nrd=n_3;       %remember we are looking at this upside down!
d=1000e-9;      %depth of NRs
radius=250e-9;  %NR radius
ff=(NR_pitch-2*radius)/NR_pitch;    %fraction of period occupied by ridge
N=21;	%order no.

%test to check N is even
% odd_test=(-1*1)^N;
% if odd_test<0
%     N=N+1;                                                                                                      
% end
m0 = (N-1)/2 + 1;
m=(-(N-1)/2:(N-1)/2)';   %mth wave (equivalent to 'i' in paper) - ith diffracted wave
%wave vectors (equations 6 & 7)
k0=2*pi/lambda;
kxi=k0.*(n_1*sin(theta)-m.*(lambda/Lambda));  %kxi is determined from Floquet condition - change i to m to avoid confusing matlab
k1_zi=conj(k0.*(n_1.^2-(kxi./k0).^2).^0.5);    %k0n1>kxi i.e. only want to look at transmitted waves (waves which are not decaying in z-axis)
k3_zi=conj(k0.*(n_3.^2-(kxi./k0).^2).^0.5);    %k0n3>kxi

h=((1-N):(N-1))';   %hth Fourier component
epsilon = (nrd^2-ngr^2) .* sin(pi .* ff .* h) ./ (pi .* h);
epsilon(N)=nrd^2*ff+(ngr^2)*(1-ff);    %equ. (2) - set this for when h=0!

%Equation 16
Kx=diag(kxi/k0);    %Kx is diagonal matrix where m,m element equal to kxi/k0
E=zeros(N);	%E is matrix formed by permittivity harmonic components: i,pth element is epsilon(i-p) - get symmetric matrix
for j=1:N       %for each row
	E(j,:)=epsilon(j:j+N-1);    
end
E = fliplr(E);
A = Kx*Kx - E;

%find eigenvalues and eigenvectors of A matrix
[W,Q]=eig(A, 'nobalance');	%W= eigenvector matrix,     %solve coupled wave equations by finding eigenvalues and eigenvectors of A (d^2Sy/dz'^2=[A][Sy]
Q=sqrt(Q);	%Q= diagonal matrix with elements being positive square root of eigenvalues (hence square root)
q=diag(Q);

%Set up matrix V
V=W*Q;  %coefficient in space harmonics for magentic tangential fields

%Equations 19,20,22,23 i.e. boundary conditions - set up matrix equation
%for cm coeffs.
X=diag(exp(-k0*q*d));		%diagonal matrix with elements exp(-k0*qm*d)
Y1=diag(k1_zi/k0);
Y3=diag(k3_zi/k0);
L=[i*Y1*W+V,(i*Y1*W-V)*X;(V-i*Y3*W)*X, -V-i*Y3*W]; %1st row eliminates Ri from (19) and (20) and second row eliminates Ti from (22) and (23)
R = [zeros(m0-1, 1); i*(k1_zi(m0)/k0 + n_1 * cos(theta)); zeros(N+m0-1, 1)]; %LHS

C = L\R;    %find cm coeffs
%C = pinv(L) * R;
% Find Ri: Eq. 21
LR1 = [eye(N); -i*Y1];
RR1 = [W, W * X; V, -V*X] * C - [zeros(m0-1, 1); 1; zeros(m0-1,1); zeros(m0-1,1); i*cos(theta)* n_1; zeros(m0-1,1)];
RI = LR1\RR1;
DEri = RI .* conj(RI) .* real(k1_zi / (k0 * n_1 * cos(theta))); %equation (25)
R = [m, DEri];

% Find Ti: Eq. 24
LT1 = [eye(N); i*Y3];
RT1 = [W*X, W; V*X, -V] * C;
TI = LT1\RT1;
DEti = TI .* conj(TI) .* real(k3_zi / (k0 * n_1 * cos(theta))); %equation (25)
T = [m, DEti];





% odd_test=(-1*1)^N;
% if odd_test<0
%     N=N+1;                                                                                                      
% end




