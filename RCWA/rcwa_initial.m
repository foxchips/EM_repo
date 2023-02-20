%Vol. 71, No. 7/July 1981/J.Opt.Soc.Am., Moharam and Gaylord
%Assume 

epsilon_1=2.4^2;  %region 1 where EM is incident and backward diffracted waves (reflected)
epsilon_3=1^2;    %region 3 containing forward diffracted waves (transmitted)

%refractive index of grating
NR_radius=140e-9;
NR_pitch=1000e-9;
NR_height=1000e-9;
tot=NR_pitch*NR_height;
Area_NR=2*(NR_radius*NR_height);
FF_NR=Area_NR/tot;
Area_gap=(NR_pitch-2*NR_radius)*NR_height;
FF_gap=Area_gap/tot;
epsilon_2=(FF_NR*epsilon_1+FF_gap*epsilon_3);   %average refractive index

Lambda=NR_pitch;    %pitch of grating
K=2*pi/Lambda;  %grating vector
lambda=450e-9;  %free space wavelength of incident light
theta=0;    %angle of incidence relative to normal to grating interface/elevation
phi=0;  %"                                                            "/azimuth - if 0 this is planar diffraction - in which can can decompose into TE and TM-polarization, conical otherwise

k_0=2*pi/lambda;    %free space wavevector
k_1=2*pi*epsilon_1^2/lambda;

no_waves=20;    %no. of plane waves to approximate field

for i=-1*no_waves/2:no_waves/2
    beta(i+)=k_1*sin(theta)-i*K*sin(phi);
end

                                                                                                                                                                                                                                                                                                                                                

