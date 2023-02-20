%rcwa_tables

clear all
close all


%Set up physical parameters of waveguide
ldnms=[460];%wavelength in nm
Tx=335;%x-period in nm (skew space, without correcting for cosd(zeta))
Ty=335;%y-period in nm (skew space, without correcting for cosd(zeta))
MS=70;%MS 
Notch=20;%Notch
%MS=[35,10];%[X GAP]=MS for cross structure
%Notch=[];%'Notch=[]' for cross instead of diamond
height_grating_nms=[80];%Height of glass grating in nm
thickness_coating_nm=0;%TiO2 coating thickness in nm
MN=[8 8]; %Number of orders to use for RCWA truncation (x,y)
zeta=30;%skew angle
%theta=30:1:88;%theta values
theta=34:88;%theta values - HYPERFOV
back_inv_structure=0;%'0' if calculating usual geometry, '1' if calculating inverted geometry (e.g. for back transmission) 




output_folder_prefix='TEST2_';%Starting character string for the output folder
output_folder_suffix='_SH_SK_0SW0';%Ending character string for the output folder
rdrive_destination_folder='E:\rcwa\results';%The destination folder in R drive, e.g. ~ 'R:\Waveguides\designs\output_gratings\akula02\hiRI'
%rdrive_destination_folder='';%The destination folder in R drive, e.g. ~ 'R:\Waveguides\designs\output_gratings\akula02\hiRI'
%pattern_input_file='AK02B_Im1D1.csv';%Provide a binary csv file that gives the shape of the structure. If not doing this, just leave it blank.
pattern_input_file='';
%resin_nk_file='C:\rcwa_versions\NKData\17BLANC.nk';
resin_nk_file='17BLANC.nk';

if back_inv_structure==1
	theta=0:1:88;%theta values
	output_folder_suffix=strcat(output_folder_suffix,'_INV');
end

%The variable phi_order takes up to six inputs. It tells the program what kind of simulation you would like to do.
%If phi_order{1}=  'tables_with_symmetry', the program runs the full tables taking into account symmetry.
%If phi_order{1}=  'tables_without_symmetry', the program runs the full tables without using any symmetry.
%If phi_order{1}=  30, the program runs the curves for phi=30 and returns the efficiencies for orders specified in phi_order{6} in an Excel sheet.

%%%%%%%%%%%%%%%%%
%{
phi_order={'tables_with_symmetry',output_folder_prefix,output_folder_suffix,rdrive_destination_folder,pattern_input_file};
phi_order={'tables_without_symmetry',output_folder_prefix,output_folder_suffix,rdrive_destination_folder,pattern_input_file};
phi_order={30,'file_prefix','file_suffix','',pattern_input_file,[-1,-1;-1 0]};
phi_order={30,'test','','',pattern_input_file,[-1 -1;-1 1]};
%}
%%%%%%%%%%%%%%%%%

%phi_order={'tables_with_symmetry',output_folder_prefix,output_folder_suffix,rdrive_destination_folder,pattern_input_file,resin_nk_file};
%phi_order={'tables_without_symmetry',output_folder_prefix,output_folder_suffix,rdrive_destination_folder,pattern_input_file};

phi_order={30,'test','','','',[0 0;-1 -1;-1 0;0 -1; 1 0; 0 1; 1 1]};

for ii=1:length(ldnms)

	ldnm=ldnms(ii)

	for jj=1:length(height_grating_nms)
		height_grating_nm=height_grating_nms(jj)
		switch back_inv_structure
			case 0
				rcwa_data_phi_theta_sym(ldnm,Tx,Ty,MS,Notch,height_grating_nm,thickness_coating_nm,MN,zeta,theta,phi_order);
				%rcwa_data_phi_theta_sym_psi(ldnm,Tx,Ty,MS,Notch,height_grating_nm,thickness_coating_nm,MN,zeta,theta,phi_order,45);
			case 1				
				rcwa_data_phi_theta_sym_invgeom(ldnm,Tx,Ty,MS,Notch,height_grating_nm,thickness_coating_nm,MN,zeta,theta,phi_order);
		end
	end
end

