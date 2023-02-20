function [DErPP,DErSS,DErSP,DErPS,order]=rcwa_SAF(ldnm,T_noskewx,T_noskewy,MS,Notch,height_grating_nm,thickness_coating_nm,MN,zeta,theta,phi)

tic
MN=[8 8];
ld=ldnm*1e-9;
k0=2*pi/ld;

n_air=1;

%x1p7=[400 459 480 486.1 546.1 587.6 643.8 656.3]*1e-9;
%n1p7=[1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5];

xn1p7=dlmread('17BLANC.nk','',2,0);
x1p7=xn1p7(:,1)*1e-9;
n1p7=complex(xn1p7(:,2),xn1p7(:,3));

n_glass=interp1(x1p7,n1p7,ld,'linear');

nxt=[n_glass,n_air];

scross=0;
if isempty(Notch)
    scross=1;
	MSd=MS(1);
	Notch=MS(2);
	MS=MSd;
	disp('RCWA for Cross Structure... ... ')
else
	disp('RCWA for Diamond-Notch Structure... ... ')
end
	

MSr=MS/100;
Notchr=Notch/100;

Tx=T_noskewx*1e-9/cosd(zeta);
Ty=T_noskewy*1e-9/cosd(zeta);

	
switch zeta
    case 0
        dphys=[Tx*tand(30) Ty]*2;
		szr=sprintf('Geometry in Orthogonal Space: Tx=%.3f nm, Ty=%.3f nm', dphys(1)*1e9, dphys(2)*1e9);
            
        [dpts, npts, tns, texture, ~]=diamond_ortho(nxt,height_grating_nm,MSr,Notchr,thickness_coating_nm,dphys,ld);

        Tnm_x=round(dphys(1)*1e9);
        Tnm_y=round(dphys(2)*1e9);
        Tnm=[Tnm_x Tnm_y];	
        [sx, sy]=retordre(2,2);
		ox=sx-sy;oy=sx+sy;

    	retnxy=[ox oy];
		iptrxy=nan(size(retnxy,1),1);
			
			
    case 30
        dphys=[Tx,Ty];
		szr=sprintf('Geometry in Skew Space: Tx=%.3f nm, Ty=%.3f nm', dphys(1)*1e9, dphys(2)*1e9);
		[dpts, npts, tns, texture, ~]=diamond_notch(nxt,height_grating_nm,MSr,Notchr,thickness_coating_nm,dphys,ld,scross);
									
        if (T_noskewx-T_noskewy)==0
            Tnm=T_noskewx;
        else

            Tnm=[T_noskewx T_noskewy];	
        end
			
        [cx, cy]=retordre(2,2);
		retnxy=[cx cy];
		iptrxy=nan(size(retnxy,1),1);	
			
end
	
disp(szr);

m=MN(1);
n=MN(2);

M=2*m+1;
N=2*n+1;

psi=[0 90];

make_sym=0;

phis=phi;
thetas=theta;

[cx, cy]=retordre(m,n);

for ii=1:length(iptrxy)
    aex=find(cx==retnxy(ii,1));
    aey=find(cy(aex)==retnxy(ii,2));
    iptrxy(ii)=min(aex)+aey-1;
end

retainord=iptrxy;
    
switch length(phis)

	case 1
		DEr=nan(length(iptrxy),length(thetas));
		DEt=DEr;
		phi=phis;
		
    for ii=1:length(thetas)
        theta=thetas(ii);	
        if abs(theta)<1e-10
            theta=1e-10;
        end
				
		kxi=k0*(nxt(1)*sind(theta)*cosd(phi) + cx.*ld/dphys(1));
		kx=kxi./k0;
		ky0=k0*nxt(1)*sind(theta)*sind(phi+zeta);
		kyi=k0*(nxt(1)*sind(theta)*sind(phi+zeta) + cy.*ld/dphys(2));
		ky=kyi/k0;			
		kxsq=kx.^2;
		kysq=ky.^2;	
		Kx=diag(kx);
		Ky=diag(ky);	

		kxky=0;
		if zeta~=0
            kxky=kx.*ky;		
		end
			
		fsum=(kxsq+kysq-2*sind(zeta)*kxky)*secd(zeta)*secd(zeta);

		k1zi=psqrt(nxt(1)^2-fsum);
		kNzi=psqrt(nxt(2)^2-fsum);
			

		[E A W1 W2 V11 layers V21 V22 q1m]=rediag_linear_conicalfull(dpts,npts,dphys,[M,N],Kx,Ky,texture,zeta);

			
			
		[Ru,Td,DEr_count,DEt_count,summe,DErPP,DErPS,DErSP,DErSS,DEtPP,DEtPS,DEtSP,DEtSS]=rezeta(W1,W2,Kx,Ky,q1m,k1zi,kNzi,tns,nxt,retainord,theta,psi,phi,zeta,k0);			
			
			
		DEr(:,ii)=DEr_count;
		DEt(:,ii)=DEt_count;
			
				
				
			
    end%ii
	
    otherwise
        DEr=nan(length(iptrxy),length(phis),length(thetas));
		DEt=DEr;
        
        DErPP=DEr;
        DErPS=DEr;
        DErSP=DEr;
        DErSS=DEr;
        
        DEtPP=DEr;
        DEtPS=DEr;
        DEtSP=DEr;
        DEtSS=DEr;
        for ii=1:length(thetas)
            theta=thetas(ii);
				
			
            for jj=1:length(phis)	
                phi=phis(jj);
				
				kxi=k0*(nxt(1)*sind(theta)*cosd(phi) + cx.*ld/dphys(1));
				kx=kxi./k0;
				ky0=k0*nxt(1)*sind(theta)*sind(phi+zeta);
				kyi=k0*(nxt(1)*sind(theta)*sind(phi+zeta) + cy.*ld/dphys(2));
				ky=kyi/k0;


				
				kxsq=kx.^2;
				kysq=ky.^2;	
				Kx=diag(kx);
				Ky=diag(ky);	

				kxky=0;
				if zeta~=0
					kxky=kx.*ky;		
				end
				
				fsum=(kxsq+kysq-2*sind(zeta)*kxky)*secd(zeta)*secd(zeta);

				k1zi=psqrt(nxt(1)^2-fsum);
				kNzi=psqrt(nxt(2)^2-fsum);
				

				[E A W1 W2 V11 signF V21 V22 q1m]=rediag_linear_conicalfull(dpts,npts,dphys,[M,N],Kx,Ky,texture,zeta);
				
				
				
				[Ru,Td,DEr_count,DEt_count,summe,DErPP0,DErPS0,DErSP0,DErSS0,DEtPP0,DEtPS0,DEtSP0,DEtSS0]=rezeta(W1,W2,Kx,Ky,q1m,k1zi,kNzi,tns,nxt,retainord,theta,psi,phi,zeta,k0);
				
				
				
				
				DEr(:,jj,ii)=DEr_count;
				DEt(:,jj,ii)=DEt_count;
                   
                DErPP(:,jj,ii)=DErPP0;
                DErPS(:,jj,ii)=DErPS0;
                DErSP(:,jj,ii)=DErSP0;
                DErSS(:,jj,ii)=DErSS0;
                    
                    
				DEtPP(:,jj,ii)=DEtPP0;
                DEtPS(:,jj,ii)=DEtPS0;
                DEtSP(:,jj,ii)=DEtSP0;
                DEtSS(:,jj,ii)=DEtSS0;
				
				
				
            end	%jj
        end%ii

end%switch length(phis)


     toc
		


switch length(phis)

	case 1
	
		pr=DEr.';
		pt=DEt.';
		
		
				
		%outfename='test1';
		
		
		tabname=strcat('phi_',num2str(phis));
		
		
		snxy=cell(1,size(retnxy,1));
		for ii=1:numel(snxy)
			snxy{ii}=strcat('[',num2str(retnxy(ii,1)),',',num2str(retnxy(ii,2)),']');
		end
			
			
			pr1=num2cell([thetas.' pr]);			
			
    
     
	
	otherwise
	
		pr=permute(DEr,[3 2 1]);
		pt=permute(DEt,[3 2 1]);
        
        prPP=permute(DErPP,[3 2 1]);
        prPS=permute(DErPS,[3 2 1]);
        prSP=permute(DErSP,[3 2 1]);
        prSS=permute(DErSS,[3 2 1]);
        
        ptPP=permute(DEtPP,[3 2 1]);
        ptPS=permute(DEtPS,[3 2 1]);
        ptSP=permute(DEtSP,[3 2 1]);
        ptSS=permute(DEtSS,[3 2 1]);
		
		%save D0x_before DEr pr DEt pt phis thetas MS Notch  retainord cx cy retnxy
		
		
		
		if make_sym==1
			qr=make_all_phi(pr,retnxy,phis,thetas);
			qt=make_all_phi(pt,retnxy,phis,thetas);
            
            qrPP=make_all_phi(prPP,retnxy,phis,thetas);
            qrPS=make_all_phi(prPS,retnxy,phis,thetas);
            qrSP=make_all_phi(prSP,retnxy,phis,thetas);
            qrSS=make_all_phi(prSS,retnxy,phis,thetas);
            
            
            qtPP=make_all_phi(ptPP,retnxy,phis,thetas);
            qtPS=make_all_phi(ptPS,retnxy,phis,thetas);
            qtSP=make_all_phi(ptSP,retnxy,phis,thetas);
            qtSS=make_all_phi(ptSS,retnxy,phis,thetas);
            
            
			phis=[-180:180];
			pr=qr;pt=qt;
            prPP=qrPP;prPS=qrPS;prSP=qrSP;prSS=qrSS;
            ptPP=qtPP;ptPS=qtPS;ptSP=qtSP;ptSS=qtSS;
		end

		
		%save D0x_after DEr pr DEt pt phis thetas MS Notch  retainord cx cy retnxy
		
		%save D06_blue_phi_theta_nosymm DEr pr DEt pt phis thetas MS Notch  retainord cx cy retnxy
	
		%save D05_nosym_ortho_test_blue DEr pr DEt pt prPP prPS prSP prSS phis thetas MS Notch  retainord cx cy retnxy
	
    for ii=1:length(retainord)
			order(ii,1)=retnxy(ii,1);
            order(ii,2)=retnxy(ii,2);
			fname=strcat(parentdir,'\',outdirld,'\',ename);
			
			           
    end
		

      
		%Convert ortho tables to skew tables if zeta=0
        if zeta==0
					

			outdir_ortho_ld=strcat(outdir_ortho,'\',num2str(ldnm));  

			mkdir(outdir_ortho_ld);			
		
			[skex, skey]=retordre(2,2);
			[qr,~]=convert_tables_ortho_to_skew(pr,pt);
            [qrPP,~]=convert_tables_ortho_to_skew(prPP,ptPP);
            [qrPS,~]=convert_tables_ortho_to_skew(prPS,ptPS);
            [qrSP,~]=convert_tables_ortho_to_skew(prSP,ptSP);
            [qrSS,~]=convert_tables_ortho_to_skew(prSS,ptSS);
            
			
				

        
        end
        
        
end

end%rcwa_data_phi_theta



function [dpts, npts, tns, texture, layers]=diamond_notch(nxt,H,MS,Notch,T,dphys,ld,cross)


if nargin<8
	cross=0;
end

n_glass=nxt(1);
n_air=nxt(2);

if T==0
	layii=1;
	layers=2+layii;
	layer_thickness=H*1e-9;

else
	layii=3;	
	layer_thickness=[T H-T T]*1e-9;
	
	xTiO2=[400 460 530 620 700]*1e-9;
	nTiO2=[2.5 2.5 2.42 2.36 2.36];
	n_TaO=interp1(xTiO2,nTiO2,ld,'linear');

end



layers=2+layii;
tns=zeros(1,layers);
tns(2:layers-1)=layer_thickness;

xpts=cell(1,layers-1);
npts=xpts;texture=xpts;

Tx=dphys(1);
Ty=dphys(2);


if cross==1
	width= MS;
	gap = Notch;
	
	obx=width*Tx;
	oby=(1-2*gap)*Ty;
	
	alp=0.5;

 
else 
	
	obx=MS*Tx;
	oby=MS*Ty;

	alp=oby/(obx+oby);
	
	onx=Notch*Tx;
	ony=Notch*Ty;

end
alpha=ones(1,layers)*alp;
show_geom_xy=0;

xpts{1}=Tx/2;ypts{1}=Ty/2;npts{1}=nxt(1);texture{1}=[nxt tns show_geom_xy];
switch cross
	case 0
		if T==0
				if Notch==0
				
					ii=2;
					xpts{ii}=Tx/2;ypts{ii}=Ty/2;npts{ii}=n_air;
					texture{ii}={[0,0,obx,oby,n_glass,1]};
				
				else
					ii=2;
					xpts{ii}=Tx/2;ypts{ii}=Ty/2;npts{ii}=n_air;
					texture{ii}={[0,0,obx,oby,n_glass,1],[-obx/2+onx/2,oby/2-ony/2,onx,ony,n_air,1],[obx/2-onx/2,-oby/2+ony/2,onx,ony,n_air,1]};
				
				
				end


		else
				if Notch==0

					ss=2;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_TaO;		
					texture{ss}={[0,0,obx,oby,n_glass,1]};
					
							
					ss=3;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
					texture{ss}={[0,0,obx,oby,n_glass,1]};
					
					
					ss=4;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
					texture{ss}={[0,0,obx,oby,n_TaO,1]};
					
				else
					
				
					ss=2;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_TaO;		
					texture{ss}={[0,0,obx,oby,n_glass,1],[-obx/2+onx/2,oby/2-ony/2,onx,ony,n_TaO,1],[obx/2-onx/2,-oby/2+ony/2,onx,ony,n_TaO,1]};
					
							
					ss=3;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
					texture{ss}={[0,0,obx,oby,n_glass,1],[-obx/2+onx/2,oby/2-ony/2,onx,ony,n_air,1],[obx/2-onx/2,-oby/2+ony/2,onx,ony,n_air,1]};
					
					
					ss=4;
					xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
					texture{ss}={[0,0,obx,oby,n_TaO,1],[-obx/2+onx/2,oby/2-ony/2,onx,ony,n_air,1],[obx/2-onx/2,-oby/2+ony/2,onx,ony,n_air,1]};
				end
				
				
				
		end
		
	case 1
		if T==0
			ii=2;
			xpts{ii}=Tx/2;ypts{ii}=Ty/2;npts{ii}=n_air;
			texture{ii}={[0,0,obx,oby,n_glass,1],[0,0,oby,obx,n_glass,1]};
		else	
			ss=2;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_TaO;		
			texture{ss}={[0,0,obx,oby,n_glass,1],[0,0,oby,obx,n_glass,1]};
			
					
			ss=3;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[0,0,obx,oby,n_glass,1],[0,0,oby,obx,n_glass,1]};
			
			
			ss=4;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[0,0,obx,oby,n_TaO,1],[0,0,oby,obx,n_TaO,1]};
		end
end

dpts.x=xpts;dpts.y=ypts;dpts.alpha=alpha;
end%diamond_notch


function [dpts, npts, tns, texture, layers]=diamond_ortho(nxt,H,MS,Notch,T,dphys,ld)

n_glass=nxt(1);
n_air=nxt(2);

if T==0
	layii=1;
	layers=2+layii;
	layer_thickness=H*1e-9;

else
	layii=3;	
	layer_thickness=[T H-T T]*1e-9;
	
	xTiO2=[400 460 530 620 700]*1e-9;
	nTiO2=[2.5 2.5 2.42 2.36 2.36];
	n_TaO=interp1(xTiO2,nTiO2,ld,'linear');

end



layers=2+layii;
tns=zeros(1,layers);
tns(2:layers-1)=layer_thickness;

xpts=cell(1,layers-1);
npts=xpts;texture=xpts;

zeta=30;

ct=cosd(zeta);st=sind(zeta);tt=tand(zeta);

Tx=dphys(1);
Ty=dphys(2);

D=dphys(2)/2;Dy=D;Dx=dphys(1)/(2*tt);

if abs(Dx-Dy)>1e-8
 disp('Period in x and y not equal!');return;
end

len=MS*Ty;
notchw=2*D/ct*(MS*st-Notch);
diaw=2*MS*D*tt;
linew=(1-MS)*D/ct;

nDp=Notch*D/ct;
nDperp=nDp*cosd(30);
	
obx=diaw;
oby=len;

alp=oby/(obx+oby);

alpha=ones(1,layers)*alp;

yoff=-len/2;

xz=[0 -diaw/2 0];
yz=[0 len/2 len]+yoff;

xzn=[0 -diaw/2+nDp/2 -notchw/2 -diaw/2+nDp/2 0];
yzn=[0 len/2-nDperp len/2 len/2+nDperp len]+yoff;



show_geom_xy=0;

xpts{1}=Tx/2;ypts{1}=Ty/2;npts{1}=nxt(1);texture{1}=[nxt tns show_geom_xy];



if T==0,
		
		
			ii=2;
			xpts{ii}=Tx/2;ypts{ii}=Ty/2;npts{ii}=n_air;
			texture{ii}={[MS,Notch,obx,oby,n_glass,2]};
		
		


else
		

			ss=2;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_TaO;		
			texture{ss}={[MS,Notch,obx,oby,n_glass,2]};
			
					
			ss=3;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[MS,Notch,obx,oby,n_glass,2]};
			
			
			ss=4;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[MS,Notch,obx,oby,n_TaO,2]};
			
			
		
end

dpts.x=xpts;dpts.y=ypts;dpts.alpha=alpha;
end%diamond_ortho

function [cx cy]=retordre(m,n)
am=-m:m;
an=-n:n;
[xam yan]=meshgrid(am,an);
cx=xam(:);
cy=yan(:);
end%retordre

function [E A W WV V11 layers V21 V22 q1m]=rediag_linear_conicalfull(dpts,npts,dphys,MN,Kx,Ky,texture,zeta)

if nargin==7
zeta=0;
end

n=prod(MN);

layers=numel(npts)+1;



q1m=nan(2*n,layers-1);
W=nan(2*n,2*n,layers-1);WV=W;

V11=nan(n,2*n,layers-1);V22=V11;V12=V11;V21=V11;


switch length(dphys)
    case 1
       if length(MN)==1
           E=nan(n,n,layers-1);A=E;
           for ii=1:layers-1
               [E(:,:,ii) A(:,:,ii)]=epsfep(dpts.x{ii},npts{ii},dphys,n);
           end
       else
           disp('Y order is defined for 1D X-grating problem!');return;
       end
       
     
              
    case 2
        if length(MN)==2
            
            txte=texture{2}{1};
            if length(txte)==6
				optn=txte(6);
			else
				optn=texture{2}{6};
            end
			
            
            switch optn
                case 1
            
                    [E A]=retoep(dpts,npts,texture,dphys,MN);	
                case {2,9}
                    [E A]=retoep2(dpts,npts,texture,dphys,MN);
				case 3
					[E A]=retoep3(dpts,npts,texture,dphys,MN);
					
            end
                
            
        else
           disp('X or Y order undefined for 2D problem!');return;
        end    
        
end



fgm=1;%'0' -- build FG matrix analytically by parts; '1' -- build FG matrix using F*G 

Id=eye(n);

for ii=1:layers-1	
	
	EIX=E(:,:,ii)\Kx;
	EIY=E(:,:,ii)\Ky;
	
	
	B=Kx*EIX - Id;
	D=Ky*EIY - Id;
	if length(MN)==1
		alpha=1;
	else
		alpha=dpts.alpha(ii);
	end
	
	Ainv=A(:,:,ii)\Id;
	
	tildeE=alpha*E(:,:,ii)+(1-alpha)*Ainv;
	barE=alpha*Ainv + (1-alpha)*E(:,:,ii);
	
	%tildeE=E(:,:,ii);barE=tildeE;
	
	
	
	OA11 = Ky*EIX;
	OA12 = Id - Ky*EIY;
	OA21 = Kx*EIX -	Id;
	OA22 = -Kx*EIY;
	
	
	F=[OA11 OA12;OA21 OA22];
	
	
	switch fgm
	
		case 0
	
			M11= Kx*Kx + D*tildeE;
			M12= Ky*(EIX*barE - Kx);
			M21= Kx*(EIY*tildeE - Ky);
			M22= Ky*Ky + B*barE;
			Ma1=[M11 M12;M21 M22];
			
			
			
		case 1
			G11=Kx*Ky;
			G12=barE-(Ky*Ky);
			G21=(Kx*Kx)-tildeE;
			G22=-Kx*Ky;
			
			if zeta~=0
			
				sinz=sind(zeta);
				cosz=cosd(zeta);
				sinz2=sinz^2;
				cosz2=cosz^2;
				OA11=OA11-sinz*Id;
				OA22=OA22+sinz*Id;	
				G11=G11-sinz*Ainv;	
				G22=G22+sinz*Ainv;
				G12=(cosz2*barE + sinz2*Ainv) - Ky*Ky;
				G21=Kx*Kx - (cosz2*tildeE + sinz2*Ainv);
				F=[OA11 OA12;OA21 OA22];
			
		
				
            end			
						
			G=[G11 G12;G21 G22];			
			Ma1=F*G;
			
	end
	
	[W1,q1]=eig(Ma1);
	
	
	W(:,:,ii)=W1;		
	
	q1m(:,ii)=sqrt(diag(q1))*secd(zeta);%heg=0 or 1,clif=0
		
	Q1=diag(q1m(:,ii));
	
	
	V=F\W1*Q1;
	V=V*cosd(zeta);	
	WV(:,:,ii)=V;
	
	
end %ii (layers)

end%rediag_linear_conicalfull


function [Ru,Td,DEr,DEt,summe,varargout]=rezeta(W1,W2,Kx,Ky,q1m,k1zi,kNzi,tns,nxt,reord,theta,psi,phi,zeta,k0)
	
layers=length(tns);

n=length(k1zi);

zerordptr=ceil(n/2);	
sigma=zeros(2*n,1);

Id=eye(n);
Id1=Id*nxt(1)^2;
IdN=Id*nxt(2)^2;

Iw1=[zeros(n) Id;Id zeros(n)];

kz1=diag(k1zi)*cosd(zeta);
kzN=diag(kNzi)*cosd(zeta);

A_mn1=(Id1-Kx*Kx)/kz1;
B_mn1=(Id1-Ky*Ky)/kz1;
C_mn1=(Kx*Ky-Id1*sind(zeta))/kz1;

A_mnN=(IdN-Kx*Kx)/kzN;
B_mnN=(IdN-Ky*Ky)/kzN;
C_mnN=(Kx*Ky-IdN*sind(zeta))/kzN;



W2(:,:,1)=([-B_mn1, -C_mn1;C_mn1, A_mn1 ])/1i; 
	
W2(:,:,layers)=([-B_mnN, -C_mnN;C_mnN, A_mnN ])/1i; 

W1(:,:,1)=Iw1;

W1(:,:,layers)=Iw1;

X1=nan(2*n,2*n,layers);
for LL=2:layers-1
		x1m=exp(-k0*q1m(:,LL)*tns(LL));
		X1(:,:,LL)=diag(x1m);		
end
X1(:,:,layers)=eye(2*n);

Id2=eye(2*n);
	
Ru=zeros(2*n);
	
Td=Id2;




for LL=layers-1:-1:1
		
	
	tilde_Ru=X1(:,:,LL+1)*Ru*X1(:,:,LL+1);
	tilde_Td=Td*X1(:,:,LL+1);		
	
	fL=W1(:,:,LL)\W1(:,:,LL+1);	
	
	gL=W2(:,:,LL)\W2(:,:,LL+1);
	

	aL=fL*(Id2+tilde_Ru);
	bL=gL*(Id2-tilde_Ru);
			
	abL=(aL+bL)\Id2;
		
	Ru=Id2-2*bL*abL;
	Td=2*tilde_Td*abL;
		
		
end%end for LL

for qq=1:length(psi)
	Ux=cosd(psi(qq))*cosd(theta)*cosd(phi)-sind(psi(qq))*sind(phi);
	Uy=cosd(psi(qq))*cosd(theta)*sind(zeta+phi)+sind(psi(qq))*cosd(zeta+phi);

	Ix=Ux*conj(Ux);
	Iy=Uy*conj(Uy);
	Ixy=Ux*conj(Uy)+Uy*conj(Ux);



	sigma0=A_mn1(zerordptr,zerordptr)*Iy+B_mn1(zerordptr,zerordptr)*Ix+C_mn1(zerordptr,zerordptr)*Ixy;


	sigma(zerordptr)=Ux/sqrt(sigma0);
	sigma(zerordptr+n)=Uy/sqrt(sigma0);



	RA=Ru*sigma;
	TA=Td*sigma;

	Rx=RA(1:n);
	Ry=RA(n+1:2*n);
	Tx=TA(1:n);
	Ty=TA(n+1:2*n);


	DErall=real(A_mn1*(Ry.*conj(Ry)) + B_mn1*(Rx.*conj(Rx)) + 2*C_mn1*real(Rx.*conj(Ry)));
	DEtall=real(A_mnN*(Ty.*conj(Ty)) + B_mnN*(Tx.*conj(Tx)) + 2*C_mnN*real(Tx.*conj(Ty)));

	
	sumR=sum(DErall);
    sumT=sum(DEtall);
    sumtot=sumR+sumT;

    summe=[sumR;sumT;sumtot];
	
	DEr=DErall(reord);
	DEt=DEtall(reord);
	
	
    psiq=psi(qq);
	switch zeta
		case 0
			input_norm=(Ux*(-sind(psiq)*sind(phi)+cosd(psiq)/cosd(theta)*cosd(phi)) ...
				+Uy*(sind(psiq)*cosd(phi)+cosd(psiq)/cosd(theta)*sind(phi)))/sqrt(sigma0);
		otherwise
		
			if phi>0
				
				input_norm=(Ux*(-sind(psiq)*sind(phi-zeta)/sind(zeta)+cosd(psiq)/cosd(theta)*cosd(zeta-phi)) ...
				+Uy*(sind(psiq)*sind(phi+zeta)/sind(zeta)+cosd(psiq)/cosd(theta)*sind(phi-zeta)*sind(phi)))/sqrt(sigma0); 
		
						
				
            else%phi<0
                
				if phi>=-30
					input_norm=(Ux*(-sind(psiq)*sind(phi-zeta)/sind(zeta)+cosd(psiq)/cosd(theta)*cosd(zeta+phi)) ...
					+Uy*(sind(psiq)*sind(phi+zeta)/sind(zeta)+cosd(psiq)/cosd(theta)*sind(phi-zeta)*sind(phi)))/sqrt(sigma0);    
				
				else
					input_norm=(Ux*(-sind(psiq)*sind(phi-zeta)/sind(zeta)+cosd(psiq)/cosd(theta)) ...
					+Uy*(sind(psiq)*sind(phi+zeta)/sind(zeta)+cosd(psiq)/cosd(theta)*sind(phi)))/sqrt(sigma0); 
								
				end
            end
	end
		
	
	    
    ky=diag(Ky);
    kx=diag(Kx);
	vphi=atan2(ky-kx*sind(zeta),kx*cosd(zeta));	
    Fc=cos(vphi);Fs=sin(vphi+zeta/180*pi);
    Rsi=(Fc.*Ry-Fs.*Rx)/input_norm;
	Tsi=(Fc.*Ty-Fs.*Tx)/input_norm;
	
	
	
	DErsiall=Rsi.*conj(Rsi).*real(k1zi/(nxt(1)*cosd(theta)));
    DEtsiall=Tsi.*conj(Tsi).*real(kNzi/(nxt(1)*cosd(theta)));
	
	sDEr=DErsiall(reord);
	sDEt=DEtsiall(reord);
	
	if length(psi)==2
		switch psi(qq)
			case 0
				DErP=DEr;
				DEtP=DEt;
				sumP=summe;
				
				DErPS=sDEr;
				DErPP=abs(DEr-sDEr);
				
				DEtPS=sDEt;
				DEtPP=abs(DEt-sDEt);
			case 90
				DErS=DEr;
				DEtS=DEt;
				sumS=summe;
				
				DErSS=sDEr;
				DErSP=abs(DEr-sDEr);
				
				DEtSS=sDEt;
				DEtSP=abs(DEt-sDEt);
		end
		if qq==2
			DEr=(DErP+DErS)/2;
			DEt=(DEtP+DEtS)/2;
			summe=[sumP sumS];
			
			
			varargout{1}=DErPP;
			varargout{2}=DErPS;
			varargout{3}=DErSP;
			varargout{4}=DErSS;
			varargout{5}=DEtPP;
			varargout{6}=DEtPS;
			varargout{7}=DEtSP;
			varargout{8}=DEtSS;
			%varargout{9}=DErP;
			%varargout{10}=DEtP;
			%varargout{11}=DErS;
			%varargout{12}=DEtS;
			
			aef=find(reord==zerordptr);
			if~isempty(aef)
				ckps=DErPS(aef)-DErSP(aef);
					
					if abs(ckps)>1e-7
						
						DErPS(aef)=DErSP(aef);							
						
					end
			end
			
		end
		
		
		
	end


	
end

end%rezeta


function qr=make_all_phi(pr,retnxy,phis,thetas)

phi=-180:180;
iptrxy=nan(size(retnxy,1),1);		
qr=nan(length(thetas),length(phi),length(iptrxy));
cx=retnxy(:,1);cy=retnxy(:,2);


 

phi_sym=phis(1);

phi_sym0=phi_sym;

switch phi_sym0

	case 30
		%[0,0] order
		ord=[0 0];
		iptrxy=0;


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		quat=90;


		ua=pr(:,:,iptrxy);

		ux=phi_sym:phis(end)+quat;
		uz=[ua fliplr(ua(:,1:quat))];

		vx=phi_sym-2*quat:ux(end);
		vz=[fliplr(uz(:,2:length(ux))) uz];

		[aef num]=min(abs(vx-phi(end)));
		nu1=1:num;
		nu2=num:length(vx)-1;
		uz=[vz(:,nu2) vz(:,nu1)];


		%figure,contourf(phi,thetas,uz,16,'linestyle','none');colorbar('location','eastoutside');
		%sr=num2str(retnxy(iptrxy,:));title(sr);

		qr(:,:,iptrxy)=uz;


		%Anti-diagonal orders
		ord=[-2 -2;2 2;-1 -1;1 1];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end
		phi_sym=phis(1);
		phimsym=phi_sym-180;


		for jj=1:2,

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phi_sym-90:phis(end);%-60 to 120
			uz1=[fliplr(ua(:,2:length(phis))) ua];

			%figure,contourf(ux1,thetas,uz1,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);


			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=phi_sym-90:phis(end);%-60 to 120
			uz2=[fliplr(ub(:,2:length(phis))) ub];

			ux1b=ux2-180;%-240 to -60
			[aef num]=min(abs(ux1b-phi(1)));
			nu1=num:length(ux2)-1;

			uz1a=uz2(:,nu1);
			uz1c=uz2(:,2:num);
			uza=[uz1a uz1 uz1c];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);

			uz2a=uz1(:,nu1);
			uz2c=uz1(:,2:num);
			uzb=[uz2a uz2 uz2c];

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;


		end



		%Diagonal orders
		ord=[1 -1;-1 1;2 -2;-2 2;];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end


		phi_sym=phis(end);
		phimsym=phi_sym-180;


		for jj=1:2,

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phis(1):phi_sym+90;%30 to 210
			uz1=[ua fliplr(ua(:,1:length(phis)-1))];

			%figure,contourf(ux1,thetas,uz1,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);


			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=phis(1):phi_sym+90;%30 to 210
			uz2=[ub fliplr(ub(:,1:length(phis)-1))];

			ux1b=ux2-180;%-150 to 30
			[aef num]=min(abs(ux1-phi(end)));
			nu1=num:length(ux1)-1;

			uz1a=uz2(:,nu1);
			uz1c=uz1(:,2:num);
			uza=[uz1a uz2 uz1c];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);

			
			
			uz2a=uz1(:,nu1);
			uz2c=uz2(:,2:num);
			uzb=[uz2a uz1 uz2c];

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;


		end

		%Orders of [0,+/-m] and [+/-m,0]
		ord=[0 -1;-1 0;0 1;1 0;0 -2;-2 0;0 2;2 0];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			udx1=120:210;
			ud2=fliplr(ud);
			[aef nu1]=min(abs(udx1-180));
			ud2b=ud2(:,2:nu1);
			ud2a=ud2(:,nu1:length(udx1)-1);
			
			uxa=[ud2a uc(:,1:length(phis)-1) fliplr(ub(:,2:length(phis))) ua ud2b]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
			qr(:,:,iptrxy(ii))=uxa;
			
			
			uxb1=fliplr(uxa);
			
			[aef nu1]=min(abs(phi-120));
			uxb1a=uxb1(:,nu1:length(phi)-1);
			uxb=[uxb1a uxb1(:,1:nu1)];
			
			%figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			[aef nu1]=min(abs(phi));
			uxc=[uxa(:,nu1:end) uxa(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
			qr(:,:,iptrxy(ii+2))=uxc;
			
			[aef nu1]=min(abs(phi));
			uxd=[uxb(:,nu1:end) uxb(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
			qr(:,:,iptrxy(ii+3))=uxd;
			
			
		   
			

		end



		ord=[-2 -1;-1 -2;2 1;1 2;1 -2;-2 1;-1 2;2 -1];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			udx1=120:210;
			ud2=fliplr(ud);
			[aef nu1]=min(abs(udx1-180));
			ud2b=ud2(:,2:nu1);
			ud2a=ud2(:,nu1:length(udx1)-1);
			
			uxa=[ud2a uc(:,1:length(phis)-1) fliplr(ub(:,2:length(phis))) ua ud2b]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
			qr(:,:,iptrxy(ii))=uxa;
			
			
			uxb1=fliplr(uxa);
			
			[aef nu1]=min(abs(phi-120));
			uxb1a=uxb1(:,nu1:length(phi)-1);
			uxb=[uxb1a uxb1(:,1:nu1)];
			
			%figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			[aef nu1]=min(abs(phi));
			uxc=[uxa(:,nu1:end) uxa(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
			qr(:,:,iptrxy(ii+2))=uxc;
			
			[aef nu1]=min(abs(phi));
			uxd=[uxb(:,nu1:end) uxb(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
			qr(:,:,iptrxy(ii+3))=uxd;
			
			
		   
			

		end
		
	case 0
			
			%[0,0] order
		ord=[0 0];
		iptrxy=0;
        

		for ii=1:length(iptrxy)
            aec=bsxfun(@minus,retnxy,ord(ii,:));
            aed=abs(aec(:,1))+abs(aec(:,2));
			iptrxy(ii)=find(~aed);
		end

		quat=90;


		ua=pr(:,:,iptrxy);

		ux=phi_sym:phis(end)+quat;
		uz=[ua fliplr(ua(:,1:quat))];

		vx=phi_sym-2*quat:ux(end);
		vz=[fliplr(uz(:,2:length(ux))) uz];

		[aef num]=min(abs(vx-phi(end)));
		nu1=1:num;
		nu2=num:length(vx)-1;
		uz=[vz(:,nu2) vz(:,nu1)];


		%figure,contourf(phi,thetas,uz,16,'linestyle','none');colorbar('location','eastoutside');
		%sr=num2str(retnxy(iptrxy,:));title(sr);
        
        %diffe=uz-bb.pr(:,:,iptrxyb(iptrxy));diffemax=max(max(abs(diffe)))

		qr(:,:,iptrxy)=uz;
        
        
        
        %Anti-diagonal orders in skew space
		ord=[0 -4;0 4;0 -2;0 2];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aec=bsxfun(@minus,retnxy,ord(ii,:));
            aed=abs(aec(:,1))+abs(aec(:,2));
			iptrxy(ii)=find(~aed);
		end
		phi_sym=phis(1);
		phimsym=phi_sym-180;


		for jj=1:2

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phi_sym:phis(end)+quat;%0 to 180
			uz1=[ua fliplr(ua(:,1:length(phis)-1))];

			%figure,contourf(ux1,thetas,uz1,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);


			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=phi_sym-180:phis(1);%-180 to 0
			uz2=[ub fliplr(ub(:,1:length(phis)-1))];

			
			nu1=1:length(ux2)-1;

			uz1a=uz2(:,nu1);
			
			uza=[uz1a uz1];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);
            
            %diffe=uza-bb.pr(:,:,iptrxyb(iptrxy(ii-1)));diffemax=max(max(abs(diffe)))
			
			uzb=fliplr(uza);

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
            
            %diffe=uzb-bb.pr(:,:,iptrxyb(iptrxy(ii)));diffemax=max(max(abs(diffe)))

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;


        end
        
        
        %Diagonal orders in skew space
		ord=[2 0;-2 0;4 0;-4 0;];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aec=bsxfun(@minus,retnxy,ord(ii,:));
            aed=abs(aec(:,1))+abs(aec(:,2));
			iptrxy(ii)=find(~aed);
		end


		phi_sym=phis(1);
		


		for jj=1:2

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phi_sym-90:phi_sym+90;%-90 to 90
			uz1=[fliplr(ua(:,2:length(phis))) ua];
            
            uz1b=ua;
			uz1c=fliplr(ua);

			
			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=ux1;
			uz2=[fliplr(ub(:,2:length(phis))) ub];			
            
            uz2b=ub;
			uz2c=fliplr(ub);

			
			uza=[uz2b(:,1:length(phis)-1) uz1 uz2c(:,2:length(phis))];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);
            %diffe=uza-bb.pr(:,:,iptrxyb(iptrxy(ii-1)));diffemax=max(max(abs(diffe)))
			
			
			uzb=[uz1b(:,1:length(phis)-1) uz2 uz1c(:,2:length(phis))];

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
            %diffe=uzb-bb.pr(:,:,iptrxyb(iptrxy(ii)));diffemax=max(max(abs(diffe)))

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;


        end
        
        
        
        %Orders of [0,+/-m] and [+/-m,0] in skew space
		ord=[1 -1;-1 -1;-1 1;1 1;2 -2;-2 -2;-2 2;2 2];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aec=bsxfun(@minus,retnxy,ord(ii,:));
            aed=abs(aec(:,1))+abs(aec(:,2));
			iptrxy(ii)=find(~aed);
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			
			ue1=uc;
            ue2=fliplr(ud);
            ue3=ua;
			ue4=fliplr(ub);
			
			uxa=[ue1 ue2(:,2:length(phis)) ue3(:,2:length(phis)-1) ue4]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
            %diffe=uxa-bb.pr(:,:,iptrxyb(iptrxy(ii)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii))=uxa;
			
			
            ue1=ud;
            ue2=fliplr(uc);
            ue3=ub;
			ue4=fliplr(ua);
			uxb=[ue1 ue2(:,2:length(phis)) ue3(:,2:length(phis)-1) ue4]; 
            
            %figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
            %diffe=uxb-bb.pr(:,:,iptrxyb(iptrxy(ii+1)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			
			uxc=fliplr(uxb);		
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
            %diffe=uxc-bb.pr(:,:,iptrxyb(iptrxy(ii+2)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+2))=uxc;			
            
            
            uxd=fliplr(uxa);			
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
            %diffe=uxd-bb.pr(:,:,iptrxyb(iptrxy(ii+3)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+3))=uxd;		   
			

        end
        
        
        
        ord=[-1 -3;1 -3;1 3;-1 3;3 -1;-3 -1;-3 1;3 1];
        %ord=[-2 -1;-1 -2;2 1;1 2;1 -2;-2 1;-1 2;2 -1]; in skew space
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aec=bsxfun(@minus,retnxy,ord(ii,:));
            aed=abs(aec(:,1))+abs(aec(:,2));
			iptrxy(ii)=find(~aed);
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			
			ue1=uc;
            ue2=fliplr(ud);
            ue3=ua;
			ue4=fliplr(ub);
			
			uxa=[ue1 ue2(:,2:length(phis)) ue3(:,2:length(phis)-1) ue4]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
            %diffe=uxa-bb.pr(:,:,iptrxyb(iptrxy(ii)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii))=uxa;
			
			
            ue1=ud;
            ue2=fliplr(uc);
            ue3=ub;
			ue4=fliplr(ua);
			uxb=[ue1 ue2(:,2:length(phis)) ue3(:,2:length(phis)-1) ue4]; 
            
            %figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
            %diffe=uxb-bb.pr(:,:,iptrxyb(iptrxy(ii+1)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			
			uxc=fliplr(uxb);		
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
            %diffe=uxc-bb.pr(:,:,iptrxyb(iptrxy(ii+2)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+2))=uxc;			
            
            
            uxd=fliplr(uxa);			
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
            %diffe=uxd-bb.pr(:,:,iptrxyb(iptrxy(ii+3)));diffemax=max(max(abs(diffe)))
			qr(:,:,iptrxy(ii+3))=uxd;		   	   
			

		end
		
	case -60
		%[0,0] order
		ord=[0 0];
		iptrxy=0;


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		quat=90;


		ua=pr(:,:,iptrxy);

		ux=phi_sym:phis(end)+quat;%-60 to 120
		uz=[ua fliplr(ua(:,1:quat))];

		vx=phi_sym-2*quat:ux(end);%-240 to 120
		vz=[fliplr(uz(:,2:length(ux))) uz];

		[aef num]=min(abs(vx-phi(1)));
		nu1=num:length(vx)-1;%-180 to 120
		nu2=1:num;
		uz=[vz(:,nu1) vz(:,nu2)];


		%figure,contourf(phi,thetas,uz,16,'linestyle','none');colorbar('location','eastoutside');
		%sr=num2str(retnxy(iptrxy,:));title(sr);

		qr(:,:,iptrxy)=uz;
		
		%Anti-diagonal orders
		ord=[-2 -2;2 2;-1 -1;1 1];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end
		phi_sym=phis(1);
		phimsym=phi_sym-180;


		for jj=1:2,

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phi_sym:phis(end)+90;%-60 to 120
			uz1=[ua fliplr(ua(:,1:length(phis)-1))];

			%figure,contourf(ux1,thetas,uz1,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);


			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=phi_sym:phis(end)+90;%-60 to 120
			uz2=[ub fliplr(ub(:,1:length(phis)-1))];

			ux1b=ux2-180;%-240 to -60
			[aef num]=min(abs(ux1b-phi(1)));
			nu1=num:length(ux2)-1;

			uz1a=uz2(:,nu1);
			uz1c=uz2(:,2:num);
			uza=[uz1a uz1 uz1c];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);

			uz2a=uz1(:,nu1);
			uz2c=uz1(:,2:num);
			uzb=[uz2a uz2 uz2c];

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;


		end

		
	
		%Diagonal orders
		ord=[1 -1;-1 1;2 -2;-2 2;];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end


		phi_sym=phis(end);
		phimsym=phi_sym-180;


		for jj=1:2,

			ii=(jj-1)*2+1;
			ua=pr(:,:,iptrxy(ii));

			ux1=phis(1)-90:phi_sym;%-150 to 30
			uz1=[fliplr(ua(:,2:length(phis))) ua];

			%figure,contourf(ux1,thetas,uz1,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);


			ii=ii+1;
			ub=pr(:,:,iptrxy(ii));
			ux2=phis(1)-90:phi_sym;%-150 to 30
			uz2=[fliplr(ub(:,2:length(phis))) ub];

			
			[aef num]=min(abs(ux1));
			nu1=num:length(ux1)-1;

			uz1a=uz2(:,nu1);
			uz1c=uz2(:,2:num);
			uza=[uz1a uz1 uz1c];

			%figure,contourf(phi,thetas,uza,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii-1),:));title(sr);

			
			
			uz2a=uz1(:,nu1);
			uz2c=uz1(:,2:num);
			uzb=[uz2a uz2 uz2c];

			%figure,contourf(phi,thetas,uzb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);

			qr(:,:,iptrxy(ii-1))=uza;
			qr(:,:,iptrxy(ii))=uzb;

			
		end
		
		
		

		%Orders of [0,+/-m] and [+/-m,0]
		ord=[0 -1;-1 0;0 1;1 0;0 -2;-2 0;0 2;2 0];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			
			[aef nu1]=min(abs(phis));
			uc3a=uc(:,2:nu1);
			uc3b=uc(:,nu1:length(ua)-1);			
					
			uxa=[uc3b fliplr(ud(:,2:length(phis))) ua fliplr(ub(:,1:length(phis)-1)) uc3a]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
			qr(:,:,iptrxy(ii))=uxa;
			
			
			uxb1=fliplr(uxa);
			
			[aef nu1]=min(abs(phi-120));
			uxb1a=uxb1(:,nu1:length(phi)-1);
			uxb=[uxb1a uxb1(:,1:nu1)];
			
			%figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			[aef nu1]=min(abs(phi));
			uxc=[uxa(:,nu1:end) uxa(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
			qr(:,:,iptrxy(ii+2))=uxc;
			
			[aef nu1]=min(abs(phi));
			uxd=[uxb(:,nu1:end) uxb(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
			qr(:,:,iptrxy(ii+3))=uxd;
			
			
		   
			

		end
	

		ord=[-2 -1;-1 -2;2 1;1 2;1 -2;-2 1;-1 2;2 -1];
		iptrxy=nan(size(ord,1),1);


		for ii=1:length(iptrxy)
			aex=find(cx==ord(ii,1));
			aey=find(cy(aex)==ord(ii,2));
			iptrxy(ii)=min(aex)+aey-1;
		end

		phi_sym=phis(1);


		for jj=1:2

			ii=(jj-1)*4+1;
			ua=pr(:,:,iptrxy(ii));
			ub=pr(:,:,iptrxy(ii+1));
			uc=pr(:,:,iptrxy(ii+2));
			ud=pr(:,:,iptrxy(ii+3));

			[aef nu1]=min(abs(phis));
			uc3a=uc(:,2:nu1);
			uc3b=uc(:,nu1:length(ua)-1);			
					
			uxa=[uc3b fliplr(ud(:,2:length(phis))) ua fliplr(ub(:,1:length(phis)-1)) uc3a]; 

			%figure,contourf(phi,thetas,uxa,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii),:));title(sr);
			qr(:,:,iptrxy(ii))=uxa;
			
			
			uxb1=fliplr(uxa);
			
			[aef nu1]=min(abs(phi-120));
			uxb1a=uxb1(:,nu1:length(phi)-1);
			uxb=[uxb1a uxb1(:,1:nu1)];
			
			%figure,contourf(phi,thetas,uxb,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+1),:));title(sr);
			qr(:,:,iptrxy(ii+1))=uxb;
			
			
			[aef nu1]=min(abs(phi));
			uxc=[uxa(:,nu1:end) uxa(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxc,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+2),:));title(sr);
			qr(:,:,iptrxy(ii+2))=uxc;
			
			[aef nu1]=min(abs(phi));
			uxd=[uxb(:,nu1:end) uxb(:,2:nu1)];
			
			%figure,contourf(phi,thetas,uxd,16,'linestyle','none');colorbar('location','eastoutside');
			%sr=num2str(retnxy(iptrxy(ii+3),:));title(sr);
			qr(:,:,iptrxy(ii+3))=uxd;			  
			

		end	
        

end

end%make_all_phi



function [E, A]=retoep(dpts,npts,texture,dphys,MN)

layers=numel(npts)+1;

nxt=texture{1}(1:2);tns=texture{1}(3:3+layers-1);do_regeom=texture{1}(end);texture{1}=[];
xpts=dpts.x;ypts=dpts.y;



M=MN(1);
N=MN(2);
m=round((M-1)/2);
n=round((N-1)/2);

E=nan(M*N,M*N,layers-1);A=E;


gxpts=cell(1,layers-1);gypts=gxpts;
gnx=gxpts;gny=gxpts;


am=-m:m;
an=-n:n;
[xam yan]=meshgrid(am,an);
cx=xam(:);
cy=yan(:);


bmtx=bsxfun(@minus,cx,cx');%=amtx
bnty=bsxfun(@minus,cy,cy');%=anty




cors={'cyan','magenta','yellow','green','aqua','red','black'};

ss=0;
for ii=1:numel(npts)
	ss=ss+length(npts{ii});
end
nq=nan(1,ss);
jj=1;
for ii=1:numel(npts)
	aef=npts{ii};
	nq(jj:jj+length(aef)-1)=aef;
	jj=jj+length(aef);
end

nq=munique([nq,nxt]);


for ii=1:layers-1
    txte=texture{ii};
    if ~isempty(txte)

        numinc=numel(txte);
        nxs=nan(1,numinc);
        
        for jj=1:numinc 
            txt=txte{jj};
            nxs(jj)=txt(5);
        end

        nq=munique([nq,nxs]);
        
    end
end

nuxx=501;nuyy=501;


plt=0;


switch plt

	case 1 

		for qq=1:layers-1
			nla=npts{qq};ela=nla.^2;
			xpts0=xpts{qq};
			ypts0=ypts{qq};
			x=munique([linspace(-dphys(1)/2,dphys(1)/2, nuxx),xpts0]);	
			y=munique([linspace(-dphys(2)/2,dphys(2)/2, nuyy),ypts0.']);
			epf=nan(length(y),length(x));
			
			if qq==1,epfl=nan(length(y),length(x),layers);end

			
			figure;
		   
			switch length(nla)
				case 1
					xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
					yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
					patch(xx,yy,char(cors(nq==nla)));%Set as background
					epf=ela*ones(size(epf));
				otherwise	
					switch length(xpts0)
						case 1

							dla=ypts0;ila=dla;
							for jj=1:length(dla)
								[bef ila(jj)]=min(abs(y-dla(jj)));
							end
							for jj=1:length(dla)
								if jj==1
									xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
									yy = [-dphys(2)/2 -dphys(2)/2 dla(jj) dla(jj)];
									patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(1:ila(jj),:)=ela(jj);
								else
									xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
									yy = [dla(jj-1) dla(jj-1) dla(jj) dla(jj)];
									patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(ila(jj-1):ila(jj),:)=ela(jj);
								end

							end

						otherwise
							dla=xpts0;ila=dla;
							for jj=1:length(dla)
								[bef ila(jj)]=min(abs(x-dla(jj)));
							end
							for jj=1:length(dla)

								if jj==1
									xx = [-dphys(1)/2 dla(jj) dla(jj) -dphys(1)/2];
									yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
									patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(:,1:ila(jj))=ela(jj);
								else
									xx = [dla(jj-1) dla(jj) dla(jj) dla(jj-1)];
									yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
									patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(:,ila(jj-1):ila(jj))=ela(jj);
								end

							end
					end
			end	
			
			txte=texture{qq};
			
				
				
				
			if~isempty(txte)
				numinc=numel(txte);
				for jj=1:numinc
					
					aef=txte{jj};
					x1=aef(1)-aef(3)/2;
					x2=aef(1)+aef(3)/2;
					y1=aef(2)-aef(4)/2;
					y2=aef(2)+aef(4)/2;
					ndx=aef(5);

					xx = [x1 x2 x2 x1];
					yy = [y1 y1 y2 y2];
					patch(xx,yy,char(cors(nq==ndx)));


					x0=[x1,x2];y0=[y1,y2];



					[aef ry1]=min(abs(y-y1));
					[aef ry2]=min(abs(y-y2));

					[aef cx1]=min(abs(x-x1));
					[aef cx2]=min(abs(x-x2));

					switch length(nla)
						case 1

							epf(ry1:ry2,cx1:cx2)=ndx.^2;

						otherwise

							switch length(xpts{ii})
								case 1
									cef=find(ila==ry2);
									if ~isempty(cef)
										epf(ry1:ry2-1,cx1:cx2)=ndx.^2;
									else
										epf(ry1:ry2,cx1:cx2)=ndx.^2;
									end	

								otherwise


									cef=find(ila==cx2);
									if ~isempty(cef)
										epf(ry1:ry2,cx1:cx2-1)=ndx.^2;
									else
										epf(ry1:ry2,cx1:cx2)=ndx.^2;
									end

							end
					end
				end



			end
			
			v=axis;
			v(1)=-dphys(1)/2;v(2)=dphys(1)/2;v(3)=-dphys(2)/2;v(4)=dphys(2)/2;
			axis(v);axis equal;
			xlabel('x (\mum)');ylabel('y (\mum)');
			sr=['Layer: #',num2str(qq)];title(sr);
			
			
			if isempty(txte) && length(nla)==1
				E(:,:,qq)=eye(M*N)*ela;
				A(:,:,qq)=eye(M*N)/ela;
			else
					
			
				epFC=ifftshift(ifft2(ifftshift(epf)));
				apFC=ifftshift(ifft2(ifftshift(1./epf)));

				offs = ceil(size(epf)/2) + 1;


				co=bmtx+offs(2);ro=bnty+offs(1);

				fr=ro;fc=co;

				fi=fr;
				for ii=1:numel(fr)
					fi(ii)=(fc(ii)-1)*size(epf,1)+fr(ii);
				end

				E1=ro;A1=ro;
				for ii=1:numel(ro)
					E1(ii)=epFC(ro(ii),co(ii));
					A1(ii)=apFC(ro(ii),co(ii));
				end
				E(:,:,qq)=E1;
				A(:,:,qq)=A1;
			end
		end%qq


	otherwise%plt~=1
		
		for qq=1:layers-1
			nla=npts{qq};ela=nla.^2;
			xpts0=xpts{qq};
			ypts0=ypts{qq};
			x=munique([linspace(-dphys(1)/2,dphys(1)/2, nuxx),xpts0]);	
			y=munique([linspace(-dphys(2)/2,dphys(2)/2, nuyy),ypts0.']);
			epf=nan(length(y),length(x));

			if qq==1,epfl=nan(length(y),length(x),layers);end
			
		   
			switch length(nla)
				case 1
					%xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
					%yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
					%patch(xx,yy,char(cors(nq==nla)));%Set as background
					epf=ela*ones(size(epf));
				otherwise	
					switch length(xpts0)
						case 1

							dla=ypts0;ila=dla;
							for jj=1:length(dla)
								[bef ila(jj)]=min(abs(y-dla(jj)));
							end
							for jj=1:length(dla)
								if jj==1
									%xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
									%yy = [-dphys(2)/2 -dphys(2)/2 dla(jj) dla(jj)];
									%patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(1:ila(jj),:)=ela(jj);
								else
									%xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
									%yy = [dla(jj-1) dla(jj-1) dla(jj) dla(jj)];
									%patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(ila(jj-1):ila(jj),:)=ela(jj);
								end

							end

						otherwise
							dla=xpts0;ila=dla;
							for jj=1:length(dla)
								[bef ila(jj)]=min(abs(x-dla(jj)));
							end
							for jj=1:length(dla)

								if jj==1
									%xx = [-dphys(1)/2 dla(jj) dla(jj) -dphys(1)/2];
									%yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
									%patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(:,1:ila(jj))=ela(jj);
								else
									%xx = [dla(jj-1) dla(jj) dla(jj) dla(jj-1)];
									%yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
									%patch(xx,yy,char(cors(nq==nla(jj))));%Set as background
									epf(:,ila(jj-1):ila(jj))=ela(jj);
								end

							end
					end
			end	
			
			txte=texture{qq};
			
				
				
				
			if~isempty(txte)
				numinc=numel(txte);
				for jj=1:numinc
					
					aef=txte{jj};
					x1=aef(1)-aef(3)/2;
					x2=aef(1)+aef(3)/2;
					y1=aef(2)-aef(4)/2;
					y2=aef(2)+aef(4)/2;
					ndx=aef(5);

					%xx = [x1 x2 x2 x1];
					%yy = [y1 y1 y2 y2];
					%patch(xx,yy,char(cors(nq==ndx)));


					x0=[x1,x2];y0=[y1,y2];



					[aef ry1]=min(abs(y-y1));
					[aef ry2]=min(abs(y-y2));

					[aef cx1]=min(abs(x-x1));
					[aef cx2]=min(abs(x-x2));

					switch length(nla)
						case 1

							epf(ry1:ry2,cx1:cx2)=ndx.^2;

						otherwise

							switch length(xpts{ii})
								case 1
									cef=find(ila==ry2);
									if ~isempty(cef)
										epf(ry1:ry2-1,cx1:cx2)=ndx.^2;
									else
										epf(ry1:ry2,cx1:cx2)=ndx.^2;
									end	

								otherwise


									cef=find(ila==cx2);
									if ~isempty(cef)
										epf(ry1:ry2,cx1:cx2-1)=ndx.^2;
									else
										epf(ry1:ry2,cx1:cx2)=ndx.^2;
									end

							end
					end
				end



			end
			
					
			
			if isempty(txte) && length(nla)==1
				E(:,:,qq)=eye(M*N)*ela;
				A(:,:,qq)=eye(M*N)/ela;
			else
					
			
				epFC=ifftshift(ifft2(ifftshift(epf)));
				apFC=ifftshift(ifft2(ifftshift(1./epf)));

				offs = ceil(size(epf)/2) + 1;


				
				co=bmtx+offs(2);ro=bnty+offs(1);

				fr=ro;fc=co;

				fi=fr;
				for ii=1:numel(fr)
					fi(ii)=(fc(ii)-1)*size(epf,1)+fr(ii);
				end

				E1=ro;A1=ro;
				for ii=1:numel(ro)
					E1(ii)=epFC(ro(ii),co(ii));
					A1(ii)=apFC(ro(ii),co(ii));
				end
				E(:,:,qq)=E1;
				A(:,:,qq)=A1;
			end
			
			
			%epfl(:,:,qq)=epf;
			
			if do_regeom==1
				gex=-dphys(1)/4;%cross-section along y at x=gex, gex=0 by default.
				gey=-dphys(2)/4;%cross-section along x at y=gey, gey=0 by default.
				
				%gex=0;%cross-section along y at x=gex, gex=0 by default.
				%gey=0;%cross-section along x at y=gey, gey=0 by default.
				
				cross_section_along_x_y=[gex gey]
				
				[hef huy]=min(abs(y-gey));
				npix=sqrt(epf(huy,:));
				dex=diff(npix);
				dex1=(~dex);
				dex2=find(dex1==0);
				
				if isempty(dex2)
					gxpts{qq}=dphys(1)/2;
					gxn{qq}=npix(1);
				else
					ggxx=nan(1,length(dex2));
					ggxn=ggxx;
					for qlo=1:length(ggxx)
						pgwx=dex2(qlo);
						ggxx(qlo)=x(pgwx+1);
						ggxn(qlo)=npix(pgwx);
					end
					gxpts{qq}=ggxx;
					gxn{qq}=ggxn;
				end
				
				[hef hux]=min(abs(x-gex));
				npiy=sqrt(epf(:,hux));
				dey=diff(npiy);
				dey1=(~dey);
				dey2=find(dey1==0);
				
				if isempty(dey2)
					gypts{qq}=dphys(2)/2;
					gyn{qq}=npiy(1);
				else
					ggyy=nan(1,length(dey2));
					ggyn=ggyy;
					for qlo=1:length(ggyy)
						pgwy=dey2(qlo);
						ggyy(qlo)=y(pgwy+1);
						ggyn(qlo)=npiy(pgwy);
					end
					gypts{qq}=ggyy;
					gyn{qq}=ggyn;	
				end
				
			end%do_regeom
			
			
			
		end%qq
		
	if do_regeom==1	
		regeom(gxpts, gxn,tns,nxt,dphys(1));
		regeoy(gypts, gyn,tns,nxt,dphys(2));
	end
	
end%switch

end%retoep

function [dpts, npts, tns, texture, layers]=diamond_ortho_inputfile(nxt,H,iptfile,T,dphys,ld)

n_glass=nxt(1);
n_air=nxt(2);

if T==0
	layii=1;
	layers=2+layii;
	layer_thickness=H*1e-9;

else
	layii=3;	
	layer_thickness=[T H-T T]*1e-9;
	
	xTiO2=[400 460 530 620 700]*1e-9;
	nTiO2=[2.5 2.5 2.42 2.36 2.36];
	n_TaO=interp1(xTiO2,nTiO2,ld,'linear');

end



layers=2+layii;
tns=zeros(1,layers);
tns(2:layers-1)=layer_thickness;

xpts=cell(1,layers-1);
npts=xpts;texture=xpts;


Tx=dphys(1);
Ty=dphys(2);


	
obx=Tx;
oby=Ty;

alp=0.64;

alpha=ones(1,layers)*alp;


show_geom_xy=0;

xpts{1}=Tx/2;ypts{1}=Ty/2;npts{1}=nxt(1);texture{1}=[nxt tns show_geom_xy];

optn=9;



if T==0,
		
			ii=2;
			xpts{ii}=Tx/2;ypts{ii}=Ty/2;npts{ii}=n_air;
			texture{ii}={[0,alp,obx,oby,n_glass,optn],[iptfile]};
		
		
		


else
		
			ss=2;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_TaO;		
			texture{ss}={[0,alp,obx,oby,n_glass,optn],[iptfile]};
			
					
			ss=3;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[0,alp,obx,oby,n_glass,optn],[iptfile]};
			
			
			ss=4;
			xpts{ss}=Tx/2;ypts{ss}=Ty/2;npts{ss}=n_air;
			texture{ss}={[0,alp,obx,oby,n_TaO,optn],[iptfile]};
			
				
		
end

dpts.x=xpts;dpts.y=ypts;dpts.alpha=alpha;

end%diamond_ortho_inputfile

function [E, A]=retoep2(dpts,npts,texture,dphys,MN)

layers=numel(npts)+1;

nxt=texture{1}(1:2);tns=texture{1}(3:3+layers-1);do_regeom=texture{1}(end);texture{1}=[];
xpts=dpts.x;ypts=dpts.y;



M=MN(1);
N=MN(2);
m=round((M-1)/2);
n=round((N-1)/2);

E=nan(M*N,M*N,layers-1);A=E;


gxpts=cell(1,layers-1);gypts=gxpts;
gnx=gxpts;gny=gxpts;


am=-m:m;
an=-n:n;
[xam yan]=meshgrid(am,an);
cx=xam(:);
cy=yan(:);


bmtx=bsxfun(@minus,cx,cx');%=amtx
bnty=bsxfun(@minus,cy,cy');%=anty




cors={'cyan','magenta','yellow','green','aqua','red','black'};

ss=0;
for ii=1:numel(npts)
	ss=ss+length(npts{ii});
end
nq=nan(1,ss);
jj=1;
for ii=1:numel(npts)
	aef=npts{ii};
	nq(jj:jj+length(aef)-1)=aef;
	jj=jj+length(aef);
end

nq=munique([nq,nxt]);


for ii=1:layers-1
    txte=texture{ii};
    if ~isempty(txte)

        numinc=numel(txte);
        nxs=nan(1,numinc);
        
        for jj=1:numinc 
            txt=txte{jj};
            nxs(jj)=txt(5);
        end

        nq=munique([nq,nxs]);
        
    end
end

txte=texture{2}{1};

optn=txte(6);


%If optn==2, let us assume all the layers have the same MS/notches.
switch optn
	case 2

		MS=txte(1); Notch=txte(2); diaw =txte(3);len= txte(4);
				
		yoff=-len/2;
		
		if Notch==0

			xz=[0 -diaw/2 0];
			yz=[0 len/2 len]+yoff;
		else
			
			D_noskew=dphys(2)/2;
			nDp=Notch*D_noskew/cosd(30);
			nDperp=nDp*cosd(30);
			notchw=dphys(2)/cosd(30)*(MS*sind(30)-Notch);

			xz=[0 -diaw/2+nDp/2 -notchw/2 -diaw/2+nDp/2 0];
			yz=[0 len/2-nDperp len/2 len/2+nDperp len]+yoff;
		end
		
		xz0=munique(xz);
		yz0=munique(yz);

		djx=logical(xz0);
		mxz0=abs(xz0(djx));

		xpt0=munique(round([-mxz0 mxz0]*1e12));
		xpt0=xpt0*1e-12;

		djy=logical(yz0);
		myz0=abs(yz0(djy));

		ypt0=munique(round([-myz0 myz0]*1e12));
		ypt0=ypt0*1e-12;
        
        
        nuxy=ceil(dphys*1e9);

        ajn=rem(nuxy,2);
        aji=find(~ajn);
        if ~isempty(aji)
            nuxy(aji)=nuxy(aji)+1;
        end
        nuxx=nuxy(1);
        nuyy=nuxy(2);

        x=munique([linspace(-dphys(1)/2,dphys(1)/2, nuxx),xpt0]);	
        y=munique([linspace(-dphys(2)/2,dphys(2)/2, nuyy),ypt0]);

        %epfl=nan(length(y),length(x),layers);
        
        numi=yz;
		for ii=1:length(numi)
			[aef numi(ii)]=min(abs(y-yz(ii)));
            
        end
      
		
		num1=numi(1);
		num2=numi(end);
		
		yqt=[y(num1:num2)].';xqt=yqt;
		count=1;	
        xytol=1e-13;
        xqt(count)=xz(1);count=count+1;
		for ii=1:length(numi)-1
			nips=numi(ii):numi(ii+1);
			xdar=[xz(ii)-xytol xz(ii+1)+xytol];
			ydar=[yz(ii)-xytol yz(ii+1)+xytol];
			for jj=2:length(nips)				
				xqt(count)=interp1(ydar,xdar,yqt(count),'linear');
				count=count+1;
			end		
		
		end
		
		
		[aef quy0]=min(abs(yqt));
		
		xdh=xqt(1:quy0);ych=yqt(1:quy0);
		xdb=xqt(quy0:end);ycb=yqt(quy0:end);
		
		xdh=xdh+dphys(1)/2;ych=ych+dphys(2)/2;
		xdb=xdb+dphys(1)/2;ycb=ycb-dphys(2)/2;
		
		xgh=-xdh;
		xgb=-xdb;
		
		%figure,plot(xqt,yqt,'k-');
		%hold on;plot(xdh,ych,'r-',xdb,ycb,'b-',xgh,ych,'r-',xgb,ycb,'b-');axis equal
		
		
        mxdh=dphys(1)/2*ones(length(ych),1);
        mxdb=dphys(1)/2*ones(length(ycb),1);
        mxgh=-dphys(1)/2*ones(length(ych),1);
        mxgb=-dphys(1)/2*ones(length(ycb),1);
        
		xytoutes={[xqt -xqt yqt],[xdh mxdh ych],[xdb mxdb ycb],[mxgh xgh ych],[mxgb xgb ycb]};
        
        
        plt=1;


        switch plt

            case 1 

                for qq=1:layers-1
                    nla=npts{qq};ela=nla.^2;

                    figure;
                    xx = [-dphys(1)/2 dphys(1)/2 dphys(1)/2 -dphys(1)/2];
                    yy = [-dphys(2)/2 -dphys(2)/2 dphys(2)/2 dphys(2)/2];
                    patch(xx,yy,char(cors(nq==nla)));%Set as background

                    epf=ela*ones(length(y),length(x));


                    txte=texture{qq};

       		

                    if~isempty(txte)
                        %xytoutes={[xqt, yqt],[xdh ych],[xdb ycb],[xgh ych],[xgb ycb]};
                        numpc=numel(xytoutes);
                        for jj=1:numpc
                            ndx=txte{1}(5);
                            xysep=xytoutes{jj};


                            xx=[xysep(:,1);flipud(xysep(:,2))];
                            yy=[xysep(:,3);flipud(xysep(:,3))];

                            patch(xx,yy,char(cors(nq==ndx)));

                            [aef num1]=min(abs(y-min(xysep(:,3))));
                            [aef num2]=min(abs(y-max(xysep(:,3))));

                            numys=num1:num2;

                            for plo=1:length(numys)
                                [aef numx1]=min(abs(x-xysep(plo,1)));
                                [aef numx2]=min(abs(x-xysep(plo,2)));
                                epf(numys(plo),numx1:numx2)=ndx.^2;  
                            end

                        end



                    end

                    v=axis;
                    v(1)=-dphys(1)/2;v(2)=dphys(1)/2;v(3)=-dphys(2)/2;v(4)=dphys(2)/2;
                    axis(v);axis equal;
                    xlabel('x (\mum)');ylabel('y (\mum)');
                    sr=['Layer: #',num2str(qq)];title(sr);

                

                    if isempty(txte) && length(nla)==1
                        E(:,:,qq)=eye(M*N)*ela;
                        A(:,:,qq)=eye(M*N)/ela;
                    else


                        epFC=ifftshift(ifft2(ifftshift(epf)));
                        apFC=ifftshift(ifft2(ifftshift(1./epf)));

                        offs = ceil(size(epf)/2) + 1;

                        co=bmtx+offs(2);ro=bnty+offs(1);

                        fr=ro;fc=co;

                        fi=fr;
                        for ii=1:numel(fr)
                            fi(ii)=(fc(ii)-1)*size(epf,1)+fr(ii);
                        end

                        E1=ro;A1=ro;
                        for ii=1:numel(ro)
                            E1(ii)=epFC(ro(ii),co(ii));
                            A1(ii)=apFC(ro(ii),co(ii));
                        end
                        E(:,:,qq)=E1;
                        A(:,:,qq)=A1;
                    end
                end%qq



        end%switch plt

		
		
	case 9
		iptfile=texture{2}{2};
		epones=csvread(iptfile);
		
		nuxy=ceil(dphys*1e9);

		ajn=rem(nuxy,2);
		aji=find(~ajn);
		if ~isempty(aji)
			nuxy(aji)=nuxy(aji)+1;
		end
		nuxx=nuxy(1);
		nuyy=nuxy(2);
		
		
		%{
		if nuxx~=size(epones,2)
			disp('No. of points along x is not consistent with period_x!');
		end
		
		if nuyy~=size(epones,1)
			disp('No. of points along y is not consistent with period_y!');
		end
		%}
		
		x=linspace(-dphys(1)/2,dphys(1)/2, nuxx);	
        y=linspace(-dphys(2)/2,dphys(2)/2, nuyy);

		
		plt=1;


        switch plt

            case 1 
				
				

                for qq=1:layers-1
                    nla=npts{qq};ela=nla.^2;
                    txte=texture{qq};       		

                    if~isempty(txte)
                                                
					   ndx=txte{1}(5);
					   epf=epones*ndx^2;
                       epf(~epf)=ela;
					else
						epf=ela*ones(length(y),length(x));
                    end
                    %{
                    figure,contourf(x*1e6,y*1e6,epf,32,'linestyle','none');colorbar('location','eastoutside')
					v=axis;
					v(1)=-dphys(1)/2;v(2)=dphys(1)/2;v(3)=-dphys(2)/2;v(4)=dphys(2)/2;v=v*1e6;
					axis(v);axis equal;
					xlabel('x (\mum)');ylabel('y (\mum)');
					sr=['Layer: #',num2str(qq)];title(sr);
                    %}
                

                    if isempty(txte) && length(nla)==1
                        E(:,:,qq)=eye(M*N)*ela;
                        A(:,:,qq)=eye(M*N)/ela;
                    else


                        epFC=ifftshift(ifft2(ifftshift(epf)));
                        apFC=ifftshift(ifft2(ifftshift(1./epf)));

                        offs = ceil(size(epf)/2) + 1;

                        co=bmtx+offs(2);ro=bnty+offs(1);

                        fr=ro;fc=co;

                        fi=fr;
                        for ii=1:numel(fr)
                            fi(ii)=(fc(ii)-1)*size(epf,1)+fr(ii);
                        end

                        E1=ro;A1=ro;
                        for ii=1:numel(ro)
                            E1(ii)=epFC(ro(ii),co(ii));
                            A1(ii)=apFC(ro(ii),co(ii));
                        end
                        E(:,:,qq)=E1;
                        A(:,:,qq)=A1;
                    end
                end%qq



        end%switch plt
		
		
	
end


end%retoep2







function [E, A]=epsfep(x0,nf,T,n)

ivft=-1;
if length(nf)==1
	E=eye(n)*nf^2;
	A=eye(n)/nf^2;

	
else

		epsilon=nf.^2;

		numepsterms=2*n-1;
		kkf=[-(n-1):n-1].';		
        		
		
		ck=kkf;
		ak=kkf;
	
		
		switch length(x0)
			case 2
				epsilon1=epsilon(1);
				epsilon2=epsilon(2);
				f=(x0(2)-x0(1))/T;
				
				ck=(epsilon2-epsilon1)*msinc(kkf*f)*f;
				ak=(1/epsilon2-1/epsilon1)*msinc(kkf*f)*f;
				ck(n)=ck(n)+epsilon1;
				ak(n)=ak(n)+1/epsilon1;
				
				shp=(sum(x0)/2/T)*2*pi;
				ck=ck.*exp(1i*shp*kkf);
				ak=ak.*exp(1i*shp*kkf);
				
			case 4
			
			
				
				ckq=nan(2,length(kkf)).';akq=ckq;
				
				
				for jj=1:2,
					epsilon1=epsilon(1);
					epsilon2=epsilon(2);
					
					switch jj
						case 1
							f=(x0(4)-x0(1))/T;
							shp=((x0(1)+x0(4))/2/T)*2*pi;
						case 2
							f=(x0(3)-x0(2))/T;
							shp=((x0(2)+x0(3))/2/T)*2*pi;
					end
					
					
					
					ckq(:,jj)=(epsilon2-epsilon1)*msinc(kkf*f)*f;
					akq(:,jj)=(1/epsilon2-1/epsilon1)*msinc(kkf*f)*f;
					ckq(n,jj)=ckq(n,jj)+epsilon1;
					akq(n,jj)=akq(n,jj)+1/epsilon1;
					ckq(:,jj)=ckq(:,jj).*exp(1i*shp*kkf);
					akq(:,jj)=akq(:,jj).*exp(1i*shp*kkf);
					
				end
				
				
				ckp=ckq(:,1)-ckq(:,2);
				akp=akq(:,1)-akq(:,2);
				
				epsilon2=epsilon(3);
				ck=(epsilon2-epsilon1)*msinc(kkf*f)*f;
				ak=(1/epsilon2-1/epsilon1)*msinc(kkf*f)*f;
				ck(n)=ck(n)+epsilon1;
				ak(n)=ak(n)+1/epsilon1;
				
				ck=ck.*exp(1i*shp*kkf);
				ak=ak.*exp(1i*shp*kkf);
				
				ak=ak+akp;
				ck=ck+ckp;
				
				
				
			otherwise
				ivft=1;
				
				nuxx=501;
				if T/(nuxx-1)>2e-9;nuxx=1001;disp('Check x resolution in epsfep');end
				xx=linspace(-T/2,T/2, nuxx);
				epf=ones(size(xx))*epsilon(1);

				for ii=1:length(x0)-1,

					[aef fu1]=min(abs(xx-x0(ii)));
					[aef fu2]=min(abs(xx-x0(ii+1)));
					epf(fu1:fu2)=epsilon(ii+1)*ones(1,fu2-fu1+1);
				end


				apf=1./epf;
				
				
				
				aa=ifftshift(ifft(ifftshift(epf)));

				mixx=ceil(nuxx/2)+1;
				ixx=fliplr(mixx-n+1:mixx+n-1);
				ck=aa(ixx).';

				aa=ifftshift(ifft(ifftshift(1./epf)));

				ak=aa(ixx).';
					
				
				
                
		
		end
		

		A=toeplitz(ak(n:end),flipud(ak(1:n)));E=toeplitz(ck(n:end),flipud(ck(1:n)));

		if ivft==-1,A=A.';E=E.';end



end
end%epsfep

function mroot=psqrt(msq)


mroot=sqrt(msq);
aef=find(imag(mroot)<0);
mroot(aef)=-mroot(aef);

end%psqrt

function q=munique(x)
y=sort(x);
dy=diff(y);
checkz=find(dy~=0);
if ~isempty(checkz),
	q=nan(1,length(checkz)+1);
	if size(y,1)>1
		q=nan(length(checkz)+1,1);
	end
	q(1)=y(1);
	
	for ii=1:length(checkz),
		q(ii+1)=y(checkz(ii)+1);
	end
	
end
end%munique

function y=msinc(x)
y=x;
kk=find(x);
A=pi*x(kk);
y(kk)=sin(A)./A;
kk=find(~x);
y(kk)=1;
end%msinc


function [pcs,ATE,ATM]=repcs(npts,DEsumTE,DEsumTM,pol)

ATE=0;
ATM=0;
pcs=1;

ss=0;
for ii=1:numel(npts),
	ss=ss+length(npts{ii});
end
nq=nan(1,ss);
jj=1;
for ii=1:numel(npts),
	aef=npts{ii};
	nq(jj:jj+length(aef)-1)=aef;
	jj=jj+length(aef);
end

nq=munique(nq);



if sum(imag(nq))==0,

	

	switch length(pol)
		case 1
			switch pol
				case 1
					pcs=mean(DEsumTE);
				case 2
					pcs=mean(DesumTM);
				case 3
					pcs=mean(DEsumTE);
				
			end
			
		case 2  
			pcs=(mean(DEsumTE)+mean(DEsumTM))/2;
	end

	
	
	if abs(pcs-1)>0.0001,
		disp('Energy conservation NOT satisfied!');
	end

else
	
	switch length(pol)
	
	
		case 1
			switch pol
				case 1
					ATE=1-DEsumTE;
				case 2
					ATM=1-DEsumTM;
				case 3
					ATE=1-DEsumTE;;
				
			end
			
		case 2  
			ATE=1-DEsumTE;
			ATM=1-DEsumTM;
	end
	
	disp('Absorption calculated from A=1-R-T');

end

end%repcs

function []=regeom(xpts, npts,d,n,T)

T=T*1e6;
d=d*1e6;

da=sum(d);

cors={'cyan','magenta','yellow','green','aqua','red','black'};
ss=0;
for ii=1:numel(npts),
	ss=ss+length(npts{ii});
end
nq=nan(1,ss);
jj=1;
for ii=1:numel(npts),
	aef=npts{ii};
	nq(jj:jj+length(aef)-1)=aef;
	jj=jj+length(aef);
end

nq=munique([n,nq]);


figure;

dmax=1.2*da;
dmin=-0.2*da;


dc=cumsum(d);
xx = [-T/2 T/2 T/2 -T/2];
yy = [dmin dmin dmax dmax];
patch(xx,yy,char(cors(nq==n(1))));%Set superstrate as background


xx = [-T/2 T/2 T/2 -T/2];
yy = [da da dmax dmax];
patch(xx,yy,char(cors(nq==n(2))));%substrate



hold on,line([-T/2 T/2],[0 0]);
hold on,line([-T/2 T/2],[1 1]*da);




for ii=2:length(d)-1,
	xp=xpts{ii}*1e6;
	np=npts{ii};
	
	switch length(xp)
		case 2
			xx = [xp(1) xp(2) xp(2) xp(1)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(2))));
			
			xx = [-T/2 xp(1) xp(1) -T/2];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			
			patch(xx,yy,char(cors(nq==np(1))));
			
			xx = [xp(2) T/2 T/2 xp(2)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];			
			patch(xx,yy,char(cors(nq==np(1))));
			
		otherwise
		
			for jj=1:length(xp)-1
				xx = [xp(jj) xp(jj+1) xp(jj+1) xp(jj)];
				yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
				patch(xx,yy,char(cors(nq==np(jj+1))));
			end	
				
			xx = [-T/2 xp(1) xp(1) -T/2];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(1))));
			
			xx = [xp(end) T/2 T/2 xp(end)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(1))));
			
	end	
	
end

axis equal
xlabel('x (\mum)');
ylabel('z (\mum)')
s=sprintf('Period T_x = %.1f nm',T*1000);
v=axis;
hold on; 
line([-T/2 -T/2],[v(3) v(4)],'linestyle','--','color',[0 0 0]);
line([T/2 T/2],[v(3) v(4)],'linestyle','--','color',[0 0 0]);
title(s);

end%regeom


function []=regeoy(xpts, npts,d,n,T)

T=T*1e6;
d=d*1e6;



da=sum(d);

cors={'cyan','magenta','yellow','green','aqua','red','black'};
ss=0;
for ii=1:numel(npts),
	ss=ss+length(npts{ii});
end
nq=nan(1,ss);
jj=1;
for ii=1:numel(npts),
	aef=npts{ii};
	nq(jj:jj+length(aef)-1)=aef;
	jj=jj+length(aef);
end

nq=munique([n,nq]);


figure;

dmax=1.2*da;
dmin=-0.2*da;



%{

dc=-cumsum(d)+da;
xx = [-T/2 T/2 T/2 -T/2];
yy = [dmin dmin dmax dmax];
patch(xx,yy,char(cors(nq==n(1))));%Set superstrate as background

xx = [-T/2 T/2 T/2 -T/2];
yy = -[0.2 0.2 0 0]*da;
patch(xx,yy,char(cors(nq==n(2))));%substrate
%}



dc=cumsum(d);
xx = [-T/2 T/2 T/2 -T/2];
yy = [dmin dmin dmax dmax];
patch(xx,yy,char(cors(nq==n(1))));%Set superstrate as background


xx = [-T/2 T/2 T/2 -T/2];
yy = [da da dmax dmax];
patch(xx,yy,char(cors(nq==n(2))));%substrate



hold on,line([-T/2 T/2],[0 0]);
hold on,line([-T/2 T/2],[1 1]*da);




for ii=2:length(d)-1,
	xp=xpts{ii}*1e6;
	np=npts{ii};
	
	switch length(xp)
	
		case 1
			doit=1;%Quick fix with doit. Need to differentiate between homogenous layer or strucutre placed without symmetry. 
			
			if doit==1
				
				
				xx = [-T/2 xp(1) xp(1) -T/2];
				yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
				
				patch(xx,yy,char(cors(nq==np(1))));
				
				xx = [xp(1) T/2 T/2 xp(1)];
				yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
				%if max(xx)>T/2,stop;end
				patch(xx,yy,char(cors(1)));
			end
		case 2
			xx = [xp(1) xp(2) xp(2) xp(1)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(2))));
			
			xx = [-T/2 xp(1) xp(1) -T/2];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			
			patch(xx,yy,char(cors(nq==np(1))));
			
			xx = [xp(2) T/2 T/2 xp(2)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			%if max(xx)>T/2,stop;end
			patch(xx,yy,char(cors(nq==np(1))));
			
		otherwise
		
			for jj=1:length(xp)-1
				xx = [xp(jj) xp(jj+1) xp(jj+1) xp(jj)];
				yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
				patch(xx,yy,char(cors(nq==np(jj+1))));
			end	
				
			xx = [-T/2 xp(1) xp(1) -T/2];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(1))));
			
			xx = [xp(end) T/2 T/2 xp(end)];
			yy = [dc(ii) dc(ii) dc(ii-1) dc(ii-1)];
			patch(xx,yy,char(cors(nq==np(1))));
			
	end	
	
end

axis equal
xlabel('y (\mum)');
ylabel('z (\mum)')
s=sprintf('Period T_y = %.1f nm',T*1000);
v=axis;
hold on; 
line([-T/2 -T/2],[v(3) v(4)],'linestyle','--','color',[0 0 0]);
line([T/2 T/2],[v(3) v(4)],'linestyle','--','color',[0 0 0]);
title(s);

end%regeoy


function [qr,qt]=convert_tables_ortho_to_skew(pr,pt)
phis=[-180:180];
%shift orthogonal orders by -60 deg to get skew orders
phi_shift0=-120;
[~, pptr]=min(abs(phis-phi_shift0));

qr=pr;qt=pt;
for ii=1:25
	v0=pr(:,:,ii);
	u1=v0(:,pptr:end);
	u2=v0(:,1:pptr);


	qr(:,:,ii)=[u1 u2(:,2:end)];

	v0=pt(:,:,ii);
	u1=v0(:,pptr:end);
	u2=v0(:,1:pptr);

	qt(:,:,ii)=[u1 u2(:,2:end)];
end
end%convert_tables_ortho_to_skew
