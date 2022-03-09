% CCP (Common Conversion Point)
% cf Zhu, EPSL, 179 -2000
% Part 2: processing
% Revised 08 Apr 2008 (velocity-models updated!)

% ___________________________________
% ---------- Parameters -------------
% velocity-model for ray-tracing:
%   1 for iasp91
%   2 for ak135
%   3 for PREM
%   structure name .mat otherwise (not yet implemented)
tic
disp('"ccp":');
Mod=1; zMoho=40; % iasp91 model with Moho down to zMoho depth
%Mod='/home/gh/EASI/Data/Bianchi2014_1Dvmod.mat';

represent=1; % to make image or not

% ray-tracing parameters
inc=1;        % step in depth
zmax=800;       % maximum depth

% determine study area (x -> perpendicular to the profile)
minx=-1200;
maxx=1200;
pasx=maxx-minx;
miny=-1150;
maxy=1600;
pasy=2;
minz=-2;
maxz=800;       % always >= zmax
pasz=2;
if exist('rms_max')==0
   rms_max=0.1; disp('rms_max autoset in ccp.m!!!!')
end


% ---------------------No changes below this line----------------
% ---------------------------------------------------------------
% ------------------------ MAIN PROGRAM -------------------------

if maxz<zmax disp('Problem: maxz < zmax !!'); beep; beep; return; end
if pasz<inc  disp('Problem: pasz < inc !!'); beep; beep; return; end

% parameters of the initial "mini" gaussian filter
% (not accounted for in ccp_smooth)
% Note: parameter unit = number of samples (NOT km!)
mu=0;sigma=1; nbm=3;
a=1;b=.5;	% ellipticity factor

% color-palette definition and loading
name_pal_col='/home/gh/Matlab/Tools/cpt_blue2red';
pal_col=load(name_pal_col);
pal_col=interp_pal(pal_col,10); 

nbtr=length(TR);

%_______________________________________
% reading and loading the velocity model
if ischar(Mod)
   isl=find(Mod=='/');
   Modtxt=Mod(isl(length(isl))+1:length(Mod));
   Vmod=4; disp(['Velocity-model is: ' Modtxt ])
elseif Mod==1
   Vmod=1;
   if exist('zMoho')
      disp(['Velocity-model is: iasp91 with Moho at ' num2str(zMoho) 'km'])
   else
      disp('Velocity-model is: iasp91')
   end
elseif Mod==2
   Vmod=2; disp('Velocity-model is: ak135')
elseif Mod==3
   Vmod=3; disp('Velocity-model is: PREM')
end

switch Vmod
case(1)
   Z=[inc:inc:zmax]; VP=pvelIASP(Z,zMoho); VS=svelIASP(Z,zMoho); Z=[0 Z];
case(2)
   Z=[inc:inc:zmax]; VP=pvel_ak135(Z); VS=svel_ak135(Z); Z=[0 Z];
case(3)
   Z=[inc:inc:zmax]; VP=pvelPREM(Z); VS=svelPREM(Z); Z=[0 Z];
case(4)
   load(Mod);
   if exist('VP0')+exist('VS0')+exist('Z0')~=3
      disp('Please define your velocity model that has Z, VP, VS!')
   else
      Z=[inc:inc:zmax];
      VP=interp1(Z0,VP0,Z,'nearest','extrap');
      VS=interp1(Z0,VS0,Z,'nearest','extrap');
      Z=[0 Z];
      clear VP0 VS0 Z0
   end
end
%_______________end of V.

for i=1:nbtr
   TR(i).Xp=[];
   TR(i).Xs=[];
   TR(i).Yp=[];
   TR(i).Ys=[];
   TR(i).Z=[];
end

%_______________________________________
% Launch rays
disp('Ray tracing')

for i=1:nbtr	
%   if mod(i,100)==0,disp([' ' num2str(i) ' / ' num2str(nbtr) ' processed']);end		
   if mod(i,100)==0
      txtout=[' ' num2str(i) ' / ' num2str(nbtr) ' processed']; lt=length(txtout);
      for iii=1:lt   fprintf('\r');   end;      fprintf(txtout);
   end

   if TR(i).prai > -1

      clear Xp;clear Xs;clear Yp;clear Ys;
      clear Tp;clear Ts;clear Td;clear Te;
	
      p=TR(i).prai/111.19;
      Xs(1)=TR(i).x0;Ys(1)=TR(i).y0;
      Xp(1)=TR(i).x0;Yp(1)=TR(i).y0;
      Tp(1)=0;Ts(1)=0;Td(1)=0;Te(1)=0;
      coslbaz=cos(TR(i).lbaz*pi/180);
      sinlbaz=sin(TR(i).lbaz*pi/180);

   %______________% Migrate with the 1-D velocity-model__________

      incidp=asin(p*VP);  % P incidence-angle matrix
      incids=asin(p*VS);  % S incidence-angle matrix

      Ss=tan(incids).*inc; % horizontal displacement
      Sp=tan(incidp).*inc;
      Xs=[Xs sinlbaz*Ss]; % position on the next layer
      Ys=[Ys coslbaz*Ss];
      Xp=[Xp sinlbaz*Sp];
      Yp=[Yp coslbaz*Sp];
      Xs=cumsum(Xs);Ys=cumsum(Ys);
      Xp=cumsum(Xp);Yp=cumsum(Yp);

      Tp=[Tp (inc./cos(incidp))./VP]; % traveltimes
      Ts=[Ts (inc./cos(incids))./VS];
      Tp=cumsum(Tp);
      Ts=cumsum(Ts);
    %____________________end of 1D migration_______

      D=sqrt((Xp-Xs).^2+(Yp-Ys).^2);
      E=sqrt((Xp-Xp(1)).^2+(Yp-Yp(1)).^2);

      Td=D.*p;
      Te=2*E.*p;

      TR(i).Z=Z+TR(i).z0;
      TR(i).Xp=Xp;
      TR(i).Yp=Yp;
      TR(i).Xs=Xs;
      TR(i).Ys=Ys;

      TR(i).amp_ps=interp1(TR(i).time',TR(i).trace,-Tp+Ts+Td);
      TR(i).amp_pps=interp1(TR(i).time',TR(i).trace,Tp+Ts+Td-Te);
      TR(i).amp_pss=interp1(TR(i).time',TR(i).trace,2*Ts+2*Td-Te);

     % theoretical traces
      tps=-Tp+Ts+Td;  tpps=Tp+Ts+Td-Te;  tpss=2*Ts+2*Td-Te;
      TR(i).amp_pps_theo=interp1(tpps,TR(i).amp_ps,tps);
      TR(i).amp_pss_theo=interp1(tpss,TR(i).amp_ps,tps);
	 
   else
      TR(i).Xp=-1;
      TR(i).Yp=-1;
      TR(i).Xs=-1;
      TR(i).Ys=-1;
      TR(i).Z=-1;
      TR(i).amp_ps=-1; TR(i).amp_pps=-1; TR(i).amp_pss=-1;
   end
end

%_____________________________________________________
% build-up 3D matrix "G"
% -------------------------------

fprintf('\n')
disp('Building the G-matrix')

xx=[minx:pasx:maxx];
yy=[miny:pasy:maxy];
zz=[minz:pasz:maxz];
G=zeros(length(xx),length(yy),length(zz));
nG=zeros(length(xx),length(yy),length(zz))+1e-8;
%[XX,YY,ZZ]=meshgrid(xx,yy,zz);	
krms=0;
for i=1:nbtr
   if TR(i).prai > -1  & rms(TR(i).trace)<=rms_max
      ix=floor((TR(i).Xs-minx)/pasx+1);
      iy=floor((TR(i).Ys-miny)/pasy+1);
      iz=floor((TR(i).Z-minz)/pasz+1);
      index=sub2ind(size(G),ix,iy,iz);	
      G(index)=G(index)+TR(i).amp_ps;			
      nG(index)=nG(index)+1;
   else
      %disp([num2str(i) ', ' num2str(rms(TR(i).trace))]);
      krms=krms+1;
   end
end

% Change to 2D
G2=squeeze(sum(G,1));
nG2=squeeze(sum(nG,1));
G2=G2./nG2;

tempnG2=nG2;tempnG2(find(nG2<1))=NaN;

% Mini gaussian filter to represent G2
C=zeros(nbm,nbm);
mm=round(nbm/2);
for i=1:nbm
   for j=1:nbm
      dist=sqrt((i-mm)*(i-mm)/(a*a)+(j-mm)*(j-mm)/(b*b));
      if dist<=mm
         C(i,j)=1/(sigma*sqrt(2*pi))*exp(-.5*(dist*dist-mu*mu)/(sigma*sigma));
      end
   end
end
C=C/sum(sum(C));
G2fil=conv2(G2,C,'same');

if represent==1
   Gp=G2fil; CL=8;
   ccp_plot;
else
  close
end 
disp(['Eliminated because of high rms-value: ' num2str(krms) ])
toc
