% CCP2D (Common Conversion Point, 2D version)
% Part 2_2D: processing
% 2-Dimensional version, based on ccp.m
% Gyorgy HETENYI, 01 Jul 2005
% Revised 08 Apr 2008 (velocity-models updated!)

% ___________________________________
% ---------- Parameters -------------
% velocity-model for ray-tracing:
%   1 for iasp91
%   2 for ak135
%   3 for PREM
%   structure name .mat otherwise
tic
disp('"ccp2d":');
Mod=1;
%Mod='/home/hetenyi/HICLIMB/Models/HC2D_6wVPVS_ALT2.mat';
Mod='/home/gh/EASI/Data/Vmodel/EPcrustEASI.mat';

represent=1; % to make image or not

% ray-tracing parameters
inc=0.5;          % step in depth
zmax=100;       % maximum depth

% determine study area (x -> perpendicular to the profile)
minx=-200;
maxx=200;
pasx=maxx-minx; % pas [FR] = step [EN]
miny=-50;
maxy=600;
pasy=1;
minz=-2;
maxz=100;       % always greater than "zmax"
pasz=0.5;
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
a=1;b=.5;       % ellipticity factor

% color-palette definition and loading
name_pal_col='/home/gh/Matlab/Tools/cpt_blue2red';
pal_col=load(name_pal_col);
pal_col=interp_pal(pal_col,10);

nbtr=length(TR);

%_______________________________________
% reading and loading the velocity model
clear M
if ischar(Mod)
   isl=find(Mod=='/');
   Modtxt=Mod(isl(length(isl))+1:length(Mod));
   Vmod=4; disp(['Velocity-model is: ' Modtxt ])
elseif Mod==1
   Vmod=1; disp('Velocity-model is: iasp91')
elseif Mod==2
   Vmod=2; disp('Velocity-model is: ak135')
elseif Mod==3
   Vmod=3; disp('Velocity-model is: PREM')
end

switch Vmod
case(1)
   Z=[0:inc:zmax]; VP=pvelIASP(Z); VS=svelIASP(Z);
case(2)
   Z=[0:inc:zmax]; VP=pvel_ak135(Z); VS=svel_ak135(Z);
case(3)
   Z=[0:inc:zmax]; VP=pvelPREM(Z); VS=svelPREM(Z);
case(4)
   load(Mod);
   NVmy=length(M.yy); NVmz=length(M.zz); % Vel.mod size
  %extend velocity-model with boundary values if the box is larger
   if miny<min(M.yy) % left
      M.yy=[miny M.yy];
      M.VP=[M.VP(:,1) M.VP];
      M.VS=[M.VS(:,1) M.VS];
      NVmy=length(M.yy);
   end
   if maxy>max(M.yy) % right
      M.yy=[M.yy maxy];
      M.VP=[M.VP M.VP(:,NVmy)];
      M.VS=[M.VS M.VS(:,NVmy)];
      NVmy=length(M.yy);
   end

  % interpolate velocities:
   M.VP0=M.VP;	M.VS0=M.VS; M.zz0=M.zz;

   M.VP=interp1(M.zz,M.VP,[0:inc:max(M.zz)]);
   M.VS=interp1(M.zz,M.VS,[0:inc:max(M.zz)]);
   M.zz=[0:inc:max(M.zz)];
   NVmz=length(M.zz);
   
   if zmax>max(M.zz) % under
      znew=M.zz(NVmz)+inc:inc:zmax;
      M.zz=[M.zz znew];
      for iy=1:length(M.yy)
         M.VP(NVmz+1:length(M.zz),iy)=pvelIASP(znew);
         M.VS(NVmz+1:length(M.zz),iy)=svelIASP(znew);
      end
      NVmz=length(M.zz);
   end
      % Note: migration is performed from the surface, z=0.
      % Velocity-model should be understood from the surface!!!
      % Topo correction should be done on the traces at reading.
      % Topo correction in ccp is performed below.
      % 2D matrix dimensions are: (z,y)
   Z=[0:inc:zmax];
end

% pseudo-2-dimensionalisation for 1-D models
if Vmod==1 | Vmod==2 | Vmod==3
   M.yy=[miny:pasy:maxy];
   M.zz=Z;
   for iv=1:length(M.yy)
      M.VP(:,iv)=VP;
      M.VS(:,iv)=VS;
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
%  if mod(i,100)==0,disp([num2str(i) ' / ' num2str(nbtr) ' processed']);end		
  if mod(i,10)==0
    txtout=[' ' num2str(i) ' / ' num2str(nbtr) ' processed']; lt=length(txtout);
    for iii=1:lt   fprintf('\r');   end;      fprintf(txtout);
  end

  if TR(i).prai > -1

    clear Xp;clear Xs;clear Yp;clear Ys;
    clear Tp;clear Ts;clear Td;clear Te;

    p=TR(i).prai/111.19;              % ray-parameter
    Xs(1)=TR(i).x0;Ys(1)=TR(i).y0;    % S-ray @ surface = @ station
    Xp(1)=TR(i).x0;Yp(1)=TR(i).y0;    % P-ray @ surface = @ station
    Tp(1)=0;Ts(1)=0;Td(1)=0;Te(1)=0;  % Arr-times @ surface (0)
    coslbaz=cos(TR(i).lbaz*pi/180);   % Local back-azimuth Y-comp.
    sinlbaz=sin(TR(i).lbaz*pi/180);   % Local back-azimuth X-comp.

   %______________% Migrate with the 2-D velocity-model__________

    incidp=asin(p*M.VP);  % P incidence-angle matrix
    incids=asin(p*M.VS);  % S incidence-angle matrix

    for iz=1:length(Z)-1  % for each layer: compute next one

     % find neighbouring indices for the y-direction
      yok=find(Yp(iz)<M.yy); yok2=yok(1); yok1=yok2-1;

     % take the distance-weighted mean of the two neighbours
      d1=( Yp(iz) - M.yy(yok1) );  % distances from grid sides
      d2=( M.yy(yok2) - Yp(iz) );

      if d1==0; w2=0; w1=1;        % if it falls on point 1
      elseif d2==0; w1=0; w2=1;    % if it falls on point 2
      else w1=1-d1/(d1+d2); w2=1-d2/(d1+d2);  % normal...
      end


     % horizontal displacement at this step
      Ss=( w1*tan(incids(iz,yok1)) + w2*tan(incids(iz,yok2)) ) *inc;

      Sp=( w1*tan(incidp(iz,yok1)) + w2*tan(incidp(iz,yok2)) ) *inc;

     % position on the next layer
      Xs(iz+1)=Xs(iz) + sinlbaz*Ss;
      Ys(iz+1)=Ys(iz) + coslbaz*Ss;
      Xp(iz+1)=Xp(iz) + sinlbaz*Sp;
      Yp(iz+1)=Yp(iz) + coslbaz*Sp;

     % following the above scheme...:
      cosincS= w1*cos(incids(iz,yok1)) + w2*cos(incids(iz,yok2)) ;

      cosincP= w1*cos(incidp(iz,yok1)) + w2*cos(incidp(iz,yok2)) ;

      VSloc= w1*M.VS(iz,yok1) + w2*M.VS(iz,yok2) ;

      VPloc= w1*M.VP(iz,yok1) + w2*M.VP(iz,yok2) ;

     % traveltimes
      Ts(iz+1)=Ts(iz) + (inc/cosincS)/VSloc;
      Tp(iz+1)=Tp(iz) + (inc/cosincP)/VPloc;
    end
    %____________________end of 2D migration_______

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
krms=0;
for i=1:nbtr
   if TR(i).prai > -1 & rms(TR(i).trace)<=rms_max
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
