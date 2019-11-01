
% 2d Movie - evolution


%% 
%FIG=1
%figure(FIG)
%close(FIG)
%figure(FIG)

%    set(FIG, 'PaperUnits', 'inches')
    
%    set(gcf, 'PaperPosition', [0 0 12.0 6.0])
%    set(FIG, 'DefaultAxesFontSize', 20)
%    set(FIG, 'DefaultTextFontSize', 20)
%    set(FIG, 'DefaultLineMarkerSize', 4)
    %set(FIG, 'DefaultLineLineWidth', 1);


%% subplot sizes and layout

wpos = [%[0.25  0.60  0.6 0.35];
        %[0.25  0.11  0.6 0.35];
        [0.09 0.60  0.34 0.30];
        [0.09 0.14  0.34 0.30];
        [0.57 0.60  0.34 0.30];
        [0.57 0.14  0.34 0.30];
 %       [0.7 0.55  0.28 0.25];
 %       [0.7 0.15  0.28 0.25];        
    ]
%wpos = [%[0.38  0.5  0.34 0.30];
        %[0.38  0.05  0.34 0.30];
%        [0.19 0.58  0.68 0.34];
%        [0.19 0.12  0.68 0.34];
%        [0.19 0.58  0.68 0.34];
%        [0.19 0.12  0.68 0.34];
        %[0.7 0.55  0.28 0.25];
        %[0.7 0.15  0.28 0.25];        
%    ]

%wpos = [[0.38  0.5  0.34 0.30];
 %       [0.38  0.05  0.34 0.30];
  %      [0.19 0.58  0.68 0.34];
   %     [0.18 0.23  0.74 0.54];
        %[0.7 0.55  0.28 0.25];
        %[0.7 0.15  0.28 0.25];        
   % ]


wM = 2
wN = 2


%%

export_cmd = {'-dpng',  '-r200' };


%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll04/';

datadir = './h5/';
%datadir = '/Volumes/Elements/PSC_DATA/shocks/derek/try19/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb/nnb001/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb_coll001/nnb004_bubble/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll01/';
%datadir = '/Volumes/Elements/PSC_DATA/try_par/nif/coll0/';


% will add .png, .mov as necessary
outfile = 'movie2';

mkdir ([datadir, '/', outfile, '.pngs']);
export_pattern = [datadir, '/', outfile, '.pngs/m%06d.png'];


%[MMi, TTe, BB0, LL0, V0, ts] = getBubbleXZRunParams(RUN)

MMi =100;
ZZ = 1;
TTe = 0.002;
n= 0.01;
LL0 = 40; %sqrt(MMi/(ZZ*n));
BB0 = 0.01;
V0 = BB0/sqrt(MMi*n);

tstart = 173000;
tstep = 1000;
ts = [tstart:tstep:173000];

particle='ele';

%V0 = V0 * sqrt(TTe/MMi);

VA = BB0 / sqrt(MMi*n);
Cs = sqrt(TTe/MMi);

cs1 = 'k';
cs2 = 'k';  % or k--

do_boundscheck = 0;

for k=1:length(ts)
  clf
  k
  
%   d = loadmats( '%s/mats/tfields_%07d.mat' \ {datadir, ts(k)}, ...
%                 '%s/mats/etmoments_%07d.mat' \ {datadir, ts(k)}, ...
%                 w'%ws/mats/itmoments_%07d.mat' \ {datadir, ts(k)} );
%  
%   
%d = h5read_all( '%s/h5/psc_%07d_sm.h5' \ {datadir, ts(k)});
%   address=strcat(datadir, particle,'_',num2str(ts(k),'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 
     
%        h5disp(address);
        
%        px=h5read(address,'/px');
%        py=h5read(address,'/py');
%        pz=h5read(address,'/pz');
%        x=h5read(address,'/x');
%        y=h5read(address,'/y');
%        z=h5read(address,'/z');
 %       dt=h5read(address,'/dt');
        
%    address0=strcat(datadir, 'ion_',num2str(0,'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 

%        h5disp(address0);

%        px0=h5read(address0,'/px');
%        py0=h5read(address0,'/py');
%        pz0=h5read(address0,'/pz');
%        x0=h5read(address0,'/x');
%        y0=h5read(address0,'/y');
%        z0=h5read(address0,'/z');



	address=strcat(datadir, 'psc_',num2str(ts(k),'%07d'),'.h5');

        NNe=h5read(address,'/NNe');
        NNi=h5read(address,'/NNi');
        dx=h5read(address,'/dx');
        dz=h5read(address,'/dz');
        dt=h5read(address,'/dt');

%        h5disp(address);

        NVxe=h5read(address,'/NVxe');
        NVye=h5read(address,'/NVye');
        NVze=h5read(address,'/NVze');

        Sxxe=h5read(address,'/Sxxe');
        Syye=h5read(address,'/Syye');
        Szze=h5read(address,'/Szze');

        xs = h5read(address,'/xs')/ sqrt(MMi/n);
        zs = h5read(address,'/zs')/ sqrt(MMi/n);

        bx = h5read(address,'/bx');
        by = h5read(address,'/by');
        bz = h5read(address,'/bz');

        jy = h5read(address,'/jy');

        ex = h5read(address,'/ex');
        ey = h5read(address,'/ey');
        ez = h5read(address,'/ez');

        vez=squeeze(NVze./NNe);         

        TTxxe = (Sxxe - NVxe.^2./NNe) ./ NNe;
        TTyye = (Syye - NVye.^2./NNe) ./ NNe;
        TTzze = (Szze - NVze.^2./NNe) ./ NNe;

        Te = .333 * (TTxxe + TTyye + TTzze);
        Te (Te < 0) = 1e-6;
       
  %     size(NN(:,1,:))
       
        
   %     h5disp(address)
      %  size(dfield);
 %   catch
 %       continue;
 %   end
    
 %   address=strcat(datadir, 'psc_',num2str(ts(k),'%07d'),'.h5')
 %   try
 %       demom = h5read(address); 
 %   catch
  %      continue;
  %  end
    
    
 %  size(dfield)
    
 %   address=strcat(datadir, 'itmoments_',num2str(ts(k),'%07d'),'.mat')
 %   try
 %       dimom = load(address,'-mat'); 
 %   catch
 %       continue;
 %   end
    
  %  dx=dfield.dx;
  %  dz=dfield.dz; 
%  xs = h5read(address,'/xs')/ sqrt(MMi);
%  zs = h5read(address,'/zs')/ sqrt(MMi);

%  bx = h5read(address,'/bx');
%  by = h5read(address,'/by');
%  bz = h5read(address,'/bz');
  
%  jy = h5read(address,'/jy');
  
%  ey = h5read(address,'/ey');
  
  
%  Sxxe = h5read(address,'/Sxxe');
  %Sxye = h5read(address,'/Sxye');
  %Sxze = h5read(address,'/Sxze');
  %Syxe = h5read(address,'/Syxe');
%  Syye = h5read(address,'/Syye');
  %Syze = h5read(address,'/Syze');
  %Szxe = h5read(address,'/Szxe');
  %Szye = h5read(address,'/Szye');
%  Szze = h5read(address,'/Szze');
  
  
  
%  NVxe = h5read(address,'/NVxe');
%  NVye = h5read(address,'/NVye');
%  NVze = h5read(address,'/NVze');
  
%  size(squeeze(NVxe))
  
  %midx = round(size(bz,2)/2);
  
  %Psi = -repmat( cumtrapz(xs, bz(:,:,:,1),2), [1 1 1 size(bz,4)] );
  %Psi = Psi + repmat( trapz(xs(1:midx), bz(:,1:midx,:,1),2), [1 size(bx,2), 1, size(bx,4)]) ;
  %Psi = Psi + cumtrapz(zs, bx(:,:,:,:), 4);
  
  %maxPsi = max(abs(Psi(:)));
  
  %cPsi = [-.95:0.1:.95]*maxPsi;
 
%  zmin = round(min(zs));
%  zmax = round(max(zs));
%  xmin = round(min(xs));
%  xmax = round(max(xs));
  
  map=zeros(300, 3);
    for i = 1:150
        map(i,1) = (i-1)/150;
        map(i,2) = (i-1)/150;
        map(i,3) = 1;%exp(log(10)*i/100)/1000;
    end
    for i = 151:300
        map(i,1) = 1;
        map(i,2) = 1-(i-151)/150;
        map(i,3) = 1-(i-151)/150;%exp(log(10)*i/100)/1000;
    end
    
    %map
    
     map2=zeros(300, 3);
    for i = 1:300
        map2(i,1) = 1-(i-1)/299;
        map2(i,2) = 1-(i-1)/299;
        map2(i,3) = 1;%exp(log(10)*i/100)/1000;
    end
%       


 % NN = demom.NNe;
    
  
 
 
 % size(squeeze(xs))
  %size(squeeze(zs))
  %size(squeeze(bx)')
  
%size(x);
%size(px);


%xedge = 50 * [0:0.01:1];
%yedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];



     h=figure;
%        FIG=1
%et(gca,'fontsize', 28);
 %   figure(FIG)
%     close(FIG)
%     figure(FIG)
%    clf

%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

%  scatter(x,pz,0.1)
set(gca,'FontSize',18)
% hgrm=histc2d(x/100,pz,xedge,yedge);
 lognne=log10(NNe/n);
% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(lognne));
axis equal
colorbar
set(gca,'YDir','normal')
 colormap(jet)
 caxis([-0.2 1.2])
  hold on
 % size(squeeze(bx))
%  xlim([0 50])
%bbcx=squeeze(bx)/BB0;
%bbcz=squeeze(bz)/BB0;

size(xs)
size(zs)
%size(bbcx)
%size(bbcz)


%quiver(xs(5:5:25600),zs(5:5:100),bbcx(5:5:100,5:5:25600),bbcz(5:5:100,5:5:25600))
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%  xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-1 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+1]);
  xlim([-30 30]);
  ylim([0 60]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%   title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )

title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));
  saveas(h,strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png')));

%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png'));


%  saveas(h,strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png'));


%      h=figure;
  
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
% lognne=log10(NNe/n);
% lognne(~isfinite(lognne))=-2;
%imagesc(xs,zs,lognne');
% colorbar
%set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
%   xlabel('x/d_i');
%    ylabel('z/d_i');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/',particle,'lognne', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);  


%      h=figure;

%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
% lognne=log10(NNi/n);
% lognne(~isfinite(lognne))=-2;
%imagesc(xs,zs,lognne');
% colorbar
%set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
%   xlabel('x/d_i');
%    ylabel('z/d_i');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/',particle,'lognni', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);


      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% hgrm=histc2d(x/100,pz,xedge,yedge);

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(bx)/BB0);
axis equal
colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-1 1])
  hold on
 %ssize(squeeze(bx))
%   xlim([0 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   xlim([-50 50])
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%  xlim([-12 12]);
%  ylim([0 40]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
  xlim([-30 30]);
  ylim([0 60]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));

  saveas(h,strcat(datadir,'movie2.pngs/','bx', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','bx', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','bx', num2str(ts(k),'%07d'),'_4.png'));





      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% hgrm=histc2d(x/100,pz,xedge,yedge);

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(by)/BB0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-10 10])
  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   xlim([-50 50]);
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]); 
%  ylim([-1 1]);
  %figure('Visible','off');
  %view(%2)
%  xlim([-12 12]);
%  ylim([0 40]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
  xlim([-30 30]);
  ylim([0 60]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));

  saveas(h,strcat(datadir,'movie2.pngs/','by', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','by', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','by', num2str(ts(k),'%07d'),'_4.png'));


%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','by', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','by', num2str(ts(k),'%07d'),'_4.png'));


      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% hgrm=histc2d(x/100,pz,xedge,yedge);

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(bz)/BB0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-1 1])
  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   xlim([-50 50]);
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%  xlim([-12 12]);
%  ylim([0 40]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
  xlim([-30 30]);
  ylim([0 60]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));

  saveas(h,strcat(datadir,'movie2.pngs/','bz', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','bz', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','bz', num2str(ts(k),'%07d'),'_4.png'));



      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(ex)/BB0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-1 1])
  hold on
 % size(squeeze(bx))
%   xlim([-50 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
%  xlim([-12 12]);
%  ylim([0 40]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));
  xlim([-30 30]);
  ylim([0 60]);

  saveas(h,strcat(datadir,'movie2.pngs/','ex', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);



      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(ey)/BB0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-1 1])
  hold on
 % size(squeeze(bx))
%   xlim([-50 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
%  xlim([-12 12]);
%  ylim([0 40%]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));
  xlim([-30 30]);
  ylim([0 60]);

  saveas(h,strcat(datadir,'movie2.pngs/','ey', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% lognne(~isfinite(lognne))=-2;
imagesc(xs,zs,squeeze(ez)/BB0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-1 1])
  hold on
 % size(squeeze(bx))
%   xlim([-50 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
%  xlim([-12 12]);
%  ylim([0 40]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
  xlim([-30 30]);
  ylim([0 60]);

  saveas(h,strcat(datadir,'movie2.pngs/','ez', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);



      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkerSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% lognne(~isfinite(lognne))=-2;w
imagesc(xs,zs,squeeze(vez)/V0);
axis equal
 colorbar
colormap(map);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([-10 10])
  hold on
 % size(squeeze(bx))
%   xlim([-50 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
%  xlim([-12 12]);
%  ylim([0 40]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));
  xlim([-30 30]);
  ylim([0 60]);

  saveas(h,strcat(datadir,'movie2.pngs/','vez', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','vez', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','vez', num2str(ts(k),'%07d'),'_4.png'));

      h=figure;
set(gca,'fontsize', 18);
%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,pz,xedge,yedge);
%    set(gcf, 'PaperPosition', [0.5 0.5 6 1])
%    set(gcf, 'DefaultAxesFontSize', 14)
%    set(gwcf, 'DefaultTextFontSize', 14)
%    set(gcf, 'DefaultLineMarkewrSize', 4)
%    set(gcf, 'DefaultLineLineWidth', 1);
%    set(gcf, 'renderer', 'painters');

% lognne(~isfinite(lognne))=-2;
log10Te=squeeze(log10(Te/TTe));
imagesc(xs,zs,log10Te);
axis equal
 colorbar
colormap(jet);
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 caxis([0 2])
  hold on
 % size(squeeze(bx))
%   xlim([-50 50]);
   xlabel('x/d_{i,up}');
    ylabel('z/d_{i,up}');
%   ylim([-1 1]);
  %figure('Visible','off');
  %view(2)
%   xlim([3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)-5 3.8+4.8*(ts(k) *(dt * BB0/MMi)-1.8)+5]);
%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
%  xlim([-12 12]);
%  ylim([0 40%%]);
%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
title(strcat('\theta_B=60^\circ, \phi=90^\circ, \Omega_{i,up}t=', num2str((ts+1000*k) *(dt * BB0/MMi),'%.2f')));
  xlim([-30 30]);
  ylim([0 60]);

  saveas(h,strcat(datadir,'movie2.pngs/','Te', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);

%u_in = mat2gray(imread(strcat(datadir,'movie2.pngs/','Te', num2str(ts(k),'%07d'),'_4.png')));
%imwrite(RemoveWhiteSpace(u_in),strcat(datadir,'movie2.pngs/','Te', num2str(ts(k),'%07d'),'_4.png'));



%      h=figure;

%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,px,xedge,yedge);
% loghgrm=log10(hgrm);
% loghgrm(~isfinite(loghgrm))=0;
%imagesc(xedge,yedge,loghgrm');
% colorbar
%set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
%   xlabel('x/d_i');
%    ylabel('p_x/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/',particle,'xpx', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);



%      h=figure;

%  scatter(x,pz,0.1)

% hgrm=histc2d(x/100,py,xedge,yedge);
% loghgrm=log10(hgrm);
% loghgrm(~isfinite(loghgrm))=0;
%imagesc(xedge,yedge,loghgrm');
% colorbar
%set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
%   xlabel('x/d_i');
%    ylabel('p_y/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/',particle,'xpy', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);  



%xedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];
%yedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];


%      h=figure;

%  scatter(x,pz,0.1)

% hgrm=histc2d(px,pz,xedge,yedge);
% loghgrm=log10(hgrm);
% loghgrm(~isfinite(loghgrm))=0;
%imagesc(xedge,yedge,loghgrm');
% colorbar
%set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%  xlim([0 50]);
%   xlabel('p_x/m_ic');
%    ylabel('p_z/m_ic');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/',particle,'pxpz', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);





  

%        h=figure;
  
%  scatter(x,pz,0.1)

% ndhist(x/100 ,px,'log');
% colorbar
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on


%    meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%   title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
 % size(squeeze(bx))
%   xlim([0 50]);
%    xlabel('x/d_i');
%    ylabel('p_x/m_i c, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2);
%  saveas(h,strcat(datadir,'movie2.pngs/','xpx', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h); 

  
%  h=figure;
  
%  scatter(x,pz,0.1)

% ndhist(px,pz,'log');
% colorbar
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
 % w tqqqwxlim([-1 1]);
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2);
%   title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%  saveas(h,strcat(datadir,'movie2.pngs/','pxpz', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);   
  
%    ElectronEnergy=100*(px.^2+py.^2+pz.^2)./2;
%    ElectronEnergy0=100*(px0.^2+py0.^2+pz0.^2)./2;
%    indices = find(abs(PhotonEnergy)>1000000000);
%    PhotonEnergyWeights(indices) = [];
%    PhotonEnergy(indices) = [];
    %PhotonEnergyWeights2=int16(PhotonEnergyWeights./1e10);
    
%    PhotonEnergy=PhotonEnergy./1000000; %to MeV
%    protonEnergy=protonEnergy./1000000;
%    ElectronEnergy=ElectronEnergy*511; %to keV
%    ElectronEnergy0=ElectronEnergy0*511; %to keV
%max(PhotonEnergy);
    
    %size(PhotonEnergyWeights2)
    
    %max(PhotonEnergyWeights2)
    %min(PhotonEnergyWeights2)
    
%    nbins=200;
    
%    EnMin=(min(ElectronEnergy)); 
    
%    EnMax=(max(ElectronEnergy)); %
    
%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy, nbins)
%    hold on
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');   

    %x = EnMin:((EnMax-EnMin)/nbins):EnMax;
    %size(x);
    %size(histw);
    %histw;
    %x(1)=[];
    
    %plot(x,histw,'LineWidth',2)
    
    %bar(x,histw,'FaceColor',[0 .5 .5],'facealpha',.7);
    %histw
    %hold on
    
    
    %EnMin=(min(protonEnergy)); 
    
    %EnMax=(max(protonEnergy)); %
    
    %size(protonEnergyWeights);
    
    %[histw, histv] = histwcv(protonEnergy, protonEnergyWeights, EnMin, EnMax, nbins);
    %h=figure;
    %x = EnMin:((EnMax-EnMin)/nbins):EnMax;
    %size(x);
    %size(histw);
    %histw;
    %x(1)=[];
    
    %plot(x,histw,'LineWidth',2)
    
    %bar(x,histw,'FaceColor',[.5 .5 0],'facealpha',.5);
    
    %hold on
    
    
   % EnMin=(min(electronEnergy)); 
    
   % EnMax=(max(electronEnergy)); %
    
   % size(electronEnergyWeights);
    
   % [histw, histv] = histwcv(electronEnergy, electronEnergyWeights, EnMin, EnMax, nbins);
   % h=figure;
   % x = EnMin:((EnMax-EnMin)/nbins):EnMax;
    %size(x);
    %size(histw);
    %histw;
   % x(1)=[];
    
   % plot(x,histw,'LineWidth',2)
    
    %bar(x,histw,'FaceColor', [.5 0 .5],'facealpha',.1);
   
    %hold on
   % legend1 = sprintf('Photons');
   % legenwwd2 = sprintf('Protons');
   % legend3 = sprintf('Electrons');
   % legend({legend1, legend2, legend3});
%       title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');w
%    xlim([EnMin EnMax]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,'ParticleEnergySpectrum_',num2str(ts(k)),'.png'),'png');
%    delete(h); 

end

quit
    %system(sprintf('/usr/local/bin/mencoder mf://%s/%s.pngs/*.png -mf fps=12:type=png -vf scale=960:720 -ovc x264  -x264encopts  qp=0  -o %s/%s.avi', datadir, outfile, datadir, outfile))

    %system('osascript -e ''tell application "QuickTime Player" to close every document whose name contains "%s.mov" '' ' \ {outfile});

    %system('rm %s/%s.mov' \ {datadir, outfile});
    %img2qt('%s/%s.mov' \ {datadir, outfile}, '1', '12', '%s/%s.pngs/*.png' \ {datadir, outfile});

    %system('open %s/%s.mov' \ {datadir, outfile} );    
