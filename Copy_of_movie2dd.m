
% 2d Movie - evolution


%% 
FIG=1
figure(FIG)
close(FIG)
figure(FIG)

    set(FIG, 'PaperUnits', 'inches')
    
    set(gcf, 'PaperPosition', [0 0 12.0 6.0])
    set(FIG, 'DefaultAxesFontSize', 20)
    set(FIG, 'DefaultTextFontSize', 20)
    set(FIG, 'DefaultLineMarkerSize', 4)
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

datadir = './';
%datadir = '/Volumes/Elements/PSC_DATA/shocks/derek/try19/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb/nnb001/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb_coll001/nnb004_bubble/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll01/';
%datadir = '/Volumes/Elements/PSC_DATA/try_par/nif/coll0/';


% will add .png, .mov as necessary
outfile = 'movie2'

mkdir ([datadir, '/', outfile, '.pngs']);
export_pattern = [datadir, '/', outfile, '.pngs/m%06d.png'];


%[MMi, TTe, BB0, LL0, V0, ts] = getBubbleXZRunParams(RUN)

MMi =100;
ZZ = 1;
TTe = 0.04;
n= 0.01;
LL0 = 40; %sqrt(MMi/(ZZ*n));
BB0 = 0.01;
V0 = BB0/sqrt(MMi);

tstart = 0;
tstep = 2000;
ts = [tstart:tstep:42000];



%V0 = V0 * sqrt(TTe/MMi);

VA = BB0 / sqrt(MMi);
Cs = sqrt(TTe/MMi);

cs1 = 'k';
cs2 = 'k';  % or k--

do_boundscheck = 0

for k=1:length(ts)
    
  clf
  k
  
%   d = loadmats( '%s/mats/tfields_%07d.mat' \ {datadir, ts(k)}, ...
%                 '%s/mats/etmoments_%07d.mat' \ {datadir, ts(k)}, ...
%                 '%s/mats/itmoments_%07d.mat' \ {datadir, ts(k)} );
%   
%   
%d = h5read_all( '%s/h5/psc_%07d_sm.h5' \ {datadir, ts(k)});
    address=strcat(datadir, 'psc_',num2str(ts(k),'%07d'),'.h5')
 %   try
        %dfield = load(address,'-mat'); 
        NN=h5read(address,'/NNe');
        dx=h5read(address,'/dx');
        dz=h5read(address,'/dz');
        dt=h5read(address,'/dt');
        
        
       size(NN(:,1,:))
       
        
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
  xs = h5read(address,'/xs')/ sqrt(MMi/n);
  zs = h5read(address,'/zs')/ sqrt(MMi/n);

  bx = h5read(address,'/bx');
  by = h5read(address,'/by');
  bz = h5read(address,'/bz');
  
  jy = h5read(address,'/jy');
  
  ey = h5read(address,'/ey');
  
          NVxe=h5read(address,'/NVxe');
        NVye=h5read(address,'/NVye');
        NVze=h5read(address,'/NVze');
  
  Sxxe = h5read(address,'/Sxxe');
  %Sxye = h5read(address,'/Sxye');
  %Sxze = h5read(address,'/Sxze');
  %Syxe = h5read(address,'/Syxe');
  Syye = h5read(address,'/Syye');
  %Syze = h5read(address,'/Syze');
  %Szxe = h5read(address,'/Szxe');
  %Szye = h5read(address,'/Szye');
  Szze = h5read(address,'/Szze');
  
  
  
  NVxe = h5read(address,'/NVxe');
  NVye = h5read(address,'/NVye');
  NVze = h5read(address,'/NVze');
  
  NVxi = h5read(address,'/NVxi');
  NVyi = h5read(address,'/NVyi');
  NVzi = h5read(address,'/NVzi');

  
  size(squeeze(NVxe))
  
  %midx = round(size(bz,2)/2);
  
  %Psi = -repmat( cumtrapz(xs, bz(:,:,:,1),2), [1 1 1 size(bz,4)] );
  %Psi = Psi + repmat( trapz(xs(1:midx), bz(:,1:midx,:,1),2), [1 size(bx,2), 1, size(bx,4)]) ;
  %Psi = Psi + cumtrapz(zs, bx(:,:,:,:), 4);
  
  %maxPsi = max(abs(Psi(:)));
  
  %cPsi = [-.95:0.1:.95]*maxPsi;
 
  zmin = round(min(zs));
  zmax = round(max(zs));
  xmin = round(min(xs));
  xmax = round(max(xs));
  
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
    
  
  BB = sqrt(bx(1,:,:,:).^2+by(1,:,:,:).^2+bz(1,:,:,:).^2);
    
  Txxe = (Sxxe(1,:,:) - NVxe(1,:,:).^2./NN(1,:,:)) ./ NN(1,:,:);
  Tyye = (Syye(1,:,:) - NVye(1,:,:).^2./NN(1,:,:)) ./ NN(1,:,:);
  Tzze = (Szze(1,:,:) - NVze(1,:,:).^2./NN(1,:,:)) ./ NN(1,:,:);
  
  Te = .333 * (Txxe + Tyye + Tzze);
  Te (Te < 0) = 1e-6;
  
  fig_title = strcat('t C_s/L = ', num2str(ts(k)*dt * Cs / (LL0*sqrt(MMi)), 3));
    
  %xs
  axis([xmin, xmax, zmin, zmax])
 
 % size(squeeze(xs))
  %size(squeeze(zs))
  %size(squeeze(bx)')
  
  xs;
 
  zs;
  
  h=figure;
  
  imagesc(xs, zs, squeeze(by)./BB0)
  colorbar
  colormap(map)
  caxis([-5 5])
  hold on
  size(squeeze(bx))    
 

   
  %figure('Visible','off');
  %view(2);
  saveas(h,strcat(datadir,'movie2.pngs/','Bx', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);   
  
  
  h=figure;
  
  size(zs)
  size(xs)
  size(Te)
  
  imagesc(xs, zs, squeeze(NN/n))
  colorbar


%  caxis([-2 0.2])
 % colormap(map)
  hold on
  size(squeeze(bx))
   
  %figure('Visible','off');
  %view(2);
  saveas(h,strcat(datadir,'movie2.pngs/','nne', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);
  
  
  
  h=figure;
  
  size(xs)
  size(zs)
  size(squeeze(NVxe))
  size(squeeze(NN))
  
  imagesc(xs, zs, squeeze(NVxe))
  colorbar


%  caxis([-2 0.2])
 % colormap(map)
  hold on
  size(squeeze(bx))
   
  %figure('Visible','off');
  %view(2);
  saveas(h,strcat(datadir,'movie2.pngs/','vxe', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);
  
  
  
    h=figure;
  
  size(xs)
  size(zs)
  size(squeeze(NVxe))
  size(squeeze(NN))
 
  imagesc(xs, zs, squeeze(NVxi))
  colorbar


%  caxis([-2 0.2])
 % colormap(map)
  hold on
  size(squeeze(bx))
   
  %figure('Visible','off');
  %view(2);
  saveas(h,strcat(datadir,'movie2.pngs/','vxi', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);
  

end

quit
    %system(sprintf('/usr/local/bin/mencoder mf://%s/%s.pngs/*.png -mf fps=12:type=png -vf scale=960:720 -ovc x264  -x264encopts  qp=0  -o %s/%s.avi', datadir, outfile, datadir, outfile))

    %system('osascript -e ''tell application "QuickTime Player" to close every document whose name contains "%s.mov" '' ' \ {outfile});

    %system('rm %s/%s.mov' \ {datadir, outfile});
    %img2qt('%s/%s.mov' \ {datadir, outfile}, '1', '12', '%s/%s.pngs/*.png' \ {datadir, outfile});

    %system('open %s/%s.mov' \ {datadir, outfile} );    
