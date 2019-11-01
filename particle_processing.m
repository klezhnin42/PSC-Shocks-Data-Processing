
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
particle='ion';
outfile = particle;

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

tstart = 0;
tstep = 1000;
ts = [tstart:tstep:100000];


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
%                 '%s/mats/itmoments_%07d.mat' \ {datadir, ts(k)} );
%  
%   
%d = h5read_all( '%s/h5/psc_%07d_sm.h5' \ {datadir, ts(k)});
    address=strcat(datadir, particle,'_',num2str(ts(k),'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 
     
%        h5disp(address);
        
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x');
        y=h5read(address,'/y');
        z=h5read(address,'/z');
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

        ey = h5read(address,'/ey');




       
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
  
size(x);
size(px);


xedge = 80 * [0:0.01:1];
yedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];


      h=figure;
  
%  scatter(x,pz,0.1)

 hgrm=histc2d(x/100,pz,xedge,yedge);
 loghgrm=log10(hgrm);
 loghgrm(~isfinite(loghgrm))=0;
imagesc(xedge,yedge,loghgrm');
 colorbar
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
  hold on
 % size(squeeze(bx))
   xlim([0 80]);
   xlabel('x/d_i');
    ylabel('p_z/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

   meanNNi=mean(NNi,1);
   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


  saveas(h,strcat(datadir,'movie2.pngs/',particle,'xpz', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);  



      h=figure;

%  scatter(x,pz,0.1)

 hgrm=histc2d(x/100,px,xedge,yedge);
 loghgrm=log10(hgrm);
 loghgrm(~isfinite(loghgrm))=0;
imagesc(xedge,yedge,loghgrm');
 colorbar
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
  hold on
 % size(squeeze(bx))
   xlim([0 80]);
   xlabel('x/d_i');
    ylabel('p_x/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

   meanNNi=mean(NNi,1);
   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


  saveas(h,strcat(datadir,'movie2.pngs/',particle,'xpx', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);



      h=figure;

%  scatter(x,pz,0.1)

 hgrm=histc2d(x/100,py,xedge,yedge);
 loghgrm=log10(hgrm);
 loghgrm(~isfinite(loghgrm))=0;
imagesc(xedge,yedge,loghgrm');
 colorbar
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
  hold on
 % size(squeeze(bx))
   xlim([0 80]);
   xlabel('x/d_i');
    ylabel('p_y/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

   meanNNi=mean(NNi,1);
   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


  saveas(h,strcat(datadir,'movie2.pngs/',particle,'xpy', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);  



xedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];
yedge = max(sqrt(px.^2+py.^2+pz.^2)) * [-1:0.01:1];


      h=figure;

%  scatter(x,pz,0.1)

 hgrm=histc2d(px,pz,xedge,yedge);
 loghgrm=log10(hgrm);
 loghgrm(~isfinite(loghgrm))=0;
imagesc(xedge,yedge,loghgrm');
 colorbar
set(gca,'YDir','normal')
 %colormap(flipud(gray))
 %caxis([0 4])
  hold on
 % size(squeeze(bx))
%  xlim([0 50]);
   xlabel('p_x/m_ic');
    ylabel('p_z/m_ic');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

%   meanNNi=mean(NNi,1);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


  saveas(h,strcat(datadir,'movie2.pngs/',particle,'pxpz', num2str(ts(k),'%07d'),'_4.png'));
  delete(h);





  

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
 %  xlim([-1 1]);
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
   % legend2 = sprintf('Protons');
   % legend3 = sprintf('Electrons');
   % legend({legend1, legend2, legend3});
%       title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
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
