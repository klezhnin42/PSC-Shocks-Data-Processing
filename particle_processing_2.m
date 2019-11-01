
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

%datadir = '/Users/klezhnin/Desktop/shock/sixty/ma30beta04/';
datadir = './h5/'
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
n= 0.05;
LL0 = 40; %sqrt(MMi/(ZZ*n));
BB0 = 0.01;
V0 = BB0/sqrt(MMi*n);

tstart = 0;
tstep = 1000;
%ts = [tstart:tstep:89000];

ts=[200000];

particle='ele';

%V0 = V0 * sqrt(TTe/MMi);

VA = BB0 / sqrt(MMi*n);
Cs = sqrt(TTe/MMi);

cs1 = 'k';
cs2 = 'k';  % or k--

do_boundscheck = 0;

h=figure;

for k=1:length(ts)
    
%  clf
  k
  
%   d = loadmats( '%s/mats/tfields_%07d.mat' \ {datadir, ts(k)}, ...
%                 '%s/mats/etmoments_%07d.mat' \ {datadir, ts(k)}, ...
%                 '%s/mats/itmoments_%07d.mat' \ {datadir, ts(k)} );
%  
%   
%d = h5read_all( '%s/h5/psc_%07d_sm.h5' \ {datadir, ts(k)});
    address=strcat(datadir,particle, '_',num2str(ts(k),'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 
      
%       h5disp(address);
        
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/sqrt(MMi/n);
        tag=h5read(address,'/tag');
 
 
        
%    address0=strcat(datadir,particle, '_',num2str(0,'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 

%       h5disp(address0);

%        px0=h5read(address0,'/px');
%        py0=h5read(address0,'/py');
%        pz0=h5read(address0,'/pz');
%        x0=h5read(address0,'/x')/sqrt(MMi/n);
%        y0=h5read(address0,'/y');
%        z0=h5read(address0,'/z');



	address=strcat(datadir, 'psc_',num2str(ts(k),'%07d'),'.h5');

        NNe=h5read(address,'/NNe');
        NNi=h5read(address,'/NNi');
        dx=h5read(address,'/dx');
        dz=h5read(address,'/dz');
        dt=h5read(address,'/dt');

    %    h5disp(address);

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


%      h=figure;
  
%  scatter(x,pz,0.1)

% ndhist(x/100,pz,'log');
% colorbar
 %colormap(flipud(gray))
 %caxis([0 4])
%  hold on
 % size(squeeze(bx))
%   xlim([0 50]);
%   xlabel('x/d_i');
%    ylabel('p_z/m_ic, 0.025 n_i/n_0');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)

   meanNNe=mean(NNe,3);
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%

%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )


%  saveas(h,strcat(datadir,'movie2.pngs/','xpz', num2str(ts(k),'%07d'),'_4.png'));
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
 %  xlim([-1 1]);
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2);
%   title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%  saveas(h,strcat(datadir,'movie2.pngs/','pxpz', num2str(ts(k),'%07d'),'_4.png'));
%  delete(h);   

    ElectronEnergy=((px).^2+py.^2+pz.^2)./2;%*4.04183354;
%    ElectronEnergy0=((px0+0.0375).^2+py0.^2+pz0.^2)./2;
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

%    xind2=6400+[0:2500];
    
%    [maxvaldens,mind] = max(meanNNe(xind2)/n);
    
%    xind3=6400+[0:mind];
    
%    xshock=mind/50;

%    xdown=xshock-1;

%    xup=xshock+1;

%    x=x/100;



%    xind2=12800+[0:8000];

%    [maxvaldens,mind] = max(meanNNe(xind2)/n);

%    xind3=12800+[0:mind];

%    xshock=mind/50;

%    xdown=xshock-1

%    densdiff=diff(squeeze(meanNNe)/n);

%    [minvaldiff,minddiff] = min(densdiff(xind2));


%    xfoot=minddiff/50;

%    xs(1)=[];

%    ffff=plot(xs,densdiff);   

%    [minvalup,mindup]=min(abs(squeeze(meanNNe(xind2))/n-1.05.*ones(length(squeeze(xind2)),1)));

%    xup=mindup/50+1

%    zindup=[-50:49];

%    plot(mind/50.*ones(length(zindup),1),zindup,'r--');
%    plot((minddiff)/50.*ones(length(zindup),1),zindup,'r--');
%    plot((mindup)/50.*ones(length(zindup),1),zindup,'r--');

%  ElectronEnergy0=ElectronEnergy;
%   indices0 = find(x0<0);
%   ElectronEnergy0(indices0) = [];


    
%    indices = find(x>xdown | x<0);
%    ElectronEnergyDown=ElectronEnergy;
%    indices1 = find(x>xshock | x<0);
%    ElectronEnergyDown(indices1) = [];
    
%    ElectronEnergyShock=ElectronEnergy;
%    indices2 = find(x>xfoot | x<xshock);
%    ElectronEnergyShock(indices2) = [];


%    ElectronEnergyFoot=ElectronEnergy;
%    indices3 = find(x>xup | x<xfoot);
%    ElectronEnergyFoot(indices3) = [];


    ElectronEnergyUp=ElectronEnergy;
    indices4 = find(z>200 | z<0);
    ElectronEnergyUp(indices4) = [];
    tag(indices4)=[];
 
%   size(ElectronEnergyUp)   
   
    ElectronEnergyUpPiston=ElectronEnergyUp;
    indices1 = find(tag==0);
    ElectronEnergyUpPiston(indices1) = [];

    ElectronEnergyUpUp=ElectronEnergyUp;
    indices2 = find(tag==1);
    ElectronEnergyUpUp(indices2) = [];








%    PhotonEnergy(indices) = [];



   % meandensdown=mean(meanNNe(xind3))/n;
    
   % meandensdownarray=meandensdown.*ones(length(xind3),1);
    
    %size(transpose(meandensdownarray))
    
    %size(xs(xind3))
    
    %zzzzz=plot(xs(xind3), transpose(meandensdownarray), 'r--');



    
    nbins=200;
   
    EnMin=0; %(min(ElectronEnergyDown)); 
    
    EnMax=2; %(max(ElectronEnergyDown)); %
    
%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold o
    en = linspace(EnMin,EnMax,200);
    hgrm=histc(ElectronEnergyUp, en);
    
    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
    plot(en,hgrm,'LineWidth',3)
    hold on

    hgrm=histc(ElectronEnergyUpUp, en);

%    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
    plot(en,hgrm,'LineWidth',3)

    hgrm=histc(ElectronEnergyUpPiston, en);

%    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
    plot(en,hgrm,'LineWidth',3)


%    hgrm2=histc(ElectronEnergy0, en);
%    plot(en,hgrm2)
   
%    plot(en,2*maxvalspec*exp(-en/(TTe*35*511))+0.001*maxvalspec*exp(-en/(TTe*150*511)),'r--')
%    plot(en,2*maxvalspec*exp(-en/(TTe*1)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*3)),'b--')
%    plot(en,maxvalspec*exp(-en/(TTe*4)),'c--')

%    legend0 = sprintf('Electron spectrum');
%    legend1 = sprintf('2* 35T_0 + 1e-3* 150T_0');
%    legend2 = sprintf('2 T_0');
%    legend3 = sprintf('3 T_0');
%    legend4 = sprintf('4 T_0');
%    legend({legend0,legend1, legend2, legend3, legend4});



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([0 EnMax]);
%    ylim([1 1e9]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumDown_',num2str(ts(k)),'.png'),'png');
%    delete(h);

%return

%    nbins=200;

%    EnMin=0; %(min(ElectronEnergyShock));

%    EnMax=2000;%(max(ElectronEnergyShock)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold o
%    en = linspace(EnMin,EnMax,200);
%    hgrm=histc(ElectronEnergyShock, en);
%    hgrm
%    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
%    plot(en,hgrm)
%    hold on

%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*100)),'b--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*300)),'c--')

%    legend0 = sprintf('Electron spectrum');
%    legend1 = sprintf('10 T_0');
%    legend2 = sprintf('30 T_0');
%    legend3 = sprintf('100 T_0');
%    legend4 = sprintf('300 T_0');
%    legend({legend0,legend1, legend2, legend3, legend4});



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([0 2000]);
%    ylim([1 1e7]);
%    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumShock_',num2str(ts(k)),'.png'),'png');
%    delete(h);


%    nbins=200;

%    EnMin=0; %(min(ElectronEnergyFoot));

%    EnMax=2000;%(max(ElectronEnergyFoot)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold o
%    en = linspace(EnMin,EnMax,200);
%    hgrm=histc(ElectronEnergyFoot, en);
%    hgrm
%    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
%    plot(en,hgrm)
%    hold on

%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*100)),'b--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*300)),'c--')

%    legend0 = sprintf('Electron spectrum');
%    legend1 = sprintf('10 T_0');
%    legend2 = sprintf('30 T_0');
%    legend3 = sprintf('100 T_0');
%    legend4 = sprintf('300 T_0');
%    legend({legend0,legend1, legend2, legend3, legend4});



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([0 2000]);
%    ylim([1 1e7]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumFoot_',num2str(ts(k)),'.png'),'png');
%    delete(h);



%    nbins=200;

%    EnMin=0; %(min(ElectronEnergyUp));

%    EnMax=100; %(max(ElectronEnergyUp)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold o
%    en = linspace(EnMin,EnMax,200);
%    hgrm=histc(ElectronEnergyUp, en);
%    hgrm
%    maxvalspec=max(hgrm);
%    en = linspace(EnMin,EnMax,200);
%   en1=en;
%    en1(1)=[];
%    size(en1)
%    size(hgrm)
%    plot(en,hgrm)
%    hold on

%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*5)),'b--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*2)),'c--')

%    legend0 = sprintf('Electron spectrum');
%    legend1 = sprintf('10 T_0');
%    legend2 = sprintf('30 T_0');
%    legend3 = sprintf('5 T_0');
%    legend4 = sprintf('2 T_0');
%    legend({legend0,legend1, legend2, legend3, legend4});



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([0 100]);
%    ylim([1 1e7]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumUp_',num2str(ts(k)),'.png'),'png');
%    delete(h);
%



%quit

%    nbins=200;

%    EnMin=(min(ElectronEnergyShock));

%    EnMax=(max(ElectronEnergyShock)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold on
%    hgrm=histogram(ElectronEnergyShock, nbins);
%    maxvalspec=max(hgrm.Values);
%    en = linspace(EnMin,EnMax,200);
%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*100)),'b--')



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([EnMin EnMax]);
%    ylim([1 1e7]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumShock_',num2str(ts(k)),'.png'),'png');
%    delete(h);




%    nbins=200;

%    EnMin=(min(ElectronEnergyFoot));

%    EnMax=(max(ElectronEnergyFoot)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold on
%    hgrm=histogram(ElectronEnergyFoot, nbins);
%    maxvalspec=max(hgrm.Values);
%    en = linspace(EnMin,EnMax,200);
%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*100)),'b--')




%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([EnMin EnMax]);
%    ylim([1 1e7]);

    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumFoot_',num2str(ts(k)),'.png'),'png');
%    delete(h);





%    nbins=200;

%    EnMin=(min(ElectronEnergyUp));

%    EnMax=(max(ElectronEnergyUp)); %

%    size(PhotonEnergyWeights);
%    h=figure;
%    histogram(ElectronEnergy0,nbins,'FaceColor','red');
%    hold on
%    hgrm=histogram(ElectronEnergyUp, nbins);
%    maxvalspec=max(hgrm.Values);
%    en = linspace(EnMin,EnMax,200);
%    plot(en,maxvalspec*exp(-en/(TTe*511*10)),'r--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*30)),'g--')
%    plot(en,maxvalspec*exp(-en/(TTe*511*100)),'b--')
    



%    title(sprintf('Simulation time, wci*t = %.3f', ts(k) *(dt * BB0/MMi)) )
%    xlabel('Particle energy, keV');
%    ylabel('Particle count');
%    xlim([EnMin EnMax]);
%    ylim([1 1e7]);
    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,particle,'EnergySpectrumUp_',num2str(ts(k)),'.png'),'png');
%    delete(h);



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
%    %set(gca, 'XScale', 'log');
%    set(gca, 'YScale', 'log');
%    saveas(h,strcat(datadir,'ParticleEnergySpectrum_',num2str(ts(k)),'.png'),'png');
%    delete(h); 
  

end

%plot(en,2*maxvalspec*exp(-en/(TTe*35*4.04183354))+1e-3*maxvalspec*exp(-en/(TTe*150*4.04183354)),'r--','LineWidth',3)
%hold on

legend0 = 'total';
legend1 = 'up';
legend2 = 'piston';
%legend3 = 't=8.2 \Omega_i^{-1}';
%legend4 = 'Fit, t=8.2 \Omega_i^{-1}';% 1.5* 14 eV + 5e-3* 889 eV');
%legend2 = sprintf('2 T_0');
%legend3 = sprintf('3 T_0');
%legend4 = sprintf('4 T_0');
legend({legend0,legend1, legend2});

%text(1,1e6,'T_e \sim 283 eV','FontSize',20)
%text(6.5,1e2,'T_e \sim 1213 eV','FontSize',20)
set(gca,'FontSize',20)
%title('Downstream electron energy spectrum','FontSize',16)
xlabel('Electron energy');
ylabel('Electron count');
xlim([0 EnMax]);
%ylim([1 2e6]);
   %set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
saveas(h,strcat(datadir,particle,'EnergySpectrumUp','.png'),'png');
delete(h);


quit
    %system(sprintf('/usr/local/bin/mencoder mf://%s/%s.pngs/*.png -mf fps=12:type=png -vf scale=960:720 -ovc x264  -x264encopts  qp=0  -o %s/%s.avi', datadir, outfile, datadir, outfile))

    %system('osascript -e ''tell application "QuickTime Player" to close every document whose name contains "%s.mov" '' ' \ {outfile});

    %system('rm %s/%s.mov' \ {datadir, outfile});
    %img2qt('%s/%s.mov' \ {datadir, outfile}, '1', '12', '%s/%s.pngs/*.png' \ {datadir, outfile});

    %system('open %s/%s.mov' \ {datadir, outfile} );    
