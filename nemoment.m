
wpos = [%[0.25  0.60  0.6 0.35];
        %[0.25  0.11  0.6 0.35];
        [0.09 0.60  0.34 0.30];
        [0.09 0.14  0.34 0.30];
        [0.57 0.60  0.34 0.30];
        [0.57 0.14  0.34 0.30];
 %       [0.7 0.55  0.28 0.25];
 %       [0.7 0.15  0.28 0.25];        
    ]
wM = 2
wN = 2


%%

export_cmd = {'-dpng',  '-r200' };


%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll04/';

%datadir = '/Users/klezhnin/Desktop/shock/sixty/ma30beta04/';
datadir = './h5_saved/'
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
TTe = 0.001;
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

xshock = [0.0, 32.0, 67.0, 104.0, 140.0];
xdown = [0.0, 42.0, 77.0, 114.0, 150.0];



h=figure;

set(h, 'DefaultAxesFontSize', 24)
set(h, 'DefaultTextFontSize', 24)

set(h, 'PaperUnits', 'inches');
set(h, 'PaperSize', [8 4]);

    address=strcat(datadir,particle, '_',num2str(0,'%07d'),'_0.h5')
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/sqrt(MMi/n);
        tag=h5read(address,'/tag');
        address=strcat(datadir, 'psc_',num2str(0,'%07d'),'.h5');

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

    ElectronEnergy=((px).^2+py.^2+pz.^2)./2*11.3676568;

    ElectronEnergyUp=ElectronEnergy;
    indices4 = find(z<xshock(1));
    ElectronEnergyUp(indices4) = [];


    ElectronEnergyPiston=ElectronEnergy;
    indices5 = find(tag==0);
    z(indices5)=[];
    ElectronEnergyPiston(indices5)=[];
    indices6 = find(z<xshock(1));
    ElectronEnergyPiston(indices6)=[];

    nbins=200;

    EnMin=0; %(min(ElectronEnergyDown)); 

    EnMax= 11.3676568; %(max(ElectronEnergyDown)); %

    enlin = linspace(EnMin,EnMax,1000);
%    en = logspace(-6.0,1.1,200);
    hgrmlin0=histc(ElectronEnergyUp, enlin);



for k=1:length(ts)
  k
    address=strcat(datadir,particle, '_',num2str(ts(k),'%07d'),'_0.h5')

        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/sqrt(MMi/n);
        tag=h5read(address,'/tag');
 


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
 
size(x);
size(px);

   meanNNe=mean(NNe,3);
   meanbx=mean(bx,3);
   meanby=mean(by,3);
   meanbz=mean(bz,3);
   bbb=sqrt(meanbx.^2+meanby.^2+meanbz.^2)/BB0;
   size(bbb)
   [Mb,Ib] = max(bbb(6800:10000));
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
   zmax=(1800+Ib)/5000*335;
   zmax=zmax;
   zmax





 
    ElectronEnergy=((px).^2+py.^2+pz.^2)./2*11.3676568;



    ElectronEnergyUp=ElectronEnergy;
    indices4 = find(z<zmax);
    ElectronEnergyUp(indices4) = [];


    ElectronEnergyPiston=ElectronEnergy;
    indices5 = find(tag==0);
    z(indices5)=[];
    ElectronEnergyPiston(indices5)=[];
    indices6 = find(z<zmax);
    ElectronEnergyPiston(indices6)=[];

  
    nbins=200;
   
    EnMin=0; %(min(ElectronEnergyDown)); 
    
    EnMax= 11.3676568; %(max(ElectronEnergyDown)); %
    enlin = linspace(EnMin,EnMax,1000);
%    en = logspace(-6.0,1.1,200);
    hgrmlin=histc(ElectronEnergyUp, enlin);

    ennth = linspace(0.2,EnMax,1000);
    hgrmnth=histc(ElectronEnergyUp, ennth);

    sum(hgrmnth)/sum(hgrmlin)   


end



quit
