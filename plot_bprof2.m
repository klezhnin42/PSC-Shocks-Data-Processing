%datadir = '/Volumes/Elements/PSC_DATA/try_collisions2/coll1/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll1/';
datadir = './h5_saved/';
%Users/klezhnin/Desktop/shock/largemmi/1d_MMI_900/theta80/try_files/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb/nnb00075/';
%datadir = '/Volumes/Elements/PSC_DATA/try_par/nif/coll1/';

% initial parameters
MMi =100;
ZZ = 1;
TTe = 0.002;
n= 0.05;
LL0 = 40.0;%sqrt(MMi/(ZZ*n));
BB0 = 0.01; %sqrt(TTe*n);
V0 = BB0/sqrt(MMi*n);
delx=2; %step in derivatives
num=1; %number of current sheet, from 1 to inf, from left to right
sizze=200;
eta0=0.0;

rates=zeros(1,sizze+1);
di=zeros(1,sizze+1);
eta=zeros(1,sizze+1);
SLund=zeros(1,sizze+1);
deltasp=zeros(1,sizze+1);

LRC=100;

%tstart = 125000;
%tstep = 25000;
%ts = [tstart:tstep:200000];

ts = [125000,138000,150000,163000,175000,188000,200000];

        FIG=1

    figure(FIG)
     close(FIG)
     figure(FIG)
    clf

    set(FIG, 'PaperPosition', [0.5 0.5 4 6])
    set(FIG, 'DefaultAxesFontSize', 14)
    set(FIG, 'DefaultTextFontSize', 14)
    set(FIG, 'DefaultLineMarkerSize', 4)
    set(FIG, 'DefaultLineLineWidth', 1);
    set(FIG, 'renderer', 'painters');

    mkdir([datadir, '/mom_balance.pngs/'])
    exportFilePattern = [datadir, '/mom_balance.pngs/eymb_%07d.png'];
    exportCmd = '-dpng';
    exportOpt = '-r300';


for p=1:length(ts)
    address=strcat(datadir, 'psc','_',num2str(ts(p),'%07d'),'.h5');
    dt=h5read(address,'/dt');
    zs=h5read(address,'/zs')/sqrt(MMi/n);
    ne=mean(h5read(address,'/NNe')/n,3);
    plot(zs-(ts(p)-125000)*dt*BB0/MMi*18.25,ne+p+3.0,'Color','black')
    hold on
end
xlim([50 200])
ylim([5 14])
%view(2);
%saveas(gcf,strcat(datadir,'neprofiles','.png'));


%save('reformation.mat','rfrmtn','-v7.3')
%address=strcat(datadir, 'psc','_',num2str(200000,'%07d'),'.h5');

%dt=h5read(address,'/dt');
%zs=h5read(address,'/zs')/sqrt(MMi/n);


%address=strcat(datadir, 'ele','_',num2str(200000,'%07d'),'_0.h5');
 %   try
        %dfield = load(address,'-mat'); 
     
%        h5disp(address);
        
%px=h5read(address,'/px');
%py=h5read(address,'/py');
%pz=h5read(address,'/pz');
%z=h5read(address,'/z')/sqrt(MMi/n);
%idLfast=h5read(address,'/idL');
%idUfast=h5read(address,'/idU');
%tag=h5read(address,'/tag');
%petot=sqrt(px.^2+py.^2+pz.^2);

%pefast=petot;

%indices1 = find(z<135);
%pefast(indices1) = [];
%z(indices1) = [];
%idLfast(indices1)=[];
%idUfast(indices1)=[];
%tag(indices1)=[];

%indices1 = find(z>145);
%pefast(indices1) = [];
%z(indices1) = [];
%idLfast(indices1)=[];
%idUfast(indices1)=[];
%tag(indices1)=[];

%indices2 = find(pefast<0.7);
%pefast(indices2) = [];
%z(indices2) = [];
%idLfast(indices2)=[];
%idUfast(indices2)=[];
%tag(indices2)=[];


%indices2 = find(tag==1);
%pefast(indices2) = [];
%z(indices2) = [];
%idLfast(indices2)=[];
%idUfast(indices2)=[];
%tag(indices2)=[];



%idLfast=idLfast(1:5);
%size(idLfast)
%idLfast=idLfast(1:20:4300);

%idLfast=idLfast(1:100:35000);



%idLfast=idLfast(1:20:4300);

idLfast=[63780156];


%size(idLfast)
%tstart = 0;
%tstep = 1000;
%ts = [tstart:tstep:200000];

%xfast=[];
%times2=[];

%for i=1:length(idLfast)
%    i
%    xfast=[];
%    times2=[];
%for p=1:length(ts)
%    address=strcat(datadir, 'ele','_',num2str(ts(p),'%07d'),'_0.h5');
%    x=h5read(address,'/z')/sqrt(MMi/n);
%    idL=h5read(address,'/idL');
%    idU=h5read(address,'/idU');
    %indices3=[];
%    indices3=find(idL==idLfast(i));
    %for i=1:length(idLfast)
    %    indices3=[indices3,find(idL==idLfast(i))];
    %end
%    if (~isempty(x(indices3)))
%        xfast=[xfast; x(indices3)];
%        times2=[times2; dt*BB0/MMi*ts(p)];
%    end
%end

%size(xfast, 1)
%size(times2,1)

%return


%times2=dt*BB0/MMi*[0:1000:300000];

%rfrmtn=load('zt_etot.mat','rfrmtn');

%save('xfast.mat','xfast','-v7.3')
%save('tfast.mat','times2','-v7.3')

xfast=load('xfast.mat','xfast');
times2=load('tfast.mat','times2');

xfast
times2


    FIG=1

%    figure(FIG)
%     close(FIG)
%     figure(FIG)
%    clf

%    set(FIG, 'PaperPosition', [0.5 2.5 6 4])
%    set(FIG, 'DefaultAxesFontSize', 14)
%    set(FIG, 'DefaultTextFontSize', 14)
%    set(FIG, 'DefaultLineMarkerSize', 4)
%    set(FIG, 'DefaultLineLineWidth', 1);
%    set(FIG, 'renderer', 'painters');

%    mkdir([datadir, '/mom_balance.pngs/'])
%    exportFilePattern = [datadir, '/mom_balance.pngs/eymb_%07d.png'];
%    exportCmd = '-dpng';
%    exportOpt = '-r300';


%times=dt*BB0/MMi*[0:10:200000];

%size(reshape(rfrmtn,[13200,12800]))
%size(times)w

%xxx=imagesc(zs,times, reshape(rfrmtn.rfrmtn,[10000,40200])');
%colorbar
%caxis([-2.5 0.0])
%colormap(jet)
%hold on

plot(xfast.xfast-(times2.times2-5.0)*18.25,times2.times2*2.0-5.0,'Color','red')

%set(gca,'YDir','normal')
    xlabel('z / d_{i, up}','FontSize',20)

    ylabel('\Omega_i t, n/n_{up}','FontSize',20)
%    ylim([0 5])
%    xlim([0 200.0])
    title('log_{10} n/n_{up}'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )

    view(2);
    saveas(gcf,strcat(datadir,'neprofiles','.png'));

%end

%idLfast

quit

%di=fliplr(di);
%eta=fwliplr(eta);
%SLund=fliplr(SLund);
%deltasp=fliplr(deltasp);
%rates=fliplr(rates);
