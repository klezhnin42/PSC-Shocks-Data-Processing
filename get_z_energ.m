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

tstart = 0;
tstep = 1000;
ts = [tstart:tstep:200000];

%save('reformation.mat','rfrmtn','-v7.3')


address=strcat(datadir, 'ele','_',num2str(200000,'%07d'),'_0.h5');
 %   try
        %dfield = load(address,'-mat'); 

%        h5disp(address);

px=h5read(address,'/px');
py=h5read(address,'/py');
pz=h5read(address,'/pz');
z=h5read(address,'/z')/sqrt(MMi/n);
idLfast=h5read(address,'/idL');
idUfast=h5read(address,'/idU');
tag=h5read(address,'/tag');

petot=sqrt(px.^2+py.^2+pz.^2);

pefast=petot;

indices1 = find(z<153);
pefast(indices1) = [];
z(indices1) = [];
idLfast(indices1)=[];
idUfast(indices1)=[];
tag(indices1)=[];

indices1 = find(z>200);
pefast(indices1) = [];
z(indices1) = [];
idLfast(indices1)=[];
idUfast(indices1)=[];
tag(indices1)=[];

indices2 = find(pefast<0.3);
pefast(indices2) = [];
z(indices2) = [];
idLfast(indices2)=[];
idUfast(indices2)=[];
tag(indices2)=[];


indices2 = find(tag==1);
pefast(indices2) = [];
z(indices2) = [];
idLfast(indices2)=[];
idUfast(indices2)=[];
tag(indices2)=[];


%pearray=[idLfast; pefast];

%zarray=[idLfast; z];

%pearray=reshape(pearray,2,size(idLfast,1));

%zarray=reshape(zarray,2,size(idLfast,1));


%size(pearray)
%size(zarray)

%pearray=sortrows(pearray.',1).';

%zarray=sortrows(zarray.',1).';

%zarray

%quit




idlarray=idLfast;

size(idlarray)

pearray=[idlarray];
pearray=sortrows(pearray.',1).';

zarray=[idlarray];
zarray=sortrows(zarray.',1).';



for p=1:length(ts)
        address=strcat(datadir, 'ele','_',num2str(ts(p),'%07d'),'_0.h5');
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        z=h5read(address,'/z')/sqrt(MMi/n);
        idL=h5read(address,'/idL');
        idU=h5read(address,'/idU');
        tag=h5read(address,'/tag');
        petot=sqrt(px.^2+py.^2+pz.^2);
 
        indices2 = find(tag==1);
        petot(indices2) = [];
        z(indices2) = [];
        idL(indices2)=[];
        idU(indices2)=[];
        tag(indices2)=[];

        Lia=ismember(idL,idlarray);
        indices3 = find(Lia==1);
        petot=petot(indices3);
        z=z(indices3);
        idL=idL(indices3);
        idU=idU(indices3);
        tag=tag(indices3);        
              

        pebuff=[idL; petot];
        pebuff=reshape(pebuff,2,size(idL,1));

        size(pebuff)

        zbuff=[idL; z];
        zbuff=reshape(zbuff,2,size(idL,1));
        
        size(zbuff)



        address=strcat(datadir, 'psc_',num2str(ts(p),'%07d'),'.h5');
        NNi=h5read(address,'/NNi');
        dt=h5read(address,'/dt');
        zs = h5read(address,'/zs')/ sqrt(MMi/n);
        by = h5read(address,'/by');




            FIG=1

        figure(FIG)
         close(FIG)
         figure(FIG)
        clf

        set(FIG, 'PaperPosition', [0.5 0.5 6 3])
        set(FIG, 'DefaultAxesFontSize', 14)
        set(FIG, 'DefaultTextFontSize', 14)
        set(FIG, 'DefaultLineMarkerSize', 4)
        set(FIG, 'DefaultLineLineWidth', 1);
        set(FIG, 'renderer', 'painters');

        mkdir([datadir, '/mom_balance.pngs/'])
        exportFilePattern = [datadir, '/mom_balance.pngs/eymb_%07d.png'];
        exportCmd = '-dpng';
        exportOpt = '-r300';



    %size(reshape(rfrmtn,[13200,12800]))
    %size(times)w

    %xxx=imagesc(xs,times, reshape(rfrmtn.rfrmtn,[75800,20200])');
    %colorbar
    %caxis([-0.2 1.2])
    %colormap(jet)
    %hold on

    xedge = 300 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

    yedge = 2.0 * [0:0.005:1];


    hgrm=histc2d(z,petot,xedge,yedge);
    loghgrm=log10(hgrm);
    loghgrm(loghgrm==-Inf)=NaN;
%   loghgrm(~isfinite(loghgrm))=0;
    size(xedge)
    size(yedge)
    size(loghgrm)
    f=pcolor(xedge,yedge,loghgrm');
    hold on
    colorbar
    caxis([0 2])
    set(gca,'YDir','normal')
    set(f,'EdgeColor', 'none')
    set(f,'HandleVisibility','off')
    colormap(jet)
    shading flat
%   hold on
%   size(squeeze(bx))w5
%   xlim([0 200]);
%   ylim([-1.5 1.5]);



%    plot(z,petot,'ko','MarkerSize',3)
%    hold on

    meanNNi=mean(NNi,3);
    size(meanNNi)
    meanBy=mean(by,3);

    yyy=plot(zs, 0.1*squeeze(meanNNi)/n,'red','LineWidth',3);%
    zzz=plot(zs, 0.1*squeeze(meanBy)/BB0,'green','LineWidth',1);%
%    plot(log10(eperp/efull(1)),times2,':','Color','red')
%    plot(log10(epar/efull(1)),times2,':','Color','blue')

%   legend('E_{k,total}','E_{k,\perp}', 'E_{k,\mid \mid}','Location','best')
%    legend boxoff
%    xlabel('log_{10} E_k/E_{k,0}')

%        ylabel('\Omega_i t','FontSize',20)
        ylim([0 2])
        xlim([0.0 300.0])
        %title('Evolution of log10(n_e/n_0) averaged along z'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )

        view(2);
        saveas(gcf,strcat(datadir,'movie2.pngs/','z_pe_',num2str(ts(p)),'.png'));









end


quit
