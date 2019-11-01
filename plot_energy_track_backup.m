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

%indices1 = find(x<60);
%pefast(indices1) = [];
%x(indices1) = [];
%idLfast(indices1)=[];
%idUfast(indices1)=[];

%indices1 = find(x>100);
%pefast(indices1) = [];
%x(indices1) = [];
%idLfast(indices1)=[];
%idUfast(indices1)=[];

%indices2 = find(pefast<1);
%pefast(indices2) = [];
%x(indices2) = [];
%idLfast(indices2)=[];
%idUfast(indices2)=[];

%idLfast=idLfast(1:5);
%size(idLfast)
%idLfast=idLfast(1:20:4300);

%idLfast=idLfast(1:100:35000);

%size(pefast)

%size(idLfast)


pearray=[idLfast; pefast];

zarray=[idLfast; z];

pearray=reshape(pearray,2,size(idLfast,1));

zarray=reshape(zarray,2,size(idLfast,1));


%size(pearray)
%size(zarray)

pearray=sortrows(pearray.',1).';

zarray=sortrows(zarray.',1).';

%zarray

%quit


idlarray=idLfast;

%idlarray=idLfast;%[63780156];
%idlarray=[25818086,49066998,62133100,117813880,179237594,183448524,186594304,206074428,224432980,224998358,236467028,239154406,251203926,255325122,266163362,320967560,334951288];
%[88285013,18881894,71621768,64284842,101333396,4946471,114788522,73341443,71547053,47966390,21424655,129288806,85862891,109775129,7354865,64284995,50497952,95803622];

for k=1:length(idlarray)
    k
%    idLfast=idlarray(k);
    
    efull=[];
    epar=[];
    eperp=[];
    muevol=[];

    xfast=[];
    times2=[];

    for p=1:length(ts)
        address=strcat(datadir, 'ele','_',num2str(ts(p),'%07d'),'_0.h5');
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        z=h5read(address,'/z')/sqrt(MMi/n);
        idL=h5read(address,'/idL');
        indices3=find(idL==idlarray(k));
        %for i=1:length(idLfast)
        %    indices3=[indices3,find(idL==idLfast(i))];
        %end
        if (~isempty(px(indices3)))
            address=strcat(datadir, 'psc','_',num2str(ts(p),'%07d'),'.h5');
            bx=h5read(address,'/bx')/BB0;
            by=h5read(address,'/by')/BB0;
            bz=h5read(address,'/bz')/BB0;
            zs = h5read(address,'/zs')/ sqrt(MMi/n);
            dt=h5read(address,'/dt');

            bxx=interp1(squeeze(zs),squeeze(mean(bx,3)),z(indices3));
            byy=interp1(squeeze(zs),squeeze(mean(by,3)),z(indices3));
            bzz=interp1(squeeze(zs),squeeze(mean(bz,3)),z(indices3));
            ppar=(px(indices3)*bxx+py(indices3)*byy+pz(indices3)*bzz)/sqrt(bxx.^2+byy.^2+bzz.^2);
            efull=[efull; px(indices3).^2+py(indices3).^2+pz(indices3).^2];
            epar=[epar; ppar*ppar];
            eperp=[eperp; px(indices3).^2+py(indices3).^2+pz(indices3).^2-ppar*ppar];
            muevol=[muevol; (px(indices3).^2+py(indices3).^2+pz(indices3).^2-ppar*ppar)/sqrt(bxx.^2+byy.^2+bzz.^2)];
            times2=[times2; dt*BB0/MMi*ts(p)];
        end
    end


            FIG=1

        figure(FIG)
         close(FIG)
         figure(FIG)
        clf

        set(FIG, 'PaperPosition', [0.5 2.5 3 5])
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

    plot(log10(efull/efull(1)),times2,'Color','black')
    hold on
    plot(log10(eperp/efull(1)),times2,':','Color','red')
    plot(log10(epar/efull(1)),times2,':','Color','blue')

    legend('E_{k,total}','E_{k,\perp}', 'E_{k,\mid \mid}','Location','best')
    legend boxoff
    xlabel('log_{10} E_k/E_{k,0}')

%        ylabel('\Omega_i t','FontSize',20)
        ylim([0 8])
        xlim([-1.0 3.0])
        %title('Evolution of log10(n_e/n_0) averaged along z'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )

        view(2);
        saveas(gcf,strcat(datadir,'tracks/','energy_track_',num2str(idlarray(k)),'.png'));
        
        
        
            FIG=1

        figure(FIG)
         close(FIG)
         figure(FIG)
        clf

        set(FIG, 'PaperPosition', [0.5 2.5 3 5])
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

    plot(log10(muevol/muevol(1)),times2,'Color','black')
    hold on

        xlabel('log_{10} \mu/\mu_0')

        %ylabel('\Omega_i t','FontSize',20)
        ylim([0 8])
        xlim([-1.0 3.0])
        %title('Evolution of log10(n_e/n_0) averaged along z'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )

        view(2);
        saveas(gcf,strcat(datadir,'tracks/','mu_',num2str(idlarray(k)),'.png'));
 
       
end


quit
