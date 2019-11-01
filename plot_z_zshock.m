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


zarray=load('zarray_full.mat','zarray');
pearray=load('pearray_full.mat','pearray');
zshock=load('zshock.mat','zshock');

z=zarray.zarray;

pe=pearray.pearray;

zsh=zshock.zshock;



for i=1:size(z,1)

    pe1=pe(i,:);
    z1=z(i,:);
    pe1(1)=[];
    z1(1)=[];
    size(z1);
    size(zsh);
           FIG=1;

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


%    xedge = 300 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

%   yedge = 2.0 * [0:0.005:1];


%    hgrm=histc2d(z1,pe1,xedge,yedge);
%    loghgrm=log10(hgrm);
%    loghgrm(loghgrm==-Inf)=NaN;
%   loghgrm(~isfinite(loghgrm))=0;
%    size(xedge)
%    size(yedge)
%    size(loghgrm)
    f=plot(smoothdata(z1-zsh'),smoothdata(pe1));
%   min(z1-zsh')
    hold on
%    colorbar
%    caxis([1e-4 2])
%    set(gca,'YDir','normal')
%    set(f,'EdgeColor', 'none')
%    set(f,'HandleVisibility','off%')
%    colormap(jet)
%    shading flat

        ylim([1e-4 2])
        xlim([min([0.0 min(z1-zsh')]) 150.0])
        set(gca, 'YScale', 'log')
        %title('Evolution of log10(n_e/n_0) averaged along z'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )
        %title(strcat('\Omega_{i,up}t =',num2str(ts(p) *(dt * BB0/MMi),'%.1f')))
        view(2);
        saveas(gcf,strcat(datadir,'movie2.pngs/','z_zshock_pe',num2str(i),'.png'));

end





quit
