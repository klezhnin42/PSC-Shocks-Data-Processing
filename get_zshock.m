
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
TTe = 0.002;
n= 0.05;
LL0 = 40; %sqrt(MMi/(ZZ*n));
BB0 = 0.01;
V0 = BB0/sqrt(MMi*n);

tstart = 0;
tstep = 1000;
ts = [tstart:tstep:200000];

%ts=[0, 50000, 100000, 150000, 200000];

particle='ele';

%V0 = V0 * sqrt(TTe/MMi);

VA = BB0 / sqrt(MMi*n);
Cs = sqrt(TTe/MMi);

cs1 = 'k';
cs2 = 'k';  % or k--

do_boundscheck = 0;

zshock=[];

for k=1:length(ts)
    
%  clf
 k;
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

    meanvz=mean(NVze./NNe,3);


    meanbx=mean(bx,3);
    meanby=mean(by,3);
    meanbz=mean(bz,3);
    bbb=sqrt(meanbx.^2+meanby.^2+meanbz.^2)/BB0;
    size(bbb);
    [Mb,Ib] = max(meanby(5100:10000));
%   yyy=plot(xs, 0.025*squeeze(meanNNi)/n,'black');%
    zmax=(100+Ib)/5000*335;
    if k>10
        if zmax>=zshock(length(zshock))
            zshock=[zshock; zmax];
        else
            zshock=[zshock; zshock(length(zshock))];
        end
    else
        zshock=[zshock; zmax];
    end
end
     zshock(1:6)=zshock(7);
     save('zshock.mat','zshock','-v7.3')

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

plot(ts,zshock)
hold on
        plot(6000,zshock(7),'ro')
        zshock(7);
        ylim([0 150])
        xlim([0 200000])
        %set(gca, 'YScale', 'log')
        %title('Evolution of log10(n_e/n_0) averaged along z'); %sprintf('wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )
        %title(strcat('\Omega_{i,up}t =',num2str(ts(p) *(dt * BB0/MMi),'%.1f')))
        view(2);
        saveas(gcf,strcat(datadir,'zshock_t','.png'));

size(zshock)

quit  
