


%datadir = '/Volumes/Elements/PSC_DATA/try_collisions2/coll1/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll1/';
%datadir = '/Users/klezhnin/Desktop/shock/flatfoil/hydrogen/2d/try_coll/trycoll300/files/';
datadir = './h5_saved/'
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
sizze=0;
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


    
for k=1:length(ts)
%  try
        address=strcat(datadir, 'psc_',num2str(ts(k),'%07d'),'.h5');

        NNe=h5read(address,'/NNe');
        NNi=h5read(address,'/NNi');
        dx=h5read(address,'/dx');
        dz=h5read(address,'/dz');
        dt=h5read(address,'/dt');
       
  %      h5disp(address);
        
        NVxe=h5read(address,'/NVxe');
        NVye=h5read(address,'/NVye');
        NVze=h5read(address,'/NVze');
        

        NVxi=h5read(address,'/NVxi');
        NVyi=h5read(address,'/NVyi');
        NVzi=h5read(address,'/NVzi');


        Sxxe=h5read(address,'/Sxxe');
        Syye=h5read(address,'/Syye');
        Szze=h5read(address,'/Szze');


        Sxxi=h5read(address,'/Sxxi');
        Syyi=h5read(address,'/Syyi');
        Szzi=h5read(address,'/Szzi');


        xs = h5read(address,'/xs')/ sqrt(MMi/n);
        zs = h5read(address,'/zs')/ sqrt(MMi/n);

        bx = h5read(address,'/bx');
        by = h5read(address,'/by');
        bz = h5read(address,'/bz');
  
        jy = h5read(address,'/jy');
 
        ey = h5read(address,'/ey');
        ez = h5read(address,'/ez');
        FIG=1

    figure(FIG)
     close(FIG)
     figure(FIG)
    clf

    set(FIG, 'PaperPosition', [0.2 0.2 8 5])
    set(FIG, 'DefaultAxesFontSize', 12)
    set(FIG, 'DefaultTextFontSize', 12)
    set(FIG, 'DefaultLineMarkerSize', 4)
    set(FIG, 'DefaultLineLineWidth', 1);
    set(FIG, 'renderer', 'painters');

    mkdir([datadir, '/mom_balance.pngs/'])
    exportFilePattern = [datadir, '/mom_balance.pngs/eymb_%07d.png'];
    exportCmd = '-dpng';
    exportOpt = '-r300';      
   
        TTxxe = (Sxxe - NVxe.^2./NNe) ./ NNe;
        TTyye = (Syye - NVye.^2./NNe) ./ NNe;
        TTzze = (Szze - NVze.^2./NNe) ./ NNe;
  
  Te = .333 * (TTxxe + TTyye + TTzze);
  Te (Te < 0) = 1e-6;
    
        TTxxi = (Sxxi - MMi* NVxi.^2./NNi) ./ NNi;
        TTyyi = (Syyi - MMi* NVyi.^2./NNi) ./ NNi;
        TTzzi = (Szzi - MMi* NVzi.^2./NNi) ./ NNi;

  Ti = .333 * (TTxxi + TTyyi + TTzzi);
  Ti (Ti < 0) = 1e-6;


    zindup=[-50:100];       
    meanby=mean(by,3);
    meanbx=mean(bx,3);
    meanbz=mean(bz,3);
    meanez=mean(ez,3);
    size(meanby);
    size(xs);
    meanTe=mean(Te,3);
    meanTi=mean(Ti,3);
  %  xs
    meanbz;
    meanNNe=mean(NNe,3);
    meanVze=mean(NVze./NNe,3);

%    xxx=plot(xs, squeeze(meanbx)/BB0);
%    hold on

%    zsdf=plot(xs, squeeze(meanby)/BB0);
    xx=plot(zs, squeeze(meanNNe)/n,'r');
    hold on
    yy=plot(zs, squeeze(meanVze)/(V0),'cyan');
    xxx=plot(zs, squeeze(meanby)/BB0);

    ylim([-1 20])
    xlim([0 200])
    xlabel('z/d_{i,up}') 
    legend({'n/n_{up}','V/V_{A,up}',' B_y/B_{up}'},'Location','northwest','FontSize',8)
    legend boxoff   
    title(strcat('\Omega_i t = ',sprintf('%.2f', ts(k) *(dt * BB0/MMi))));


    view(2);
    saveas(gcf,strcat(datadir,'1devolgif_', num2str(ts(k),'%07d'),'_6.png'));
%  end 
  
end;

quit
