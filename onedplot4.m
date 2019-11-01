


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

%tstart = 40000;
%tstep = 2000;
ts = 200000;% = [tstart:tstep:60000];

 
    
for k=0:sizze
%  try
        address=strcat(datadir, 'psc_',num2str(ts+1000*k,'%07d'),'.h5');

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

    set(FIG, 'PaperPosition', [0.2 0.2 5 8])
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

%  TTxxi2=mean(TTxxi,1);
%  TTyyi2=mean(TTyyi,1);
%  TTzzi2=mean(TTzzi,1);

%  mean(TTxxi2(10500:13000))
%  mean(TTyyi2(10500:13000))
%  mean(TTzzi2(10500:13000))

    Nu = eta0 * BB0 .* NNe .* ((TTe./Te).^(1.5));
    Xi = bx ./ Nu;
    D = Xi.^4 + 14.79*Xi.^2 + 3.77;
    F1 = 1 - (6.42*Xi.^2 + 1.84) ./ D;
    Eta_perp = Nu ./ NNe;
    etaJ=NNe .* Eta_perp .* F1 .* jy;
    zindup=[-50:100];       
    meanby=mean(by,3);
    meanbx=mean(bx,3);
    meanbz=mean(bz,3);
    meanez=mean(ez,3);
    size(meanby)
    size(xs)
    meanTe=mean(Te,3);
    meanTi=mean(Ti,3);
  %  xs
    meanbz
    meanNNe=mean(NNe,3);
    meanVze=mean(NVze./NNe,3);

%    xxx=plot(xs, squeeze(meanbx)/BB0);
%    hold on

%    zsdf=plot(xs, squeeze(meanby)/BB0);
    subplot(4,1,1)
    xx=plot(zs, squeeze(meanNNe)/n,'r');
    hold on
    yy=plot(zs, squeeze(meanVze)/(V0),'cyan');
    xxx=plot(zs, squeeze(meanby)/BB0);

    fffff=plot((1567.0/14.91).*ones(length(zindup)),zindup,'k--','LineWidth',1)
    ggg=plot((2089.0/14.91).*ones(length(zindup)),zindup,'k--','LineWidth',1)
    gggg=plot((153.0).*ones(length(zindup)),zindup,'k--','LineWidth',1)
%    xxx=plot(zs, squeeze(meanby)/BB0);
    ylim([-1 20])
    xlim([0 200]) 
    legend({'n/n_{up}','V/V_{A,up}',' B_y/B_{up}'},'Location','northwest','FontSize',8)
    legend boxoff   

    text(30,5,'Piston','Rotation',0,'Color','black','FontSize',18)
    text(115,4,'Downstream','Rotation',70,'Color','black','FontSize',11)
    text(155,3,'Upstream','Color','black','FontSize',12 )
    text(146,4,'Shock layer','Rotation',90,'Color','black','FontSize',11)
    text(-15,-4,'(a)','FontSize',10)
    text(-15,-32,'(b)','FontSize',10)
    text(-15,-60,'(c)','FontSize',10)
    text(-15,-88,'(d)','FontSize',10)
    text(-15,-115,'(e)','FontSize',10)
    %set(gca,'XTick',[]);

    subplot(4,1,2)
    address=strcat(datadir, 'ion','_',num2str(ts,'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 

%        h5disp(address);
%    return   
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/ sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/ sqrt(MMi/n);
%       dt=h5read(address,'/dt');
        tag=h5read(address,'/tag');



    indices2 = find(tag==2);
%    pefast(indices2) = [];
    z(indices2) = [];
    px(indices2) = [];
    py(indices2) = [];
    pz(indices2) = [];
    
    xedge = 200 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

    yedge = 0.25 * [-1:0.005:1];


    hgrm=histc2d(z,pz,xedge,yedge);
    loghgrm=log10(hgrm);
    loghgrm(loghgrm==-Inf)=NaN;
%   loghgrm(~isfinite(loghgrm))=0;
    f=pcolor(xedge,yedge,loghgrm');
%    colorbar
    caxis([0 5])
    set(gca,'YDir','normal')
    set(f,'EdgeColor', 'none')
    set(f,'HandleVisibility','off')
    colormap(jet)
    shading flat
    
    hold on
 % size(squeeze(bx))w5
    xlim([0 200]);
    ylim([-0.15 0.2]);

   meanNNi=mean(NNi,3);
   size(meanNNi)
   yyy=plot(zs, 0.025*squeeze(meanNNi)/n,'red')%'LineWidth',2);%

   zindup=[-50:49];

   asdf=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
   fsad=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
   sdaf=plot((153).*ones(length(zindup)),zindup,'k--')

%   legend({'\color{black}0.025 n_i/n_{i,up}'},'Location','southwest','FontSize',10,'Color','black')
%   legend boxoff
%   text(30,0.15,'z-p_{iz} phase plot','Rotation',0,'Color','black','FontSize',10)
    ylabel('z-p_{iz}, 0.025n_i/n_{i,up}','FontSize',12)

 subplot(4,1,3)
    address=strcat(datadir, 'ele','_',num2str(ts,'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 

%        h5disp(address);
%    return   
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/ sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/ sqrt(MMi/n);
%       dt=h5read(address,'/dt');
        tag=h5read(address,'/tag');



    indices2 = find(tag==2);
%    pefast(indices2) = [];
    z(indices2) = [];
    px(indices2) = [];
    py(indices2) = [];
    pz(indices2) = [];

    xedge = 200 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

    yedge = 1.5 * [-1:0.005:1];


    hgrm=histc2d(z,pz,xedge,yedge);
    loghgrm=log10(hgrm);
    loghgrm(loghgrm==-Inf)=NaN;
%   loghgrm = loghgrm/sum(sum(loghgrm));
%   loghgrm(~isfinite(loghgrm))=0;
    f=pcolor(xedge,yedge,loghgrm');
%    colorbar
    caxis([0 5])
    set(gca,'YDir','normal')
    set(f,'EdgeColor', 'none')
    set(f,'HandleVisibility','off')
    colormap(jet)
    shading flat
    hold on
 % size(squeeze(bx))w5
    xlim([0 200]);
    ylim([-1.5 1.5]);
%    xlabel('z/d_{i,up}');
%   xlabel('z / d_{i,up}','FontSize',12)


   meanNNi=mean(NNi,3);
   size(meanNNi)
   yyy=plot(zs, 0.2*squeeze(meanNNi)/n,'red')%,'LineWidth',2);%

   zindup=[-50:49];

   asdf=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
   fsad=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
   sdaf=plot((153).*ones(length(zindup)),zindup,'k--')

%   legend({'\color{black}0.2 n_i/n_{i,up}'},'Location','southwest','FontSize',10,'Color','black')
%   legend boxoff
%   text(30,1.2,'z-p_{ez} phase plot','Rotation',0,'Color','black','FontSize',10)
    ylabel('z-p_{ez}, 0.2n_i/n_{i,up}','FontSize',12)

 subplot(4,1,4)
    address=strcat('../../try2/h5_saved/', 'ele','_',num2str(ts,'%07d'),'_0.h5')
 %   try
        %dfield = load(address,'-mat'); 

%        h5disp(address);
%    return   
        px=h5read(address,'/px');
        py=h5read(address,'/py');
        pz=h5read(address,'/pz');
        x=h5read(address,'/x')/ sqrt(MMi/n);
        y=h5read(address,'/y');
        z=h5read(address,'/z')/ sqrt(MMi/n);
%       dt=h5read(address,'/dt');
        tag=h5read(address,'/tag');



        address=strcat('../../try2/h5_saved/', 'psc_',num2str(ts+1000*k,'%07d'),'.h5');

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

    meanNNe=mean(NNe,3);


    indices2 = find(tag==2);
%    pefast(indices2) = [];
    z(indices2) = [];
    px(indices2) = [];
    py(indices2) = [];
    pz(indices2) = [];

    xedge = 200 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

    yedge = 1.5 * [-1:0.005:1];

    hgrm=histc2d(z,pz,xedge,yedge);
    loghgrm=log10(hgrm);
    loghgrm(loghgrm==-Inf)=NaN;
%    loghgrm = loghgrm/sum(sum(loghgrm));
%   loghgrm(~isfinite(loghgrm))=0;
    f=pcolor(xedge,yedge,loghgrm');
%    colorbar
    caxis([0.2 5.2])
    set(gca,'YDir','normal')
    set(f,'EdgeColor', 'none')
    set(f,'HandleVisibility','off')
    colormap(jet)
    shading flat
    hold on
 % size(squeeze(bx))w5
    xlim([0 200]);
    ylim([-1.5 1.5]);
%    xlabel('z/d_{i,up}');
    xlabel('z/d_{i,up}','FontSize',12)
    ylabel('z-p_{ez}, 0.2n_i/n_{i,up}','FontSize',12)

   meanNNi=mean(NNi,3);
   size(meanNNi)
   yyy=plot(zs, 0.2*squeeze(meanNNi)/n,'red')%,'LineWidth',2);%

   zindup=[-50:49];

%  asdf=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
%   fsad=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
%   sdaf=plot((2223/14.91).*ones(length(zindup)),zindup,'k--')

%   legend({'\color{white}0.2 n_i/n_{i,up}'},'Location','southwest','FontSize',10,'Color','white')
%   legend boxoff
%  text(30,1.2,'z-p_{ez} phase plot','Rotation',0,'Color','white','FontSize',10)
   text(10,0.0,'\color{white}Collisionless, T_{heating}=8 \Omega_i^{-1}','Rotation',0,'Color','white','FontSize',10)   

%    title(sprintf('Simulation time, wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )
%  title(strcat('Ablation',', \Omega_{i,up}t =',num2str((ts+1000*k) *(dt * BB0/MMi),'%.1f')))
    view(2);
    saveas(gcf,strcat(datadir,'1d_', num2str(ts+1000*k,'%07d'),'_6.png'));
%  end 
   
end;

quit

