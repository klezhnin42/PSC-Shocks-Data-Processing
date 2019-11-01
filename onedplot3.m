


%datadir = '/Volumes/Elements/PSC_DATA/try_collisions2/coll1/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll1/';
datadir = './h5_saved/';
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

    set(FIG, 'PaperPosition', [0.3 0.5 5 10])
    set(FIG, 'DefaultAxesFontSize', 14)
    set(FIG, 'DefaultTextFontSize', 14)
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
    subplot(5,1,1)
    xx=plot(zs, squeeze(meanNNe)/n,'r');
    hold on
    yy=plot(zs, squeeze(meanVze)/(V0),'cyan');
    fffff=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
    ggg=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
    gggg=plot((2223/14.91).*ones(length(zindup)),zindup,'k--')
    ylim([0 20])
    xlim([0 200]) 
    legend({'n/n_{up}','V/V_{A,up}'},'Location','northwest','FontSize',8)
    legend boxoff   

    text(30,5,'Piston','Rotation',0,'Color','black','FontSize',18)
    text(113,3.5,'Downstream','Rotation',70,'Color','black','FontSize',12)
    text(155,3,'Upstream','Color','black','FontSize',12 )
    text(144,4,'Shock layer','Rotation',90,'Color','black','FontSize',12)


    subplot(5,1,2)
    xxx=plot(zs, squeeze(meanby)/BB0);
    hold on
    yyy=plot(zs, 10*squeeze(meanez)/BB0);
    fffff=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
    ggg=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
    gggg=plot((2223/14.91).*ones(length(zindup)),zindup,'k--')
    legend({'B_y/B_{up}','10c E_z/B_{up}'},'Location','northwest','FontSize',8)    
    legend boxoff   

    ylim([-1.1 6])
    xlim([0 200 ])

    subplot(5,1,3)  
    
%    yyy=plot(zs, squeeze(meanNNe)/n,'r');%
%    hold on

%    zzz=plot(zs, squeeze(meanVze)/(V0),'cyan');
    fff=plot(zs,squeeze(meanTe)/TTe,'g');
    hold on
    ffff=plot(zs,squeeze(meanTi)/TTe,'m');
    fffff=plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
    ggg=plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
    gggg=plot((2223/14.91).*ones(length(zindup)),zindup,'k--')
    ylim([0 70])
    legend({'T_e/T_{e,up}','T_i/T_{i,up}'},'Location','northwest','FontSize',8)
    legend boxoff
%    set(gca, 'YScale', 'log')   
%    ffff=plot(xs,squeeze(meanTi)/TTe);
%    meanTe=mean(Te,1);
  
%    zzzz=plot(xs, squeeze(meanTe)/%TTe);
%    xind =[1567:2089];  
%    xind2=5000+xind;

%    [maxvaldens,mind] = max(meanNNe(xind2)/n);

%    xind3=12800+[0:mind];

%    xshock=mind/50;

%    xfoot=xshock-1;

%    densdiff=diff(squeeze(meanNNe)/n);

%    [minvaldiff,minddiff] = min(densdiff(xind2));

%    xs(1)=[];

%    ffff=plot(xs,densdiff);   

%    [minvalup,mindup]=min(abs(squeeze(meanNNe(xind2))/n-1.2.*ones(length(squeeze(xind2)),1)));

%    mindup; 

%    zindup=[-50:100];
  
%    plot(mind/50.*ones(length(zindup),1),zindup,'r--');
%    plot((minddiff)/50.*ones(length(zindup),1),zindup,'g--');
%    plot((mindup)/50.*ones(length(zindup),1),zindup,'b--');
%   asdf=plot(xind/(14.91),mean(meanNNe(xind2))/n.*ones(length(xind2)),'r--')
%   scw=plot(xind/(14.91),0.5*mean(meanTe(xind2))/TTe.*ones(length(xind2)),'g--')
%   dfgh=plot(xind/(14.91),0.1*mean(meanTi(xind2))/TTe.*ones(length(xind2)),'m--')
   
%   plot((1567/14.91).*ones(length(zindup)),zindup,'k--')
%   plot((2089/14.91).*ones(length(zindup)),zindup,'k--')
%   plot((2223/14.91).*ones(length(zindup)),zindup,'k--')

%   text(10,15,'Piston','Rotation', 0)
%   text(110,15,'Downstream','Rotation', 80)
%   text(160,4,'Upstream')
%   text(147 ,17,'Shock layer','Rotation',90)
%  text(-20,-230,'(a)','FontSize',20)
%mean(meanNNe(xind3))/n

%mean(meanTe(xind3)/TTe)

%mean(meanTi(xind3)/TTe)
   
   % rhs  = + Vze(zind,ixX) .* bx(zind,ixX)...
   %     - Vxe(zind,ixX) .* bz(zind,ixX)...
   %     +  (Tyze(zind+delx,ixX) - Tyze(zind-delx,ixX))  / (2*delx*dfield.dz*dfield.decZ) ...
   %     +  (Txye(zind,ixX+delx) - Txye(zind,ixX-delx))  / (2*delx*dfield.dx*dfield.decX) ...
   %     +  (Pyer(zind,ixX) - Pyel(zind, ixX) ) / (dfield.dt*(4000))...
   %     -  etaJ(zind,ixX);
   % plot(zs(zind)-zX, -rhs / KK) 
    
%   set(gca,'fontsize',20,'LineWidth',2)
%    set(xxx,'LineWidth',2)
%    set(yyy,'LineWidth',2)
%    set(zzz,'LineWidth',2)
%    set(zzzz,'LineWidth',2)
%    
 
%    legend({'B_y/B_{up}','n/n_{up}','V/V_{A,up}', '0.5T_e/T_{up}', '0.1T_i/T_{up}'},'FontSize',10,'Location','northeast')
 
%    xlabel('z / d_{i,up}','FontSize',20)
 
%    ylabel('{B_x [3.0 B0], n_e [0.6n_{\rm max}], v_{ez} [0.42 C_{s}] }','FontSize',20)
%    ylim([0 30])
%   xlim([0.0-107.5 250.0-107.5])
    xlim([0 200])

    subplot(5,1,4)
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
    loghgrm(~isfinite(loghgrm))=0;
    imagesc(xedge,yedge,loghgrm');
%   colorbar
    set(gca,'YDir','normal')
 %colormap(flipud(gray))
    caxis([0 5])
    colormap(jet)
    hold on
 % size(squeeze(bx))w5
    xlim([0 200]);
    ylim([-0.2 0.2]);
%    xlabel('z/d_{i,up}');
%    xlabel('z / d_{i,up}','FontSize',20)
%    ylabel('p_{iz}/m_ic, 0.025 n_i/n_{up}');
 %  ylim([-1 1]);
  %figure('Visible','off');
  %view(2)ww

   meanNNi=mean(NNi,3);
   size(meanNNi)
   yyy=plot(zs, 0.025*squeeze(meanNNi)/n,'red')%'LineWidth',2);%

   zindup=[-50:49];

%   legend({'0.025 n_i/n_{i,up}'},'Location','northwest','FontSize',8,'Color','white')
%   legend boxoff
%   text(3,0.0,'Ion z-p_{iz}','Rotation',0,'Color','white','FontSize',18)

%  plot((1530/14.91-97.1).*ones(length(zindup)),zindup,'w--')
%   plot((1960/14.91-97.1).*ones(length(zindup)),zindup,'w--')
%   plot((2100/14.91-97.1).*ones(length(zindup)),zindup,'w--')

%   text(3,0.12,'Piston','Rotation',90,'Color','white','FontSize',18)
%   text(10,0.23,'Downstream','Color','white','FontSize',18)
%   text(47,0.04,'Upstream','Color','white','FontSize',18)
%   text(37,0.13,'Shock layer','Rotation',70,'Color','white','FontSize',18)

   asdf=plot((1567/14.91).*ones(length(zindup)),zindup,'w--')
   fsad=plot((2089/14.91).*ones(length(zindup)),zindup,'w--')
   sdaf=plot((2223/14.91).*ones(length(zindup)),zindup,'w--')

   legend({'\color{white}0.025 n_i/n_{i,up}'},'Location','southwest','FontSize',10,'Color','white')
   legend boxoff
   text(30,0.15,'z-p_{iz} phase plot','Rotation',0,'Color','white','FontSize',10)


 subplot(5,1,5)
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
    loghgrm(~isfinite(loghgrm))=0;
    imagesc(xedge,yedge,loghgrm');
%   colorbar
    set(gca,'YDir','normal')
 %colormap(flipud(gray))
    caxis([0 5])
    colormap(jet)
    hold on
 % size(squeeze(bx))w5
    xlim([0 200]);
    ylim([-1.5 1.5]);
%    xlabel('z/d_{i,up}');
    xlabel('z / d_{i,up}','FontSize',20)


   meanNNi=mean(NNi,3);
   size(meanNNi)
   yyy=plot(zs, 0.2*squeeze(meanNNi)/n,'red')%,'LineWidth',2);%

   zindup=[-50:49];

%  plot((1530/14.91-97.1).*ones(length(zindup)),zindup,'w--')
%   plot((1960/14.91-97.1).*ones(length(zindup)),zindup,'w--')
%   plot((2100/14.91-97.1).*ones(length(zindup)),zindup,'w--')

%   text(3,0.12,'Piston','Rotation',90,'Color','white','FontSize',18)
%   text(10,0.23,'Downstream','Color','white','FontSize',18)
%   text(47,0.04,'Upstream','Color','white','FontSize',18)
%   text(37,0.13,'Shock layer','Rotation',70,'Color','white','FontSize',18)

%   legend({'0.25 n_i/n_{i,up}'},'Location','northwest','FontSize',8,'Color','white')
%   legend boxoff
%   text(3,0.0,'Electron z-p_{iz}','Rotation',0,'Color','white','FontSize',18)

   asdf=plot((1567/14.91).*ones(length(zindup)),zindup,'w--')
   fsad=plot((2089/14.91).*ones(length(zindup)),zindup,'w--')
   sdaf=plot((2223/14.91).*ones(length(zindup)),zindup,'w--')

   legend({'\color{white}0.2 n_i/n_{i,up}'},'Location','southwest','FontSize',10,'Color','white')
   legend boxoff
   text(30,1.2,'z-p_{ez} phase plot','Rotation',0,'Color','white','FontSize',10)


%    title(sprintf('Simulation time, wci*t = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )
%  title(strcat('Ablation',', \Omega_{i,up}t =',num2str((ts+1000*k) *(dt * BB0/MMi),'%.1f')))
    view(2);
    saveas(gcf,strcat(datadir,'1d_', num2str(ts+1000*k,'%07d'),'_6.png'));
%  end 
   imout=RemoveWhiteSpace(gcf)
   saveas(imout,strcat(datadir,'1d_', num2str(ts+1000*k,'%07d'),'_6.png'));
end;


quit

%di=fliplr(di);
%eta=fliplr(eta);
%SLund=fliplr(SLund);
%deltasp=fliplr(deltasp);
%rates=fliplr(rates);



%FIG=1
%figure(FIG)
%close(FIG)
%figure(FIG)

%x = linspace((ts-sizze*2000)*(dfield.dt * Cs / (LL0*sqrt(MMi))), ts*(dfield.dt * Cs / (LL0*sqrt(MMi))), sizze+1);

%plot(x,di/max(di))
%hold on

%plot(x,eta/max(eta))
%hold on

%plot(x,SLund/max(SLund))
%hold on

%plot(x,deltasp/max(deltasp))
%hold on

%plot(x,rates/max(rates))
%hold on


%legend(strcat('d_i, max(d_i)=',num2str(max(di))),strcat('\eta, max(\eta)=',num2str(max(eta))),strcat('S, max(S)=',num2str(max(SLund))),strcat('\delta_{\rm SP}, max(\delta_{\rm SP})=',num2str(max(deltasp))),strcat('rate, max(rate)=',num2str(max(rates))),'Location','southwest');

%xlabel('timestep number, t C_s/L');
%ylabel('Various parameters, normalized by their maximum over the simulation time');
%title('d_i, \eta, S, \delta_{\rm SP} and E_y / v_A B_0 evolution');

%view(2);
%saveas(gcf,strcat(datadir,'/mom_balance.pngs/param_point',num2str(num),'.png'));
