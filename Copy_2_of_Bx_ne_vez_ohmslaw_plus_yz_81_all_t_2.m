


%datadir = '/Volumes/Elements/PSC_DATA/try_collisions2/coll1/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nif/coll1/';
datadir = './h5/';
%datadir = '/Volumes/Elements/PSC_DATA/try_nnb/nnb00075/';
%datadir = '/Volumes/Elements/PSC_DATA/try_par/nif/coll1/';

% initial parameters
MMi =100;
ZZ = 1;
TTe = 0.002;
n= 0.01;
LL0 = 40.0;%sqrt(MMi/(ZZ*n));
BB0 = 0.01; %sqrt(TTe*n);
V0 = BB0/sqrt(MMi*n);
delx=2; %step in derivatives
num=1; %number of current sheet, from 1 to inf, from left to right
sizze=65;
eta0=0.0;

rates=zeros(1,sizze+1);
di=zeros(1,sizze+1);
eta=zeros(1,sizze+1);
SLund=zeros(1,sizze+1);
deltasp=zeros(1,sizze+1);

LRC=100;

%tstart = 40000;
%tstep = 2000;
ts = 10000;% = [tstart:tstep:60000];

 
    
for k=0:sizze

        address=strcat(datadir, 'psc_',num2str(ts+1000*k,'%07d'),'.h5');

        NNe=h5read(address,'/NNe');
        NNi=h5read(address,'/NNi');
        dx=h5read(address,'/dx');
        dz=h5read(address,'/dz');
        dt=h5read(address,'/dt');
        
        h5disp(address);
        
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
        FIG=1

    figure(FIG)
     close(FIG)
     figure(FIG)
    clf

    set(FIG, 'PaperPosition', [0.5 2.5 6 4])
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
    
    Nu = eta0 * BB0 .* NNe .* ((TTe./Te).^(1.5));
    Xi = bx ./ Nu;
    D = Xi.^4 + 14.79*Xi.^2 + 3.77;
    F1 = 1 - (6.42*Xi.^2 + 1.84) ./ D;
    Eta_perp = Nu ./ NNe;
    etaJ=NNe .* Eta_perp .* F1 .* jy;
       
    meanby=mean(by,1);
    
    size(meanby)
    size(xs)
    meanTe=mean(Te,1);
  %  xs
    
    meanNNe=mean(NNe,1);
    meanVxe=mean(NVxe./NNe,1);
   
     xxx=plot(xs, squeeze(meanby)/BB0);
    hold on
   
    yyy=plot(xs, squeeze(meanNNe)/n);%
    zzz=plot(xs, squeeze(meanVxe)/(V0));
    fff=plot(xs,squeeze(meanTe)/TTe);
%    meanTe=mean(Te,1);
    
%    zzzz=plot(xs, squeeze(meanTe)/TTe);
    
    xind2=6400+[0:2500];

    [maxvaldens,mind] = max(meanNNe(xind2)/n);

    xind3=6400+[0:mind];

    xshock=mind/50;

    xdown=xshock-1;

    xup=xshock+1;

    zindup=[-50:49];
    
    plot(mind/50.*ones(length(zindup),1),zindup,'r--');
    plot((mind-50)/50.*ones(length(zindup),1),zindup,'r--');
    plot((mind+50)/50.*ones(length(zindup),1),zindup,'r--');



%    x=x/100;
   
   % rhs  = + Vze(zind,ixX) .* bx(zind,ixX)...
   %     - Vxe(zind,ixX) .* bz(zind,ixX)...
   %     +  (Tyze(zind+delx,ixX) - Tyze(zind-delx,ixX))  / (2*delx*dfield.dz*dfield.decZ) ...
   %     +  (Txye(zind,ixX+delx) - Txye(zind,ixX-delx))  / (2*delx*dfield.dx*dfield.decX) ...
   %     +  (Pyer(zind,ixX) - Pyel(zind, ixX) ) / (dfield.dt*(4000))...
   %     -  etaJ(zind,ixX);
   % plot(zs(zind)-zX, -rhs / KK) 
    
   set(gca,'fontsize',20,'LineWidth',2)
    set(xxx,'LineWidth',2)
%    set(yyy,'LineWidth',2)
%    set(zzz,'LineWidth',2)
%    set(zzzz,'LineWidth',2)
%    
    
    legend({'B_y','n','v_{ex}', 'T_{e}'},'FontSize',15,'Location','northeast')
    
%    xlabel('z / d_{i0}','FontSize',20)
    
%    ylabel('{B_x [3.0 B0], n_e [0.6n_{\rm max}], v_{ez} [0.42 C_{s}] }','FontSize',20)
    ylim([-15.0 30.0])
    xlim([0.0 50.0])
    title(sprintf('Profile in X point, t C_s/L = %.3f', (ts+1000*k) *(dt * BB0/MMi)) )

    view(2);
    saveas(gcf,strcat(datadir,'1d_', num2str(ts+1000*k,'%07d'),'_2.png'));
   
   

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
