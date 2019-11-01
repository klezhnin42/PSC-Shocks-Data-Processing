FIG=1
figure(FIG)
close(FIG)
figure(FIG)

set(FIG, 'PaperUnits', 'inches')
set(gcf, 'PaperPosition', [0 0 9.6 4])
set(FIG, 'DefaultAxesFontSize', 12)
set(FIG, 'DefaultTextFontSize', 12)
set(FIG, 'DefaultLineMarkerSize', 4)



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
   

%datadir = '/Users/klezhnin/Desktop/shock/eighty/try3/';
datadir = './h5_saved/';
% initial parameters
MMi =100;
ZZ = 1;
TTe = 0.002;
n= 0.05;
BB0 = 0.01; %sqrt(TTe*n);
V0 = BB0/sqrt(MMi*n);
sizze=200;

%tstart = 40000;ww
%tstep = 2000;
ts = 0;% = [tstart:tstep:60000];

for k=0:sizze

        k
        address=strcat(datadir, 'psc_',num2str(ts+k*1000,'%07d'),'.h5');

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

        ex = h5read(address,'/ex'); 
        ey = h5read(address,'/ey');
        ez = h5read(address,'/ez');

        
        size(ex);
        size(ey);
        size(ez);
       
lz=10000;
lx=10;
        
  kx = (2*pi/(dx*lx))*[0:(lx/2-1) (-lx/2):(-1)]; % Vector of wavenumbers
  %ky = (2*pi/(dx*ly))*[0:(ly/2-1) (-ly/2):(-1)]; % Vector of wavenumbers
  kz = (2*pi/(dx*lz))*[0:(lz/2-1) (-lz/2):(-1)];
  
  [KX, KZ]  = meshgrid(kx,kz); % Matrix of (x,y) wavenumbers corresponding

  %% This will be the matrix to convert from rho to Phi in fourier space
  delsq = (KX.^2+ KZ.^2);
  
  %% Global term will be necessarily set to smallw
  delsq(1,1,1) = 100000000;
        
        
[Emx,Emy,Emz, Esx,Esy,Esz] = EM_ES_2D(ex,ey,ez,KX, KZ, delsq );

   meanNNe=mean(NNe,3);
%   xind2=12800+[0:5000];
   size(meanNNe);
%   [maxvaldens,mind] = max(meanNNe(xind2));
   meanex=mean(ex,3);
   size(meanex);
   size(Esx);
   meanEsx=mean(Esx,2);
   size(meanEsx);
   meanEsy=mean(Esy,2);
   meanEsz=mean(Esz,2);
   size(meanEsz);
   meanbtot=mean(sqrt(bx.^2+by.^2+bz.^2),3);
   meanEstot=sqrt(meanEsx.^2+meanEsy.^2+meanEsz.^2);
%   max(meanEmx)
%   max(meanEmy)
%   max(meanEmz)
   meanEmx=mean(Emx,2);
   size(meanEmx);
   meanEmy=mean(Emy,2);
   meanEmz=mean(Emz,2);
   meanEmtot=sqrt(meanEmx.^2+meanEmy.^2+meanEmz.^2);
   max(meanEmx)
   max(meanEmy)
   max(meanEmz)   
   xfast=load('xfast.mat','xfast');
   ptot=load('ptot.mat','ptot');
   zind=linspace(-10.0,10.0,200);
   figure(1)
   clf

   plot(zs,meanNNe/n)
   hold on
   plot(zs,meanbtot/BB0)
   plot(zs,10.0*(meanEmx/BB0))
   plot(zs,10.0*(meanEmy/BB0))
   plot(zs,10.0*(meanEsz/BB0))
%   plot(zs,meanEmtot/BB0)
   plot(xfast.xfast(k+1)*ones(length(zind),1),zind,'r--')
   plot(xfast.xfast(k+1),ptot.ptot(k+1),'ro')
%   plot()
%   plot(zs,meanEsy/BB0)
%   plot(zs,meanEsz/BB0)
   xlim([xfast.xfast(k+1)-15 xfast.xfast(k+1)+5])
   ylim([-10 10])
%   set(gca, 'YScale', 'log')
   legend('n','Btot','10Emx', '10Emy','10Esz')
   title(strcat('lin n, Btot, Emx, Emy, Esz,',' \Omega_i t=',num2str(BB0*dt*(ts+1000*k)/MMi,'%.2f')))
   view(2);
   saveas(gcf,strcat(datadir,'1dems_comp_lin', num2str(ts+1000*k,'%07d'),'_4.png'));
%xshock=mind/50;


%figure(1)
%clf




%subplot(2,2,1)

%axes('position', [0 0 1 1], 'visible', 'off')
%imagesc(xs, zs, squeeze(ex./BB0))
%xlim([xshock xshock+5])
%axis equal
%xlim([-0.2 0.2])
%ylim([0 300])
%title('E_x/B_0')
%set(gca,'YDir','normal')
%plot(xs,ex/BB0)
%caxis([-1 1])
%colormap(gca,map)
%colorbar
%xlabel('x/di')
%ylabel('z/di')




%subplot(2,2,2)
%axes('position', [0 0 1 1], 'visible', 'off')
%imagesc(xs, zs, log10(squeeze(abs(Esx/BB0))))
%xlim([xshock xshock+5])
%axis equal
%xlim([-0.2 0.2])
%ylim([0 300])
%set(gca,'YDir','normal')
%title('log_{10}|E_{sx}/B_0|')
%colormap(gca,jet);
%caxis([-3 1])
%colorbar;


%subplot(2,2,3)
%axes('position', [0 0 1 1], 'visible', 'off')
%imagesc(xs, zs, squeeze(ez/BB0))
%set(gca,'YDir','normal')
%title('E_z/B_0')
%caxis([-1 1])
%colormap(map)w
%caxis([-1 1])
%colormap(gca,map)
%colorbar
%plot(xs,ey/BB0)
%shading interpw
%colorbar;
%caxis([-1 1])
%xlabel('x/di')
%ylabel('z/di')
%xlim([xshock xshock+5])
%axis equal
%ylim([0 300])
%xlim([-0.2 0.2])

%subplot(2,2,4)

%axes('position', [0 0 1 1], 'visible', 'off')
%imagesc(xs, zs, log10(squeeze(abs(Esz/BB0))))
%set(gca,'YDir','normal')
%title('log_{10}|E_{sz}/B_0|')
%colorbar;
%colormap(jet);
%caxis([-3 1])
%plot(xs,Es_y)
%shading interp
%caxis([-1 1])
%xlabel('x/di')
%ylabel('z/di')
%xlim([xshock xshock+5])
%axis equal
%ylim([0 300])
%xlim([-0.2 0.2])
%   view(2);
%    saveas(gcf,strcat(datadir,'es_fields', num2str(ts+1000*k,'%07d'),'_4.png'));

end

quit
