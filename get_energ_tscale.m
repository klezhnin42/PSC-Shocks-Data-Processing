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

z=zarray.zarray;

pe=pearray.pearray;

tenerg=[];

for i=1:size(z,1)

    pe1=pe(i,:);
    z1=z(i,:);
    pe1(1)=[];
    z1(1)=[];
    

%    xedge = 300 * [0:0.005:1];

%xedge = max(sqrt(px.^2)) * [-1:0.005:1];

%   yedge = 2.0 * [0:0.005:1];


%    hgrm=histc2d(z1,pe1,xedge,yedge);
%    loghgrm=log10(hgrm);
%    loghgrm(loghgrm==-Inf)=NaN;
%   loghgrm(~isfinite(loghgrm))=0;
%    size(xedge)
%    size(yedge)
%    size(loghgrm) s
    smoothpe=smoothdata(pe1);
    [Mmax,Imax]=max(smoothpe);
    [Mefold,Iefold]=min(abs(smoothpe-max(smoothpe)/exp(1.0)));
    tenerg=[tenerg; (Imax-Iefold)/length(ts)*8.0];
%    colorbar
%    caxis([1e-4 2])
%    set(gca,'YDir','normal')
%    set(f,'EdgeColor', 'none')
%    set(f,'HandleVisibility','off%')
%    colormap(jet)
%    shading flat


end

    sum(tenerg)/size(z,1)



quit
