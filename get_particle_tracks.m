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

pearray=idlarray;
pearray=sortrows(pearray.',1).';

zarray=idlarray;
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

        pebuff=zeros(1,size(idlarray,1));
        zbuff=zeros(1,size(idlarray,1));

        for i=1:length(idlarray)
	    indadd=find(idL==idlarray(i));
	    if isempty(indadd)
                indtopush=find(idlarray==idlarray(i));
                pebuff(indtopush)=NaN;
                zbuff(indtopush)=NaN;
            else
                indtopush=find(idlarray==idlarray(i));
                pebuff(indtopush)=petot(indadd);
                zbuff(indtopush)=z(indadd);
            end
        end

%        pearray
%        pebuff

        pearray=[pearray; pebuff'];
%        pearray=reshape(pearray,1+p,size(idlarray,1));

%        size(pebuff);

        zarray=[zarray; zbuff'];
%       zarray=reshape(zarray,1+p,size(idlarray,1));

%        size(pearray)
%        size(zarray)
           


end
 
    zarray=reshape(zarray,size(idlarray,1),size(ts,2)+1);
    pearray=reshape(pearray,size(idlarray,1),size(ts,2)+1);


    save('zarray_full.mat','zarray','-v7.3')

    save('pearray_full.mat','pearray','-v7.3')



quit
