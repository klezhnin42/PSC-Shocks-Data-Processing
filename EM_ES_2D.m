function [ Emx,Emy,Emz, Esx,Esy,Esz ] = EM_ES_2D( ex,ey,ez,KX, KZ, delsq )
  FEx = fftn(squeeze(ex));
  FEy = fftn(squeeze(ey));
  FEz = fftn(squeeze(ez));
  FEx(1,1) = 0;
  FEy(1,1) = 0;
  FEz(1,1) = 0;
 
  size(FEx);
  size(FEy);
  size(FEz);
  
  size(KX);
  size(KZ);
%  size(FEz)

Frho =  1j*KX.*FEx + 1j*KZ.*FEz;

%% Calculation of Fourier transfrom of Phi from Frho
FPhi = Frho./delsq;

%% Calculate E_static Fourier components from F Phi
FEz_s =  -1j*KZ.*FPhi;
FEy_s =  0.*FPhi;
FEx_s =  -1j*KX.*FPhi;

Es_x = real(ifftn(FEx_s ));
Es_y = real(ifftn(FEy_s ));
Es_z = real(ifftn(FEz_s));

%% Compute real E_electromagnetic components
Em_x = real(ifftn(FEx -FEx_s ));
Em_y = real(ifftn(FEy -FEy_s ));
Em_z = real(ifftn(FEz -FEz_s ));

%% Compute real E_electromagnetic components
Esx = Es_x;
Esy = Es_y;
Esz = Es_z;

Emx = Em_x;
Emy = Em_y;
Emz = Em_z;

end

