function [CoupMatLPab] = fmf6CoupMat_DEN_MULT(f0,L,dz,XTavg_km)
% INPUTS
%    f0        - Center Frequency
%    L         - Fiber length (m)
%    dz        - Step length (m) (typical value: 100)
%
% OUTPUTS
%    CoupMatLPab - coupling matrix per step
%
% REFERENCES
%
% REVISION
%	created by Filipe Ferreira @ 2022

load fNmodes6.mat
load pc6.mat
load Disp6.mat
load S6.mat
load KuvSurf6.mat
load rho_6.mat
load XT1d_rho_max_6.mat
xt = 10*log10(mean(10.^(XT1d_rho_max/10),2));

disp('Calculating the linear mode coupling operator...')

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Propagations constants load or calculation
nModes = 6;

%% parameters and global variables
w0     = 2*pi*f0;
nsteps = ceil(L/dz);

%% deltaBeta1 for LPuv,ab,xy
for k1 = 1:nModes
    for k2 = 1:nModes
        deltaBeta1(k1,k2) = ( pc(k1) - pc(k2) ) * w0 / c0 / 2;
    end
end
dB12 = deltaBeta1(1,2); dB13 = deltaBeta1(1,3); dB14 = deltaBeta1(1,4); dB15 = deltaBeta1(1,5); dB16 = deltaBeta1(1,6);
dB23 = deltaBeta1(2,3); dB24 = deltaBeta1(2,4); dB25 = deltaBeta1(2,5); dB26 = deltaBeta1(2,6);
dB34 = deltaBeta1(3,4); dB35 = deltaBeta1(3,5); dB36 = deltaBeta1(3,6);
dB45 = deltaBeta1(4,5); dB46 = deltaBeta1(4,6);
dB56 = deltaBeta1(5,6);

%% Offset
if ~isinf(XTavg_km)
    XTavg_dz = XTavg_km + 10*log10(dz/1000);
    rhoExt   = rho(2):1/10000:0.3;
    xtExt    = interp1(rho(2:end),xt(2:end),rhoExt);
    [~,Ind]  = min(abs(xtExt -XTavg_dz));
    %rhoExt(Ind)
end

%% Coupling Matrix between LPuv,ab (ignoring polarizations)
coupMat = zeros(nModes,nModes,1e4);
Aout    = zeros(nModes,nModes,1e4);

XT1d_temp = zeros(nModes,1e4);

h = tic;
for kx = 1:nsteps
    if rem(kx,100) == 0
        tt = toc(h);
        disp([num2str(kx),' out of ',num2str(nsteps),' - tte = ',num2str(tt*(nsteps-kx)/60/100),' min'])
        h = tic;
    end
    
    while 1
        % Random Fiber Core Displacement
        rInd = Ind + round(2*(rand-0.5));
        if rInd < 1; rInd = 1; end
        drho = rhoExt( rInd );
        dphi = 2*pi * rand(1);
        
        dx = drho .* cos(dphi);
        dy = drho .* sin(dphi);
        
        varC = 1;
        
        k12 = KuvSurf(1,2).s(dx , dy)/varC; k13 = KuvSurf(1,3).s(dx , dy)/varC; k14 = KuvSurf(1,4).s(dx , dy)/varC; k15 = KuvSurf(1,5).s(dx , dy)/varC; k16 = KuvSurf(1,6).s(dx , dy)/varC;
        k23 = KuvSurf(2,3).s(dx , dy)/varC; k24 = KuvSurf(2,4).s(dx , dy)/varC; k25 = KuvSurf(2,5).s(dx , dy)/varC; k26 = KuvSurf(2,6).s(dx , dy)/varC;
        k34 = KuvSurf(3,4).s(dx , dy)/varC; k35 = KuvSurf(3,5).s(dx , dy)/varC; k36 = KuvSurf(3,6).s(dx , dy)/varC;
        k45 = KuvSurf(4,5).s(dx , dy)/varC; k46 = KuvSurf(4,6).s(dx , dy)/varC;
        k56 = KuvSurf(5,6).s(dx , dy)/varC;
        
        kuvX = [0   k12 k13 k14 k15 k16
            k12 0   k23 k24 k25 k26
            k13 k23 0   k34 k35 k36
            k14 k24 k34 0   k45 k46
            k15 k25 k35 k45 0   k56
            k16 k26 k36 k46 k56 0  ];
        
        [coupMat] = coupMat6x6_mat(dB12, dB13, dB14, dB15, dB16, dB23, dB24, dB25, dB26, dB34, dB35, dB36, dB45, dB46, dB56, k12, k13, k14, k15, k16, k23, k24, k25, k26, k34, k35, k36, k45, k46, k56);
        
        PoutX(:,:,:) = abs(coupMat).^2;
        for nM = 1:6
            if nM == 3 || nM == 4
                PoutZ(:,:) = abs(coupMat(:,nM,:)).^2;
                XT1d_temp(nM,:) = 10*log10((sum(PoutZ,1)-sum(PoutZ(3:4,:)))./sum(PoutZ(3:4,:)));
            elseif nM == 5 || nM == 6
                PoutZ(:,:) = abs(coupMat(:,nM,:)).^2;
                XT1d_temp(nM,:) = 10*log10((sum(PoutZ,1)-sum(PoutZ(5:6,:)))./sum(PoutZ(5:6,:)));
            else
                PoutZ(:,:) = abs(coupMat(:,nM,:)).^2;
                XT1d_temp(nM,:) = 10*log10((sum(PoutZ,1)-PoutZ(nM,:))./PoutZ(nM,:));
            end
            
            XT1d_temp(isinf(XT1d_temp)) = -100;
        end
        
        [~,IndX] = min( abs( XT1d_temp(2,:) - XTavg_dz ) );
        XT = mean(XT1d_temp(2,IndX(1)));
        
        aux = coupMat(:,:,IndX);
        if sum(abs(aux(:)).^2) < 6.001
            break;
        else
            123;
        end
    end
    
    CoupMatLPab(kx).c = aux;
end