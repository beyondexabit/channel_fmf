function [CoupMatLPab] = fmf6CoupMat_expm(f0,L,dz,XTavg_km)
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

for k1 = 1:nModes
    B(k1) = ( pc(k1) ) * w0 / c0 / 2;
end

%% Offset
if ~isinf(XTavg_km)
    XTavg_dz = XTavg_km + 10*log10(dz/1000);
    rhoExt   = [0:1e-8:1e-6 2e-6:1e-6:rho(end)];
    xtExt    = interp1(rho(1:end),xt(1:end),rhoExt);
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
        
        for k1 = 1:6
            for k2 = 1:6
                kuvX(k1,k2) = KuvSurf(k1,k2).s(dx , dy)/varC;
            end
        end
        aux = expm(1i*diag(B)+1i*kuvX);      
        
        if sum(abs(aux(:)).^2) < 6.001
            break;
        else
            123;
        end
    end
    
    CoupMatLPab(kx).c = aux;
end