function [CoupMatLPabxy] = calculateCoupMat_2pol(f0,Lspan,dz,XTavg,method,seed)
% REVISION
% created by Filipe Ferreira @ 2022

s = RandStream('swb2712','Seed',seed); RandStream.setGlobalStream(s);
nsteps = ceil(Lspan/dz);
nModes = 6;
nPols  = 2;

if strcmp(method,'unitaryLc')
    Lc = 10^(-XTavg/10)*1000;
    zc = round(Lc/dz)*dz;
    
    for z = 1:nsteps
        polRot1 = eye(nPols*nModes); polRot2 = eye(nPols*nModes); polRot3 = eye(nPols*nModes);
        polRot4 = eye(nPols*nModes); polRot5 = eye(nPols*nModes); polRot6 = eye(nPols*nModes);
        [polRot1( 1:1: 2, 1:1: 2),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot2( 3:1: 4, 3:1: 4),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot3( 5:1: 6, 5:1: 6),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot4( 7:1: 8, 7:1: 8),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot5( 9:1: 10,9:1:10),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot6(11:1:12,11:1:12),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        polRot = polRot1*polRot2*polRot3*polRot4*polRot5*polRot6;
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = polRot;
        
        realZ = z*dz;
        if rem(realZ,zc) == 0
            auxCoup = sqrt(1/2) * (randn(nPols*nModes,nPols*nModes) + 1i*randn(nPols*nModes,nPols*nModes));
            [Q,R] = qr(auxCoup);
            CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = Q;
            disp('.')
        end
    end
elseif strcmp(method,'justPolRot')
    for z = 1:nsteps
        polRot1 = eye(nPols*nModes); polRot2 = eye(nPols*nModes); polRot3 = eye(nPols*nModes);
        polRot4 = eye(nPols*nModes); polRot5 = eye(nPols*nModes); polRot6 = eye(nPols*nModes);
        [polRot1( 1:1: 2, 1:1: 2),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot2( 3:1: 4, 3:1: 4),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot3( 5:1: 6, 5:1: 6),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot4( 7:1: 8, 7:1: 8),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot5( 9:1: 10,9:1:10),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot6(11:1:12,11:1:12),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        polRot = polRot1*polRot2*polRot3*polRot4*polRot5*polRot6;
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = polRot;
    end
    return;
elseif strcmp(method,'fullCoupAll') || XTavg == Inf
    for z = 1:nsteps
        auxCoup = sqrt(1/2) * (randn(nPols*nModes,nPols*nModes) + 1i*randn(nPols*nModes,nPols*nModes));
        [Q,~] = qr(auxCoup);
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = Q;
    end
elseif strcmp(method,'noCoup') || XTavg == -Inf
    for z = 1:nsteps
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = eye(nPols*nModes);
    end
else
    if contains(method,'expm')
        [CoupMatLPab       ] = fmf6CoupMat_expm(f0,Lspan,dz,XTavg);
    elseif contains(method,'sa')
        [CoupMatLPab       ] = fmf6CoupMat_DEN_MULT(f0,Lspan,dz,XTavg);
    else 
        error('unknow method')
    end
    
    
    for p1 = 1:length(CoupMatLPab)
        d1(p1) = sum(abs(CoupMatLPab(p1).c(:)).^2);
    end
    if max(d1) > 6.1
        warning('fibre matrix is not be unitary')
    end
    
    for z = 1:nsteps
        %% Modal coupling LP??ab,xy
        polRot1 = eye(nPols*nModes); polRot2 = eye(nPols*nModes); polRot3 = eye(nPols*nModes);
        polRot4 = eye(nPols*nModes); polRot5 = eye(nPols*nModes); polRot6 = eye(nPols*nModes);
        [polRot1( 1:1: 2, 1:1: 2),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot2( 3:1: 4, 3:1: 4),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot3( 5:1: 6, 5:1: 6),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot4( 7:1: 8, 7:1: 8),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot5( 9:1: 10,9:1:10),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot6(11:1:12,11:1:12),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        polRot = polRot1*polRot2*polRot3*polRot4*polRot5*polRot6;
        
        [Q,~]             = qr(CoupMatLPab(z).c(1:1:nModes,1:1:nModes));
        CoupMatLPabxyTemp = zeros(nPols*nModes);
        CoupMatLPabxyTemp(1:2:nPols*nModes,1:2:nPols*nModes) = Q;
        CoupMatLPabxyTemp(2:2:nPols*nModes,2:2:nPols*nModes) = Q;
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = CoupMatLPabxyTemp*polRot;
    end
end