clear
close all

addpath('data_and_exec')

%% transmission signal
M = 16;      % Modulation order
k = log2(M); % Bits/symbol
n = M*k*32;  % Transmitted bits
nSamp = 4;   % Samples per symbol

txfilter  = comm.RaisedCosineTransmitFilter('RolloffFactor',0.1, 'FilterSpanInSymbols',16,'OutputSamplesPerSymbol',nSamp);
rxfilter  = comm.RaisedCosineReceiveFilter('RolloffFactor',0.1, 'FilterSpanInSymbols',16,'InputSamplesPerSymbol',nSamp,'DecimationFactor',nSamp);
filtDelay = k*16;
x      = randi([0 1],n,6);
modSig = qammod(x,M,'InputType','bit');

for k1 = 1:6
    txSig(:,2*k1-1) = 6*1e-2*txfilter(modSig(:,k1));
    txSig(:,2*k1  ) = 6*1e-2*txfilter(modSig(:,k1));
end
%scatterplot(txSig(filtDelay+1:4:end-filtDelay))

%% fibre span parameters
L        = 10e3; % meters
dz       = 100;  % meters
XTavg    = -25;  % average fibre XT dB/km
methodXT = 'sa';     %  'sa' semi-analytical, 'expm' exp matrix sol.,...
                       %  'unitaryLc', 'justPolRot', 'fullCoupAll', 'noCoup'
                       % note that Kuv was calculated for method 'sa' so
                       % when running with 'expm' one may get less XT
                       % accumulated
lambda0   = 1550e-9;
f0        = physconst('LightSpeed')/lambda0;

load fNmodes6.mat
load pc6.mat
load Disp6.mat
load S6.mat
load vg6.mat
load nlCoef6.mat

FiberParameters.D           = Disp;
FiberParameters.S           = S;
FiberParameters.nlCoef      = nlCoef;
FiberParameters.lossCoef    = 0.2;
FiberParameters.ModeDelay   = (1./vg-1/vg(1)); %s/m
FiberParameters.dgd_pol     = 1e-16;

%% get step-by-step coupling matrices
seed  = randi(100); % for random fibre perturbations
[CoupMatLPabVect] = calculateCoupMat_2pol(f0,L,dz,XTavg,methodXT,seed); % this takes time worth to save when re-running is necessary

% check accumulated XT
coupMat = eye(12);
for k1 = 1:L/dz; coupMat = coupMat*CoupMatLPabVect(k1).c; end
figure(); imagesc(10*log10(abs(coupMat).^2),[-40 0]); colorbar; title('span coupling matrix [dB]'); xlabel('mode index'); ylabel('mode index')
p = abs(coupMat).^2;
M = eye(12); inds = find(M); M(inds(1:2:end-1)+1) = 1; M(inds(2:2:end  )-1) = 1;
fprintf(['actual accumulated XT [dB]: ',num2str(10*log10(sum(p(find(~M)))/sum(p(find(M))))),'\n'])
fprintf(['target accumulated XT [dB]: ',num2str(XTavg + 10*log10(L/1e3)    ),'\n'])
fprintf(['seed: ',num2str(seed),'\n'])

M = eye(12); inds = find(M); M(inds(1:2:end-1)+1) = 1; M(inds(2:2:end  )-1) = 1;
for k1 = 1:L/dz
    coupMat = abs(CoupMatLPabVect(k1).c).^2;    
	xt(k1) = 10*log10(sum(coupMat(find(~M)))/sum(coupMat(find(M))));
end
figure()
plot(xt)
title('xt per step')

%% fibre transmission
nlInd     = 1;      % 0 no nonlinearities - 1 with nonlinearities
minStep   = 1;      maxStep   = 1e3;
delta_tol = 1e-5;   % local error for adapatative step
fineStep  = 0;      % 1 - using adaptative step, 0 - just fixed coarse step

Sin.E(1:12,:) = txSig.';
Sin.T         = length(txSig);
Sin.dt        = 1/nSamp*1/1e9;
Sout = FMF_transmission_6Modes_2pol('stochastic',FiberParameters,CoupMatLPabVect,L,dz,nlInd,minStep,maxStep,delta_tol,fineStep,Sin);
Aout = Sout.E;

%% Output constellation with matched filter only
figure('Name','constellations - 12 pol modes')
for k1 = 1:12
    subplot(3,4,k1)
    rxSig = rxfilter(Aout(k1,:).');
    scatter(real(rxSig(filtDelay:end-filtDelay)),imag(rxSig(filtDelay:end-filtDelay)),'.')
end
