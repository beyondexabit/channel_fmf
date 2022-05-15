function varargout = FMF_transmission_6Modes_2pol(methodSSMF,FiberParameters,CoupMatLPab,L,dz,nlInd,minStep,maxStep,delta_tol,fineStep,varargin)
% REVISION
% created by Filipe Ferreira @ 2022

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Source parameters
lambda0  = 1550e-9;
w0       = 2*pi*c0/lambda0;
f0       = w0/(2*pi);

% Calling from Matlab. Extract signals and create matrix Ain
Ain   = varargin{1}.E;

Disp(1:2:12)    = FiberParameters.D; Disp(2:2:12) = FiberParameters.D;
S(1:2:12)       = FiberParameters.S; S(2:2:12) = FiberParameters.S;
lossCoef        = FiberParameters.lossCoef;
DMD(1:2:12)     = FiberParameters.ModeDelay-FiberParameters.dgd_pol/2; DMD(2:2:12) = FiberParameters.ModeDelay+FiberParameters.dgd_pol/2;

% Define frequency vector f
TimeWindow     	= varargin{1}.T*varargin{1}.dt;             % Duration of signal in seconds
Fs            	= 1/TimeWindow;                             % Frequency spacing
f             	= (-varargin{1}.T/2:varargin{1}.T/2-1).'*Fs;

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Fiber Model General Parameters
nlCoefOrig = FiberParameters.nlCoef;
nlCoefOrig2pols(1:2:12,1:2:12) = nlCoefOrig; nlCoefOrig2pols(2:2:12,2:2:12) = nlCoefOrig; nlCoefOrig2pols(1:2:12,2:2:12) = nlCoefOrig; nlCoefOrig2pols(2:2:12,1:2:12) = nlCoefOrig;

switch methodSSMF
    case 'stochastic'
        for k1 = 1:6
            nlCoefWeigthed(2*k1-1,1:2:12) = 2/1;
            nlCoefWeigthed(2*k1  ,2:2:12) = 2/1;
            
            nlCoefWeigthed(2*k1-1,2:2:12) = 2/3;
            nlCoefWeigthed(2*k1  ,1:2:12) = 2/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 1/1;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 1/1;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 2/3; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 2/3; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
    case 'wcoupled'
        for k1 = 1:6
            nlCoefWeigthed(2*k1-1,1:2:12) = 4/3;
            nlCoefWeigthed(2*k1  ,1:2:12) = 4/3;
            
            nlCoefWeigthed(2*k1-1,2:2:12) = 4/3;
            nlCoefWeigthed(2*k1  ,2:2:12) = 4/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 8/9;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 8/9;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 8/9; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 8/9; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
    case 'scoupled'
        nlCoef2pols = 4/3*12/13*mean(mean(nlCoefOrig2pols))*ones(12,12);
    case 'unitaryLc'
        for k1 = 1:6
            nlCoefWeigthed(2*k1-1,1:2:12) = 2/1;
            nlCoefWeigthed(2*k1  ,2:2:12) = 2/1;
            
            nlCoefWeigthed(2*k1-1,2:2:12) = 2/3;
            
            nlCoefWeigthed(2*k1  ,1:2:12) = 2/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 1/1;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 1/1;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 2/3; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 2/3; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
end
nlCoef2pols = nlInd*nlCoef2pols;

%% Transmission
[Aout,~       ] = fmf_NL_xM_2pol(Ain,f,f0,DMD, Disp,  S, nlCoef2pols,lossCoef,L,dz,dz,CoupMatLPab,minStep,maxStep,delta_tol,fineStep);

varargout{1}       = varargin{1};
varargout{1}.E     = Aout;

