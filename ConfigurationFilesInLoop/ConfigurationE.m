%% Main file - Configuration E
%
%Case Study scenario:
%Damage growth is induced in the initial state of the frame representing
%a localized deterioration scenario (first floor links only).
%Here, the initial/healthy state is represented the frame with activated 
%links and the following set-up::
% beta=16.32,gama= -16.32, a= 0.10, k= 5e7, AmpBW=1.0, deltah=deltav=0
%The damage/deterioration scenario is represented by the hysteretic Bouc-Wen 
%links simulated with the following parametric set-up (first floor links only):
% beta=16.32,gama= -16.32, a= 0.10, k= 5e7, AmpBW=1.0, deltan=4.85,
% deltav=2.31
%
%For each input signal both the initial(MODELH) and the damaged (MODELD) 
%model are simulated. Their displacement output is stored on the respective
%workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
%each simulation is characterized as hard/medium (errors: <10/10-20%/)
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.
%and
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14

%% Clear workspace and load folders
clc;clear;
addpath InputFiles
addpath core

%% Define Input

%Select Model
MODELHl = InputFileLinksFirst(); ndim=6*numel(MODEL.nodes(:,1)); 

%Configuration E 
%Localized deterioration scenario
%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dn
%bw_k=k
InputD.bw_a = 0.10; InputD.bw_k = 5e7; InputD.Alpha=1.0;
InputD.N=1; InputD.Beta=16.32; InputD.Gamma= -16.32; 
InputD.deltav = 2.31;  InputD.deltan=4.85;
% Additional amplitude parameter for BW.
%Multiplies only hysteretic term to magnify influence
InputD.AmpBW=1.0;

MODELDm = InputFileLinksFirst_LocalizeE(InputD);

%Initial conditions
InputD.u0=zeros(ndim,1); InputD.v0=zeros(ndim,1); InputD.a0=zeros(ndim,1);

%Damping parameters
InputD.zeta = [0.02 0.02]; InputD.OmegaIndexes = [1 2]; 
% Zeta ration for the first two modes.

% Time integration parameters.
fs = 150;              % working sampling frequency.
upsamp = 1;            % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;     % integration sampling frequency.
InputD.dt = 1/fsint;    % integration time step.

%% Excitation signal design. 

P = 5;                  % number of excitation periods.
N = 12000;               % number of points per period.
Nint = N*upsamp;        % number of points per period during integration.
fmin = 2;               % excitation bandwidth.
fmax = 50;
Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;
fres = fsint/Nint;

Psis = rand(10,1);
As = zeros(10,1);
for w=1:10
    As(w)= (w-1)*0.10+ rand*0.10;
end

Parameters = zeros(10,1);

filename = 'ConfigurationE1_';
RefAmp = 5e3;
for p=1:10
    A = RefAmp + As(p)*(1.25*RefAmp-RefAmp); % excitation amplitude.
    Input.Angle= pi/4;

    Q = zeros(Nint,1);      % definition of the multisine excitation.
    exclines = 1:ceil(fmax/fres);
    exclines(exclines < floor(fmin/fres)) = [];

    %Set random seed
    Q(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));

    f = 2*real(ifft(Q));
    f = A*f/std(f);
    f = repmat(f,[P 1]);
    InputD.SynthesizedAccelerogram=f; 

    %Define initial state
    InputH=InputD;
    InputH.deltav = 0.0; InputH.deltan=0.0;

    %% Evaluate model
    MODELD = BoucWenRun(InputH,MODELDm);
    MODELH = BoucWenRun(InputH,MODELHl);
    
    if sum(sum(isnan(MODELD.U)))==0 && sum(sum(isnan(MODELH.U)))==0
        Error = CheckErrorStruct(MODELH,MODELD);

        Results.OutputU = sparse(MODELD.U); Results.OutputUH = sparse(MODELH.U);
        Input.SynthesizedAccelerogram(2:end+1) = Input.SynthesizedAccelerogram(1:end) ;
        Input.SynthesizedAccelerogram(1)=0; 
        Input.SynthesizedAccelerogram = sparse(Input.SynthesizedAccelerogram);
        Results.Input = Input;
        Results.A = A; Results.fmin = fmin; Results.fmax = fmax;
        Parameters(p,1)=A; Parameters(p,2)=Input.Angle;
        Parameters(p,3)=fmin; Parameters(p,4)=Error.Norm; 
        filenamesave=strcat(filename, 'SimulationNo_',int2str(p));
        save(filenamesave,'Results','Parameters','-v7.3')
    else
        Parameters(p,1)=A; Parameters(p,2)=Input.Angle; Parameters(p,3)=fmin; Parameters(p,4)=0;          
    end
end
%save(filename,'-v7.3')
