%% Main file - Configuration C 
%
%Case Study scenario:
%Model calibration scenario during operation stages
%Here, the initial/healthy state is represented by the frame with activated
%links and the following set-up:
% beta= 4.67,gama= 1.21, a= 0.10, k= 1e8, AmpBW=1.0, deltan=deltav=0
%The model calibration part requires estimating the parameters of the model
%for a changed state, represented by a new hysteretic links set-up:
% beta= 30.87,gama= -40.21, a= 0.10, k= 1e8, AmpBW=1.0, deltan=deltav=0
%
%For each input signal both the healthy(MODELH) and the calibrated (MODELD) 
%model are simulated. Their displacement output is stored on the respective
%workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
%each simulation is characterized as hard/medium/easy (errors: <10/10-20%/>20%)
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.
%and
% J. Noel and M. Schoukens, 
% Hysteretic benchmark with a dynamic nonlinearity,
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7–14
%
%% Clear workspace and load folders
clc;clear;
addpath InputFiles
addpath core

%% Define Input

%Select Model
MODEL = InputFileLinksFirst(); ndim=6*numel(MODEL.nodes(:,1));

%Configuration C
%Model calibration scenario
%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.bw_a = 0.10; Input.bw_k = 1e8; Input.Alpha=0.75;
Input.N=1; Input.Beta=+30.87; Input.Gamma=-40.21; Input.deltav = 0; Input.deltan=0;

%Additional amplitude parameter for BW.
%Multiplies only hysteretic term to magnify influence
Input.AmpBW=1.0;

%Initial conditions
Input.u0=zeros(ndim,1); Input.v0=zeros(ndim,1); Input.a0=zeros(ndim,1);

%Damping parameters
Input.zeta = [0.02 0.02]; Input.OmegaIndexes = [1 2]; 
% Zeta ration for the first two modes.

% Time integration parameters.
fs = 150;              % working sampling frequency.
upsamp = 1;            % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;     % integration sampling frequency.
Input.dt = 1/fsint;    % integration time step.

%% Excitation signal design. 

P = 5;                  % number of excitation periods.
N = 12000;               % number of points per period.
Nint = N*upsamp;        % number of points per period during integration.

fmin = 5;               % excitation bandwidth.
fmax = 50;
Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;
fres = fsint/Nint;

Psis = rand(10,1);
As = rand(10,1);
%As = zeros(10,1);
%for w=1:10
%    As(w)= (w-1)*0.10+ rand*0.10;
%end

Parameters = zeros(10,1);
filename = 'ConfigurationC1_';
RefAmp = 2e3;
for w=1:10
    A = RefAmp + As(w)*(1.50*RefAmp-RefAmp); % excitation amplitude.
    Input.Angle= pi/4; 
    fmin=w; rng(w*77);

    Q = zeros(Nint,1);      % definition of the multisine excitation.
    exclines = 1:ceil(fmax/fres);
    exclines(exclines < floor(fmin/fres)) = [];

    %Set random seed
    Q(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));

    f = 2*real(ifft(Q));
    f = A*f/std(f);
    f = repmat(f,[P 1]);
    Input.SynthesizedAccelerogram=f; 

    if length(f)>72000
       continue; 
    end

    %Define initial state
    InputH=Input;
    InputH.Beta=4.67; InputH.Gamma=1.21;


    %% Evaluate model
    MODELH = BoucWenRun(InputH,MODEL);
    MODELD = BoucWenRun(Input,MODEL);
    
    if sum(sum(isnan(MODELD.U)))==0 && sum(sum(isnan(MODELH.U)))==0
        Error = CheckErrorStruct(MODELH,MODELD);

        Results.OutputU = sparse(MODELD.U); Results.OutputUH = sparse(MODELH.U);
        Input.SynthesizedAccelerogram(2:end+1) = Input.SynthesizedAccelerogram(1:end);
        Input.SynthesizedAccelerogram(1)=0; 
        Input.SynthesizedAccelerogram = sparse(Input.SynthesizedAccelerogram);
        Results.Input = Input;
        Results.A = A; Results.fmin = fmin; Results.fmax = fmax;
        Parameters(w,1)=A; Parameters(w,2)=Input.Angle;
        Parameters(w,3)=fmin; Parameters(w,4)=Error.Norm; 
        filenamesave=strcat(filename, 'SimulationNo_',int2str(w));
        save(filenamesave,'Results','Parameters','-v7.3')
    else
        Parameters(w,1)=A; Parameters(w,2)=Input.Angle; Parameters(w,3)=fmin; Parameters(w,4)=0;          
    end
end
%save(filename,'-v7.3')
