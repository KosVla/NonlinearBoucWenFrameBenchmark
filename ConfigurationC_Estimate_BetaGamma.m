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
%each simulation is characterized as easy/medium/hard (errors: <10/10-20%/>20%)
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
MODEL = InputFileLinksFirst(); ndim=6*numel(MODEL.nodes(:,1));

%Configuration C
%Model calibration scenario
%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.bw_a = 0.10; Input.bw_k = 1e8; Input.Alpha=1.0;
Input.N=1; Input.Beta=+30.87; Input.Gamma=-40.21; Input.deltav = 0; Input.deltan=0;

%Additional amplitude parameter for BW.
%Multiplies only hysteretic term to magnify influence
Input.AmpBW=1.0;

%Initial conditions
Input.u0=zeros(ndim,1); Input.v0=zeros(ndim,1); Input.a0=zeros(ndim,1);

%Damping parameters
Input.zeta = [0.02 0.02]; Input.OmegaIndexes = [1 2]; 
% Zeta ratio for the first two modes.

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

A =2e4;                 % excitation amplitude.

Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;
F = zeros(Nint,1);      % definition of the multisine excitation.
fres = fsint/Nint;
exclines = 1:ceil(fmax/fres);
exclines(exclines < floor(fmin/fres)) = [];
        
%Set random seed
% rng(106)
F(exclines+1) = exp(complex(0,2*pi*rand(size(exclines))));

f = 2*real(ifft(F));
f = A*f/std(f);
f = repmat(f,[P 1]);
Input.SynthesizedAccelerogram=f; Input.Angle=pi/4;

%Define initial state
InputH=Input;
InputH.Beta=4.67; InputH.Gamma=1.21;

%% Evaluate model
MODELD = BoucWenRun(Input,MODEL);
MODELH = BoucWenRun(InputH,MODEL);
%% If upsampling is performed downsampling is also needed following the 
% the template code of the ExampleMainFile*
if upsamp>1
    MODELD.Uups = MODELD.U;
    MODELD.U = Downsampling(MODELD.U,upsamp,P+1,N);
    MODELD.Vups = MODELD.V;
    MODELD.V = Downsampling(MODELD.V,upsamp,P+1,N);
    MODELH.Uups = MODELH.U;
    MODELH.U = Downsampling(MODELH.U,upsamp,P+1,N);
    MODELH.Vups = MODELH.V;
    MODELH.V = Downsampling(MODELH.V,upsamp,P+1,N);
    f = downsample(Input.SynthesizedAccelerogram,upsamp);
    % Removal of the last simulated period to eliminate the edge effects
    %due to the low-pass filter.
    f = f(1:(P+1-1)*N,:);
    Input.SynthesizedAccelerogramUps=f;
    Input.SynthesizedAccelerogram=f;
    P = P+1-1;
end
%% Results visualization

%Plot deformed state of frame
% scale = 3*1/max(abs(MODEL.u([1:6:end,2:6:end,3:6:end])));
% plot_model( MODEL, scale )

%Plot time history
% dofsXY = [1:6:ndim 2:6:ndim]; 
% [~,pos]=max(max(abs(MODELD.U(dofsXY',:)),[],2));
% ndof = dofsXY(pos);
% plot(MODELD.U(ndof,:));
% hold on
% plot(MODELH.U(ndof,:));
Error = CheckErrorStruct(MODELH,MODELD)


%Plot hysterisis loop on one of the nonlinear links
[~,ndof]=max(max(abs(MODELD.HistU),[],2));
figure
plot(MODELD.HistU(ndof, 1:(end-1)),MODELD.HistR(ndof,1:(end-1)),'-b');
% plot(MODELH.HistU(ndof, 1:(end-1)),MODELH.HistR(ndof,1:(end-1)),'-b');
