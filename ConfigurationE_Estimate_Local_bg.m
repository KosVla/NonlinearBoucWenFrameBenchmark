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
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7�14

%% Clear workspace and load folders
clc;clear;
addpath InputFiles
addpath core

%% Define Input

%Select Model
MODELHl = InputFileLinksFirst(); ndim=6*numel(MODELHl.nodes(:,1)); 

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
% Zeta ratio for the first two modes.

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
fmax = 70;
A = 1e4;                 % excitation amplitude.

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
InputD.SynthesizedAccelerogram=f; InputD.Angle=pi/4;

%Define initial state 
InputH=InputD;
InputH.deltav = 0.0; InputH.deltan=0.0;

%% Evaluate model
MODELD = BoucWenRun(InputH,MODELDm);
MODELH = BoucWenRun(InputH,MODELHl);
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
