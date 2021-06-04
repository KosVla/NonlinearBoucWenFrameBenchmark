%% Main file - Configuration B
%
%Case Study scenario:
%Damage growth is induced in the initial state of the frame representing
%a deterioration scenario. Here, the initial/healthy state is represented
%the frame with activated links and the following set-up:
% beta= 3.54,gama= 6.39, a= 0.10, k= 6e7, AmpBW=1.0, deltan=deltav=0
%The damage/deterioration scenario is represented by the hysteretic Bouc-Wen 
%links simulated with the following parametric set-up:
% beta= 3.54,gama= 6.39, a= 0.10, k= 6e7, AmpBW=1.0, (same as the initial)
% deltav = 2.57, deltan=1.42.
%
%For each input signal both the initial(MODELH) and the damaged (MODELD) 
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
% in Workshop on non-linear system identification benchmarks, 2016, pp. 7�14

%% Clear workspace and load folders
clc;clear;
addpath InputFiles
addpath core

%% Define Input

%Select Model
MODEL = InputFileLinksFirst(); ndim=6*numel(MODEL.nodes(:,1));

%Configuration B
%Deterioration scenarion on existing hysteretic frame
%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
Input.bw_a = 0.10; Input.bw_k = 6e7; Input.Alpha=1.0;
Input.N=1; Input.Beta=3.54; Input.Gamma=6.39; Input.deltav = 2.57; Input.deltan=1.42;

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

fmin = 4;               % excitation bandwidth.
fmax = 50;
A = 1e3;                 % excitation amplitude.

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
InputH.deltav = 0; Input.deltan=0;

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