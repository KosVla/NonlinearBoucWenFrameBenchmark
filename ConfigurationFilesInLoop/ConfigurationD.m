%% Main file - Configuration D
%
%Case Study scenario:
%Localized damage concentrated on all links of the first floor (beams+columns)
%Here, the initial/healthy state is represented by the frame with activated
%links and the following set-up:
% beta=21.82,gama=9.71, a= 0.6235, k= 1e8, AmpBW=1.0, deltan=deltav=0
%The induced mechanical damage is represented by the hysteretic Bouc-Wen 
%links (first floor only) simulated with the following parametric set-up:
% beta=21.82,gama=9.71, a=0.1236, k=7.7234e7, AmpBW=1.0, deltan=deltav=0
%
%For each input signal both the initial(MODELH) and the calibrated (MODELD) 
%model are simulated. Their displacement output is stored on the respective
%workspace. Based on the relative error norm(MODELD.U-MODELH.U)/norm(MODELH.U)
%each simulation is characterized as hard/medium (errors: <10/10-20%)
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
MODELHl = InputFileLinksFirst(); ndim=6*numel(MODELHl.nodes(:,1)); 

%Configuration D 
%Localized mechanical damage changes bw_a,bw_k
%Bouc-Wen hysteresis model parameters - Notation follows description pdf
%bw_a=a / Alpha=A / N=w / Beta=b / Gamma=g / deltav=dv / deltan=dnu
%bw_k=k
InputD.bw_a = 0.1236; InputD.bw_k = 7.7234e7; InputD.Alpha=1.0;
InputD.N=1; InputD.Beta=21.82; InputD.Gamma=9.71; InputD.deltav = 0; InputD.deltan=0;
%Additional amplitude parameter for BW.
%Multiplies only hysteretic term to magnify influence
InputD.AmpBW=1.0;

MODELDm = InputFileLinksFirst_Localize(InputD);

%Initial conditions
InputD.u0=zeros(ndim,1); InputD.v0=zeros(ndim,1); InputD.a0=zeros(ndim,1);

%Damping parameters
InputD.zeta = [0.02 0.02]; InputD.OmegaIndexes = [1 2]; 
% Zeta ration for the first two modes.

% Time integration parameters.
fs = 150;              % working sampling frequency.
upsamp = 1;            % upsampling factor to ensure integration accuracy.
fsint = fs*upsamp;     % integration sampling frequency.
Input.dt = 1/fsint;    % integration time step.

%% Excitation signal design. 

P = 5;                  % number of excitation periods.
N = 4000;               % number of points per period.
Nint = N*upsamp;        % number of points per period during integration.
Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;
fres = fsint/Nint;

%Amplitude range
Ampls = [1.0e3 2.0e3];

%Fmin range
fmins = [2 5];
%Fmax range
fmaxs = [25 50];

rng(96)
As = rand(10,1); fs = rand(10,1); fsM = rand(10,1);

Input.Angle= pi/4;

Parameters = zeros(10,1);
filename = 'ConfigurationD1_';

countHard=0; countEasy=0; countMedium=0;

for w=1:10
    % Excitation amplitude
    A = Ampls(1) + As(w)*(Ampls(2) - Ampls(1)); 

    %Excitation frequencies
    fmin = fmins(1) + fs(w)*(fmins(2) - fmins(1)); 
    fmax = fmaxs(1) + fsM(w)*(fmaxs(2) - fmaxs(1));

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
    InputH.bw_a = 0.6235; InputH.bw_k = 1e8;


    %% Evaluate model
    MODELD = BoucWenRun(InputH,MODELDm);
    MODELH = BoucWenRun(InputH,MODELHl);
        
    if sum(sum(isnan(MODELD.U)))==0 && sum(sum(isnan(MODELH.U)))==0
        Error = CheckErrorStruct(MODELH,MODELD);

        Results.OutputU = sparse(MODELD.U);
        Results.OutputUH = sparse(MODELH.U);
        Input.SynthesizedAccelerogram = sparse(Input.SynthesizedAccelerogram);
        Results.Input = Input;
        Results.A = A; Results.fmin = fmin; Results.fmax = fmax;
        
        if Error.Norm<0.10
            countHard=countHard+1;
            filenamesave=strcat(filename, 'HardTask_No',int2str(countHard));            
        elseif Error.Norm<0.20
            countMedium=countMedium+1;
            filenamesave=strcat(filename, 'MediumTask_No',int2str(countMedium));                
        else
            countEasy=countEasy+1;
            filenamesave=strcat(filename, 'EasyTask_No',int2str(countEasy));                
        end
        
%         filenamesave=strcat(filename, 'Simulation_No_',int2str(w));
        save(filenamesave,'Results','-v7.3')
        
        Parameters(w,1)=A; Parameters(w,2)=Input.Angle;
        Parameters(w,3)=fmin; Parameters(w,4)=Error.Norm; 
    end
end
save('Parameters1.mat','Parameters','-v7.3')