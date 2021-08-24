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

%% Clear workspace and load folders
clc;clear;
addpath InputFiles
addpath core

%% Define Input

%Select Model
MODEL = InputFileLinksFirst(); ndim=6*numel(MODEL.nodes(:,1));

%Configuration B
%Deterioration scenarion on existing hysteretic frame
%Configuration A: 
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
Pfilter = 1;            % extra period to avoid edge effects during low-pass filtering (see line 59).
P = P + Pfilter;
fres = fsint/Nint;

%Amplitude range
Ampls = [1.0e3 2.0e3];

%Fmin range
fmins = [2 5];
%Fmax range
fmaxs = [25 50];

rng(232)
As = rand(10,1); fs = rand(10,1); fsM = rand(10,1);

Input.Angle= pi/4;

Parameters = zeros(10,1);
filename = 'ConfigurationB1_';

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
    Input.SynthesizedAccelerogram=f; 
    
    if length(f)>72000
       continue; 
    end

    %Define initial state
    InputH=Input;
    InputH.deltav = 0; InputH.deltan=0;


    %% Evaluate model
    MODELH = BoucWenRun(InputH,MODEL);
    MODELD = BoucWenRun(Input,MODEL);
        
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

