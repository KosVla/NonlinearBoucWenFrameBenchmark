function [MODEL] = GetHystereticParams(MODEL,InputParams)
%Initializes and assembles on the MODEL struct all matrices needed 
% for time integration
%
%Input:
%  InputParams : struct / Contains all input parameters relevant to 
%                         Bouc-Wen, time integration, Boundary conditions etc.
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. Here, the vectors and matrices
%                   needed for the integration are assembled.
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties,
%                   parameters and matrices. 
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

%Properties of the nl links
bw_a= InputParams.bw_a; 
Alpha=InputParams.Alpha; 
Amp = InputParams.AmpBW; N = InputParams.N;
deltav = InputParams.deltav; deltan=InputParams.deltan;
Gamma=InputParams.Gamma; Beta=InputParams.Beta;
bw_k = InputParams.bw_k; 

MODEL.nl_link_bw_properties = [bw_a bw_k Alpha Beta Gamma N Amp deltav deltan];

%Assembly of matrix of BW links
%Each row corresponds to a degree of freedom hysteretic link
%Here an individual degree of freedom properties can be changed by changing the
%respective row.
%The elements contained in MODEL.nl_links_alternate have been treated and
%the different BW parameters are saved on the activated dofs
MODEL = GetDofsLoop(MODEL);

%Initial conditions and matrices initialization
nt = MODEL.dyn.nt;
HIST.R = zeros(size(MODEL.DofsLoop,1),1); HIST.Um = HIST.R;
HIST.E = HIST.R; HIST.Zeta = HIST.R;

if (isfield(MODEL,'DofsLoopLinear'))
    HISTL.R = zeros(size(MODEL.DofsLoopLinear,1),1); HISTL.Um = HISTL.R;
    HISTL.E = HISTL.R; HISTL.Zeta = HISTL.R;
    MODEL.HistBWL=HISTL;
end

MODEL.HistBW=HIST;
MODEL.HistR = zeros(size(MODEL.DofsLoop,1),nt);
MODEL.HistU = MODEL.HistR;

MODEL.nt=1; MODEL= assemble_nlBW( MODEL );
MODEL= apply_bc_nl( MODEL );

%Apply numerical correction for singular M
test = diag(MODEL.M);
m=min(test(test>0))/1000;
for i=1:size(MODEL.M,1)
    if MODEL.M(i,i)==0
        MODEL.M(i,i)=m;
    end
end

if (isfield(InputParams,'zeta')~=0)
    Kstiff = MODEL.K(MODEL.freedofs,MODEL.freedofs);
    Mmass = MODEL.M(MODEL.freedofs,MODEL.freedofs);
    [MODEL.dyn.a,MODEL.dyn.b,~] = GetRayleighDamping(Kstiff,Mmass,InputParams.zeta,InputParams.OmegaIndexes);
end

neq = numel(MODEL.u);
MODEL.U = zeros(neq,nt); MODEL.V = zeros(neq,nt);
MODEL.A = zeros(neq,nt);

MODEL.U(:,1)=InputParams.u0;
MODEL.V(:,1)=InputParams.v0;
MODEL.A(:,1)=InputParams.a0;

%Assembly of external excitation matrices
[S,F] = AssembleSandZ(InputParams,MODEL);
Fmatrix = S*F'; 
MODEL.Rmatrix = Fmatrix;
 
[~,pos]=max(F);
MODEL.pos=pos;


end

