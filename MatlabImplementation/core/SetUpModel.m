function [MODEL] = SetUpModel(MODEL, InputParams)

%Properties of the Bouc-Wen links
bw_a= MODEL.BW.bw_a; Alpha=MODEL.BW.Alpha; 
N = MODEL.BW.N; bw_k = MODEL.BW.bw_k; 
deltav = MODEL.BW.deltav; deltan=MODEL.BW.deltan;
Gamma=MODEL.BW.Gamma; Beta=MODEL.BW.Beta;

MODEL.nl_link_bw_properties =...
    [bw_a bw_k Alpha Beta Gamma N deltav deltan];

%Assembly of matrix of Bouc-Wen links
MODEL.BW.HystLinks = AssembleHystereticLinks(MODEL);

%Initial conditions and matrices initialization
HIST.R = zeros(size(MODEL.BW.HystLinks,1),1); HIST.Um = HIST.R;
HIST.E = HIST.R; HIST.Zeta = HIST.R;
MODEL.BW.HistBW=HIST;

MODEL.BW.HistR = zeros(size(MODEL.BW.HystLinks,1),MODEL.dyn.nt);
MODEL.BW.HistU = MODEL.BW.HistR;

MODEL.nt=1; MODEL= assemble_nlBW( MODEL );
MODEL.Mall = MODEL.M; 
%Apply numerical correction for singular M
if MODEL.correctM
    test = diag(MODEL.M);
    m=min(test(test>0))/1000;
    test(test==0)=m;
    MODEL.M(1:MODEL.ndim+1:end) = test;
end

MODEL= apply_bc_nl( MODEL );

Kstiff = MODEL.K(MODEL.freedofs,MODEL.freedofs);
Mmass = MODEL.M(MODEL.freedofs,MODEL.freedofs);
if (isfield(InputParams,'zeta')~=0)
    [MODEL.dyn.a,MODEL.dyn.b] = GetRayleighDamping(Kstiff,Mmass,InputParams.zeta,InputParams.OmegaIndexes);
end
MODEL.C = MODEL.dyn.a *  MODEL.M + MODEL.dyn.b * MODEL.K;


nt = MODEL.dyn.nt; neq = numel(MODEL.u);
MODEL.U = zeros(neq,nt); MODEL.V = zeros(neq,nt); MODEL.A = zeros(neq,nt);

MODEL.U(:,1)=InputParams.u0; MODEL.V(:,1)=InputParams.v0;
% MODEL.A(:,1)=InputParams.a0;

% Earthquake forcing: F_eq = -M·ẍ_g (inertia forces opposite to ground acceleration)
[S,F] = AssembleSandZ(InputParams,MODEL);
Fmatrix = -S*F'; 
MODEL.Rmatrix = Fmatrix;
 
[~,pos]=max(F);
MODEL.pos=pos;

end

