function [ r,k, hist ] = bw_1d( u, hist, par )
%Evaluates the Bouc-Wen hysteretic forcing of a single link
%
%Input:
%  u : vector / Displacements on current timestep 
%  hist: struct / Struct containing the values of the Bouc-Wen parameters 
%                   at the previous timestep (e.g. zeta, Energy etc.)
%  par: vector / Properties of the Bouc-Wen links
%
%Returns:
%  r : vector / Updated restoring forces
%  k : matrix / Stiffness terms of the hysteretic link
%  hist: struct / Updated history struct for the Bouc-Wen model
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

%Compute current stiffness and restoring force
%Incremental displacement
dxi = u - hist.uj;

%Strength deterioration
Energy = hist.E;
ve = 1 + par.deltav*Energy;
%Stiffness degradation
ne = 1 + par.deltan*Energy;

%Stiffness terms
zetaprior = hist.zeta;
kbw = (par.Alpha - ve*(par.Beta*sign(zetaprior)*sign(dxi) + par.Gamma) * abs(zetaprior)^(par.N))/ne;

%Force terms
LinearTerm = par.a*par.k*dxi;
dzeta = kbw*dxi;
HystereticTerm = (1-par.a)*par.k*dzeta;
k = par.a*par.k + (1-par.a)*par.k*kbw;

%Restoring force
r = hist.rj + par.Amp*HystereticTerm + LinearTerm;

%Update history

hist.rj = r;
hist.uj = hist.uj + dxi;
hist.zeta = zetaprior+dzeta;
hist.E = hist.E+ dzeta*dxi;
end