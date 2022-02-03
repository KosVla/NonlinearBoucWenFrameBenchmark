function [ r,k, HIST ] = bw_Ndall( U, DofsLoop, HIST, link_properties )
%Evaluates the Bouc-Wen hysteretic forcing for all nonlinear links
%
%Input:
%  U : vector / Displacements on current timestep on dofs of 
%               all hysteretic links (activated or not)
%  DofsLoop: matrix / Properties of the models of the hysteretic links.
%                     Each row corresponds to a degree of freedom, each
%                     link has 6 dofs. The first two columns of each row
%                     correspond to the degrees of freedom of links.
%                     First column is the starting dof, second column 
%                     is the end dof of the hysteretic link. 
%  HIST: struct / Struct containing the values of the Bouc-Wen parameters
%                   at the previous timestep (e.g. zeta, Energy etc.)
%  link_properties: matrix / Properties of the models of the hysteretic 
%                           links. Each row corresponds to a degree of 
%                           freedom, each link has 6 dofs. The columns 
%                           contain the Bouc-Wen parameters, following
%                           the notation on the GetDofsLoop function:
%                       [bw_a bw_k Alpha Beta Gamma N Amp deltav deltan]
%Returns:
%  r : vector / Updated restoring forces
%  k : matrix / Stiffness terms of the hysteretic links
%  HIST: struct / Updated history struct for the Bouc-Wen model
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

Amp = link_properties(:,7); a = link_properties(:,1);
k = link_properties(:,2); Alpha=link_properties(:,3);
Beta=link_properties(:,4); Gamma=link_properties(:,5);
N=link_properties(:,6); deltav=link_properties(:,8); 
deltan=link_properties(:,9); 

%Incremental displacement
Um = U(DofsLoop(:,1))-U(DofsLoop(:,2));
Dxi = Um - HIST.Um;

%Strength deterioration
Energy = HIST.E;
ve = 1 + deltav.*Energy;
%Stiffness degradation
ne = 1 + deltan.*Energy;

%Stiffness terms
zetaprior = HIST.Zeta; 
kbw = (Alpha - ve.*(Beta.*sign(zetaprior).*sign(Dxi) + Gamma).* abs(zetaprior).^(N))./ne;

%Force terms
LinearTerm = a.*k.*Dxi;
dzeta = kbw.*Dxi;
HystereticTerm = (1-a).*k.*dzeta;
k = a.*k + (1-a).*k.*kbw;

%Restoring force
r = HIST.R + Amp.*HystereticTerm + LinearTerm;

%Account for nonactivated links
r(a==1) = link_properties(a==1,end).*(U(DofsLoop(a==1,1))-U(DofsLoop(a==1,2))); 
k(a==1) = link_properties(a==1,end);

%Update history
HIST.R = r;
HIST.Um = HIST.Um + Dxi;
HIST.Zeta = zetaprior+dzeta;
HIST.E = HIST.E+ zetaprior.*Dxi;

end