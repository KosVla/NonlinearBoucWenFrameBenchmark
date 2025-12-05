function [ r,k, HIST ] = bw_RK4( U, DofsLoop, HIST, link_properties )
%Evaluates the Bouc-Wen hysteretic forcing using 4th-order Runge-Kutta
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
%  k : matrix / Stiffness terms of the hysteretic links (at final state)
%  HIST: struct / Updated history struct for the Bouc-Wen model
%
%Integration scheme: Classical RK4 (4th order accurate)

a = link_properties(:,1);k = link_properties(:,2); 
Alpha=link_properties(:,3); Beta=link_properties(:,4);
Gamma=link_properties(:,5); N=link_properties(:,6);
deltav=link_properties(:,7); deltan=link_properties(:,8);

%Incremental displacement
Um = U(DofsLoop(:,1))-U(DofsLoop(:,2));
Dxi = Um - HIST.Um;

%Strength deterioration
Energy = HIST.E;
ve = 1 + deltav.*Energy;
%Stiffness degradation
ne = 1 + deltan.*Energy;

%RK4 Integration - Classical 4th order method
zetaprior = HIST.Zeta;

% Stage 1: Evaluate at z_n
kbw1 = compute_kbw(zetaprior, Dxi, Alpha, Beta, Gamma, N, ve, ne);
k1 = kbw1.*Dxi;

% Stage 2: Evaluate at z_n + k1/2
zeta2 = zetaprior + k1/2;
kbw2 = compute_kbw(zeta2, Dxi, Alpha, Beta, Gamma, N, ve, ne);
k2 = kbw2.*Dxi;

% Stage 3: Evaluate at z_n + k2/2
zeta3 = zetaprior + k2/2;
kbw3 = compute_kbw(zeta3, Dxi, Alpha, Beta, Gamma, N, ve, ne);
k3 = kbw3.*Dxi;

% Stage 4: Evaluate at z_n + k3
zeta4 = zetaprior + k3;
kbw4 = compute_kbw(zeta4, Dxi, Alpha, Beta, Gamma, N, ve, ne);
k4 = kbw4.*Dxi;

% Weighted average (RK4 formula)
dzeta = (k1 + 2*k2 + 2*k3 + k4) / 6;

%Force terms
LinearTerm = a.*k.*Dxi;
HystereticTerm = (1-a).*k.*dzeta;

%Restoring force
r = HIST.R + HystereticTerm + LinearTerm;

%Tangent stiffness - evaluate at final state
zeta_final = zetaprior + dzeta;
kbw_final = compute_kbw(zeta_final, Dxi, Alpha, Beta, Gamma, N, ve, ne);
k = a.*k + (1-a).*k.*kbw_final;

%Update history
HIST.R = r;
HIST.Um = HIST.Um + Dxi;
HIST.Zeta = zeta_final;
%Energy update - Trapezoidal rule
HIST.E = HIST.E+ (zetaprior + HIST.Zeta)/2 .* Dxi;

end