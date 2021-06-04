function [MODEL] = InputFileLinksFirst_LocalizeE(HistDamaged)
%Input file - Assembles the MODEL struct
%Hysteretic links are located in all horizontal beams and in the columns of
%the basement. In the columns the link is located at the top.
%
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. 
%                   For example: Nodal coordinates in MODEL.nodes,
%                   Connectivity in MODEL.beam_elements and so on
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.


%Define input geometry of single frame panel
lx=7.50; ly=5.00; lz=3.2;

% Node definition
% x,y,z coordinates of the nodes are stored in nodes matrix
MODEL.nodes = [0  0  0
               lx 0  0
               2*lx 0 0
               2*lx ly 0
               lx ly 0
               0  ly 0
               0  0  lz
               lx 0  lz
               2*lx 0 lz
               2*lx ly lz
               lx ly lz
               0  ly lz
               0  0  2*lz
               lx 0  2*lz
               2*lx 0 2*lz
               2*lx ly 2*lz
               lx ly 2*lz
               0  ly 2*lz
               0  0  0
               lx 0  0
               2*lx 0 0
               2*lx ly 0
               lx ly 0
               0  ly 0
               0  0  lz
               0  0  lz
               0  0  lz
               lx 0  lz
               lx 0  lz
               lx 0  lz
               lx 0  lz
               2*lx 0 lz
               2*lx 0 lz
               2*lx 0 lz
               2*lx ly lz
               2*lx ly lz
               2*lx ly lz
               lx ly lz
               lx ly lz
               lx ly lz
               lx ly lz
               0  ly lz
               0  ly lz
               0  ly lz
               0  0  2*lz
               0  0  2*lz
               lx 0  2*lz
               lx 0  2*lz
               lx 0  2*lz
               2*lx 0  2*lz
               2*lx 0  2*lz
               2*lx ly 2*lz
               2*lx ly 2*lz
               lx ly 2*lz
               lx ly 2*lz
               lx ly 2*lz
               0  ly 2*lz
               0  ly 2*lz              
               ];

% Definition of beam elements
% The first element of each row corresponds to the beam element type:
% 0 -> normal beam element (the only one provided in this template)
% The second and third elements correspond to the numbers of the first and
% second node of the beam

MODEL.beam_elements = [0 1 27
                       0 2 31
                       0 3 34
                       0 4 37
                       0 5 40
                       0 6 44
                       0 26 28
                       0 30 32
                       0 33 35
                       0 36 41
                       0 39 42
                       0 43 25
                       0 29 38
                       0 7 13
                       0 8 14
                       0 9 15
                       0 10 16
                       0 11 17
                       0 12 18
                       0 46 47
                       0 49 50
                       0 51 52
                       0 53 54
                       0 56 57
                       0 58 45
                       0 48 55];

% Definition of hysteretic links
% The first and second elements correspond to the numbers of the first and
% second node of the link. The nodes are the virual copies assembled in 
%MODEL.nodes           

MODEL.nl_link_elements = [19 1 %Ground links
                          20 2
                          21 3
                          22 4
                          23 5
                          24 6
                          25 7 %First floor links
                          7 26 
                          27 7
                          28 8 
                          29 8 %
                          8 30 %
                          31 8 
                          32 9 %
                          9 33 %
                          34 9
                          35 10 %
                          10 36 %
                          37 10
                          38 11 %
                          11 39
                          40 11
                          41 11 %
                          42 12
                          12 43
                          44 12
                          45 13
                          13 46
                          47 14
                          14 48
                          14 49
                          50 15
                          15 51
                          52 16
                          16 53
                          54 17
                          55 17
                          17 56
                          57 18
                          18 58
                          ];

                    
% Material properties definition
% Each row corresponds to a material
% The first element of each row corresponds to the elastic modulus, second
% to Poisson's ratio, third to density

%Reinforced concrete
% E = 30e09;
% nee = 0.15;
% rho = 2400;

%Steel
E = 210e09;
nee = 0.30;
rho = 8000;

MODEL.material_properties = [E nee rho];


% Definition of the material properties assigned to each beam
% Each row corresponds to an element, the entry of each row corresponds to
% the material assigned to the element

b=size(MODEL.beam_elements,1);
% All beam elements assume the same cross section
MODEL.beam_material_properties = ones(b,1);

%Define which Bouc-Wen links should have different properties
%The first term should indicate the hysteretic link, the next 9 terms the
%BW parameters of the link following the notation:
%[bw_a bw_k Alpha Beta Gamma N Amp deltav deltan]
bw_a= HistDamaged.bw_a; 
Alpha=HistDamaged.Alpha; 
Amp = HistDamaged.AmpBW; N = HistDamaged.N;
deltav = HistDamaged.deltav; deltan=HistDamaged.deltan;
Gamma=HistDamaged.Gamma; Beta=HistDamaged.Beta;
bw_k = HistDamaged.bw_k; 


nl_links_alternate=[11 12 14 15 17 18 20 23];
MODEL.nl_links_alternate =zeros(length(nl_links_alternate),10);
MODEL.nl_links_alternate(:,1)=nl_links_alternate;
MODEL.nl_links_alternate(:,2:end)=repmat(...
    [bw_a bw_k Alpha Beta Gamma N Amp deltav deltan],length(nl_links_alternate),1);

%Dofs of the hysteretic links where the nl links are activated
%This applies to elements with different BW properties as well
a=size(MODEL.nl_link_elements,1);
MODEL.nl_link_flags = ones(a,6);
% MODEL.nl_link_flags(:,4)=0;
%For example
%MODEL.nl_link_flags(:,4) deactivates the 4th degree of freedom of all
%links, whereas MODEL.nl_link_flags(5,4) only for the 5th link.

% Cross section properties definition

%Reinforced concrete cross section
% a1 = 0.4;
% b1 = 0.4;
% 
% A1  = a1*b1*ones(1,1);
% I11 = 0.141*a1*b1^3*ones(1,1);
% I21 = a1*b1^3/12*ones(1,1);
% I31 = b1*a1^3/12*ones(1,1);

%Steel HEA 200 cross section
A1  = 5383*1e-06;
I11 = 204.3*10^3*1e-12;
I21 = 36.92*10^6*1e-12;
I31 = 13.36*10^6*1e-12;


a2 = 0.3;
b2 = 0.6;

A2  = a2*b2*ones(1,1);
I12 = 0.141*a2*b2^3*ones(1,1);
I22 = a2*b2^3/12*ones(1,1);
I32 = b2*a2^3/12*ones(1,1);

% Definition of the different cross section, the first entry corresponds to
% the cross sectional area, the rest of the entries correspond to the
% moments of inertia in the three local directions

MODEL.cross_sections = [A1 I11 I21 I31
                        A2 I12 I22 I32];

% Definition of the cross sections assigned to each beam
% Each row corrsponds to an element, the first entry of each row corresponds to
% the cross section properties assigned to the left end, the second entry
% corresponds to the cross section assigned to the right end

MODEL.beam_cross_sections = ones(b,2);

% Distributed loads definition
% Each row corresponds to an element
% The first three entries in each row correspond to the x,y,z, componets of
% the distributed load at the left end of the beam
% The last three entries in each row correspond to the x,y,z, componets of
% the distributed load at the right end of the beam

p11 = 0;
p12 = 0;
p13 = -1e3;

p21 = 0;
p22 = 0;
p23 = -1e3;

% MODEL.beam_loads = [zeros(4,6);
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23
%                     zeros(4,6)
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23
%                     p11 p12 p13 p21 p22 p23];
                
MODEL.beam_loads = zeros(b,6);


% Nodal loads definition
% Each row corresponds to a node, each column corresponds to the component
% of the load in each direction, moment loads are also possible

p=2*15e3;

MODEL.nodal_loads = [ ];

% MODEL.nodal_loads = [12 p 0 0 0 0 0
%                      5 p 0 0 0 0 0
%                      8 p 0 0 0 0 0
%                      9 p 0 0 0 0 0];
%MODEL.nodal_loads = [3 0 0 p 0 0 0];
% MODEL.nodal_loads = [];


% Nodal displacement definition
% Each row corresponds to a node
% The first six columns are flags to indicate whether a constraint is
% applied in the corresponding dof
% The last six columns are the values of the applied displacements
                 
MODEL.nodal_displacements = [19 1 1 1 1 1 1 0 0 0 0 0 0
                             20 1 1 1 1 1 1 0 0 0 0 0 0
                             21 1 1 1 1 1 1 0 0 0 0 0 0
                             22 1 1 1 1 1 1 0 0 0 0 0 0
                             23 1 1 1 1 1 1 0 0 0 0 0 0
                             24 1 1 1 1 1 1 0 0 0 0 0 0];

% MODEL.nodal_displacements = [];

%Nodal springs definition
% Each row corresponds to a node
% The first entry corresponds to the number of the node.
% The rest entries to the spring coefficient in each degree of freedom
% The same notation applies for lumped masses


% mult = 1;
%                          
% kx = mult*2.6614e+07;
% ky = mult*2.7283e+07;
% kz = mult*3.6472e+07;
% kxx = mult*1.9085e+07;
% kyy = mult*2.9274e+07;
% kzz = mult*2.8611e+07;
%                          
% MODEL.springs = [1 kx ky kz kxx kyy kzz
%                  2 kx ky kz kxx kyy kzz
%                  3 kx ky kz kxx kyy kzz
%                  4 kx ky kz kxx kyy kzz];

MODEL.springs = [];

MODEL.masses = [];

%Definition of plates
% Each row corresponds to a plate element, the entries are the
% corresponding nodes
% MODEL.plate_elements = [7 8 11 12];

% Plate material properties
% Same notation as beam elements [E, nee, rho]
% Each row corresponds to a material
% The first element of each row corresponds to the elastic modulus, second
% to Poisson's ratio, third to density

% MODEL.plate_material_properties = 2*ones(size(MODEL.plate_elements,1));
% MODEL.plate_thickness=0.2*ones(2,4);


%Dynamic analysis parameters

MODEL.dyn.dt = 0.01; %time step
MODEL.dyn.nt = 3050; %duration of the analysis
MODEL.dyn.a = 0; %Rayleigh damping coefficient alpha - Mass
MODEL.dyn.b = 0.05; %Rayleigh damping coefficient beta - Stiffness

end

