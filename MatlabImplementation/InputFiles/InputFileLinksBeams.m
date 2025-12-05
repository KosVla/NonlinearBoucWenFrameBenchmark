function [MODEL] = InputFileLinksBeams()
% Assembles system struct for a two-story frame with hysteretic links in
% in all horizontal beams
%
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. 
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

%Define input geometry of single frame panel
lx=7.50; ly=5.00; lz=3.2;

% Node definition
% x,y,z coordinates of the nodes are stored in nodes matrix
MODEL.nodes = [0  0  0      %1
               lx 0  0      %2
               2*lx 0 0     %3
               2*lx ly 0    %4
               lx ly 0      %5
               0  ly 0      %6
               0  0  lz     %7
               lx 0  lz     %8
               2*lx 0 lz    %9
               2*lx ly lz   %10
               lx ly lz     %11
               0  ly lz     %12
               0  0  2*lz   %13
               lx 0  2*lz   %14
               2*lx 0 2*lz  %15
               2*lx ly 2*lz %16
               lx ly 2*lz   %17
               0  ly 2*lz   %18
               0  0  0      %19
               lx 0  0      %20
               2*lx 0 0     %21
               2*lx ly 0    %22
               lx ly 0      %23
               0  ly 0      %24
               0  0  lz     %25
               0  0  lz     %26
               lx 0  lz     %27
               lx 0  lz     %28
               lx 0  lz     %29
               2*lx 0 lz    %30
               2*lx 0 lz    %31
               2*lx ly lz   %32
               2*lx ly lz   %33
               lx ly lz     %34
               lx ly lz     %35
               lx ly lz     %36
               0  ly lz     %37
               0  ly lz     %38
               0  0  2*lz   %39
               0  0  2*lz   %40
               lx 0  2*lz   %41
               lx 0  2*lz   %42
               lx 0  2*lz   %43
               2*lx 0  2*lz %44
               2*lx 0  2*lz %45
               2*lx ly 2*lz %46
               2*lx ly 2*lz %47
               lx ly 2*lz   %48
               lx ly 2*lz   %49
               lx ly 2*lz   %50
               0  ly 2*lz   %51
               0  ly 2*lz   %52      
               ];     

% Definition of beam elements
% The first element of each row corresponds to the beam element type:
% 0 -> normal beam element (the only one provided in this template)
% The second and third elements correspond to the numbers of the first and
% second node of the beam

MODEL.beam_elements = [0 1 7
                       0 2 8
                       0 3 9
                       0 4 10
                       0 5 11
                       0 6 12
                       0 7 13
                       0 8 14
                       0 9 15
                       0 10 16
                       0 11 17
                       0 12 18
                       0 38 25
                       0 26 27
                       0 29 30
                       0 28 35
                       0 31 32
                       0 33 34
                       0 36 37
                       0 52 39
                       0 40 41
                       0 43 44
                       0 42 49
                       0 45 46
                       0 47 48
                       0 50 51];

% Definition of hysteretic links
% The first and second elements correspond to the numbers of the first and
% second node of the link. The nodes are the virual copies assembled in 
%MODEL.nodes           

MODEL.nl_link_elements = [19 1
                          20 2
                          21 3
                          22 4
                          23 5
                          24 6
                          25 7
                          7 26
                          27 8
                          8 28
                          8 29
                          30 9
                          9 31
                          32 10
                          10 33
                          34 11
                          35 11
                          11 36
                          37 12
                          12 38
                          39 13
                          13 40
                          41 14
                          14 43
                          14 42
                          44 15
                          15 45
                          46 16
                          16 47
                          48 17
                          49 17
                          17 50
                          51 18
                          18 52                          
                          ];
            
% Material properties definition
% Each row corresponds to a material
% The first element of each row corresponds to the elastic modulus, second
% to Poisson's ratio, third to density

%Steel
E = 210e09; nee = 0.30; rho = 8000;

%Reinforced concrete
Ec = 30e09; neec = 0.15; rhoc = 2400;

MODEL.material_properties = [E nee rho;
                             Ec neec rhoc];

% Definition of the material properties assigned to each beam
% Each row corrsponds to an element, the entry of each row corresponds to
% the material assigned to the element

b=size(MODEL.beam_elements,1);
% All beam elements assume the same cross section
MODEL.beam_material_properties = ones(b,1);

%Define which Bouc-Wen links should have different properties
%The first term should indicate the hysteretic link, the next 9 terms the
%BW parameters of the link following the notation:
%[bw_a bw_k Alpha Beta Gamma N Amp deltav deltan]
% MODEL.nl_links_alternate = [5 0.4 2e8 1.0 3.0 2.0 1.0 1.0 0 0];

%Dofs of the hysteretic links where the nl links are activated
%This applies to elements with different BW properties as well
a=size(MODEL.nl_link_elements,1);
MODEL.nl_link_flags = ones(a,6);
% MODEL.nl_link_flags(:,4)=0;
%For example
%MODEL.nl_link_flags(:,4) deactivates the 4th degree of freedom of all
%links, whereas MODEL.nl_link_flags(5,4) only for the 5th link.

% Cross section properties definition
%Steel HEA 200 cross section
A1  = 5383*1e-06; I11 = 204.3*10^3*1e-12;
I21 = 36.92*10^6*1e-12; I31 = 13.36*10^6*1e-12;

a2 = 0.3; b2 = 0.6;
A2  = a2*b2*ones(1,1); I12 = 0.141*a2*b2^3*ones(1,1);
I22 = a2*b2^3/12*ones(1,1); I32 = b2*a2^3/12*ones(1,1);

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

% p11 = 0; p12 = 0; p13 = -1e3;
% 
% p21 = 0; p22 = 0; p23 = -1e3;
% 
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
MODEL.nodal_loads = [ ];

% p=2*15e3;
% MODEL.nodal_loads = [12 p 0 0 0 0 0
%                      5 p 0 0 0 0 0
%                      8 p 0 0 0 0 0
%                      9 p 0 0 0 0 0];

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
MODEL.dyn.a = 0; %Rayleigh damping coefficient alpha
MODEL.dyn.b = 0.05; %Rayleigh damping coefficient beta
MODEL.ndim=6*numel(MODEL.nodes(:,1)); 
MODEL.correctM = false;

end

