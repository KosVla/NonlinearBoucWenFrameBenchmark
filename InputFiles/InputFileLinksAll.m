function [MODEL] = InputFileLinksAll()
%Input file - Assembles the MODEL struct
%Hysteretic links are located in all horizontal beams and columns of the
%two-story frame. In the columns there is a link at the bottom and one at
%the top.
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
MODEL.nodes = [0  0  0          %1
               lx 0  0          %2
               2*lx 0 0         %3
               2*lx ly 0        %4
               lx ly 0          %5
               0  ly 0          %6
               0  0  lz         %7
               lx 0  lz         %8
               2*lx 0 lz        %9
               2*lx ly lz       %10
               lx ly lz         %11
               0  ly lz         %12
               0  0  2*lz       %13
               lx 0  2*lz       %14
               2*lx 0 2*lz      %15
               2*lx ly 2*lz     %16
               lx ly 2*lz       %17
               0  ly 2*lz       %18 - End of "real" nodes
               0  0  0          %19 - Start of hysteretic links - Ground
               lx 0  0          %20
               2*lx 0 0         %21
               2*lx ly 0        %22
               lx ly 0          %23
               0  ly 0          %24 - End of ground links
               0  0  lz         %25(-7) Links at columns - Basement
               lx 0  lz         %26(-8)
               2*lx 0 lz        %27(-9)
               2*lx ly lz       %28(-10)
               lx ly lz         %29(-11)
               0  ly lz         %30(-12)
               0  0  lz         %(7-)31 Links at the base of the first floor columns
               lx 0  lz         %(8-)32
               2*lx 0 lz        %(9-)33
               2*lx ly lz       %(10-)34
               lx ly lz         %(11-)35
               0  ly lz         %(12-)36
               0  0  2*lz       %37(-13) Links at the top of the first floor columns
               lx 0  2*lz       %38(-14)
               2*lx 0 2*lz      %39(-15)
               2*lx ly 2*lz     %40(-16)
               lx ly 2*lz       %41(-17)
               0  ly 2*lz       %42(-18)
               0  0  lz         %43 Links at the lx-beams of the 1st floor
               lx 0  lz         %44
               lx 0 lz          %45
               2*lx 0 lz        %46
               2*lx ly lz       %47
               lx  ly lz        %48
               lx  ly lz        %49
               0  ly lz         %50
               0  0  2*lz       %51 Links at the lx-beams of the 2nd floor
               lx 0  2*lz       %52
               lx 0 2*lz        %53
               2*lx 0 2*lz      %54
               2*lx ly 2*lz     %55
               lx  ly 2*lz      %56
               lx  ly 2*lz      %57
               0  ly 2*lz       %58
               2*lx  0  lz      %59 Links at the ly-beams of the 1st floor
               2*lx ly  lz      %60
               lx 0 lz          %61
               lx ly lz         %62
               0 ly lz          %63
               0 0 lz           %64
               2*lx  0  2*lz    %65 Links at the ly-beams of the 2nd floor
               2*lx ly  2*lz    %66
               lx 0 2*lz        %67
               lx ly 2*lz       %68
               0 ly 2*lz        %69
               0 0 2*lz         %70
               ];

% Definition of beam elements
% The first element of each row corresponds to the beam element type:
% 0 -> normal beam element (the only one provided in this template)
% The second and third elements correspond to the numbers of the first and
% second node of the beam

MODEL.beam_elements = [0 1 25 %Basement columns
                       0 2 26
                       0 3 27
                       0 4 28
                       0 5 29
                       0 6 30
                       0 31 37 %First floor columns
                       0 32 38
                       0 33 39
                       0 34 40
                       0 35 41
                       0 36 42                       
                       0 43 44 %Start of frame lx elements - 1st floor
                       0 45 46
                       0 47 48
                       0 49 50
                       0 51 52 %Start of frame lx elements - 2nd floor
                       0 53 54
                       0 55 56
                       0 57 58
                       0 59 60 %Start of frame ly elements - 1st floor
                       0 61 62
                       0 63 64
                       0 65 66 %Start of frame ly elements - 2nd floor
                       0 67 68
                       0 69 70
                       ];

% Definition of hysteretic links
% The first and second elements correspond to the numbers of the first and
% second node of the link. The nodes are the virual copies assembled in 
%MODEL.nodes           

MODEL.nl_link_elements = [19 1  %Start of hysteretic links - Ground
                          20 2
                          21 3
                          22 4
                          23 5
                          24 6
                          25 7  %Basement column links 
                          26 8
                          27 9
                          28 10
                          29 11
                          30 12
                          7 31  %Links at the basis of the first floor columns
                          8 32
                          9 33
                          10 34
                          11 35
                          12 36
                          37 13 %Links at the top of the first floor columns
                          38 14
                          39 15
                          40 16
                          41 17
                          42 18
                          7 43 %Links at the lx-beams - 1st floor
                          44 8
                          8 45
                          46 9
                          10 47
                          48 11
                          11 49
                          50 12
                          13 51 %Links at the lx-beams - 2nd floor
                          52 14
                          14 53
                          54 15
                          16 55
                          56 17
                          17 57
                          58 18
                          9 59 %Links at the ly-beams - 1st floor
                          60 10
                          8 61
                          62 11
                          12 63
                          64 7
                          15 65 %Links at the ly-beams - 2nd floor
                          66 16
                          14 67
                          68 17
                          18 69
                          70 13
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

