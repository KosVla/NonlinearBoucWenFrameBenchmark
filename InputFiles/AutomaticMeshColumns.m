function [MODEL] = AutomaticMeshColumns(...
    NoOfFloors,NoOfFramesx,NoOfFramesy,dimensions)
%UNDER DEVELOPMENT
%Input file - Assembles the MODEL struct
%Hysteretic links are located in all columns. In the 
%columns there is a link at the bottom and one at the top.
%
%Inputs:
%   NoOfFloors: int / Number of stories for the frame (excluding the basement, so additional to the basement)
%   NoOfFramesx: int / Number of frames in the x-direction
%   NoOfFramesy: int /Number of frames in the y-direction
%   dimensions: vector / [lx, ly, lz] dimensions of the frame
%                           lx: Length of each frame in x-direction
%                           ly: Length of each frame in y-direction
%                           lz: Height of each floor in z
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

lx = dimensions(1);
ly = dimensions(2);
lz = dimensions(3);

%Initialize
Nodes = zeros((NoOfFloors+2)*(NoOfFramesx+1)*(NoOfFramesy+1),3);

count=1;
for i=1:(NoOfFramesy+1)
    for k=1:(NoOfFramesx+1)
        Nodes(count,:) = [(k-1)*lx (i-1)*ly 0]; 
        count=count+1;
    end
end

for r = 1:(NoOfFloors+1)
    Nodes(r*(count-1)+1:(r+1)*(count-1),:) = Nodes(1:count-1,:);
    Nodes(r*(count-1)+1:(r+1)*(count-1),3) = lz*r;
end

planenodes = (NoOfFramesx+1)*(NoOfFramesy+1);

%Define hysteretic links - Ground
k = size(Nodes,1);
Nodes(k+1:k+planenodes,:) = Nodes(1:planenodes,:);
NodeNumbersLinkGrounds = k+1:k+planenodes;

%Links at columns - Basement
k = size(Nodes,1);
Nodes(k+1:k+planenodes,:) = Nodes(1:planenodes,:);
Nodes(k+1:k+planenodes,3) = lz;
NodeNumbersLinkBasementColumns = k+1:k+planenodes;

NodeNumbersBasis=cell(NoOfFloors,1);
NodeNumbersTop=cell(NoOfFloors,1);

%Links at columns - Floors
for r = 1:NoOfFloors
    %Basis of column
    k = size(Nodes,1);
    Nodes(k+1:k+planenodes,:) = Nodes(1:planenodes,:);
    Nodes(k+1:k+planenodes,3) = r*lz;
    NodeNumbersBasis{r} = k+1:k+planenodes;
    %Top of column
    k = size(Nodes,1);
    Nodes(k+1:k+planenodes,:) = Nodes(1:planenodes,:);
    Nodes(k+1:k+planenodes,3) = (r+1)*lz;
    NodeNumbersTop{r} = k+1:k+planenodes;
end

%Beam connectivity

%Basement columns
BeamElements = [(1:planenodes)' NodeNumbersLinkBasementColumns'];

%Floor columns
for r = 1:NoOfFloors
    k = size(BeamElements,1);
    BeamElements(k+1:k+planenodes,:) = [NodeNumbersBasis{r}' NodeNumbersTop{r}'];
end

%Lx beams
for r=1:(NoOfFloors+1)
    corner = (NoOfFramesx+1)*(NoOfFramesy+1)+1+(r-1)*planenodes;
    for m=1:(NoOfFramesy+1)
        startn = corner + (m-1)*(NoOfFramesx+1) ;
        for t=1:NoOfFramesx
            k = size(BeamElements,1);
            if m==1
                BeamElements(k+1,:) =[startn+(t-1) startn+t];
            else
                BeamElements(k+1,:) =[startn+t startn+(t-1)];
            end
        end
    end
end

%Ly beams
for r=1:(NoOfFloors+1)
    corner = (NoOfFramesx+1)*(NoOfFramesy+1)+1+(r-1)*planenodes;
    for m=1:(NoOfFramesx+1)
        startn = corner + (m-1);
        for t=1:NoOfFramesy
            k = size(BeamElements,1);
            if m==1
                BeamElements(k+1,:) =...
                    [startn+(NoOfFramesx+1)*t startn+(NoOfFramesx+1)*(t-1)];
            else
                BeamElements(k+1,:) =...
                    [startn+(NoOfFramesx+1)*(t-1) startn+(NoOfFramesx+1)*t];
            end
        end
    end
end

%Nonlinear links connectivity
%Start of hysteretic links - Ground
nl_link_elements = [NodeNumbersLinkGrounds' (1:planenodes)'];

%Basement column links 
k = size(nl_link_elements,1);
nl_link_elements(k+1:k+planenodes,:)=...
    [NodeNumbersLinkBasementColumns' (planenodes+1:2*planenodes)'];

%Links at the basis and the top of the floor columns
for r = 1:(NoOfFloors)
    k = size(nl_link_elements,1);
    nl_link_elements(k+1:k+planenodes,:) =...
        [(r*planenodes+1:(r+1)*planenodes)' NodeNumbersBasis{r}'];
    nl_link_elements(k+planenodes+1:k+2*planenodes,:) =...
        [NodeNumbersTop{r}' ((r+1)*planenodes+1:(r+2)*planenodes)'];
end

MODEL.nodes=Nodes;
MODEL.beam_elements = [zeros(length(BeamElements),1) BeamElements];
MODEL.nl_link_elements = nl_link_elements;

%Boundary conditions
nodal_displacements = zeros(length(NodeNumbersLinkGrounds), 13);
nodal_displacements(:,1) = NodeNumbersLinkGrounds;
nodal_displacements(:,2:7)=1;
nodal_displacementsz = zeros(planenodes, 13);
nodal_displacementsz(:,1) = 1:planenodes;
nodal_displacementsz(:,2:4)=1;
MODEL.nodal_displacements = [nodal_displacementsz; nodal_displacements];
                         
%Steel
E = 210e09; nee = 0.30; rho = 8000;
MODEL.material_properties = [E nee rho];


b=size(MODEL.beam_elements,1);
% All beam elements assume the same cross section
MODEL.beam_material_properties = ones(b,1);

%Dofs of the hysteretic links where the nl links are activated
%This applies to elements with different BW properties as well
a=size(MODEL.nl_link_elements,1);
MODEL.nl_link_flags = zeros(a,6);
MODEL.nl_link_flags(:,[5 6])=1;

%Steel HEA 200 cross section
A1  = 5383*1e-06;
I11 = 204.3*10^3*1e-12;
I21 = 36.92*10^6*1e-12;
I31 = 13.36*10^6*1e-12;

% Definition of the different cross section, the first entry corresponds to
% the cross sectional area, the rest of the entries correspond to the
% moments of inertia in the three local directions

MODEL.cross_sections = [A1 I11 I21 I31];

MODEL.beam_cross_sections = ones(b,2);
MODEL.beam_loads = zeros(b,6);              
MODEL.springs = [];

MODEL.masses = [];
MODEL.nodal_loads = [ ];

%Definition of plates
% Each row corresponds to a plate element, the entries are the
% corresponding nodes
% MODEL.plate_elements = [7 8 11 12];


%Dynamic analysis parameters
MODEL.dyn.dt = 0.01; %time step
MODEL.dyn.nt = 3050; %duration of the analysis
MODEL.dyn.a = 0; %Rayleigh damping coefficient alpha - Mass
MODEL.dyn.b = 0.05; %Rayleigh damping coefficient beta - Stiffness

end

