function [ MODEL ] = assemble_nlBW_SDOF( MODEL )
%Assembles the updated system matrices based on the displacement of the 
%   current timestep
%
%Input:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. Here the MODEL.u represents 
%                   the displacement of the current timestep and is used
%       			as an input variable to update the internal forces and
%                   compute the Bouc-Wen hysteretic forces of the links.
%                   The Bouc-Wen variables of the previous step are located
%                   in MODEL.HistBW. Additional properties stored in the 
%                   MODEL struct are also used.
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties,
%                   parameters and matrices. The stiffness matrix and the 
%                   internal forces vector are updated, along with the
%       			respective history variables of the Bouc-Wen model 
%                   (zeta, x etc -> MODEL.HistBW)
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

ndof = 6*numel(MODEL.nodes(:,1));

nt = MODEL.nt;
MODEL.K = sparse(ndof,ndof);
MODEL.M = sparse(ndof,ndof);
MODEL.fint = zeros(ndof,1);
MODEL.f = zeros(ndof,1);

condition=false;
if (isfield(MODEL,'u')==0)
    MODEL.u = zeros(ndof,1);
    condition=true;
end

% Assembly of the beam elements

n_beams =size(MODEL.beam_elements,1);

for e=1:n_beams
    if condition
        %beam_type = MODEL.beam_elements(e,1);
        beam_nodes = MODEL.beam_elements(e,2:end);
        beam_COORDS = MODEL.nodes(beam_nodes,:);
        beam_material_properties = MODEL.material_properties(MODEL.beam_material_properties(e),:);
        beam_cross_section_properties = MODEL.cross_sections(MODEL.beam_cross_sections(e,:)',:);
        beam_loads = MODEL.beam_loads(e,:);

        [ beam_dofs ] = element_global_dofs( beam_nodes );
        [Ke,Me,fe] = beam_mass_stiffness_rhs(beam_COORDS, beam_material_properties, beam_cross_section_properties, beam_loads);
        MODEL.ElementsK{e}=Ke;
        MODEL.ElementsDofs{e}=beam_dofs;
        MODEL.K(beam_dofs,beam_dofs) = MODEL.K(beam_dofs,beam_dofs) + Ke;
        MODEL.M(beam_dofs,beam_dofs) = MODEL.M(beam_dofs,beam_dofs) + Me;
        MODEL.f(beam_dofs) = MODEL.f(beam_dofs) + fe;
    else
        beam_dofs=MODEL.ElementsDofs{e};
        Ke = MODEL.ElementsK{e};
        MODEL.K(beam_dofs,beam_dofs) = MODEL.K(beam_dofs,beam_dofs) + Ke;
        MODEL.M=MODEL.Mbeam;
        MODEL.f=MODEL.fbeam;
    end
    
    ue = MODEL.u(beam_dofs);
    MODEL.fint(beam_dofs) = MODEL.fint(beam_dofs) + Ke*ue;
end

if condition
    MODEL.Kbeam=MODEL.K;
    MODEL.Mbeam=MODEL.M;
    MODEL.fbeam=MODEL.f;
end
 
% Assembly of the plate elements

if (isfield(MODEL,'plate_elements'))

    n_plates = numel(MODEL.plate_elements(:,1));
else
    n_plates = 0;
end

for e=1:n_plates
    plate_nodes = MODEL.plate_elements(e,1:end);
    plate_COORDS = MODEL.nodes(plate_nodes,:);
    plate_material_properties = MODEL.material_properties(MODEL.plate_material_properties(e),:);
    plate_thickness = MODEL.plate_thickness(e,:);
    plate_loads = MODEL.plate_loads(e,:);
    
    [Ke,Me,fe] = shell_mass_stiffness_rhs( plate_COORDS, plate_material_properties, plate_thickness, plate_loads);
    [ plate_dofs ] = element_global_dofs( plate_nodes );
    ue = MODEL.u(plate_dofs);

    MODEL.K(plate_dofs,plate_dofs) = MODEL.K(plate_dofs,plate_dofs) + Ke;
    MODEL.M(plate_dofs,plate_dofs) = MODEL.M(plate_dofs,plate_dofs) + Me;
    MODEL.fint(plate_dofs) = MODEL.fint(plate_dofs) + Ke*ue;
    MODEL.f(plate_dofs) = MODEL.f(plate_dofs) + fe;
end

% Assembly of nonlinear links

if (isfield(MODEL,'nl_link_elements'))
    n_links = numel(MODEL.nl_link_elements(:,1));
%     ks = max(max(diag(MODEL.K)))*1e3; %Coefficient for dofs of hysteretic links not activated
    if (~isfield(MODEL,'nl_link_hist'))
        MODEL.nl_link_hist = cell(n_links,1);
        MODEL.nl_link_hist(:) ={struct('uj',{0,0,0,0,0,0},'rj',{0,0,0,0,0,0},'zeta',{0,0,0,0,0,0},'E',{0,0,0,0,0,0})};
    end
else
    n_links = 0;
end

for e=1:n_links
    
    link_nodes = MODEL.nl_link_elements(e,1:end);
    link_properties = MODEL.DofsLoop((e-1)*6+1:e*6,3:end);
    link_hist = MODEL.nl_link_hist{e,:};
    link_dofs= element_global_dofs( link_nodes );
    
    ue = MODEL.u(link_dofs);
    
    [finte, Ke, link_hist] = link_residual_stiffness(ue, link_properties, link_hist);
    
    MODEL.nl_link_hist{e,:} = link_hist;
    MODEL.HistR((e-1)*6+1:e*6,nt)=finte(1:6);
    MODEL.HistU((e-1)*6+1:e*6,nt)=ue(1:6)-ue(7:12);

    MODEL.K(link_dofs,link_dofs) = MODEL.K(link_dofs,link_dofs) + Ke;
    MODEL.fint(link_dofs) = MODEL.fint(link_dofs) + finte;
end

% Assembly of the springs

if (isfield(MODEL,'springs')&&(numel(MODEL.springs)>0))

    n_springs = numel(MODEL.springs(:,1));
else
    n_springs = 0;
end

for s=1:n_springs
    node = MODEL.springs(s,1);
    spring_constants = MODEL.springs(s,2:end);
    [ dofs ] = node_global_dofs( node );
    MODEL.K(dofs,dofs) = MODEL.K(dofs,dofs) + diag(eye(numel(dofs))*spring_constants');
end

% Assembly of masses

if (isfield(MODEL,'masses')&&(numel(MODEL.masses)>0))

    n_masses = numel(MODEL.masses(:,1));
else
    n_masses = 0;
end

for m=1:n_masses
    node = MODEL.masses(m,1);
    mass = MODEL.masses(m,2:end);
    [ dofs ] = node_global_dofs( node );
    MODEL.M(dofs,dofs) = MODEL.M(dofs,dofs) + diag(eye(numel(dofs))*mass');
end