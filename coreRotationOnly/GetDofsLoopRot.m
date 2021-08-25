function [MODEL] = GetDofsLoopRot(MODEL)
%Assembles the pairs of degrees of freedom  and the properties
% of all hysteretic links in a single matrix
%
%Input:
%  MODEL : struct / The MODEL struct contains all system properties,
%                   parameters and matrices. Here, the geometry 
%                   configuration is used to assemble the matrix 
%                   of dof pairs of the links
%Returns:
%  MODEL : struct / The MODEL struct contains all system properties, 
%                   parameters and matrices. The DofsLoop matrix is 
%                   assembled and stored in the MODEL struct. 
%                   It is a 13-column matrix containing the start and end
%                   node of all hysteretic links and their BW parameters,
%                   following the notation:
%                   [bw_a bw_k Alpha Beta Gamma N Amp deltav deltan length EI coeff]
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

DofsLoopHist=zeros(size(MODEL.nl_link_elements,1)*6,13);
count=1;

for e=1:size(MODEL.nl_link_elements,1)
    link_nodes = MODEL.nl_link_elements(e,1:end);
    link_dofs = element_global_dofs( link_nodes );
    link_flags = MODEL.nl_link_flags(e,:);
    
    beamelement = FindEquivalentBeam(link_nodes, MODEL.beam_elements);
    material_properties = MODEL.material_properties(MODEL.beam_material_properties(beamelement),:);
    cross_section_properties = MODEL.cross_sections(MODEL.beam_cross_sections(beamelement,:)',:);          
    beam_nodes = MODEL.beam_elements(beamelement,2:3);
    beam_COORDS = MODEL.nodes(beam_nodes,:);
    
    L = norm(beam_COORDS(2,:)-beam_COORDS(1,:));
    E   = material_properties(1);
    nee = material_properties(2);
    A   = cross_section_properties(1,1);
    I1  = cross_section_properties(1,2); %Ix
    I2  = cross_section_properties(1,3); %Iy
    I3  = cross_section_properties(1,4); %Iz
    
    check_alt=false;
    if (isfield(MODEL,'nl_links_alternate'))
        [check_alt,check_pos] = ismember(e,MODEL.nl_links_alternate(:,1)); 
    end
    
    for d=1:6
        DofsLoopHist(count,1:2)=[link_dofs(d) link_dofs(d+6)];
        %[linkstart linkend bw_a bw_k Alpha Beta Gamma N Amp deltav deltan]
        DofsLoopHist(count,3:11)=MODEL.nl_link_bw_properties;
        DofsLoopHist(count,12) = L;
        
        if d==1
            DofsLoopHist(count,end)= E*A;     
        elseif d==2 || d==6
            DofsLoopHist(count,end)= E*I3;              
        elseif d==3 || d==5
            DofsLoopHist(count,end)= E*I2;              
        elseif d==4
            DofsLoopHist(count,end)= E/(2*(1+nee))*I1;                         
        end
        
        if (link_flags(d)==0)
            %Change parameters of links to force not activation (a=0.0, AmpBW=0.0)
            %This implies that the link will have a constant bw_k equal to the max
            %of the K diagonal
            DofsLoopHist(count,3)=1; DofsLoopHist(count,9)=0.0;
        elseif check_alt 
            %Change parameters of links if it is included in the different
            %Bouc-Wen parameters matrix
            
            DofsLoopHist(count,3:11)=MODEL.nl_links_alternate(check_pos,2:end);
        end
        count=count+1;
    end
end

MODEL.DofsLoop=DofsLoopHist;

end

