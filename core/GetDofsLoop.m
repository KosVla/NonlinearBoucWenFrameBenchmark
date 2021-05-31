function [MODEL] = GetDofsLoop(MODEL)
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
%                   It is a 11-column matrix containing the start and end
%                   node of all hysteretic links and their BW parameters,
%                   following the notation:
%                   [bw_a bw_k Alpha Beta Gamma N Amp deltav deltan]
%
%Please cite as:
% K. Vlachas, K. Tatsis, K. Agathos, A. Brink, and E. Chatzi,
% A local basis approximation approach for nonlinearparametric model order reduction,
% Journal of Sound and Vibration, vol. 502, p. 116055, 2021.

DofsLoopHist=zeros(size(MODEL.nl_link_elements,1)*6,11);
count=1;

for e=1:size(MODEL.nl_link_elements,1)
    link_nodes = MODEL.nl_link_elements(e,1:end);
    link_dofs = element_global_dofs( link_nodes );
    link_flags = MODEL.nl_link_flags(e,:);
    
    check_alt=false;
    if (isfield(MODEL,'nl_links_alternate'))
        [check_alt,check_pos] = ismember(e,MODEL.nl_links_alternate(:,1)); 
    end
    
    for d=1:6
        DofsLoopHist(count,1:2)=[link_dofs(d) link_dofs(d+6)];
        DofsLoopHist(count,3:end)=MODEL.nl_link_bw_properties;
        
        if (link_flags(d)==0)
            %Change parameters of links to force linear behavior if degree
            %of freedom is not activeted (a=1.0, AmpBW=0.0)
            DofsLoopHist(count,3)=1.0;
            DofsLoopHist(count,9)=0.0;
        elseif check_alt 
            %Change parameters of links if it is included in the different
            %Bouc-Wen parameters matrix
            DofsLoopHist(count,3:end)=MODEL.nl_links_alternate(check_pos,2:end);
        end
        count=count+1;
    end
end

MODEL.DofsLoop=DofsLoopHist;

end

