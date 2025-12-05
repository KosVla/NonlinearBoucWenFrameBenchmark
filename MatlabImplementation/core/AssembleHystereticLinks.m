function HystLinks = AssembleHystereticLinks(MODEL)
%Assembles the pairs of degrees of freedom  and the properties
% of all hysteretic links in a single matrix
%Each row corresponds to a degree of freedom hysteretic link
%Here an individual degree of freedom properties can be changed by changing the
%respective row.
%The elements contained in MODEL.nl_links_alternate have been treated and
%the different BW parameters are saved on the activated dofs
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

HystLinks=zeros(size(MODEL.nl_link_elements,1)*6,10);
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
        HystLinks(count,1:2)=[link_dofs(d) link_dofs(d+6)];
        HystLinks(count,3:end)=MODEL.nl_link_bw_properties;
        
        if (link_flags(d)==0)
            %If link is not activated, set a=1.0
            HystLinks(count,3)=1.0;
        elseif check_alt 
            %Set link parameters
            HystLinks(count,3:end)=MODEL.nl_links_alternate(check_pos,2:end);
        end
        count=count+1;
    end
end

end

