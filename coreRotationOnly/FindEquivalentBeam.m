function [beamelement] = FindEquivalentBeam(link_nodes, beam_elements)

[loc,pos] = ismember(link_nodes(1), beam_elements(:,2));

if ~loc
    [loc,pos] = ismember(link_nodes(1), beam_elements(:,3));
end
if ~loc
    [loc,pos] = ismember(link_nodes(2), beam_elements(:,2));
end
if ~loc
    [loc,pos] = ismember(link_nodes(2), beam_elements(:,3));
end

if loc
    beamelement = pos;
else
    beamelement = Nan;
end

end

