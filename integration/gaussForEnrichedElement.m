function [W,Q] = gaussForEnrichedElement(e,noGPs,sctrV,levelSets,xTip,...
                                         split_elem, splitElems,tip_elem,...
                                         itip_nodes,iTip_nodes)

global element node
            
sctr   = element(e,:);      

if    ( (ismember(e,split_elem)) || (ismember(e,splitElems)) ) && ...% split element by cracks
        (ismember(e,tip_elem) ~= 1)
    %[W,Q] = discontQ4quad(7,levelSets(1,sctrV,1));
    [W,Q] = quadrature(20,'GAUSS',2);
elseif (ismember(e,tip_elem))   % tip element
    %[W,Q] = disTipQ4quad(7,levelSets(1,sctrV),node(sctrV,:),xTip);
    [W,Q] = quadrature(20,'GAUSS',2);
elseif (any(intersect(itip_nodes,sctr)) ~= 0) || ...
       (any(intersect(iTip_nodes,sctr)) ~= 0) % having tip enriched nodes
    [W,Q] = quadrature(16,'GAUSS',2);
    %     elseif (ismember(e,splitElems))     % split element by mat. interface
    %         [W,Q] = discontQ4quad(3,chi(1,sctrV));
else
    [W,Q] = quadrature(noGPs,'GAUSS',2);
end