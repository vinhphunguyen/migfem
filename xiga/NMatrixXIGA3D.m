function [Nv] = NMatrixXIGA3D(e,enrich_node,N,pt)

% compute the vector of shape functions of element 'e'
% for extended IGA.
% Nv = [N | enrichment functions]
% Input:
%   e: element under consideration
%   enrich_node: enriched nodes
%   N: NURBS shape functions 
%   dNdxi: derivatives of N w.r.t xi
%
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

global controlPts element elementV xCr xTip levelSets levelSetsB8 node

sctr = element(e,:);
nn   = length(sctr);

% GP in global coord, used

pts        = controlPts(sctr,:);
Gpt        = N * pts;                  
Phi         = levelSets(sctr,:);

% B8 element for level sets

% [N8,dN8dxi] = lagrange_basis('B8',pt);
% sctrB8      = elementV(e,:);
% Phi         = levelSetsB8(sctrB8,:);



% Nfem is always included

%nn3 = 3*nn;

%Nfem = zeros(1,nn);
Nfem = N;

% Switch between non-enriched and enriched elements

if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    Nv = Nfem;
else                               % Enriched elements
    Nxfem = [] ;
    
    % loop on nodes, check node is enriched ...
    
    for in = 1 : nn
        nodeId  = sctr(in);
        enrNoId = enrich_node(nodeId);
        Nin     = N(in);
         
        if ( enrNoId == 1)     % H(x) enriched node
            
            % Enrichment function, H(x) at global Gauss point
            
            dist = Gpt(1,3)-xCr(1,2);
            Hgp  = heaviside(dist);
            
            % Enrichment function, H(x) at node "in"
            
            dist = controlPts(nodeId,3)-xCr(1,2);
            Hi   = heaviside(dist); Hi=0;
            HH   = Hgp - Hi;
            
            aa   = Nin*HH;
            
            % Add to the total Bxfem
            Nxfem = [Nxfem aa];
            
        elseif ( enrNoId == 2) % B(x) enriched node   
            % branch functions in terms of B8 level sets
            
            %[Br] = branchFunctions3d(Phi,N8);                       
            
            % branch functions in terms of NURBS level sets

            [Br] = branchFunctions3d(Phi,N');
                                   
            % compute branch functions at node "in"
            
            %[BrI] = branchNode(levelSets(nodeId,1),levelSets(nodeId,2));
            BrI = [0 0 0 0];
            % components of Nxfem vector
                                                         
            aa1 = Nin*(Br(1)-BrI(1));
            aa2 = Nin*(Br(2)-BrI(2));
            aa3 = Nin*(Br(3)-BrI(3));
            aa4 = Nin*(Br(4)-BrI(4));
          
            Nxfem = [Nxfem aa1 aa2 aa3 aa4];
        end
    end          % end of loop on nodes
    
    % N vector
    
    Nv = [Nfem Nxfem];
    clear Nfem; 
    clear Nxfem;
end  % end of switch between enriched and non-enriched elements

