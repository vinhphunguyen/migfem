function [fu,k,en] = penaltyBoundaryCondition(elRangeV,elConnV,noElemsV,...
                       bndMesh,bndPoints,E,nu,stressState,...
                       sigmato,xTip,seg,cracklength,pos,mode)

global p q controlPts weights xCrack split_nodes
global uKnot vKnot noDofs levelSets noGPs1

noCtrPts = size(controlPts,1);

fu = zeros(noDofs,1);
k  = zeros(noDofs,noDofs);

[W1,Q1] = quadrature(noGPs1, 'GAUSS', 1 );

S = [1 0 ; 0 1];

for e=1:noElemsV
    xiE   = elRangeV(e,:); % [xi_i,xi_i+1]
    conn  = elConnV(e,:);
    noFns = length(conn);
    sctr  = bndMesh(e,:);
    le    = length(sctr);
    en    = zeros(1,2*le);
    force = zeros(1,2*le);
    
    %     if (length(intersect(sctr,split_nodes))==le)
    %         enrichedElem = 1;
    %         en    = zeros(1,4*le);
    %         force = zeros(1,4*le);
    %     else
    %         enrichedElem = 0;
    %     end
    
    % loop over Gauss points
    for gp=1:size(W1,1)
        xi      = Q1(gp,:);
        wt      = W1(gp);
        Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * xi + xiE(2) + xiE(1));
        J2      = 0.5 * ( xiE(2) - xiE(1) );
        
        [N dNdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights);
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        x        = N    *bndPoints(conn,:); % global coord of GP
        jacob1   = dNdxi*bndPoints(conn,:);
        J1       = norm (jacob1);

        if mode == 1            
            [ux,uy] = exactDispModeI(x,E,nu,stressState,...
                                     sigmato,xTip,seg,cracklength);
        else
            
            [ux,uy] = exactDispModeII(x,E,nu,stressState,...
                                      sigmato,xTip,seg,cracklength);
        end
        
        for j = 1 : le
            en(2*j-1) = 2*sctr(j)-1;
            en(2*j  ) = 2*sctr(j)  ;
            
            force(2*j-1) = S(1,1)*N(j)*ux + S(1,2)*N(j)*uy;
            force(2*j  ) = S(2,1)*N(j)*ux + S(2,2)*N(j)*uy;
        end
        
        fu(en(1:2*le)) = fu(en(1:2*le)) + J1 * J2 * wt * force(1:2*le)';
        
        %         if (enrichedElem==1)
        %             for in = 1 : le
        %                 nodeId = sctr(in)
        %                 Ni     = N(in);
        %                 xCr    = reshape(xCrack(1,:,:),2,2);
        %
        %                 % Enrichment function, H(x) at global Gauss point
        %                 dist = signed_distance(xCr,x);
        %                 Hgp  = heaviside(dist);
        %
        %                 % Enrichment function, H(x) at node "in"
        %                 dist = signed_distance(xCr,controlPts(nodeId,:));
        %                 Hi   = heaviside(dist);
        %
        %                 en(2*(in+le)-1) = 2*pos(nodeId)-1;
        %                 en(2*(in+le)  ) = 2*pos(nodeId)  ;
        %
        %                 force(2*(in+le)-1) = S(1,1)*Ni*(Hgp-Hi)*ux + ...
        %                     S(1,2)*Ni*(Hgp-Hi)*uy;
        %                 force(2*(in+le)  ) = S(2,1)*Ni*(Hgp-Hi)*ux + ...
        %                     S(2,2)*Ni*(Hgp-Hi)*uy;
        %             end
        %
        %             fu(en(2*le+1:4*le)) = fu(en(2*le+1:4*le)) + ...
        %                                 J1 * J2 * wt * force(2*le+1:4*le)';
        %
        %             for i = 1:le
        %                 nodeI = sctr(i);
        %                 dist = signed_distance(xCr,controlPts(nodeI,:));
        %                 Hi   = heaviside(dist);
        %                 Hi   = Hgp - Hi;
        %
        %                 id=[2*pos(nodeI)-1:2*pos(nodeI)];
        %
        %                 for j=1:le
        %                     nodeJ = sctr(j);
        %                     dist = signed_distance(xCr,controlPts(nodeJ,:));
        %                     Hj   = heaviside(dist);
        %                     Hj   = Hgp - Hj;
        %
        %                     jd=[2*pos(nodeJ)-1:2*pos(nodeJ)];
        %
        %                     k(id,jd) = k(id,jd) + J1 * J2 * wt * N(i)*Hi*N(j)*Hj*S ;
        %                 end
        %             end
        %
        %         end
        
        for i = 1:le
            id=[2*sctr(i)-1:2*sctr(i)];
            for j=1:le
                jd=[2*sctr(j)-1:2*sctr(j)];
                k(id,jd) = k(id,jd) + J1 * J2 * wt * N(i)*N(j)*S ;
            end
        end
    end % end of loop over Gauss points
end % end of loop over boundary elements
