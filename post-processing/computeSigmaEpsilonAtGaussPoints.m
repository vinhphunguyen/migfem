% stresses and strains visualization (for XIGA with material interfaces)
% stresses are computed at Gauss points used in integration of
% the stiffness matrix. 
% For visualization, a triangulation of the GPs in global coordinate system
% is performed. 

gpi = 1;
stress = [];
strain = [];

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector of NURBS mesh
    sctrV  = elementV(e,:);         %  element scatter vector of Q4 FE mesh
    nn     = length(sctr);
    levelS = CHI(1,sctr);
    pts    = controlPts(sctr,:);
   
    elemDisp  = element_disp(e,pos,enrich_node,U);
    
    [W,Q] = gaussForEnrichedElement(e,noGPs,sctrV,levelSets,xTip,...
                                    split_elem, splitElems,tip_elem,...
                                    itip_nodes,iTip_nodes);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt = Q(gp,:);
        
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        
        [R dRdxi dRdeta] = BSPLINE2DBasisDers([Xi; Eta],p,q,uKnot,vKnot);
        
        [B,w] = BMatrixXIGA(Xi,Eta,e,enrich_node,R,dRdxi,dRdeta);
        
        [N,dNdxi] = lagrange_basis('Q4',pt);
        %[B,J1] = BMatrixXIGAGen(Xi,Eta,e,enrich_node,R,dRdxi,dRdeta,N,dNdxi);
        %[exN] = NMatrixXIGA(e,enrich_node,R);
        
        x         = R * pts; % global coord of GP
        levelset  = x(2) - y0;
        
        if (levelset >= 0)
            C = Cm;
        else
            C = Ci;
        end
        
        strainGp        = B*elemDisp;
        sigma           = C*strainGp;
        stress          = [stress sigma];
        strain          = [strain strainGp];
        
        gpi = gpi + 1;
    end
end

tri = delaunay(gps(:,1),gps(:,2));

