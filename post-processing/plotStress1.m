
buildVisualizationMesh;
stress = zeros(noElems,size(elementV,2),3);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector    
    nn     = length(sctr);
    
    sctrB(1,1:2:2*nn) = 2*sctr-1;    
    sctrB(1,2:2:2*nn) = 2*sctr  ;
    
    B      = zeros(3,2*nn);
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
    
    % loop over Gauss points
    
    gp = 1;
    for iv=1:2
        if (iv==2)
            xiE = sort(xiE,'descend');
        end
        for iu=1:2            
            Xi  = xiE(iu);
            Eta = etaE(iv);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob      = [dRdxi; dRdeta] * pts;
            
            if (det(jacob) <= 1e-8)
                [N dRdxi dRdeta] = NURBS2DBasisDers([Xi-0.1; ...
                    Eta],p,q,uKnot,vKnot,weights');
                jacob      = [dRdxi; dRdeta] * pts;
            end
            
            % Jacobian inverse and spatial derivatives
            
            invJacob   = inv(jacob);
            dRdx       = invJacob * [dRdxi; dRdeta];
            
            % B matrix
            
            B(1,1:2:2*nn)  = dRdx(1,:);
            B(2,2:2:2*nn)  = dRdx(2,:);
            B(3,1:2:2*nn)  = dRdx(2,:);
            B(3,2:2:2*nn)  = dRdx(1,:);
            
            strain          = B*U(sctrB);
            stress(e,gp,:)  = C*strain;
            disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr)];
            
            gp = gp +1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaXY(nid) = stress(e,in,3);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
    end
end

VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
             [sigmaXX sigmaYY sigmaXY],[dispX dispY]);



