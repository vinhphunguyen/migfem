function writeVTK(vtuFile,U,node,elementV)

global element p q r index elRangeU elRangeV elRangeW uKnot vKnot wKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP vIdGrd vSprGK...
    noPtsX noPtsY noPtsZ C

noElems = size(element,1);
stress  = zeros(noElems,size(elementV,2),6);
disp    = zeros(noElems,size(elementV,2),3);

Ux = U(1:3:noDofs);
Uy = U(2:3:noDofs);
Uz = U(3:3:noDofs);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
    
    sctr   = element(e,:);          %  element scatter vector
    
    nn              = length(sctr);    
    sctrB           = zeros(1,3*nn);
    sctrB(1:3:3*nn) = 3*sctr-2;
    sctrB(2:3:3*nn) = 3*sctr-1;
    sctrB(3:3:3*nn) = 3*sctr;
    
    B      = zeros(6,3*nn);
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1),  uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1), vKnot);
    wspan = FindSpan(noPtsZ-1,r,zetaE(1),wKnot);
    
    elemDisp = [Ux(sctr) Uy(sctr) Uz(sctr)];
    
    % loop over Gauss points
    
    gp = 1;
    
    for iw=1:2
        Zeta = zetaE(iw);
        for iv=1:2
            Eta  = etaE(iv);
            for iu=1:2                
                Xi   = xiE(iu);                                                
                [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDersSpecial([Xi;Eta;Zeta],...
                    p,q,r,uKnot,vKnot,wKnot,weights',[uspan;vspan;wspan]);
                
                % compute the jacobian of physical and parameter domain mapping
                % then the derivative w.r.t spatial physical coordinates
                
                jacob  = pts' * [dRdxi' dRdeta' dRdzeta'];
                             
                %             if (det(jacob) <= 1e-8)
                %                 [N dRdxi dRdeta] = NURBS2DBasisDers([Xi-0.01; ...
                %                     Eta],p,q,uKnot,vKnot,weights');
                %                 jacob      = pts' * [dRdxi' dRdeta'];
                %             end
                
                % Jacobian inverse and spatial derivatives
                
                invJacob   = inv(jacob);
                dRdx       = [dRdxi' dRdeta' dRdzeta'] * invJacob;
                
                % B matrix
                
                B = strainDispMatrix3d(nn,dRdx);
                
                strain          = B*U(sctrB);
                stress(e,gp,:)  = C*strain;
                disp  (e,gp,:)  = N*elemDisp;               
                gp = gp +1;
            end
        end
    end % end of gp loops
    
    % disp stored in IGA element connectivity
    % change positions according to standard FE connectivity
    
    col3 = disp(e,3,:);
    col4 = disp(e,4,:);
    col7 = disp(e,7,:);
    col8 = disp(e,8,:);
    
    disp(e,3,:) = col4;
    disp(e,4,:) = col3;
    disp(e,7,:) = col8;
    disp(e,8,:) = col7;
    
end % end of element loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

numNode = size(node,1);

% normal stresses
sigmaXX = zeros(numNode,1);
sigmaYY = zeros(numNode,1);
sigmaZZ = zeros(numNode,1);

% shear stresses
sigmaXY = zeros(numNode,1);
sigmaYZ = zeros(numNode,1);
sigmaZX = zeros(numNode,1);

% displacements
dispX   = zeros(numNode,1);
dispY   = zeros(numNode,1);
dispZ   = zeros(numNode,1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:8
        nid = connect(in);
        sigmaXX(nid) = stress(e,in,1);
        sigmaYY(nid) = stress(e,in,2);
        sigmaZZ(nid) = stress(e,in,3);
        sigmaXY(nid) = stress(e,in,4);
        sigmaYZ(nid) = stress(e,in,5);
        sigmaZX(nid) = stress(e,in,6);
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
        dispZ(nid) = disp(e,in,3);
    end
end

% write to VTS (structured grid)
%
      
mshToVTK (node, [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],...
          [dispX dispY dispZ], vtuFile)
      
% write to VTU file      
 VTKPostProcess3d(node,elementV,'B8',vtuFile,...
     [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);



