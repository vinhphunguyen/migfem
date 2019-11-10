function ok  = plotStress3DForPatch(data,ip,vtuFile,U,damage,materials)
%
% compute stresses and displacements at nodes for patch ip
% Averaged stresses are computed using nodal averaging technique.
%
% U: (noNodes,3) matrix of displacements.
%
% VP Nguyen
% Cardiff University, Wales, UK

% build visualization B8 mesh

vmesh      = data.vmesh{ip};
elementV   = vmesh.element;
node       = vmesh.node;
index      = data.mesh{ip}.index;
elRangeU   = data.mesh{ip}.rangeU;
elRangeV   = data.mesh{ip}.rangeV;
elRangeW   = data.mesh{ip}.rangeW;
globElems  = data.mesh{ip}.globElems;
locElems   = data.mesh{ip}.locElems;
uKnot      = data.mesh{ip}.uKnot;
vKnot      = data.mesh{ip}.vKnot;
wKnot      = data.mesh{ip}.wKnot;
controlPts = data.mesh{ip}.controlPts;
noPtsX     = data.mesh{ip}.noPtsX;
noPtsY     = data.mesh{ip}.noPtsY;
noPtsZ     = data.mesh{ip}.noPtsZ;
p          = data.mesh{ip}.p;
q          = data.mesh{ip}.q;
r          = data.mesh{ip}.r;
weights    = data.mesh{ip}.weights;
matMap      = data.matMap{ip};

noElems  = size(elementV,1);

stress = zeros(noElems,size(elementV,2),8);
disp   = zeros(noElems,size(elementV,2),3);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
    
    sctrg   = globElems(e,:);         %  global element scatter vector
    sctrl   = locElems (e,:);         %  local element scatter vector
    nn      = length(sctrg);
    
    sctrB(1:3:3*nn)    = 3*sctrg-2;
    sctrB(2:3:3*nn)    = 3*sctrg-1;
    sctrB(3:3:3*nn)    = 3*sctrg-0;
    
    B      = zeros(6,3*nn);
    pts    = controlPts(sctrl,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1),  uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1), vKnot);
    wspan = FindSpan(noPtsZ-1,r,zetaE(1),wKnot);
    
    elemDisp = [U(sctrg,1) U(sctrg,2) U(sctrg,3)];
    
    mat = materials{matMap(e)};

    
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
                
                if (abs(det(jacob)) <= 1e-6)
                        [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDersSpecial([Xi+0.01;Eta+0.01;Zeta],...
                     p,q,r,uKnot,vKnot,wKnot,weights',[uspan;vspan;wspan]);
                                        jacob  = pts' * [dRdxi' dRdeta' dRdzeta'];
                    det(jacob);
                end
                
                % Jacobian inverse and spatial derivatives
                
                invJacob   = inv(jacob);
                dRdx       = [dRdxi' dRdeta' dRdzeta'] * invJacob;
                
                % B matrix                                                
                B(1,1:nn)       = dRdx(:,1)';
                B(2,nn+1:2*nn)  = dRdx(:,2)';
                B(3,2*nn+1:end) = dRdx(:,3)';
                
                B(4,1:nn)       = dRdx(:,2)';
                B(4,nn+1:nn*2)  = dRdx(:,1)';
                
                B(5,2*nn+1:end) = dRdx(:,2)';
                B(5,nn+1:nn*2)  = dRdx(:,3)';
                
                B(6,1:nn)       = dRdx(:,3)';
                B(6,2*nn+1:end) = dRdx(:,1)';
                
                strain          = B*[U(sctrg,1); U(sctrg,2); U(sctrg,3)];
                sigma           = mat.stiffMat*strain;
                stress(e,gp,1:6)= sigma;
                stress(e,gp,7)  = N*damage(sctrg);
                stress(e,gp,8)  = sqrt(sigma(1)^2+sigma(2)^2+sigma(3)^2-...
                    sigma(1)*sigma(2)-sigma(2)*sigma(3)-sigma(3)*sigma(1)+...
                    3*(sigma(4)^2+sigma(5)^2+sigma(6)^2));
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

% displacements
dispX   = zeros(numNode,1);
dispY   = zeros(numNode,1);
dispZ   = zeros(numNode,1);

% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaZZ = zeros(numNode,2);

% shear stresses
sigmaXY = zeros(numNode,2);
sigmaYZ = zeros(numNode,2);
sigmaZX = zeros(numNode,2);
% von Mises stress
sigmaVM = zeros(numNode,2);
damage  = zeros(numNode,2);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:8
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaZZ(nid,:) = sigmaZZ(nid,:) + [stress(e,in,3) 1];
        sigmaXY(nid,:) = sigmaXY(nid,:) + [stress(e,in,4) 1];
        sigmaYZ(nid,:) = sigmaYZ(nid,:) + [stress(e,in,5) 1];
        sigmaZX(nid,:) = sigmaZX(nid,:) + [stress(e,in,6) 1];
        damage(nid,:)  = damage(nid,:) + [stress(e,in,7) 1];
        sigmaVM(nid,:) = sigmaVM(nid,:) + [stress(e,in,8) 1];
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
        dispZ(nid) = disp(e,in,3);
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaZZ(:,1) = sigmaZZ(:,1)./sigmaZZ(:,2); sigmaZZ(:,2) = [];
sigmaXY(:,1) = sigmaXY(:,1)./sigmaXY(:,2); sigmaXY(:,2) = [];
sigmaYZ(:,1) = sigmaYZ(:,1)./sigmaYZ(:,2); sigmaYZ(:,2) = [];
sigmaZX(:,1) = sigmaZX(:,1)./sigmaZX(:,2); sigmaZX(:,2) = [];
sigmaVM(:,1) = sigmaVM(:,1)./sigmaVM(:,2); sigmaVM(:,2) = [];
damage(:,1) = damage(:,1)./damage(:,2); damage(:,2) = [];

VTKPostProcess3d(node,elementV,'B8',vtuFile,...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX sigmaVM damage],[dispX dispY dispZ]);

ok = 1;
