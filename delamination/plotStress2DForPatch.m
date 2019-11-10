function ok  = plotStress2DForPatch(data,ip,vtuFile,U,materials)
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
globElems  = data.mesh{ip}.globElems;
locElems   = data.mesh{ip}.locElems;
uKnot      = data.mesh{ip}.uKnot;
vKnot      = data.mesh{ip}.vKnot;
controlPts = data.mesh{ip}.controlPts;
noPtsX     = data.mesh{ip}.noPtsX;
noPtsY     = data.mesh{ip}.noPtsY;
p          = data.mesh{ip}.p;
q          = data.mesh{ip}.q;
weights    = data.mesh{ip}.weights;
matMap     = data.matMap{ip};

noElems  = size(elementV,1);

stress = zeros(noElems,size(elementV,2),3);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctrg   = globElems(e,:);         %  global element scatter vector
    sctrl   = locElems (e,:);         %  local element scatter vector
    nn      = length(sctrg);
    
    sctrB(1:2:2*nn)    = 2*sctrg-1;
    sctrB(2:2:2*nn)    = 2*sctrg-0;
    
    B      = zeros(3,2*nn);
    pts    = controlPts(sctrl,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1),  uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1), vKnot);
    
    elemDispS = [U(sctrg,1); U(sctrg,2)];
    elemDispU = [U(sctrg,1) U(sctrg,2)];
    
    mat = materials{matMap(e)};
    
    % loop over Gauss points
    
    gp = 1;
    
    for iv=1:2
        Eta  = etaE(iv);
        for iu=1:2
            Xi   = xiE(iu);
            [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi;Eta],...
                p,q,uKnot,vKnot,weights',[uspan;vspan]);
            
            % compute the jacobian of physical and parameter domain mapping
            % then the derivative w.r.t spatial physical coordinates
            
            jacob  = pts' * [dRdxi' dRdeta'];
            %
            if (abs(det(jacob)) <= 1e-6)
                [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi+0.01;Eta],...
                    p,q,uKnot,vKnot,weights',[uspan;vspan]);
                jacob  = pts' * [dRdxi' dRdeta'];
                det(jacob);
            end
            
            % Jacobian inverse and spatial derivatives
            
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            
            % B matrix
            B(1,1:nn)       = dRdx(:,1)';
            B(2,nn+1:2*nn)  = dRdx(:,2)';
            B(3,1:nn)       = dRdx(:,2)';
            B(3,nn+1:2*nn)  = dRdx(:,1)';
            
            Nm(1,1:2:nn*2)  = N';
            Nm(2,2:2:nn*2)  = N';
            
            strain          = B*elemDispS;
            sigma           = mat.stiffMat*strain;
            stress(e,gp,:)  = sigma;
            disp  (e,gp,:)  = N*elemDispU;
            gp = gp +1;
        end
    end % end of gp loops
    
    % disp stored in IGA element connectivity
    % change positions according to standard FE connectivity
    
    col3 = disp(e,3,:);
    col4 = disp(e,4,:);
    
    disp(e,3,:) = col4;
    disp(e,4,:) = col3;
    
end % end of element loop



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

numNode = size(node,1);

% displacements
dispX   = zeros(numNode,1);
dispY   = zeros(numNode,1);


% normal stresses
sigmaXX = zeros(numNode,2);
sigmaYY = zeros(numNode,2);
sigmaZZ = zeros(numNode,2);


for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nid = connect(in);
        sigmaXX(nid,:) = sigmaXX(nid,:) + [stress(e,in,1) 1];
        sigmaYY(nid,:) = sigmaYY(nid,:) + [stress(e,in,2) 1];
        sigmaZZ(nid,:) = sigmaZZ(nid,:) + [stress(e,in,3) 1];
        
        
        dispX(nid) = disp(e,in,1);
        dispY(nid) = disp(e,in,2);
    end
end

% Average nodal stress values (learned from Mathiew Pais XFEM code)
sigmaXX(:,1) = sigmaXX(:,1)./sigmaXX(:,2); sigmaXX(:,2) = [];
sigmaYY(:,1) = sigmaYY(:,1)./sigmaYY(:,2); sigmaYY(:,2) = [];
sigmaZZ(:,1) = sigmaZZ(:,1)./sigmaZZ(:,2); sigmaZZ(:,2) = [];

activeElems = setdiff(1:size(elementV,1),data.voids);

VTKPostProcess(node,elementV(activeElems,:),3,'Quad4',vtuFile,...
    [sigmaXX sigmaYY sigmaZZ],[dispX dispY dispY]);

ok = 1;
