% compute stresses and displacements at nodes of
% the visualization mesh for 3D extended IGA.
% Also export this mesh together with stresses
% and displacements to VTK file which can be processed by Paraview.
% Vinh Phu Nguyen
% Delft University of Technology, The Netherlands

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stresses and displacements

stress = zeros(noElems,size(elementV,2),6);
disp   = zeros(noElems,size(elementV,2),3);

xx     = [-1 1];

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    idw    = index(e,3);
    
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    zetaE  = elRangeW(idw,:); % [zeta_k,zeta_k+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    
    pts    = controlPts(sctr,:);
    
    uspan = FindSpan(noPtsX-1,p,xiE(1),  uKnot);
    vspan = FindSpan(noPtsY-1,q,etaE(1), vKnot);
    wspan = FindSpan(noPtsZ-1,r,zetaE(1),wKnot);
    
    % dofs of element 'e' including enriched dofs
    % if 'e' is an enriched element
    
    elemDisp  = elementDisp3d(e,pos,enrich_node,U);
        
    % loop over Gauss points
    
    gp = 1;
    for iw=1:2
        Zeta = zetaE(iw);
        for iv=1:2
            Eta  = etaE(iv);
            for iu=1:2
                Xi   = xiE(iu);
                
                pt   = [xx(iu) xx(iv) xx(iw)];
                
                [N dRdxi dRdeta dRdzeta] = NURBS3DBasisDersSpecial([Xi;Eta;Zeta],...
                    p,q,r,uKnot,vKnot,wKnot,weights',[uspan;vspan;wspan]);
                
                [B,w] = BMatrixXIGA3D(e,enrich_node,N,dRdxi,dRdeta,dRdzeta,pt);
                [exN] = NMatrixXIGA3D(e,enrich_node,N,pt);
               
                strain          = B*elemDisp;
                stress(e,gp,:)  = C*strain;
                
                % the next is wrong
                %disp(e,gp,:)    = N*[Ux(sctr) Uy(sctr) Uz(sctr)];
                
                disp(e,gp,:)    = exN*[elemDisp(1:3:end)...
                                       elemDisp(2:3:end)...
                                       elemDisp(3:3:end)];                
                gp = gp +1;
            end
        end
    end
    
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
    
end

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

for e=size(elementV,1):-1:1
%for e=1:size(elementV,1)
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

VTKPostProcess3d(node,elementV,'B8',vtuFile,...
    [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);



