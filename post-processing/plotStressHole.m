
buildVisualizationMesh;
stress = zeros(noElems,size(elementV,2),4);
disp   = zeros(noElems,size(elementV,2),2);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    sctrV  = elementV(e,:);        %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
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
            
            jacob      = pts' * [dRdxi' dRdeta'];
            
            % Jacobian inverse and spatial derivatives
            
            invJacob   = inv(jacob);
            dRdx       = [dRdxi' dRdeta'] * invJacob;
            x          = N  * pts; % global coord of GP
            
            if (CHI(sctrV(gp)) >= 0)
                B(1,1:nn)       = dRdx(:,1)';
                B(2,nn+1:2*nn)  = dRdx(:,2)';
                B(3,1:nn)       = dRdx(:,2)';
                B(3,nn+1:2*nn)  = dRdx(:,1)';
                
                strain          = B*U(sctrB);
                sigma           = C*strain;
            else
                sigma           = zeros(3,1);
            end
            
            stress(e,gp,1:3)= sigma;
            
            % von Mises stress
            stress(e,gp,4)  = sqrt(sigma(1)^2+sigma(2)^2-...
                sigma(1)*sigma(2)+3*sigma(3)^2);
            disp(e,gp,:) = N*[Ux(sctr) Uy(sctr)];
            gp           = gp +1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export to VTK format to plot in Mayavi or Paraview

nNodes = size(node,1);

Sxx = zeros(nNodes,2); Sxy = Sxx; Syy = Sxx; Svm = Sxx;

dispX = zeros(size(node,1),1);
dispY = zeros(size(node,1),1);

for e=1:size(elementV,1)
    connect = elementV(e,:);
    for in=1:4
        nNode = connect(in);
                
        Sxx(nNode,:) = Sxx(nNode,:) + [stress(e,in,1) 1];
        Sxy(nNode,:) = Sxy(nNode,:) + [stress(e,in,3) 1];
        Syy(nNode,:) = Syy(nNode,:) + [stress(e,in,2) 1];
        Svm(nNode,:) = Svm(nNode,:) + [stress(e,in,4) 1];
                
        dispX(nNode) = disp(e,in,1);
        dispY(nNode) = disp(e,in,2);
    end
end

% Average nodal stress values
Sxx(:,1) = Sxx(:,1)./Sxx(:,2); Sxx(:,2) = [];
Sxy(:,1) = Sxy(:,1)./Sxy(:,2); Sxy(:,2) = [];
Syy(:,1) = Syy(:,1)./Syy(:,2); Syy(:,2) = [];
Svm(:,1) = Svm(:,1)./Svm(:,2); Svm(:,2) = [];


VTKPostProcess(node,elementV,2,'Quad4',vtuFile,...
    [Sxx(:,1) Syy(:,1) Sxy(:,1) Svm(:,1)],[dispX dispY]);

stressComp=4;
figure
clf
plot_field(node,elementV,'Q4',Svm(:,1));
hold on
colorbar
title('von Mises stress')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

% figure
% clf
% plot_field(node,elementV,'Q4',disp(:,:,2));
% hold on
% colorbar
% title('Displacement in x direction')
% axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)

