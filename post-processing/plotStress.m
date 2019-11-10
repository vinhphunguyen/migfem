% compute stresses at integration points
% extrapolate these values to nodal points 
% todo: nodal averaging
%
% The code does not use this script but plotStress1.m.
% Theory: lecture node IFEM, Prof. Felipa at Colorado.


sqrt3Inv = 1/sqrt(3);
sqrt3    = sqrt(3);
stressPoints=sqrt3Inv*[-1 -1;1 -1;1 1;-1 1];

extrapolateMatrix = [1+0.5*sqrt3 -0.5 1-0.5*sqrt3 -0.5;
                     -0.5 1+0.5*sqrt3 -0.5 1-0.5*sqrt3;
                     1-0.5*sqrt3 -0.5 1+0.5*sqrt3 -0.5;
                     -0.5 1-0.5*sqrt3 -0.5 1+0.5*sqrt3];

buildVisualizationMesh;
stress=zeros(noElems,size(elementV,2),3);

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(3,2*nn);
    pts    = controlPts(sctr,:);
    
    % loop over Gauss points
    for gp=1:4
    %gp = 1;
    %for iv=1:2
     %   if (iv==2)
      %      xiE = sort(xiE,'descend');
       % end
        %for iu=1:2
        
        pt=stressPoints(gp,:);     
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        
        %Xi  = xiE(iu);
        %Eta = etaE(iv);
        
        [dRdxi dRdeta] = NURBS2Dders([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob      = pts' * [dRdxi' dRdeta'];
        xieta=[Xi Eta]
        %N
        %N*pts
        
%         if (det(jacob) <= 1e-8)
%           [N dRdxi dRdeta] = NURBS2DBasisDers([Xi-0.01; ...
%               Eta],p,q,uKnot,vKnot,weights')
%           jacob      = pts' * [dRdxi' dRdeta'];
%         end
       
        % Jacobian inverse and spatial derivatives
        
        invJacob   = inv(jacob);
        dRdx       = [dRdxi' dRdeta'] * invJacob;
        
        % B matrix
        
        B(1,1:nn)       = dRdx(:,1)';
        B(2,nn+1:2*nn)  = dRdx(:,2)';
        B(3,1:nn)       = dRdx(:,2)';
        B(3,nn+1:2*nn)  = dRdx(:,1)';
        
        strain          = B*U(sctrB);
        stress(e,gp,:)  = C*strain; stress(e,:,1);
        
        %gp = gp +1;
        %end
    end
    
    stress(e,:,1)  = extrapolateMatrix*stress(e,:,1)';
    stress(e,:,2)  = extrapolateMatrix*stress(e,:,2)';
    stress(e,:,3)  = extrapolateMatrix*stress(e,:,3)';
    
end

stressComp=1;
figure
clf
plot_field(node,elementV,'Q4',stress(:,:,stressComp));
hold on
colorbar
title('Stress in x direction')
axis off
%plot_mesh(node,elementV,'Q4','g.-');

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)
