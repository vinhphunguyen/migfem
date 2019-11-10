addpath ../fem_util/;
addpath ../nurbs-geopdes/inst/
addpath ../nurbs-util/
addpath ../meshing/
addpath ../fem-functions/
addpath ../post-processing/
addpath ../iga/
addpath ../nurbs-geopdes/inst/
addpath ../C_files/

clear all
clc

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)

global p q uKnot vKnot weights

controlPts=zeros(4,2,2);
controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [1;0];
controlPts(1:2,1,2) = [0;1];
controlPts(1:2,2,2) = [1;1];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

%% build NURBS object
solid = nrbmak(controlPts,{uKnot,vKnot});

figure
nrbctrlplot(solid);
view([0 90])

% evaluate order 

solid = nrbdegelev(solid,[2 2]); 

uKnot     = cell2mat(solid.knots(1));
vKnot     = cell2mat(solid.knots(2));


refineCountX = 3;
for i=1:refineCountX
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
   
    newKnotsY = [];
    
    % new knots along two directions (uniform)
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    
    newKnots  = {newKnotsX newKnotsY};
    
    % h-refinement
    
    solid     = nrbkntins(solid,newKnots);
    
    uKnot      = cell2mat(solid.knots(1));
    vKnot      = cell2mat(solid.knots(2));
end

figure
nrbctrlplot(solid);
view([0 90])

convert2DNurbs
generateIGA2DMesh
buildVisualizationMesh;

figure
hold on
plot_mesh(node,elementV,'Q4','r-',1.4);
n5 = plot(controlPts(:,1),controlPts(:,2),'go');
set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
axis off  
set(gcf, 'color', 'white');

% Gauss quadrature rule
[W,Q]=quadrature(  6, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

x= [];
jacob=[];

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector
    nn     = length(sctr);
    
    pts    = controlPts(sctr,:);
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        [R,dRdx,J]  = getShapeGrads2D(pt,xiE,etaE,pts);
        x0          = R * pts;
        jacob = [jacob;J];
        x     = [x; x0];
    end
end

jac = zeros(noElems,size(elementV,2),1);


for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);         %  element scatter vector    
    nn     = length(sctr);
   
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
            
           
            jac(e,gp,1)  = det(jacob);                        
            gp = gp +1;
        end
    end
end

figure
plot_field(node,elementV,'Q4',jac);
colorbar
axis off

%plot_mesh(x,tri,'T3','b-',0.9);

%%
%nonlinear parametrization

% controlPts1=zeros(4,3,3);
% controlPts1(1:2,1,1) = [0;0];
% controlPts1(1:2,2,1) = [0.5;0];
% controlPts1(1:2,3,1) = [1;0];
% 
% controlPts1(1:2,1,2) = [0;0.5];
% controlPts1(1:2,2,2) = [0.5;0.5];
% controlPts1(1:2,3,2) = [1;1/2];
% 
% controlPts1(1:2,1,3) = [0;1];
% controlPts1(1:2,2,3) = [1/2;1];
% controlPts1(1:2,3,3) = [1;1];
% 
% controlPts1(4,:,:)   = 1;
% 
% uKnot = [0 0 0 1 1 1];
% vKnot = [0 0 0 1 1 1];
% 
% %% build NURBS object
% solid = nrbmak(controlPts1,{uKnot,vKnot});
% 
% figure
% nrbctrlplot(solid);
% view([0 90])
% 
% % evaluate order 
% 
% solid = nrbdegelev(solid,[1 1]); 
% 
% uKnot     = cell2mat(solid.knots(1));
% vKnot     = cell2mat(solid.knots(2));
% 
% 
% refineCountX = 0;
% for i=1:refineCountX
%     uKnotVectorU = unique(uKnot);
%     uKnotVectorV = unique(vKnot);
%    
%     newKnotsY = [];
%     
%     % new knots along two directions (uniform)
%     
%     newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
%     newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
%     
%     newKnots  = {newKnotsX newKnotsY};
%     
%     % h-refinement
%     
%     solid     = nrbkntins(solid,newKnots);
%     
%     uKnot      = cell2mat(solid.knots(1));
%     vKnot      = cell2mat(solid.knots(2));
% end
% 
% figure
% nrbctrlplot(solid);
% view([0 90])
% 
% figure
% hold on
% plot(controlPts(:,1),controlPts(:,2),'rs');
% plot(controlPts1(:,1),controlPts1(:,2),'r*');

% p = 3;
% q = 3;
% 
% noPtsX = 10;
% noPtsY = 10;
% 
% gcoord=meshRectangularCoord(1,1,noPtsX-1,noPtsY-1);
% controlPts=gcoord;
% 
% weights = ones(1,noPtsX*noPtsY)';
% 
% knotUTemp = linspace(0,1,noPtsX-p+1);
% knotVTemp = linspace(0,1,noPtsY-q+1);
% 
% uKnot = [0 0 0 knotUTemp 1 1 1];
% vKnot = [0 0 0 knotVTemp 1 1 1];
% 
% noCtrPts       = noPtsX * noPtsY;
% noDofs         = noCtrPts * 2;
% 
% % sometimes h-refinement process gives NAN new control pts
% % simply remove them with the following 
% 
% controlPts  = controlPts(1:noCtrPts,:);
% weights     = weights(1:noCtrPts);
% 
% % convert2DNurbs
% generateIGA2DMesh
% buildVisualizationMesh;
% 
% figure
% hold on
% plot_mesh(node,elementV,'Q4','r-',1.4);
% n5 = plot(controlPts(:,1),controlPts(:,2),'go');
% set(n5,'MarkerSize',6, 'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',0.9);
% axis off  
% set(gcf, 'color', 'white');
% 
% % Gauss quadrature rule
% [W,Q]=quadrature(  16, 'GAUSS', 2 ); % noGPs x noGPs point quadrature
% 
% x= [];
% jacob=[];
% 
% for e=1:noElems
%     idu    = index(e,1);
%     idv    = index(e,2);
%     xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
%     etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
%     
%     sctr   = element(e,:);         %  element scatter vector
%     nn     = length(sctr);
%     
%     pts    = controlPts(sctr,:);
%     
%     for gp=1:size(W,1)
%         pt      = Q(gp,:);
%         wt      = W(gp);
%         
%         [R,dRdx,J]  = getShapeGrads2D(pt,xiE,etaE,pts);
%         x0          = R * pts;
%         jacob = [jacob;J];
%         x     = [x; x0];
%     end
% end
% 
% 
% tri = delaunay(x);
% figure
% plot_field(x,tri ,'T3',jacob'); axis('equal');
% colorbar('vert'); axis off;
%plot_mesh(x,tri,'T3','b-',0.9);