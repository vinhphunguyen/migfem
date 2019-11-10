% Plot 2D NURBS mesh
% VP Nguyen
% Cardiff University, UK
% nvinhphu@gmail.com
%
% The NURBS mesh is the image of the knot lines
% in the physical space.
%
% Usage
% plotMesh (controlPts,weights,uKnot,vKnot,p,q,90,'r-','try.eps');

function plotMesh (controlPts,weights, uKnot,vKnot,...
                   p,q,resolution, se, fileName)
               
              
noPtsX      = length(uKnot)-p-1;
noPtsY      = length(vKnot)-q-1;

controlPtsX = controlPts(:,1);
controlPtsY = controlPts(:,2);

% discretize the xi and eta directions

noPts   = resolution;
xiVec   = linspace(0,max(uKnot),noPts);
etaVec  = linspace(0,max(vKnot),noPts);

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec); %noKnotsV=1;

x1  = zeros(noKnotsU,noPts);
y1  = zeros(noKnotsU,noPts);

x2  = zeros(noKnotsV,noPts);
y2  = zeros(noKnotsV,noPts);

% NURBS curves of knot lines corresponding to
% xi direction

projcoord = nurb2proj(noPtsX*noPtsY, controlPts, weights);


dim=size(projcoord,2);

for uk=1:noKnotsU
    xi = uKnotVec(uk);
    for i=1:noPts
        eta = etaVec(i);
        tem = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
            projcoord,dim,xi,eta);
        x1(uk,i) = tem(1)/tem(3);
        y1(uk,i) = tem(2)/tem(3);
    end
end

for vk=1:noKnotsV
    eta = vKnotVec(vk);
    
    for i=1:noPts
        xi  = xiVec(i);
        tem =  SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
            projcoord,dim,xi,eta);
        
        x2(vk,i) = tem(1)/tem(3);
        y2(vk,i) = tem(2)/tem(3);
    end
end

%figure
hold on
plot(x1',y1',se,'LineWidth',1.2);
plot(x2',y2',se,'LineWidth',1.2);

for i=1:noPtsY
   idx = (i-1)*noPtsX+1:i*noPtsX; 
   plot(controlPtsX(idx), controlPtsY(idx),'bo',...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',9,'LineWidth',1.2);
end
% 
% for i=1:noPtsX
%    idx = i:noPtsX:i+noPtsX*(noPtsY-1);
%    plot(controlPtsX(idx), controlPtsY(idx),'bo',...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',9,'LineWidth',1.2);
% end

axis off 
axis equal
axis tight
set(gcf,'color','white')

% opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
% exportfig(gcf,fileName,opts)




