% Plot 3D NURBS mesh
% VP Nguyen
% Cardiff University, UK
% nvinhphu@gmail.com
%
% The NURBS mesh is the image of the knot lines
% in the physical space.
%
% Usage
% plotMesh (controlPts,weights,uKnot,vKnot,p,q,90,'r-','try.eps');
%
% It seems that this is obsolete, using NURBS toolbox instead
% nrbkntplot (solid)

function plotMesh3 (controlPts,weights, uKnot,vKnot,wKnot,...
                   p,q,r,resolution,se,fileName)
               
              
noPtsX      = length(uKnot)-p-1;
noPtsY      = length(vKnot)-q-1;
noPtsZ      = length(wKnot)-r-1;

controlPtsX = controlPts(:,1);
controlPtsY = controlPts(:,2);
controlPtsZ = controlPts(:,3);

% discretize the xi, eta and zeta directions

noPts   = resolution;

xiVec    = linspace(0,max(uKnot),noPts);
etaVec   = linspace(0,max(vKnot),noPts);
zetaVec  = linspace(0,max(wKnot),noPts);

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);
wKnotVec = unique(wKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec); 
noKnotsW = length(wKnotVec); 

x1  = zeros(noKnotsU*noKnotsV,noPts);
y1  = zeros(noKnotsU*noKnotsV,noPts);
z1  = zeros(noKnotsU*noKnotsV,noPts);

x2  = zeros(noKnotsV*noKnotsW,noPts);
y2  = zeros(noKnotsV*noKnotsW,noPts);
z2  = zeros(noKnotsV*noKnotsW,noPts);

x3  = zeros(noKnotsW*noKnotsU,noPts);
y3  = zeros(noKnotsW*noKnotsU,noPts);
z3  = zeros(noKnotsW*noKnotsU,noPts);

% NURBS curves of knot lines corresponding to
% xi direction

projcoord = nurb2proj(noPtsX*noPtsY*noPtsZ, controlPts, weights);


dim = 4;

count = 0;
for uk=1:noKnotsU    
    xi = uKnotVec(uk);
    for vk=1:noKnotsV
        count = (uk-1)*noKnotsU + vk;
        eta = vKnotVec(vk);
        for i=1:noPts
            zeta = zetaVec(i);
            tem  = SolidPoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                              noPtsZ-1,r,wKnot,projcoord,dim,xi,eta,zeta);
            x1(count,i) = tem(1)/tem(4);
            y1(count,i) = tem(2)/tem(4);
            z1(count,i) = tem(3)/tem(4);
        end
    end
end

count = 0;
for wk=1:noKnotsW   
    zeta = wKnotVec(wk);
    for vk=1:noKnotsV
        count = (wk-1)*noKnotsW + vk;
        eta = vKnotVec(vk);
        for i=1:noPts
            xi   = xiVec(i);
            tem  = SolidPoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                              noPtsZ-1,r,wKnot,projcoord,dim,xi,eta,zeta);
            x2(count,i) = tem(1)/tem(4);
            y2(count,i) = tem(2)/tem(4);
            z2(count,i) = tem(3)/tem(4);
        end
    end
end

count = 0;
for wk=1:noKnotsW
    zeta = wKnotVec(wk);
    for uk=1:noKnotsU
        count = (wk-1)*noKnotsW + uk;
        xi = uKnotVec(uk);
        for i=1:noPts
            eta  = etaVec(i);
            tem  = SolidPoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                              noPtsZ-1,r,wKnot,projcoord,dim,xi,eta,zeta);
            x3(count,i) = tem(1)/tem(4);
            y3(count,i) = tem(2)/tem(4);
            z3(count,i) = tem(3)/tem(4);
        end
    end
end


hold on
plot3(x1',y1',z1',se,'LineWidth',1.1);
plot3(x2',y2',z2',se,'LineWidth',1.1);
plot3(x3',y3',z3',se,'LineWidth',1.1);
view(3)

% plot3(controlPtsX, controlPtsY, controlPtsZ,'ko',...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',7,'LineWidth',1.2);


axis off 
axis equal
axis tight

%opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)




