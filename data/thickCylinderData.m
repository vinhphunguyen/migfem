% Data for the thick cylinder problem
% Quadratic NURBS for all three directions
% Vinh Phu Nguyen
% Johns Hopkins University

R      = 10;
L      = 15;
t      = 2;
Ri     = R -t;
Rii    = 0.5*(R+Ri);
hL     = 0.5*L;

controlPts=[0 Ri 0;-Ri Ri 0;-Ri 0 0;-Ri -Ri 0;0 -Ri 0;
            0 Rii 0;-Rii Rii 0;-Rii 0 0;-Rii -Rii 0;0 -Rii 0;
            0 R 0;-R R 0;-R 0 0;-R -R 0;0 -R 0;
            ...
            0 Ri hL;-Ri Ri hL;-Ri 0 hL;-Ri -Ri hL;0 -Ri hL;
            0 Rii hL;-Rii Rii hL;-Rii 0 hL;-Rii -Rii hL;0 -Rii hL;
            0 R hL;-R R hL;-R 0 hL;-R -R hL;0 -R hL;
            ...
            0 Ri L;-Ri Ri L;-Ri 0 L;-Ri -Ri L;0 -Ri L;
            0 Rii L;-Rii Rii L;-Rii 0 L;-Rii -Rii L;0 -Rii L;
            0 R L;-R R L;-R 0 L;-R -R L;0 -R L;
           ];

p = 2;
q = 2;
r = 2;

% only one quadratic element along the thickness (q)
% uKnot will be refined using knot insertion

uKnot = [0 0 0 0.5 0.5 1 1 1];
vKnot = [0 0 0 1 1 1];
wKnot = [0 0 0 1 1 1];

noPtsX = length(uKnot)-p-1;
noPtsY = length(vKnot)-q-1;
noPtsZ = length(wKnot)-r-1;

weights = ones(1,noPtsX*noPtsY*noPtsZ)';

fac = 1/sqrt(2);

weights([2,4,7,9,12,14,17,19,22,24,27,29,32,34,37,39,42,44]) = fac;

% refinement along p direction

refineCount = 4;

for i=1:refineCount
    oldControlPts = controlPts;
    
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    
    % new knots along three directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    
    %% h-refinement (NURBS) in x-direction
    dim = size(controlPts,2);
    
    nonewkX      = size(newKnotsX,2);
    newprojcoord = zeros(noPtsX*noPtsY+nonewkX*noPtsY,dim+1);
    
    rstart = 1;
    wstart = 1;
    
    for j=1:noPtsY
        rstop = rstart + noPtsX-1;
        wstop = wstart + noPtsX-1 + nonewkX;
        
        locCP        = controlPts(rstart:rstop,:);
        locweights   = weights   (rstart:rstop);
        locprojcoord = nurb2proj(noPtsX, locCP, locweights);
        
        % refinement of x
        [tempknotVectorX,tempControlPts] = ...
            RefineKnotVectCurve(noPtsX-1,p,uKnot,locprojcoord,newKnotsX,nonewkX-1);
        
        newprojcoord(wstart:wstop,:)=tempControlPts;
        wstart = wstop+1;
        rstart = rstop+1;
    end
    
    uKnot                 = tempknotVectorX;
    [controlPts, weights] = proj2nurbs(newprojcoord);
    noPtsX                = noPtsX+nonewkX;
end

% refinement along r direction

noElemsZ = 7;

knotWTemp = linspace(0,1,noElemsZ+2);
wKnot     = [0 0 knotWTemp 1 1];
noPtsZ = length(wKnot)-r-1;

for i=1:noPtsZ-1
    a = i*L/(noPtsZ-1);
    controlPts=[controlPts;
        [controlPts(1:noPtsX,1:2)          a*ones(noPtsX,1)];
        [controlPts(noPtsX+1:2*noPtsX,1:2) a*ones(noPtsX,1)];
        [controlPts(2*noPtsX+1:3*noPtsX,1:2) a*ones(noPtsX,1)];
        ];
    
    weights=[weights;
        weights(1:noPtsX);weights(1:noPtsX);weights(1:noPtsX);        
        ];
end
%%%

buildVisualization3dMesh

figure
hold on
plot_mesh(node,elementV,'B8','b.-',1.1);
view(3)
% plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');
% 
% numNode = size(node,1);
% 
% % normal stresses
% sigmaXX = zeros(numNode,1);
% sigmaYY = zeros(numNode,1);
% sigmaZZ = zeros(numNode,1);
% 
% % shear stresses
% sigmaXY = zeros(numNode,1);
% sigmaYZ = zeros(numNode,1);
% sigmaZX = zeros(numNode,1);
% 
% % displacements
% dispX   = zeros(numNode,1);
% dispY   = zeros(numNode,1);
% dispZ   = zeros(numNode,1);
% 
% 
% VTKPostProcess3d(node,elementV,'B8',...
%     [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);
